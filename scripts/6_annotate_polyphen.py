"""Annotate nearby/distant mutation strings with PolyPhen-2 predictions.

Queries the myvariant.info dbNSFP data to add a (PP:D/P/B) tag to each
mutation entry in the proximity database:

  D = Probably Damaging
  P = Possibly Damaging
  B = Benign

Mutations already tagged or not found in dbNSFP are left unchanged.
Results are cached in data/cache/polyphen.tsv so subsequent runs skip
already-queried (gene, mutation) pairs.
"""
from __future__ import annotations

import re
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import pandas as pd
import requests
from tqdm import tqdm

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent
PROXIMITY_DB = PROJECT_ROOT / "Output" / "ptm_mutation_proximity_db.tsv"
CACHE_FILE = PROJECT_ROOT / "data" / "cache" / "polyphen.tsv"

_API_URL = "https://myvariant.info/v1/query"
_MAX_WORKERS = 10

# Columns that contain formatted mutation strings to be re-tagged
_MUTATION_COLS = ["mutations_within_5_positions", "mutations_more_than_5_positions"]

# Parses one entry: base mutation + existing tags + distance/PAE info
# e.g. "R175H(isoform?)-3.52Å(PAE:2.1)"  →  groups: "R175H", "(isoform?)", "-3.52Å(PAE:2.1)"
_ENTRY_RE = re.compile(r"^([A-Z]\d+[A-Z*])((?:\([^)]+\))*)(-.+)$")

_MUT_RE = re.compile(r"^([A-Z])(\d+)([A-Z*])$")

# Severity order for picking the worst prediction when multiple transcripts exist
_SEVERITY = {"D": 2, "P": 1, "B": 0}


# ── Cache helpers ─────────────────────────────────────────────────────────────

def _load_cache() -> dict[tuple[str, str], tuple[str, str]]:
    """Return {(gene, mutation): (pred, score)} from the TSV cache."""
    if not CACHE_FILE.exists():
        return {}
    df = pd.read_csv(CACHE_FILE, sep="\t", dtype=str, keep_default_na=False)
    return {(row["gene"], row["mutation"]): (row["pred"], row["score"])
            for _, row in df.iterrows()}


def _save_cache(cache: dict[tuple[str, str], tuple[str, str]]) -> None:
    CACHE_FILE.parent.mkdir(parents=True, exist_ok=True)
    rows = [{"gene": g, "mutation": m, "pred": p, "score": s}
            for (g, m), (p, s) in cache.items()]
    pd.DataFrame(rows, columns=["gene", "mutation", "pred", "score"]).to_csv(
        CACHE_FILE, sep="\t", index=False
    )


# ── API fetch ─────────────────────────────────────────────────────────────────

def _best_prediction(hits: list[dict]) -> tuple[str, str]:
    """Return the most severe HDIV (pred, score) across all transcript hits, or ('', '')."""
    best_pred = ""
    best_score = -1.0

    for hit in hits:
        hdiv = hit.get("dbnsfp", {}).get("polyphen2", {}).get("hdiv", {})
        preds = hdiv.get("pred", [])
        scores = hdiv.get("score", [])

        if isinstance(preds, str):
            preds = [preds]
        if isinstance(scores, (int, float)):
            scores = [scores]

        for pred, score in zip(preds, scores):
            if pred not in _SEVERITY:
                continue
            if _SEVERITY[pred] > _SEVERITY.get(best_pred, -1):
                best_pred = pred
                best_score = float(score) if isinstance(score, (int, float)) else -1.0

    if not best_pred:
        return "", ""
    score_str = f"{best_score:.3f}" if best_score >= 0 else ""
    return best_pred, score_str


def fetch_polyphen(gene: str, mutation: str) -> tuple[str, str]:
    """Query myvariant.info for PolyPhen-2 HDIV prediction.

    Returns (pred, score) where pred is 'D', 'P', 'B', or '' if not found.
    Empty strings are also cached to avoid re-querying missing variants.
    """
    m = _MUT_RE.match(mutation)
    if not m:
        return "", ""
    ref, pos, alt = m.group(1), m.group(2), m.group(3)
    if alt == "*":
        return "", ""

    q = (f"dbnsfp.genename:{gene} AND dbnsfp.aa.ref:{ref} "
         f"AND dbnsfp.aa.alt:{alt} AND dbnsfp.aa.pos:{pos}")
    try:
        resp = requests.get(
            _API_URL,
            params={"q": q, "fields": "dbnsfp.polyphen2,dbnsfp.aa", "size": 10},
            timeout=15,
        )
        resp.raise_for_status()
        return _best_prediction(resp.json().get("hits", []))
    except Exception:
        return "", ""


# ── String annotation ─────────────────────────────────────────────────────────

def annotate_mutation_string(
    mutation_str: str,
    gene: str,
    cache: dict[tuple[str, str], tuple[str, str]],
) -> str:
    """Insert (PP:X) tags into a formatted mutation string for one row."""
    if not mutation_str or not mutation_str.strip():
        return mutation_str

    parts = []
    for entry in mutation_str.split(", "):
        m = _ENTRY_RE.match(entry.strip())
        if not m:
            parts.append(entry)
            continue

        base = m.group(1)        # e.g. R175H
        existing = m.group(2)    # e.g. (isoform?) or ""
        rest = m.group(3)        # e.g. -3.52Å(PAE:2.1)

        if "(PP:" in existing:   # already tagged, leave alone
            parts.append(entry)
            continue

        pred, score = cache.get((gene, base), ("", ""))
        pp_tag = f"(PP:{pred},{score})" if pred else ""
        parts.append(f"{base}{existing}{pp_tag}{rest}")

    return ", ".join(parts)


# ── Main ──────────────────────────────────────────────────────────────────────

def main() -> None:
    print(f"Reading proximity DB: {PROXIMITY_DB}")
    df = pd.read_csv(PROXIMITY_DB, sep="\t", encoding="utf-16", dtype=str,
                     keep_default_na=False)

    cache = _load_cache()

    # Collect all unique (gene, mutation) pairs that need a lookup
    needed: set[tuple[str, str]] = set()
    for _, row in df.iterrows():
        gene = row.get("gene", "")
        if not gene:
            continue
        for col in _MUTATION_COLS:
            cell = row.get(col, "")
            if not cell:
                continue
            for entry in cell.split(", "):
                m = _ENTRY_RE.match(entry.strip())
                if not m:
                    continue
                base = m.group(1)
                if "(PP:" not in m.group(2) and (gene, base) not in cache:
                    needed.add((gene, base))

    already_cached = sum(1 for k in needed if k in cache)
    to_fetch = [(g, mut) for g, mut in needed if (g, mut) not in cache]
    print(f"{len(needed)} unique (gene, mutation) pairs — "
          f"{len(cache) - already_cached + len(to_fetch)} new to fetch, "
          f"{len(cache)} already cached")

    if to_fetch:
        with ThreadPoolExecutor(max_workers=_MAX_WORKERS) as pool:
            futures = {pool.submit(fetch_polyphen, g, mut): (g, mut)
                       for g, mut in to_fetch}
            for future in tqdm(as_completed(futures), total=len(futures),
                               desc="Fetching PolyPhen-2 scores"):
                g, mut = futures[future]
                pred, score = future.result()
                cache[(g, mut)] = (pred, score)

        _save_cache(cache)
        print("Cache saved.")

    # Re-annotate mutation strings
    for col in _MUTATION_COLS:
        if col not in df.columns:
            continue
        df[col] = [
            annotate_mutation_string(row[col], row.get("gene", ""), cache)
            for _, row in df.iterrows()
        ]

    tagged = sum(
        1 for _, row in df.iterrows()
        for col in _MUTATION_COLS
        if "(PP:" in str(row.get(col, ""))
    )
    print(f"Tagged {tagged} mutation entries with PolyPhen-2 predictions.")

    df.to_csv(PROXIMITY_DB, sep="\t", index=False, encoding="utf-16")
    print(f"Updated proximity DB written to: {PROXIMITY_DB}")


if __name__ == "__main__":
    main()
