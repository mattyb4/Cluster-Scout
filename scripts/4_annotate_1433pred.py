"""Annotate the PTM proximity database with 14-3-3-Pred binding-site predictions.

For each Ser/Thr PTM site in the output, queries the 14-3-3-Pred API
(https://www.compbio.dundee.ac.uk/1433pred/) and adds two columns:

  1433pred_binding_site  — "Yes" if Consensus > 0, "No" if <= 0, blank if
                           the site is not Ser/Thr or has no prediction.
  1433pred_consensus     — The raw consensus score (float string), blank
                           where no prediction is available.

API responses are cached per UniProt ID under data/cache/1433pred/ so
subsequent runs do not re-query the server.
"""
from __future__ import annotations

import json
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import pandas as pd
import requests
from tqdm import tqdm

sys.path.insert(0, str(Path(__file__).resolve().parent))
from pipeline_utils import project_root, SITE_RE  # noqa: E402

PROJECT_ROOT = project_root(__file__)
PROXIMITY_DB = PROJECT_ROOT / "Output" / "ptm_mutation_proximity_db.tsv"
CACHE_DIR = PROJECT_ROOT / "data" / "cache" / "1433pred"
CONFIRMED_SITES_FILE = PROJECT_ROOT / "data" / "14-3-3 interactors with known P sites.xlsx"

_API_URL = "https://www.compbio.dundee.ac.uk/1433pred/pid={uid}&out=json"
_MAX_WORKERS = 5


# ── API / cache helpers ───────────────────────────────────────────────────────

def fetch_1433pred(uniprot_id: str) -> list[dict] | None:
    """Return the raw 14-3-3-Pred JSON for *uniprot_id*, hitting the cache first."""
    cache_file = CACHE_DIR / f"{uniprot_id}.json"
    if cache_file.exists():
        with cache_file.open(encoding="utf-8") as f:
            return json.load(f)

    try:
        resp = requests.get(_API_URL.format(uid=uniprot_id), timeout=30)
    except requests.RequestException:
        return None

    if resp.status_code != 200:
        return None

    try:
        data = resp.json()
    except ValueError:
        return None

    if not isinstance(data, list):
        return None

    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    with cache_file.open("w", encoding="utf-8") as f:
        json.dump(data, f)

    return data


def build_site_score_map(data: list[dict]) -> dict[int, float]:
    """Convert a 14-3-3-Pred response into a {position: consensus_score} dict."""
    scores: dict[int, float] = {}
    for entry in data:
        try:
            pos = int(entry["Site"])
            score = float(entry["Consensus"])
            scores[pos] = score
        except (KeyError, ValueError, TypeError):
            continue
    return scores


# ── Annotation logic ──────────────────────────────────────────────────────────

def annotate_row(ptm_site: str, score_map: dict[int, float]) -> tuple[str, str]:
    """Return (1433pred_binding_site, 1433pred_consensus) for one PTM-site row.

    Both values are empty strings when the site is not Ser/Thr or has no
    prediction available for that position.
    """
    if not isinstance(ptm_site, str):
        return "", ""

    m = SITE_RE.match(ptm_site.strip())
    if not m:
        return "", ""

    residue, position = m.group(1), int(m.group(2))

    if residue not in ("S", "T"):
        return "", ""

    if position not in score_map:
        return "", ""

    score = score_map[position]
    binding = "Yes" if score > 0 else "No"
    return binding, str(round(score, 3))


# ── Confirmed-site lookup ─────────────────────────────────────────────────────

def load_confirmed_sites(path: Path) -> dict[tuple[str, int], str]:
    """Load the known 14-3-3 interactor dataset.

    Returns {(uniprot_id, position): pmid} for all valid S/T entries.
    """
    if not path.exists():
        print(f"Warning: confirmed-sites file not found: {path}", file=__import__("sys").stderr)
        return {}

    df = pd.read_excel(path, dtype=str)
    df.columns = [c.strip() for c in df.columns]
    df["Residue"] = df["Residue"].str.strip()
    df["PMID"] = df["PMID"].str.strip().str.lstrip("\xa0")

    confirmed: dict[tuple[str, int], str] = {}
    for _, row in df.iterrows():
        uid = str(row.get("Uniprot ID", "")).strip()
        residue = str(row.get("Residue", "")).strip()
        pmid = str(row.get("PMID", "")).strip()
        try:
            pos = int(float(str(row.get("Site", "")).strip()))
        except (ValueError, TypeError):
            continue
        if residue not in ("S", "T"):
            continue
        confirmed[(uid, pos)] = pmid

    return confirmed


def annotate_confirmed(uid: str, ptm_site: str, confirmed: dict[tuple[str, int], str]) -> tuple[str, str]:
    """Return (1433_confirmed_site, 1433_confirmed_pmid) for one row."""
    m = SITE_RE.match(ptm_site.strip()) if isinstance(ptm_site, str) else None
    if not m:
        return "", ""
    pos = int(m.group(2))
    pmid = confirmed.get((uid, pos), "")
    if pmid:
        return "Yes", pmid
    return "", ""


# ── Main ──────────────────────────────────────────────────────────────────────

def main() -> None:
    """Annotate the proximity DB with 14-3-3 predicted and confirmed binding-site columns."""
    confirmed_sites = load_confirmed_sites(CONFIRMED_SITES_FILE)
    print(f"Loaded {len(confirmed_sites)} confirmed 14-3-3 binding sites")

    print(f"Reading proximity DB: {PROXIMITY_DB}")
    df = pd.read_csv(PROXIMITY_DB, sep="\t", encoding="utf-16", dtype=str,
                     keep_default_na=False)

    unique_uniprots = df["UniProt"].dropna().unique().tolist()

    # Fetch predictions — only hit the API for IDs not already cached.
    already_cached = sum(
        1 for uid in unique_uniprots if (CACHE_DIR / f"{uid}.json").exists()
    )
    print(
        f"{already_cached}/{len(unique_uniprots)} UniProt IDs already cached; "
        f"fetching {len(unique_uniprots) - already_cached} new..."
    )

    score_maps: dict[str, dict[int, float]] = {}
    with ThreadPoolExecutor(max_workers=_MAX_WORKERS) as pool:
        futures = {pool.submit(fetch_1433pred, uid): uid for uid in unique_uniprots}
        for future in tqdm(as_completed(futures), total=len(futures), desc="Fetching 14-3-3-Pred data"):
            uid = futures[future]
            data = future.result()
            if data is not None:
                score_maps[uid] = build_site_score_map(data)

    binding_sites, consensus_scores = [], []
    confirmed_col, confirmed_pmid_col = [], []
    for _, row in df.iterrows():
        uid = row.get("UniProt", "")
        ptm_site = row.get("ptm_site", "")

        smap = score_maps.get(uid, {})
        binding, score = annotate_row(ptm_site, smap)
        binding_sites.append(binding)
        consensus_scores.append(score)

        conf, pmid = annotate_confirmed(uid, ptm_site, confirmed_sites)
        confirmed_col.append(conf)
        confirmed_pmid_col.append(pmid)

    df["1433pred_binding_site"] = binding_sites
    df["1433pred_consensus"] = consensus_scores
    df["1433_confirmed_site"] = confirmed_col
    df["1433_confirmed_pmid"] = confirmed_pmid_col

    predicted = sum(1 for b in binding_sites if b == "Yes")
    not_predicted = sum(1 for b in binding_sites if b == "No")
    n_confirmed = sum(1 for c in confirmed_col if c == "Yes")
    print(
        f"Annotated {len(df)} rows: "
        f"{predicted} predicted binding sites, {not_predicted} not predicted, "
        f"{len(df) - predicted - not_predicted} not applicable (non-Ser/Thr or no data).\n"
        f"{n_confirmed} experimentally confirmed 14-3-3 binding sites."
    )

    df.to_csv(PROXIMITY_DB, sep="\t", index=False, encoding="utf-16")
    print(f"Updated proximity DB written to: {PROXIMITY_DB}")


if __name__ == "__main__":
    main()
