"""Annotate the proximity database with 14-3-3 binding predictions, PolyPhen-2
pathogenicity scores, and predicted upstream kinases.

Reads the proximity DB once, runs three annotation phases, and writes back once:

  Phase 1 — 14-3-3: Queries the 14-3-3-Pred API for predicted binding-site
            scores and cross-references experimentally confirmed interactors.
  Phase 2 — PolyPhen-2: Queries myvariant.info for pathogenicity predictions
            and tags each mutation with a (PP:D/P/B,score) label.
  Phase 3 — Kinases: Uses the Kinase Library to predict the top 5 upstream
            kinases for each phosphorylation site based on the ±7 residue
            sequence window from AlphaFold CIF structures.
"""
from __future__ import annotations

import json
import re
import sys
import warnings
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import pandas as pd
import requests
from tqdm import tqdm

sys.path.insert(0, str(Path(__file__).resolve().parent))
from pipeline_utils import (  # noqa: E402
    project_root, AA3TO1, MUT_RE, SITE_RE,
    find_canonical_cif, load_first_chain,
    input_dir, resolve_input_file, INTERACTORS_1433_INPUT_DIR,
)

PROJECT_ROOT = project_root(__file__)
MODELS_ROOT = PROJECT_ROOT / "cif_models"
DEFAULT_OUTPUT_DIR = PROJECT_ROOT / "Output"

_NUM_PHASES = 3


def _emit_progress(phase: int, phase_pct: float, desc: str) -> None:
    """Print overall progress for the app to parse. phase is 0-indexed, phase_pct is 0-100."""
    overall = int((phase * 100 + phase_pct) / _NUM_PHASES)
    print(f"\r##PROGRESS## {overall} {desc}", end="", flush=True)

# ═════════════════════════════════════════════════════════════════════════════
# Phase 1: 14-3-3-Pred binding-site predictions + confirmed interactors
# ═════════════════════════════════════════════════════════════════════════════

_1433_CACHE_DIR = PROJECT_ROOT / "data" / "cache" / "1433pred"
_1433_API_URL = "https://www.compbio.dundee.ac.uk/1433pred/pid={uid}&out=json"
_1433_MAX_WORKERS = 5


def fetch_1433pred(uniprot_id: str) -> list[dict] | None:
    """Return the raw 14-3-3-Pred JSON for *uniprot_id*, hitting the cache first."""
    cache_file = _1433_CACHE_DIR / f"{uniprot_id}.json"
    if cache_file.exists():
        with cache_file.open(encoding="utf-8") as f:
            return json.load(f)
    try:
        resp = requests.get(_1433_API_URL.format(uid=uniprot_id), timeout=30)
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
    _1433_CACHE_DIR.mkdir(parents=True, exist_ok=True)
    with cache_file.open("w", encoding="utf-8") as f:
        json.dump(data, f)
    return data


def build_site_score_map(data: list[dict]) -> dict[int, float]:
    """Convert a 14-3-3-Pred response into a {position: consensus_score} dict."""
    scores: dict[int, float] = {}
    for entry in data:
        try:
            scores[int(entry["Site"])] = float(entry["Consensus"])
        except (KeyError, ValueError, TypeError):
            continue
    return scores


def annotate_1433_row(ptm_site: str, score_map: dict[int, float]) -> tuple[str, str]:
    """Return (binding_site, consensus) for one PTM site row."""
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
    return ("Yes" if score > 0 else "No"), str(round(score, 3))


def load_confirmed_sites(path: Path) -> dict[tuple[str, int], str]:
    """Load known 14-3-3 interactors. Returns {(uniprot_id, position): pmid}."""
    if not path.exists():
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


def annotate_confirmed(uid: str, ptm_site: str,
                       confirmed: dict[tuple[str, int], str]) -> tuple[str, str]:
    """Return (confirmed_site, pmid) for one row."""
    m = SITE_RE.match(ptm_site.strip()) if isinstance(ptm_site, str) else None
    if not m:
        return "", ""
    pos = int(m.group(2))
    pmid = confirmed.get((uid, pos), "")
    return ("Yes", pmid) if pmid else ("", "")


def run_1433_phase(df: pd.DataFrame) -> None:
    """Phase 1: fetch 14-3-3 predictions and confirmed sites, add columns to df."""
    print("\n── Phase 1: 14-3-3 binding-site predictions ──")

    interactors_dir = input_dir(PROJECT_ROOT, INTERACTORS_1433_INPUT_DIR)
    try:
        confirmed_file = resolve_input_file(interactors_dir, (".xlsx", ".xls"))
        confirmed_sites = load_confirmed_sites(confirmed_file)
        print(f"Loaded {len(confirmed_sites)} confirmed 14-3-3 binding sites")
    except FileNotFoundError:
        confirmed_sites = {}
        print("No 14-3-3 interactors file found — skipping confirmed-site annotation")

    unique_uniprots = df["UniProt"].dropna().unique().tolist()
    already_cached = sum(
        1 for uid in unique_uniprots if (_1433_CACHE_DIR / f"{uid}.json").exists()
    )
    print(f"{already_cached}/{len(unique_uniprots)} UniProt IDs already cached; "
          f"fetching {len(unique_uniprots) - already_cached} new...")

    score_maps: dict[str, dict[int, float]] = {}
    with ThreadPoolExecutor(max_workers=_1433_MAX_WORKERS) as pool:
        futures = {pool.submit(fetch_1433pred, uid): uid for uid in unique_uniprots}
        done = 0
        total = len(futures)
        for future in tqdm(as_completed(futures), total=total,
                           desc="Fetching 14-3-3-Pred data"):
            uid = futures[future]
            data = future.result()
            if data is not None:
                score_maps[uid] = build_site_score_map(data)
            done += 1
            _emit_progress(0, done / total * 100, f"14-3-3 predictions: {done}/{total}")

    binding_sites, consensus_scores = [], []
    confirmed_col, confirmed_pmid_col = [], []
    for _, row in df.iterrows():
        uid = row.get("UniProt", "")
        ptm_site = row.get("ptm_site", "")
        binding, score = annotate_1433_row(ptm_site, score_maps.get(uid, {}))
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
    n_confirmed = sum(1 for c in confirmed_col if c == "Yes")
    print(f"  {predicted} predicted binding sites, {n_confirmed} experimentally confirmed")


# ═════════════════════════════════════════════════════════════════════════════
# Phase 2: PolyPhen-2 mutation pathogenicity scores
# ═════════════════════════════════════════════════════════════════════════════

_PP_CACHE_FILE = PROJECT_ROOT / "data" / "cache" / "polyphen.tsv"
_PP_API_URL = "https://myvariant.info/v1/query"
_PP_MAX_WORKERS = 10
_PP_SEVERITY = {"D": 2, "P": 1, "B": 0}

_MUTATION_COLS = ["mutations_within_5_positions", "mutations_more_than_5_positions"]
_ENTRY_RE = re.compile(r"^([A-Z]\d+[A-Z*])((?:\([^)]+\))*)(-.+)$")


def _pp_load_cache() -> dict[tuple[str, str], tuple[str, str]]:
    """Return {(gene, mutation): (pred, score)} from the TSV cache."""
    if not _PP_CACHE_FILE.exists():
        return {}
    df = pd.read_csv(_PP_CACHE_FILE, sep="\t", dtype=str, keep_default_na=False)
    return {(row["gene"], row["mutation"]): (row["pred"], row["score"])
            for _, row in df.iterrows()}


def _pp_save_cache(cache: dict[tuple[str, str], tuple[str, str]]) -> None:
    """Persist the PolyPhen cache to disk."""
    _PP_CACHE_FILE.parent.mkdir(parents=True, exist_ok=True)
    rows = [{"gene": g, "mutation": m, "pred": p, "score": s}
            for (g, m), (p, s) in cache.items()]
    pd.DataFrame(rows, columns=["gene", "mutation", "pred", "score"]).to_csv(
        _PP_CACHE_FILE, sep="\t", index=False)


def _pp_best_prediction(hits: list[dict]) -> tuple[str, str]:
    """Return the most severe HDIV (pred, score) across all hits."""
    best_pred, best_score = "", -1.0
    for hit in hits:
        dbnsfp = hit.get("dbnsfp", {})
        polyphen2 = dbnsfp.get("polyphen2", {})
        if isinstance(polyphen2, list):
            pp2_entries = polyphen2
        elif isinstance(polyphen2, dict):
            pp2_entries = [polyphen2]
        else:
            continue
        for pp2 in pp2_entries:
            hdiv = pp2.get("hdiv", {}) if isinstance(pp2, dict) else {}
            preds = hdiv.get("pred", [])
            scores = hdiv.get("score", [])
            if isinstance(preds, str):
                preds = [preds]
            if isinstance(scores, (int, float)):
                scores = [scores]
            for pred, score in zip(preds, scores):
                if pred not in _PP_SEVERITY:
                    continue
                if _PP_SEVERITY[pred] > _PP_SEVERITY.get(best_pred, -1):
                    best_pred = pred
                    best_score = float(score) if isinstance(score, (int, float)) else -1.0
    if not best_pred:
        return "", ""
    return best_pred, (f"{best_score:.3f}" if best_score >= 0 else "")


def fetch_polyphen(gene: str, mutation: str) -> tuple[str, str]:
    """Query myvariant.info for PolyPhen-2 HDIV prediction."""
    m = MUT_RE.match(mutation)
    if not m:
        return "", ""
    ref, pos, alt = m.group(1), m.group(2), m.group(3)
    if alt == "*":
        return "", ""
    q = (f"dbnsfp.genename:{gene} AND dbnsfp.aa.ref:{ref} "
         f"AND dbnsfp.aa.alt:{alt} AND dbnsfp.aa.pos:{pos}")
    try:
        resp = requests.get(
            _PP_API_URL,
            params={"q": q, "fields": "dbnsfp.polyphen2,dbnsfp.aa", "size": 10},
            timeout=15,
        )
        resp.raise_for_status()
        return _pp_best_prediction(resp.json().get("hits", []))
    except Exception:
        return "", ""


def annotate_mutation_string(
    mutation_str: str, gene: str,
    cache: dict[tuple[str, str], tuple[str, str]],
) -> str:
    """Insert (PP:X,score) tags into a formatted mutation string."""
    if not mutation_str or not mutation_str.strip():
        return mutation_str
    parts = []
    for entry in mutation_str.split(", "):
        m = _ENTRY_RE.match(entry.strip())
        if not m:
            parts.append(entry)
            continue
        base, existing, rest = m.group(1), m.group(2), m.group(3)
        if "(PP:" in existing:
            parts.append(entry)
            continue
        pred, score = cache.get((gene, base), ("", ""))
        pp_tag = f"(PP:{pred},{score})" if pred else ""
        parts.append(f"{base}{existing}{pp_tag}{rest}")
    return ", ".join(parts)


def run_polyphen_phase(df: pd.DataFrame) -> None:
    """Phase 2: fetch PolyPhen-2 scores and tag mutation strings."""
    print("\n── Phase 2: PolyPhen-2 pathogenicity scores ──")

    cache = _pp_load_cache()

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

    to_fetch = list(needed)
    print(f"{len(to_fetch)} unique (gene, mutation) pairs to fetch "
          f"({len(cache)} already in cache from prior runs)")

    if to_fetch:
        with ThreadPoolExecutor(max_workers=_PP_MAX_WORKERS) as pool:
            futures = {pool.submit(fetch_polyphen, g, mut): (g, mut)
                       for g, mut in to_fetch}
            done = 0
            total = len(futures)
            for future in tqdm(as_completed(futures), total=total,
                               desc="Fetching PolyPhen-2 scores"):
                g, mut = futures[future]
                pred, score = future.result()
                cache[(g, mut)] = (pred, score)
                done += 1
                _emit_progress(1, done / total * 100, f"PolyPhen-2 scores: {done}/{total}")
        _pp_save_cache(cache)
    else:
        _emit_progress(1, 100, "PolyPhen-2 scores: all cached")

    tagged = 0
    for col in _MUTATION_COLS:
        if col not in df.columns:
            continue
        new_values = []
        for _, row in df.iterrows():
            annotated = annotate_mutation_string(row[col], row.get("gene", ""), cache)
            if "(PP:" in annotated:
                tagged += 1
            new_values.append(annotated)
        df[col] = new_values

    print(f"  Tagged {tagged} mutation entries with PolyPhen-2 predictions")


# ═════════════════════════════════════════════════════════════════════════════
# Phase 3: Kinase Library predictions for phosphorylation sites
# ═════════════════════════════════════════════════════════════════════════════

_KIN_TOP_K = 5
_KIN_WINDOW = 7
_KIN_CACHE_FILE = PROJECT_ROOT / "data" / "cache" / "kinase_predictions.tsv"


def _kin_load_cache() -> dict[str, str]:
    """Return {window_15mer: formatted_prediction} from the TSV cache."""
    if not _KIN_CACHE_FILE.exists():
        return {}
    df = pd.read_csv(_KIN_CACHE_FILE, sep="\t", dtype=str, keep_default_na=False)
    return {row["window"]: row["prediction"] for _, row in df.iterrows()}


def _kin_save_cache(cache: dict[str, str]) -> None:
    """Persist the kinase prediction cache to disk."""
    _KIN_CACHE_FILE.parent.mkdir(parents=True, exist_ok=True)
    rows = [{"window": w, "prediction": p} for w, p in cache.items()]
    pd.DataFrame(rows, columns=["window", "prediction"]).to_csv(
        _KIN_CACHE_FILE, sep="\t", index=False)


def extract_sequence(chain) -> dict[int, str]:
    """Build a {position: one_letter_aa} dict from a biotite chain's CA atoms."""
    ca_mask = chain.atom_name == "CA"
    ca_atoms = chain[ca_mask]
    return {
        int(ca_atoms.res_id[i]): AA3TO1.get(str(ca_atoms.res_name[i]), "X")
        for i in range(len(ca_atoms))
    }


def build_kinase_window(pos_to_aa: dict[int, str], site_pos: int) -> str | None:
    """Build a 15-mer sequence window centered on site_pos with lowercase phosphosite."""
    residue = pos_to_aa.get(site_pos)
    if not residue or residue not in ("S", "T", "Y"):
        return None
    chars = []
    for offset in range(-_KIN_WINDOW, _KIN_WINDOW + 1):
        p = site_pos + offset
        chars.append(residue.lower() if offset == 0 else pos_to_aa.get(p, "_"))
    return "".join(chars)


def predict_kinases(window: str) -> str:
    """Run the Kinase Library on a 15-mer window and return a formatted top-5 string."""
    import kinase_library as kl
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        sub = kl.Substrate(window)
        result = sub.predict()
    top = result.head(_KIN_TOP_K)
    parts = []
    for kinase, row in top.iterrows():
        parts.append(f"{kinase}({row['Score']:.2f},{row['Percentile']:.1f}%)")
    return "; ".join(parts)


def run_kinase_phase(df: pd.DataFrame) -> None:
    """Phase 3: predict upstream kinases for phosphorylation sites using CIF sequences."""
    print("\n── Phase 3: Kinase predictions ──")

    unique_uniprots = df["UniProt"].unique().tolist()
    seq_maps: dict[str, dict[int, str]] = {}
    skipped = 0
    for uid in unique_uniprots:
        uniprot_dir = MODELS_ROOT / uid
        cif_file = find_canonical_cif(uniprot_dir) if uniprot_dir.is_dir() else None
        if cif_file is None:
            skipped += 1
            continue
        chain = load_first_chain(cif_file)
        if chain is None:
            skipped += 1
            continue
        seq_maps[uid] = extract_sequence(chain)

    if skipped:
        print(f"  Skipped {skipped} proteins (no CIF or unparseable)")
    print(f"  Loaded sequences for {len(seq_maps)} proteins")

    cache = _kin_load_cache()
    cache_hits = 0
    new_predictions = 0

    predictions: list[str] = []
    annotated = 0
    total_rows = len(df)
    for row_idx, (_, row) in enumerate(tqdm(df.iterrows(), total=total_rows,
                                            desc="Predicting kinases")):
        uid = row.get("UniProt", "")
        ptm_site = row.get("ptm_site", "")
        ptm_type = row.get("ptm_type", "")

        _emit_progress(2, (row_idx + 1) / total_rows * 100,
                       f"Kinase predictions: {row_idx + 1}/{total_rows}")

        if "phosphorylation" not in ptm_type.lower():
            predictions.append("")
            continue

        m = SITE_RE.match(ptm_site.strip()) if ptm_site else None
        if not m:
            predictions.append("")
            continue

        pos = int(m.group(2))
        pos_to_aa = seq_maps.get(uid)
        if pos_to_aa is None:
            predictions.append("")
            continue

        window = build_kinase_window(pos_to_aa, pos)
        if window is None:
            predictions.append("")
            continue

        if window in cache:
            predictions.append(cache[window])
            cache_hits += 1
            annotated += 1
            continue

        try:
            pred = predict_kinases(window)
            cache[window] = pred
            predictions.append(pred)
            annotated += 1
            new_predictions += 1
        except Exception:
            predictions.append("")

    if new_predictions:
        _kin_save_cache(cache)

    df["kinase_predictions"] = predictions
    print(f"  Annotated {annotated}/{len(df)} rows "
          f"({cache_hits} cached, {new_predictions} newly predicted)")


# ═════════════════════════════════════════════════════════════════════════════
# Main
# ═════════════════════════════════════════════════════════════════════════════

def main() -> None:
    """Run all three annotation phases on the proximity database."""
    import argparse
    parser = argparse.ArgumentParser(description="Annotate the proximity database.")
    parser.add_argument(
        "--output-dir",
        default=str(DEFAULT_OUTPUT_DIR),
        help="Directory containing the proximity DB (default: Output/)",
    )
    args = parser.parse_args()

    proximity_db = Path(args.output_dir) / "ptm_mutation_proximity_db.tsv"
    print(f"Reading proximity DB: {proximity_db}")
    df = pd.read_csv(proximity_db, sep="\t", encoding="utf-16", dtype=str,
                     keep_default_na=False)
    print(f"{len(df)} rows, {df['UniProt'].nunique()} unique proteins\n")

    run_1433_phase(df)
    run_polyphen_phase(df)
    run_kinase_phase(df)

    df.to_csv(proximity_db, sep="\t", index=False, encoding="utf-16")
    print(f"\nUpdated proximity DB written to: {proximity_db}")


if __name__ == "__main__":
    main()
