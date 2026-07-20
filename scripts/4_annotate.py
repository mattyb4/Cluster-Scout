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
import time
import warnings
from concurrent.futures import ThreadPoolExecutor, as_completed
from functools import lru_cache
from pathlib import Path

import pandas as pd
import requests
from tqdm import tqdm

sys.path.insert(0, str(Path(__file__).resolve().parent))
from pipeline_utils import (  # noqa: E402
    project_root, AA3TO1, MUT_RE, SITE_RE, fmt_time,
    find_canonical_cif, load_first_chain,
    input_dir, resolve_input_file, INTERACTORS_1433_INPUT_DIR,
)

PROJECT_ROOT = project_root(__file__)
MODELS_ROOT = PROJECT_ROOT / "cif_models"
DEFAULT_OUTPUT_DIR = PROJECT_ROOT / "Output"

_NUM_PHASES = 5


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


def run_1433_phase(df: pd.DataFrame) -> tuple[dict, dict]:
    """Phase 1: fetch 14-3-3 predictions and confirmed sites, add columns to df.
    Returns (score_maps, confirmed_sites) for use in long-format annotation."""
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
    return score_maps, confirmed_sites


# ═════════════════════════════════════════════════════════════════════════════
# Phase 2: PolyPhen-2 mutation pathogenicity scores
# ═════════════════════════════════════════════════════════════════════════════

_PP_CACHE_FILE = PROJECT_ROOT / "data" / "cache" / "polyphen.tsv"
_PP_API_URL = "https://myvariant.info/v1/query"
_PP_MAX_WORKERS = 30
_pp_session = requests.Session()
_PP_SEVERITY = {"D": 2, "P": 1, "B": 0}
_PP_CLASS = {"D": "probably_damaging", "P": "possibly_damaging", "B": "benign"}
_PP_CODE_MAP = {"benign": "B", "possibly_damaging": "P", "probably_damaging": "D"}
_PP_TAG_RE = re.compile(r"\(PP:([DPB]),")
_MUT_POS_RE = re.compile(r"[A-Z](\d+)[A-Z*]")

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
    """Query myvariant.info for PolyPhen-2 HDIV prediction using a shared session."""
    m = MUT_RE.match(mutation)
    if not m:
        return "", ""
    ref, pos, alt = m.group(1), m.group(2), m.group(3)
    if alt == "*":
        return "", ""
    q = (f"dbnsfp.genename:{gene} AND dbnsfp.aa.ref:{ref} "
         f"AND dbnsfp.aa.alt:{alt} AND dbnsfp.aa.pos:{pos}")
    try:
        resp = _pp_session.get(
            _PP_API_URL,
            params={"q": q, "fields": "dbnsfp.polyphen2", "size": 10},
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


def run_polyphen_phase(df: pd.DataFrame) -> dict:
    """Phase 2: fetch PolyPhen-2 scores and tag mutation strings. Returns the full cache."""
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
    return cache


# ═════════════════════════════════════════════════════════════════════════════
# Phase 3: Kinase Library predictions for phosphorylation sites
# ═════════════════════════════════════════════════════════════════════════════

_KIN_TOP_K = 5
_KIN_WINDOW = 7
_KIN_CACHE_FILE = PROJECT_ROOT / "data" / "cache" / "kinase_predictions.tsv"
_KIN_MAX_WORKERS = 6


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


def _speed_up_kinase_library() -> None:
    """Memoize kinase_library's reference-data loaders for this process.

    Profiling showed Substrate.predict() re-reads the same kinome/matrix
    reference files from disk on every single call (~50 pd.read_csv calls per
    prediction) rather than caching them, even though those files never change
    during a run. get_kinase_list/get_kinome_info are the two most-called
    loaders and, in every call path reachable from predict(), are only ever
    invoked with hashable scalar arguments (kin_type, non_canonical) — so
    wrapping them in an unbounded in-memory cache is safe and eliminates the
    bulk of that redundant I/O. Idempotent: safe to call more than once.
    """
    import kinase_library.modules.data as kl_data
    if getattr(kl_data, "_cluster_scout_cached", False):
        return
    kl_data.get_kinase_list = lru_cache(maxsize=None)(kl_data.get_kinase_list)
    kl_data.get_kinome_info = lru_cache(maxsize=None)(kl_data.get_kinome_info)
    kl_data._cluster_scout_cached = True


def predict_kinases(window: str) -> str:
    """Run the Kinase Library on a 15-mer window and return a formatted top-5 string."""
    import kinase_library as kl
    _speed_up_kinase_library()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        sub = kl.Substrate(window)
        result = sub.predict()
    top = result.head(_KIN_TOP_K)
    parts = []
    for kinase, row in top.iterrows():
        parts.append(f"{kinase}({row['Score']:.2f},{row['Percentile']:.1f}%)")
    return "; ".join(parts)


def run_kinase_phase(df: pd.DataFrame) -> tuple[dict, dict]:
    """Phase 3: predict upstream kinases for phosphorylation sites using CIF sequences.
    Returns (seq_maps, kin_cache) for use in long-format annotation."""
    print("\n── Phase 3: Kinase predictions ──")

    unique_uniprots = df["UniProt"].unique().tolist()
    print(f"  Loading sequences for {len(unique_uniprots)} proteins...")
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

    # First pass: resolve each row's kinase window (or None if not applicable),
    # without predicting anything yet, so identical windows shared across rows
    # (a common case — many PTM sites recur across proteins) are only ever
    # predicted once instead of once per row.
    row_windows: list[str | None] = []
    for _, row in df.iterrows():
        uid = row.get("UniProt", "")
        ptm_site = row.get("ptm_site", "")
        ptm_type = row.get("ptm_type", "")

        window = None
        if "phosphorylation" in ptm_type.lower():
            m = SITE_RE.match(ptm_site.strip()) if ptm_site else None
            if m:
                pos_to_aa = seq_maps.get(uid)
                if pos_to_aa is not None:
                    window = build_kinase_window(pos_to_aa, int(m.group(2)))
        row_windows.append(window)

    all_windows = {w for w in row_windows if w is not None}
    to_predict = sorted(all_windows - cache.keys())
    print(f"  {len(all_windows) - len(to_predict)}/{len(all_windows)} windows cached; "
          f"predicting {len(to_predict)} new...")

    failed_windows: set[str] = set()
    if to_predict:
        with ThreadPoolExecutor(max_workers=_KIN_MAX_WORKERS) as pool:
            futures = {pool.submit(predict_kinases, w): w for w in to_predict}
            done = 0
            total = len(futures)
            for future in tqdm(as_completed(futures), total=total, desc="Predicting kinases"):
                window = futures[future]
                try:
                    cache[window] = future.result()
                except Exception:
                    failed_windows.add(window)
                done += 1
                _emit_progress(2, done / total * 100, f"Kinase predictions: {done}/{total}")
        new_count = len(to_predict) - len(failed_windows)
        if new_count:
            _kin_save_cache(cache)
    else:
        new_count = 0
        _emit_progress(2, 100, "Kinase predictions: all cached")

    # Second pass: assemble the final per-row predictions from the (now fully
    # populated) cache, in the original row order.
    predictions = [
        "" if w is None or w in failed_windows else cache.get(w, "")
        for w in row_windows
    ]
    annotated = sum(1 for p in predictions if p)

    df["kinase_predictions"] = predictions
    print(f"  Annotated {annotated}/{len(df)} rows "
          f"({len(all_windows) - len(to_predict)} windows cached, "
          f"{new_count} windows newly predicted)")
    return seq_maps, cache


# ═════════════════════════════════════════════════════════════════════════════
# PolyPhen class filter
# ═════════════════════════════════════════════════════════════════════════════


def _filter_mut_str(mutation_str: str, exclude_codes: set[str]) -> str:
    """Remove mutation entries whose PP code is in exclude_codes; keep unscored entries."""
    if not mutation_str or not exclude_codes:
        return mutation_str
    kept = []
    for entry in mutation_str.split(", "):
        entry = entry.strip()
        if not entry:
            continue
        m = _PP_TAG_RE.search(entry)
        code = m.group(1) if m else ""
        if code and code in exclude_codes:
            continue
        kept.append(entry)
    return ", ".join(kept)


def _positions_from_str(mutation_str: str) -> list[int]:
    """Return unique residue positions from a formatted mutation string, in order."""
    seen: set[int] = set()
    result = []
    for entry in (mutation_str or "").split(", "):
        m = _MUT_POS_RE.search(entry)
        if m:
            pos = int(m.group(1))
            if pos not in seen:
                seen.add(pos)
                result.append(pos)
    return result


def apply_polyphen_filter(df: pd.DataFrame, exclude_classes: list[str]) -> pd.DataFrame:
    """Remove mutations of excluded PP classes from the wide-format proximity DB.

    Filters mutation strings, recomputes count and distance columns, and drops PTM
    rows where no qualifying mutations remain.  *_total_patient_count columns are
    left unchanged (they reflect pre-filter totals and cannot be recomputed here).
    """
    exclude_codes = {_PP_CODE_MAP[c] for c in exclude_classes if c in _PP_CODE_MAP}
    if not exclude_codes:
        return df

    print(f"\nApplying PolyPhen filter — excluding: {', '.join(exclude_classes)}")
    print("  Note: *_total_patient_count columns retain pre-filter totals")

    within_col = "mutations_within_5_positions"
    beyond_col = "mutations_more_than_5_positions"
    disrupting_col = "confirmed_disrupting_mutations"

    for col in (within_col, beyond_col, disrupting_col):
        if col in df.columns:
            df[col] = df[col].fillna("").apply(
                lambda s: _filter_mut_str(s, exclude_codes)
            )

    def _count(s: str) -> int:
        return len([e for e in s.split(", ") if e.strip()]) if s else 0

    def _uniq_pos(s: str) -> int:
        return len(set(_positions_from_str(s)))

    for prefix, col in [("within_5", within_col), ("more_than_5", beyond_col)]:
        df[f"mutation_count_{prefix}_positions"] = df[col].apply(_count)
        df[f"unique_mutation_position_count_{prefix}_positions"] = df[col].apply(_uniq_pos)

    # Recompute morethan5_linear_distance from filtered beyond-5 string
    if "morethan5_linear_distance" in df.columns and "ptm_site" in df.columns:
        def _linear_dists(row) -> str:
            ptm_m = re.search(r"(\d+)", str(row.get("ptm_site", "")))
            if not ptm_m:
                return ""
            ptm_pos = int(ptm_m.group(1))
            return ",".join(str(abs(p - ptm_pos)) for p in _positions_from_str(row[beyond_col]))
        df["morethan5_linear_distance"] = df.apply(_linear_dists, axis=1)

    # Recompute mutation_at_ptm_site from filtered within-5 string
    if "mutation_at_ptm_site" in df.columns and "ptm_site" in df.columns:
        def _at_ptm(row) -> str:
            ptm_m = re.search(r"(\d+)", str(row.get("ptm_site", "")))
            if not ptm_m:
                return "no"
            return "yes" if int(ptm_m.group(1)) in _positions_from_str(row[within_col]) else "no"
        df["mutation_at_ptm_site"] = df.apply(_at_ptm, axis=1)

    before = len(df)
    df = df[
        (df[within_col].fillna("").str.len() > 0) |
        (df[beyond_col].fillna("").str.len() > 0)
    ].reset_index(drop=True)
    removed = before - len(df)
    print(f"  Removed {removed} PTM rows with no qualifying mutations; "
          f"{len(df)} rows remaining")
    return df


# ═════════════════════════════════════════════════════════════════════════════
# Phase 4: AIUPred disorder scores
# ═════════════════════════════════════════════════════════════════════════════

_AIUPRED_API_URL = "https://aiupred.elte.hu/rest_api"
_AIUPRED_CACHE_FILE = PROJECT_ROOT / "data" / "cache" / "aiupred_disorder.tsv"
_aiupred_session = requests.Session()
# Modest concurrency: aiupred.elte.hu is a single-instance academic server, not a
# scalable production API (unlike myvariant.info in Phase 2) — stay polite to it.
_AIUPRED_MAX_WORKERS = 5

# Maps response JSON keys to canonical type names stored in the cache
_AIUPRED_KEY_MAP = {
    "AIUPred": "general",
    "AIUPred-binding": "binding",
    "AIUPred-linker": "linker",
}

# One binding call yields both general disorder and binding-region scores
_AIUPRED_CALLS = [
    ("binding", ["general", "binding"]),
]


def _aiupred_load_cache() -> dict[tuple[str, str], dict[int, float]]:
    """Load cached AIUPred scores. Returns {(uniprot_id, analysis_type): {pos: score}}."""
    if not _AIUPRED_CACHE_FILE.exists():
        return {}
    try:
        df = pd.read_csv(_AIUPRED_CACHE_FILE, sep="\t", dtype=str, keep_default_na=False)
    except Exception:
        return {}
    cache: dict[tuple[str, str], dict[int, float]] = {}
    for _, row in df.iterrows():
        uid = str(row.get("uniprot_id", ""))
        atype = str(row.get("analysis_type", ""))
        try:
            scores = {int(k): float(v) for k, v in json.loads(row.get("scores_json", "{}")).items()}
        except Exception:
            scores = {}
        if uid and atype:
            cache[(uid, atype)] = scores
    return cache


def _aiupred_save_cache(cache: dict[tuple[str, str], dict[int, float]]) -> None:
    _AIUPRED_CACHE_FILE.parent.mkdir(parents=True, exist_ok=True)
    rows = [
        {"uniprot_id": uid, "analysis_type": atype,
         "scores_json": json.dumps({str(k): v for k, v in scores.items()})}
        for (uid, atype), scores in cache.items()
    ]
    pd.DataFrame(rows, columns=["uniprot_id", "analysis_type", "scores_json"]).to_csv(
        _AIUPRED_CACHE_FILE, sep="\t", index=False
    )


def fetch_aiupred_all(api_type: str, uniprot_id: str) -> dict[str, dict[int, float]]:
    """Fetch AIUPred scores for one protein, returning all score arrays in the response.

    Returns {type_name: {1-based position: score}}.  A binding call returns both
    'general' and 'binding' scores; a linker call returns only 'linker'.
    """
    params: dict[str, str] = {"accession": uniprot_id, "smoothing": "default",
                               "analysis_type": api_type}
    try:
        resp = _aiupred_session.get(_AIUPRED_API_URL, params=params, timeout=60)
        resp.raise_for_status()
        data = resp.json()
    except Exception:
        return {}
    if not isinstance(data, dict):
        return {}
    result: dict[str, dict[int, float]] = {}
    for key, val in data.items():
        type_name = _AIUPRED_KEY_MAP.get(key)
        if type_name is None:
            continue
        if isinstance(val, list) and val and isinstance(val[0], (int, float)):
            result[type_name] = {i + 1: float(s) for i, s in enumerate(val)}
    return result


def run_aiupred_phase(df: pd.DataFrame) -> dict[str, dict[str, dict[int, float]]]:
    """Phase 4: fetch general and binding disorder scores for each protein.

    Makes one API call per protein (binding), which yields both general and
    binding scores. Returns {type_name: {uniprot_id: {position: score}}} for
    types general/binding.
    """
    print("\n── Phase 4: AIUPred disorder scores (general + binding) ──")
    cache = _aiupred_load_cache()
    uniprots = [u for u in df["UniProt"].unique() if u]

    needs = [
        (api_type, [u for u in uniprots if any((u, t) not in cache for t in produced)])
        for api_type, produced in _AIUPRED_CALLS
    ]
    total_fetches = sum(len(n) for _, n in needs)
    done = 0
    any_new = False

    for api_type, need in needs:
        if not need:
            print(f"  {api_type}: all {len(uniprots)} proteins cached")
            continue
        print(f"  Fetching {len(need)} proteins via {api_type} call")
        with ThreadPoolExecutor(max_workers=_AIUPRED_MAX_WORKERS) as pool:
            futures = {pool.submit(fetch_aiupred_all, api_type, uid): uid for uid in need}
            for future in tqdm(as_completed(futures), total=len(futures), desc=f"  AIUPred/{api_type}"):
                uid = futures[future]
                for t, scores in future.result().items():
                    cache[(uid, t)] = scores
                done += 1
                _emit_progress(3, done / max(total_fetches, 1) * 100,
                               f"AIUPred {uid}: {done}/{total_fetches}")
        any_new = True

    if any_new:
        _aiupred_save_cache(cache)

    return {
        t: {uid: cache.get((uid, t), {}) for uid in uniprots}
        for t in ("general", "binding")
    }


# ═════════════════════════════════════════════════════════════════════════════
# Phase 5: InterPro functional domains
# ═════════════════════════════════════════════════════════════════════════════

_INTERPRO_CACHE_FILE = PROJECT_ROOT / "data" / "cache" / "interpro_domains.tsv"
_INTERPRO_API_URL = "https://www.ebi.ac.uk/interpro/api/entry/interpro/protein/uniprot/{uid}/"
_INTERPRO_MAX_WORKERS = 5


def _interpro_load_cache() -> dict[str, list[dict]]:
    """Load cached InterPro entries. Returns {uniprot_id: [{"name","type","start","end"}, ...]}."""
    if not _INTERPRO_CACHE_FILE.exists():
        return {}
    try:
        df = pd.read_csv(_INTERPRO_CACHE_FILE, sep="\t", dtype=str, keep_default_na=False)
    except Exception:
        return {}
    cache: dict[str, list[dict]] = {}
    for _, row in df.iterrows():
        uid = str(row.get("uniprot_id", ""))
        if not uid:
            continue
        try:
            cache[uid] = json.loads(row.get("entries_json", "[]"))
        except Exception:
            cache[uid] = []
    return cache


def _interpro_save_cache(cache: dict[str, list[dict]]) -> None:
    _INTERPRO_CACHE_FILE.parent.mkdir(parents=True, exist_ok=True)
    rows = [
        {"uniprot_id": uid, "entries_json": json.dumps(entries)}
        for uid, entries in cache.items()
    ]
    pd.DataFrame(rows, columns=["uniprot_id", "entries_json"]).to_csv(
        _INTERPRO_CACHE_FILE, sep="\t", index=False
    )


def fetch_interpro_domains(uniprot_id: str) -> list[dict]:
    """Fetch curated InterPro entries (domains/families/sites/etc.) for one protein.

    Queries /entry/interpro/ (curated InterPro-integrated entries only, not
    every individual member-database signature under /entry/all/ — the
    latter returns many overlapping near-duplicates, e.g. Pfam, CDD, and
    PROSITE each separately flagging essentially the same domain).

    Returns a list of {"name", "type", "start", "end"} dicts, one per
    (entry, fragment) pair — an entry can have multiple discontinuous
    fragments, each treated as its own range for containment checks.
    """
    entries: list[dict] = []
    url = _INTERPRO_API_URL.format(uid=uniprot_id)
    try:
        while url:
            resp = requests.get(url, timeout=30)
            resp.raise_for_status()
            data = resp.json()
            for result in data.get("results", []):
                md = result.get("metadata", {})
                name = md.get("name", "")
                etype = md.get("type", "")
                for protein in result.get("proteins", []):
                    for loc in protein.get("entry_protein_locations", []):
                        for frag in loc.get("fragments", []):
                            start, end = frag.get("start"), frag.get("end")
                            if start is not None and end is not None:
                                entries.append({"name": name, "type": etype, "start": start, "end": end})
            url = data.get("next")
    except Exception:
        return entries
    return entries


def run_interpro_phase(df: pd.DataFrame) -> dict[str, list[dict]]:
    """Phase 5: fetch InterPro functional-domain entries for each protein."""
    print("\n── Phase 5: InterPro functional domains ──")
    cache = _interpro_load_cache()
    uniprots = [u for u in df["UniProt"].unique() if u]
    need = [u for u in uniprots if u not in cache]

    if not need:
        print(f"  all {len(uniprots)} proteins cached")
    else:
        print(f"  Fetching {len(need)} proteins")
        with ThreadPoolExecutor(max_workers=_INTERPRO_MAX_WORKERS) as pool:
            futures = {pool.submit(fetch_interpro_domains, uid): uid for uid in need}
            for i, future in enumerate(tqdm(as_completed(futures), total=len(futures), desc="  InterPro"), 1):
                uid = futures[future]
                cache[uid] = future.result()
                _emit_progress(4, i / max(len(need), 1) * 100, f"InterPro {uid}: {i}/{len(need)}")
        _interpro_save_cache(cache)

    return {uid: cache.get(uid, []) for uid in uniprots}


def find_domain_at_position(entries: list[dict], position: int | None) -> str:
    """Return a "name (type, start-end)" string for every InterPro entry
    containing *position*, semicolon-joined (a residue can legitimately fall
    inside more than one entry, e.g. a specific domain nested inside a
    broader homologous-superfamily call), or "" if none / position is None.
    """
    if position is None:
        return ""
    hits = [e for e in entries if e["start"] <= position <= e["end"]]
    return "; ".join(f"{e['name']} ({e['type']}, {e['start']}-{e['end']})" for e in hits)


# ═════════════════════════════════════════════════════════════════════════════
# Long-format annotation
# ═════════════════════════════════════════════════════════════════════════════


def annotate_long_format(
    df: pd.DataFrame,
    score_maps: dict,
    confirmed_sites: dict,
    pp_cache: dict,
    seq_maps: dict,
    kin_cache: dict,
    disorder_maps: dict | None = None,
    domain_maps: dict | None = None,
) -> None:
    """Fill annotation columns in the long-format PTM/mutation table in-place."""
    pred_1433, consensus_1433, confirmed_1433 = [], [], []
    pp_scores, pp_classes = [], []
    kinase_preds = []
    _ATYPES = ("general", "binding")
    ptm_disorder: dict[str, list] = {t: [] for t in _ATYPES}
    mut_disorder: dict[str, list] = {t: [] for t in _ATYPES}
    ptm_domains, mut_domains = [], []

    for _, row in df.iterrows():
        uid = str(row.get("uniprot_id", "") or "")
        ptm_site = str(row.get("ptm_position", "") or "")
        ptm_type = str(row.get("ptm_type", "") or "")
        gene = str(row.get("gene", "") or "")
        mutation = str(row.get("mutation", "") or "")
        clean_mut = mutation.replace("(isoform?)", "")

        # 14-3-3 — only applicable to S/T residues
        m_site = SITE_RE.match(ptm_site.strip()) if ptm_site else None
        is_st = bool(m_site and m_site.group(1) in ("S", "T"))
        if is_st:
            binding, score = annotate_1433_row(ptm_site, score_maps.get(uid, {}))
            pred_1433.append(binding.lower() if binding else "no")
            consensus_1433.append(score)
            conf, _ = annotate_confirmed(uid, ptm_site, confirmed_sites)
            confirmed_1433.append("yes" if conf == "Yes" else "no")
        else:
            pred_1433.append("")
            consensus_1433.append("")
            confirmed_1433.append("")

        # PolyPhen-2
        pred, pp_score = pp_cache.get((gene, clean_mut), ("", ""))
        pp_scores.append(pp_score)
        pp_classes.append(_PP_CLASS.get(pred, ""))

        # Kinase Library — phosphorylation sites only
        if "phosphorylation" in ptm_type.lower() and m_site:
            pos = int(m_site.group(2))
            pos_to_aa = seq_maps.get(uid)
            window = build_kinase_window(pos_to_aa, pos) if pos_to_aa else None
            kinase_preds.append(kin_cache.get(window, "") if window else "")
        else:
            kinase_preds.append("")

        # AIUPred disorder — compute PTM/mut positions once, look up all three types
        ptm_pos_m = SITE_RE.match(ptm_site.strip()) if ptm_site else None
        ptm_pos = int(ptm_pos_m.group(2)) if ptm_pos_m else None
        mut_pos_m = MUT_RE.search(clean_mut)
        mut_pos = int(mut_pos_m.group(2)) if mut_pos_m else None
        for t in _ATYPES:
            if disorder_maps is not None:
                ps = disorder_maps[t].get(uid, {})
                ptm_disorder[t].append(f"{ps[ptm_pos]:.3f}" if ptm_pos in ps else "")
                mut_disorder[t].append(f"{ps[mut_pos]:.3f}" if mut_pos in ps else "")
            else:
                ptm_disorder[t].append("")
                mut_disorder[t].append("")

        # InterPro functional domains — same ptm_pos/mut_pos computed above
        if domain_maps is not None:
            entries = domain_maps.get(uid, [])
            ptm_domains.append(find_domain_at_position(entries, ptm_pos))
            mut_domains.append(find_domain_at_position(entries, mut_pos))
        else:
            ptm_domains.append("")
            mut_domains.append("")

    df["polyphen_score"] = pp_scores
    df["polyphen_class"] = pp_classes
    df["1433_predicted"] = pred_1433
    df["1433_predicted_consensus"] = consensus_1433
    df["1433_confirmed"] = confirmed_1433
    df["kinase_predictions"] = kinase_preds
    for t in _ATYPES:
        df[f"ptm_aiupred_{t}"] = ptm_disorder[t]
        df[f"mut_aiupred_{t}"] = mut_disorder[t]
    df["ptm_is_disordered"] = [
        "yes" if v and float(v) > 0.5 else "no" for v in ptm_disorder["general"]
    ]
    df["ptm_is_binding"] = [
        "yes" if v and float(v) > 0.5 else "no" for v in ptm_disorder["binding"]
    ]
    df["mut_is_disordered"] = [
        "yes" if v and float(v) > 0.5 else "no" for v in mut_disorder["general"]
    ]
    df["mut_is_binding"] = [
        "yes" if v and float(v) > 0.5 else "no" for v in mut_disorder["binding"]
    ]
    df["ptm_domain"] = ptm_domains
    df["mutation_domain"] = mut_domains


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
    parser.add_argument(
        "--pp-exclude",
        nargs="+",
        choices=["benign", "possibly_damaging", "probably_damaging"],
        default=[],
        metavar="CLASS",
        help="Exclude mutations with these PolyPhen-2 classes from the output",
    )
    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    proximity_db = output_dir / "ptm_mutation_proximity_db.tsv"
    print(f"Reading proximity DB: {proximity_db}")
    df = pd.read_csv(proximity_db, sep="\t", encoding="utf-16", dtype=str,
                     keep_default_na=False)
    print(f"{len(df)} rows, {df['UniProt'].nunique()} unique proteins\n")

    t0 = time.time()
    score_maps, confirmed_sites = run_1433_phase(df)
    print(f"  Phase 1 (14-3-3) completed in {fmt_time(time.time() - t0)}")

    t0 = time.time()
    pp_cache = run_polyphen_phase(df)
    print(f"  Phase 2 (PolyPhen-2) completed in {fmt_time(time.time() - t0)}")

    t0 = time.time()
    seq_maps, kin_cache = run_kinase_phase(df)
    print(f"  Phase 3 (Kinase) completed in {fmt_time(time.time() - t0)}")

    t0 = time.time()
    disorder_maps = run_aiupred_phase(df)
    print(f"  Phase 4 (AIUPred) completed in {fmt_time(time.time() - t0)}")
    for atype in ("general", "binding"):
        col_vals = []
        for _, row in df.iterrows():
            uid = str(row.get("UniProt", "") or "")
            ptm_site = str(row.get("ptm_site", "") or "")
            m = SITE_RE.match(ptm_site.strip()) if ptm_site else None
            pos = int(m.group(2)) if m else None
            ps = disorder_maps[atype].get(uid, {})
            col_vals.append(f"{ps[pos]:.3f}" if pos in ps else "")
        df[f"ptm_aiupred_{atype}"] = col_vals
    df["ptm_is_disordered"] = df["ptm_aiupred_general"].apply(
        lambda v: "yes" if v and float(v) > 0.5 else "no"
    )
    df["ptm_is_binding"] = df["ptm_aiupred_binding"].apply(
        lambda v: "yes" if v and float(v) > 0.5 else "no"
    )

    t0 = time.time()
    domain_maps = run_interpro_phase(df)
    print(f"  Phase 5 (InterPro) completed in {fmt_time(time.time() - t0)}")
    col_vals = []
    for _, row in df.iterrows():
        uid = str(row.get("UniProt", "") or "")
        ptm_site = str(row.get("ptm_site", "") or "")
        m = SITE_RE.match(ptm_site.strip()) if ptm_site else None
        pos = int(m.group(2)) if m else None
        col_vals.append(find_domain_at_position(domain_maps.get(uid, []), pos))
    df["ptm_domain"] = col_vals

    if args.pp_exclude:
        df = apply_polyphen_filter(df, args.pp_exclude)

    df.to_csv(proximity_db, sep="\t", index=False, encoding="utf-16")
    print(f"\nUpdated proximity DB written to: {proximity_db}")

    long_db = output_dir / "ptm_mutation_proximity_long.tsv"
    if long_db.exists():
        print(f"\nAnnotating long-format DB: {long_db}")
        df_long = pd.read_csv(long_db, sep="\t", encoding="utf-16", dtype=str,
                              keep_default_na=False)
        print(f"{len(df_long)} rows")
        annotate_long_format(
            df_long, score_maps, confirmed_sites, pp_cache, seq_maps, kin_cache,
            disorder_maps, domain_maps,
        )
        if args.pp_exclude:
            before = len(df_long)
            df_long = df_long[
                ~df_long["polyphen_class"].isin(args.pp_exclude)
            ].reset_index(drop=True)
            print(f"  Long format: removed {before - len(df_long)} rows by PolyPhen filter")
        df_long.to_csv(long_db, sep="\t", index=False, encoding="utf-16")
        print(f"Updated long-format DB written to: {long_db}")
    else:
        print(f"\nNo long-format DB found at {long_db} — skipping")


if __name__ == "__main__":
    main()
