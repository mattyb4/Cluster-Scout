"""Annotate the pipeline's output with 14-3-3 binding predictions, PolyPhen-2
pathogenicity scores, predicted upstream kinases, AIUPred disorder, and
InterPro functional domains.

  Phase 1 — 14-3-3: Queries the 14-3-3-Pred API for predicted binding-site
            scores and cross-references experimentally confirmed interactors.
  Phase 2 — PolyPhen-2: Queries myvariant.info for pathogenicity predictions
            and tags each mutation with a (PP:D/P/B,score) label.
  Phase 3 — Kinases: Uses the Kinase Library to predict the top 5 upstream
            kinases for each phosphorylation site based on the ±7 residue
            sequence window from AlphaFold CIF structures.
  Phase 4 — AIUPred: Predicted intrinsic disorder / binding-region scores.
  Phase 5 — InterPro: Curated functional-domain lookup per position.

--mode ptm-proximity (default) runs all 5 phases against
ptm_mutation_proximity_db.tsv/_long.tsv. --mode mutation-clustering runs only
Phases 2/4/5 (PolyPhen/AIUPred/InterPro are mutation/position-level; 14-3-3
and Kinase require a curated PTM site, which mutation-clustering mode has no
concept of) against mutation_cluster_db.tsv/_long.tsv. --ptm-source psp
annotates the PhosphoSitePlus-based psp_ptm_mutation_proximity_db.tsv/_long.tsv
instead.
"""
from __future__ import annotations

import json
import os
import re
import shutil
import subprocess
import sys
import threading
import time
from collections.abc import Callable
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Any

import pandas as pd
import requests
from tqdm import tqdm

sys.path.insert(0, str(Path(__file__).resolve().parent))
from pipeline_utils import (  # noqa: E402
    AA3TO1,
    DEFAULT_PTM_SOURCE,
    INTERACTORS_1433_INPUT_DIR,
    MUT_RE,
    PTM_SOURCES,
    SITE_RE,
    find_canonical_cif,
    fmt_time,
    hotspots_tsv_path,
    input_dir,
    load_first_chain,
    project_root,
    ptm_output_paths,
    resolve_input_file,
)

PROJECT_ROOT = project_root(__file__)
MODELS_ROOT = PROJECT_ROOT / "cif_models"
DEFAULT_OUTPUT_DIR = PROJECT_ROOT / "Output"

_NUM_PHASES = 5


# Each phase's own progress (0-100), by phase index. Phases run concurrently
# (see _run_phases), so overall progress is their average rather than
# "phases before this one are done".
_phase_progress: dict[int, float] = {}
_progress_lock = threading.Lock()


def _emit_progress(phase: int, phase_pct: float, desc: str, num_phases: int = _NUM_PHASES) -> None:
    """Print overall progress for the app to parse. phase is 0-indexed, phase_pct is 0-100."""
    with _progress_lock:
        _phase_progress[phase] = phase_pct
        overall = int(sum(_phase_progress.get(i, 0.0) for i in range(num_phases)) / num_phases)
        print(f"\r##PROGRESS## {overall} {desc}", end="", flush=True)


def _run_phases(phases: dict[str, tuple[int, int, Callable[[], Any]]]) -> dict[str, Any]:
    """Run annotation phases at the same time and return {name: result}.

    *phases* maps a name to (phase index, number of phases, callable). The
    phases wait on different web services (or, for kinases, local helper
    processes), so running them together makes step 4 about as long as its
    slowest phase instead of the sum -- without sending any one service more
    requests than before. Each callable must work on its own copy of the data
    it needs; the caller merges results back in a fixed order, so output is
    the same as running them one after another.
    """
    _phase_progress.clear()
    results: dict[str, Any] = {}
    started = time.time()
    with ThreadPoolExecutor(max_workers=len(phases)) as pool:
        futures = {pool.submit(fn): (name, idx, total) for name, (idx, total, fn) in phases.items()}
        for future in as_completed(futures):
            name, idx, total = futures[future]
            results[name] = future.result()
            _emit_progress(idx, 100, f"{name} done", total)
            print(f"\n  {name} finished after {fmt_time(time.time() - started)}")
    return results

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
    # Write-then-rename, so a fetch stopped mid-write (the prefetch that runs
    # during steps 2-3 is stopped when step 4 starts) never leaves a partial file
    tmp_file = cache_file.with_name(f"{cache_file.name}.{threading.get_ident()}.tmp")
    with tmp_file.open("w", encoding="utf-8") as f:
        json.dump(data, f)
    os.replace(tmp_file, cache_file)
    return data


def proteins_with_st_sites(uniprots, sites) -> list[str]:
    """The proteins (in first-seen order) that have at least one Ser/Thr PTM
    site -- the only residues 14-3-3-Pred scores, so the only proteins whose
    14-3-3 request can change the output."""
    wanted: dict[str, None] = {}
    for uid, site in zip(uniprots, sites):
        m = SITE_RE.match(str(site or "").strip())
        if uid and m and m.group(1) in ("S", "T"):
            wanted.setdefault(uid, None)
    return list(wanted)


def fetch_1433_for(uniprots: list[str], progress=None) -> dict[str, dict[int, float]]:
    """{uniprot: {position: consensus score}} for *uniprots*, fetched with
    _1433_MAX_WORKERS concurrent requests (cached ones read from disk)."""
    score_maps: dict[str, dict[int, float]] = {}
    with ThreadPoolExecutor(max_workers=_1433_MAX_WORKERS) as pool:
        futures = {pool.submit(fetch_1433pred, uid): uid for uid in uniprots}
        done = 0
        for future in tqdm(as_completed(futures), total=len(futures), desc="Fetching 14-3-3-Pred data"):
            uid = futures[future]
            data = future.result()
            if data is not None:
                score_maps[uid] = build_site_score_map(data)
            done += 1
            if progress:
                progress(done, len(futures))
    return score_maps


def prefetch_1433(step1_file: Path) -> None:
    """Fetch 14-3-3-Pred results into the cache for a step-1 table's
    proteins with Ser/Thr PTM sites. Run in the background during steps 2-3
    (the protein list is known after step 1), so step 4 finds most of them
    cached; stopping it partway is safe."""
    table = pd.read_csv(step1_file, sep="\t", dtype=str, keep_default_na=False)
    uniprots, sites = [], []
    for uid, field in zip(table["uniprot_id"], table["ptms_on_protein"]):
        for token in field.split(";"):
            uniprots.append(uid)
            sites.append(token.strip().split(":", 1)[0])
    wanted = proteins_with_st_sites(uniprots, sites)
    missing = [u for u in wanted if not (_1433_CACHE_DIR / f"{u}.json").exists()]
    print(f"14-3-3 prefetch: {len(wanted) - len(missing)}/{len(wanted)} proteins with Ser/Thr sites cached; "
          f"fetching {len(missing)}", flush=True)
    fetch_1433_for(missing)
    print("14-3-3 prefetch finished", flush=True)


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

    # Only proteins with a Ser/Thr site get 14-3-3 values, so only they need a request
    unique_uniprots = proteins_with_st_sites(df["UniProt"], df["ptm_site"])
    already_cached = sum(
        1 for uid in unique_uniprots if (_1433_CACHE_DIR / f"{uid}.json").exists()
    )
    print(f"{already_cached}/{len(unique_uniprots)} proteins with Ser/Thr sites already cached; "
          f"fetching {len(unique_uniprots) - already_cached} new...")

    score_maps = fetch_1433_for(
        unique_uniprots,
        progress=lambda done, total: _emit_progress(0, done / total * 100, f"14-3-3 predictions: {done}/{total}"),
    )

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
# myvariant.info rate-limits by burst/concurrency, not steady request rate -- 30
# concurrent workers triggered ~90% HTTP 429s within a second in testing, while
# sequential requests paced this closely succeeded ~100% of the time. So this
# phase fetches sequentially with a small delay instead of via a thread pool --
# but asks for many variants per request (see fetch_polyphen_batch), which cut
# the time per variant from ~250 ms to ~17 ms in testing.
_PP_REQUEST_DELAY = 0.1
_PP_BATCH_POSITIONS = 100      # positions per batched request
_PP_BATCH_MAX_HITS = 1000      # myvariant.info's largest page; bigger chunks are split
_PP_BATCH_FIELDS = "dbnsfp.polyphen2,dbnsfp.aa.ref,dbnsfp.aa.alt,dbnsfp.aa.pos"
_PP_RATE_LIMIT_RETRIES = 3
_PP_RATE_LIMIT_BACKOFF = 1.0
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


def _as_list(value) -> list:
    return value if isinstance(value, list) else [value]


def _pp_best_prediction(hits: list[dict]) -> tuple[str, str]:
    """Return the most severe HDIV (pred, score) across all hits and
    transcripts -- and within that class, the highest score, so the answer
    doesn't depend on the order myvariant.info returns records in."""
    best_pred, best_score = "", -1.0
    for hit in hits:
        for dbnsfp in _as_list(hit.get("dbnsfp") or {}):
            polyphen2 = dbnsfp.get("polyphen2", {}) if isinstance(dbnsfp, dict) else {}
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
                    value = float(score) if isinstance(score, (int, float)) else -1.0
                    if (_PP_SEVERITY[pred], value) > (_PP_SEVERITY.get(best_pred, -1), best_score):
                        best_pred, best_score = pred, value
    if not best_pred:
        return "", ""
    return best_pred, (f"{best_score:.3f}" if best_score >= 0 else "")


def fetch_polyphen(gene: str, mutation: str) -> tuple[str, str] | None:
    """Query myvariant.info for PolyPhen-2 HDIV prediction using a shared session.

    Returns ("", "") for a successful query that genuinely found no HDIV score --
    a real, cacheable answer. Returns None if the request itself failed (network
    error, timeout, non-200 after retries): the caller must NOT cache this, so a
    transient failure gets retried on the next run instead of being permanently
    mistaken for "no data".

    Retries on HTTP 429 with backoff as a safety net -- the caller is expected to
    already be pacing requests sequentially (see _PP_REQUEST_DELAY), so 429s here
    should be rare, not the normal case.
    """
    m = MUT_RE.match(mutation)
    if not m:
        return "", ""
    ref, pos, alt = m.group(1), m.group(2), m.group(3)
    if alt == "*":
        return "", ""
    q = (f"dbnsfp.genename:{gene} AND dbnsfp.aa.ref:{ref} "
         f"AND dbnsfp.aa.alt:{alt} AND dbnsfp.aa.pos:{pos}")
    backoff = _PP_RATE_LIMIT_BACKOFF
    for attempt in range(_PP_RATE_LIMIT_RETRIES):
        try:
            resp = _pp_session.get(
                _PP_API_URL,
                params={"q": q, "fields": "dbnsfp.polyphen2", "size": 10},
                timeout=15,
            )
        except Exception:
            return None
        if resp.status_code == 429:
            if attempt < _PP_RATE_LIMIT_RETRIES - 1:
                time.sleep(backoff)
                backoff *= 2
                continue
            return None
        try:
            resp.raise_for_status()
        except Exception:
            return None
        return _pp_best_prediction(resp.json().get("hits", []))
    return None


def _pp_get(params: dict) -> dict | None:
    """One myvariant.info query, retrying HTTP 429 with backoff; None if it failed."""
    backoff = _PP_RATE_LIMIT_BACKOFF
    for attempt in range(_PP_RATE_LIMIT_RETRIES):
        try:
            resp = _pp_session.get(_PP_API_URL, params=params, timeout=60)
        except Exception:
            return None
        if resp.status_code == 429:
            if attempt < _PP_RATE_LIMIT_RETRIES - 1:
                time.sleep(backoff)
                backoff *= 2
                continue
            return None
        try:
            resp.raise_for_status()
            return resp.json()
        except Exception:
            return None
    return None


def _hit_aa_values(hit: dict) -> tuple[set[str], set[str], set[int]]:
    """(reference residues, alternate residues, protein positions) a dbNSFP
    record lists, pooled across its transcripts -- the same values the
    single-variant query's field matches search."""
    refs: set[str] = set()
    alts: set[str] = set()
    positions: set[int] = set()
    for dbnsfp in _as_list(hit.get("dbnsfp") or {}):
        if not isinstance(dbnsfp, dict):
            continue
        for aa in _as_list(dbnsfp.get("aa") or {}):
            if not isinstance(aa, dict):
                continue
            refs.update(str(r) for r in _as_list(aa.get("ref")) if r is not None)
            alts.update(str(a) for a in _as_list(aa.get("alt")) if a is not None)
            for pos in _as_list(aa.get("pos")):
                try:
                    positions.add(int(pos))
                except (TypeError, ValueError):
                    pass
    return refs, alts, positions


def fetch_polyphen_batch(gene: str, mutations: list[str],
                         request_delay: float = _PP_REQUEST_DELAY) -> dict[str, tuple[str, str] | None]:
    """PolyPhen-2 HDIV predictions for many of one gene's mutations at once.

    Asks myvariant.info for every dbNSFP record of *gene* at up to
    _PP_BATCH_POSITIONS positions per request, then gives each mutation the
    records whose reference residue, alternate residue and (any transcript's)
    position match it -- the same records fetch_polyphen's one-variant query
    finds -- scored by _pp_best_prediction. A chunk over myvariant.info's
    1000-record page is split in half.

    Returns {mutation: (pred, score)}, or None for a mutation whose request
    failed (so the caller doesn't cache it). Stop codons and unparseable
    labels get ("", "") without a request, as in fetch_polyphen.
    """
    results: dict[str, tuple[str, str] | None] = {}
    by_position: dict[int, list[tuple[str, str, str]]] = {}
    for mutation in mutations:
        m = MUT_RE.match(mutation)
        if not m or m.group(3) == "*":
            results[mutation] = ("", "")
            continue
        by_position.setdefault(int(m.group(2)), []).append((mutation, m.group(1), m.group(3)))

    positions = sorted(by_position)
    queue = [positions[i:i + _PP_BATCH_POSITIONS] for i in range(0, len(positions), _PP_BATCH_POSITIONS)]
    while queue:
        chunk = queue.pop(0)
        query = f"dbnsfp.genename:{gene} AND dbnsfp.aa.pos:({' OR '.join(map(str, chunk))})"
        body = _pp_get({"q": query, "fields": _PP_BATCH_FIELDS, "size": _PP_BATCH_MAX_HITS})
        if request_delay:
            time.sleep(request_delay)
        if body is None:
            for pos in chunk:
                for mutation, _ref, _alt in by_position[pos]:
                    results[mutation] = None
            continue
        if body.get("total", 0) > _PP_BATCH_MAX_HITS:
            if len(chunk) > 1:
                half = len(chunk) // 2
                queue[:0] = [chunk[:half], chunk[half:]]
            else:
                # One position with more records than a page -- ask per variant
                for mutation, _ref, _alt in by_position[chunk[0]]:
                    results[mutation] = fetch_polyphen(gene, mutation)
            continue

        matched: dict[str, list[dict]] = {m: [] for pos in chunk for m, _r, _a in by_position[pos]}
        wanted = set(chunk)
        for hit in body.get("hits", []):
            refs, alts, hit_positions = _hit_aa_values(hit)
            for pos in hit_positions & wanted:
                for mutation, ref, alt in by_position[pos]:
                    if ref in refs and alt in alts:
                        matched[mutation].append(hit)
        for mutation, hits in matched.items():
            results[mutation] = _pp_best_prediction(hits)
    return results


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


def run_polyphen_phase(
    df: pd.DataFrame, mutation_cols: list[str] = _MUTATION_COLS,
    bare_mutation_cols: tuple[str, ...] = (),
    phase_idx: int = 1, num_phases: int = _NUM_PHASES,
    request_delay: float = _PP_REQUEST_DELAY,
) -> dict:
    """Phase 2: fetch PolyPhen-2 scores and tag mutation strings. Returns the full cache.

    *bare_mutation_cols* are columns holding a single bare mutation label with
    no distance suffix (e.g. anchor_mutation) -- scanned for fetching but not
    tagged inline (no entry format to tag); read back afterward via
    _pp_lookup_single once the cache is populated.

    Fetches gene by gene, many variants per request (fetch_polyphen_batch),
    sequentially with *request_delay* between calls rather than concurrently
    -- see the comment above _PP_REQUEST_DELAY.
    """
    print("\n── Phase 2: PolyPhen-2 pathogenicity scores ──")

    cache = _pp_load_cache()

    needed: set[tuple[str, str]] = set()
    for _, row in df.iterrows():
        gene = row.get("gene", "")
        if not gene:
            continue
        for col in mutation_cols:
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
        for col in bare_mutation_cols:
            bare = str(row.get(col, "") or "").replace("(isoform?)", "").strip()
            if bare and MUT_RE.match(bare) and (gene, bare) not in cache:
                needed.add((gene, bare))

    to_fetch = list(needed)
    print(f"{len(to_fetch)} unique (gene, mutation) pairs to fetch "
          f"({len(cache)} already in cache from prior runs)")

    if to_fetch:
        failed = 0
        total = len(to_fetch)
        by_gene: dict[str, list[str]] = {}
        for g, mut in to_fetch:
            by_gene.setdefault(g, []).append(mut)
        done = 0
        with tqdm(total=total, desc="Fetching PolyPhen-2 scores") as bar:
            for n_genes, (g, muts) in enumerate(sorted(by_gene.items()), 1):
                for mut, result in fetch_polyphen_batch(g, muts, request_delay).items():
                    if result is not None:
                        cache[(g, mut)] = result
                    else:
                        failed += 1
                done += len(muts)
                bar.update(len(muts))
                _emit_progress(phase_idx, done / total * 100, f"PolyPhen-2 scores: {done}/{total}", num_phases)
                if n_genes % 25 == 0:
                    _pp_save_cache(cache)  # keep progress if the run is interrupted
        if failed:
            print(f"  {failed} lookup(s) failed and will be retried on the next run (not cached)")
        _pp_save_cache(cache)
    else:
        _emit_progress(phase_idx, 100, "PolyPhen-2 scores: all cached", num_phases)

    tagged = 0
    for col in mutation_cols:
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


def _pp_lookup_single(mutation: str, gene: str, cache: dict[tuple[str, str], tuple[str, str]]) -> tuple[str, str]:
    """Look up PolyPhen (pred, score) for a single bare mutation label with no distance suffix."""
    clean = (mutation or "").replace("(isoform?)", "").strip()
    return cache.get((gene, clean), ("", ""))


# ═════════════════════════════════════════════════════════════════════════════
# Phase 3: Kinase Library predictions for phosphorylation sites
# ═════════════════════════════════════════════════════════════════════════════

_KIN_TOP_K = 5
_KIN_WINDOW = 7
_KIN_CACHE_FILE = PROJECT_ROOT / "data" / "cache" / "kinase_predictions.tsv"
_KIN_MAX_WORKERS = 6          # helper processes run in parallel
_KIN_CHUNK_SIZE = 250         # windows per helper process
# Kinase Library predictions run in a separate uv-managed environment -- see
# scripts/kinase_predictor.py for why (its numpy/pandas pins conflict with ours)
_KIN_HELPER = Path(__file__).resolve().parent / "kinase_predictor.py"
_KIN_HELPER_TIMEOUT = 3600    # seconds; the first call also builds the environment


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


def _kinase_helper_command() -> list[str] | None:
    uv = shutil.which("uv")
    return [uv, "run", "--quiet", "--script", str(_KIN_HELPER)] if uv else None


def _kinase_helper_env() -> dict[str, str]:
    env = dict(os.environ)
    env.pop("VIRTUAL_ENV", None)  # the project's env -- irrelevant to the helper's own
    env["PYTHONIOENCODING"] = "utf-8"
    return env


def check_kinase_helper() -> str | None:
    """Build (first time) and verify the Kinase Library helper's environment.
    Returns None if it works, else a description of what went wrong."""
    command = _kinase_helper_command()
    if command is None:
        return ("uv isn't on PATH -- kinase predictions run in their own uv-managed environment "
                "(see scripts/kinase_predictor.py)")
    try:
        proc = subprocess.run(command + ["--check"], capture_output=True, text=True, encoding="utf-8",
                              env=_kinase_helper_env(), timeout=_KIN_HELPER_TIMEOUT)
    except (OSError, subprocess.TimeoutExpired) as exc:
        return f"couldn't start the Kinase Library helper ({exc})"
    if proc.returncode != 0 or '"ok": true' not in proc.stdout:
        last = (proc.stderr.strip().splitlines() or ["no error output"])[-1]
        return f"the Kinase Library helper's environment failed to set up ({last})"
    return None


def _predict_chunk(windows: list[str]) -> tuple[dict[str, str], dict[str, str]]:
    """One helper process over *windows*: ({window: prediction}, {window: error})."""
    command = _kinase_helper_command()
    predictions: dict[str, str] = {}
    errors: dict[str, str] = {}
    if command is None:
        return predictions, {w: "uv isn't on PATH" for w in windows}
    try:
        proc = subprocess.run(command, input="\n".join(windows) + "\n", capture_output=True, text=True,
                              encoding="utf-8", env=_kinase_helper_env(), timeout=_KIN_HELPER_TIMEOUT)
        output, stderr = proc.stdout, proc.stderr
    except (OSError, subprocess.TimeoutExpired) as exc:
        output, stderr = "", str(exc)
    for line in output.splitlines():
        try:
            result = json.loads(line)
        except ValueError:
            continue
        if "prediction" in result:
            predictions[result["window"]] = result["prediction"]
        elif "error" in result:
            errors[result["window"]] = result["error"]
    last = (stderr.strip().splitlines() or ["no output"])[-1]
    for window in windows:
        if window not in predictions and window not in errors:
            errors[window] = f"the helper produced no result ({last})"
    return predictions, errors


def predict_kinase_windows(windows: list[str], progress=None) -> tuple[dict[str, str], dict[str, str]]:
    """Kinase Library top-5 predictions for many windows, spread over
    _KIN_MAX_WORKERS helper processes. Returns ({window: prediction},
    {window: error}); *progress(done, total)* is called as chunks finish."""
    chunks = [windows[i:i + _KIN_CHUNK_SIZE] for i in range(0, len(windows), _KIN_CHUNK_SIZE)]
    predictions: dict[str, str] = {}
    errors: dict[str, str] = {}
    done = 0
    with ThreadPoolExecutor(max_workers=_KIN_MAX_WORKERS) as pool:
        futures = {pool.submit(_predict_chunk, chunk): chunk for chunk in chunks}
        for future in tqdm(as_completed(futures), total=len(futures), desc="Predicting kinases (chunks)"):
            chunk_predictions, chunk_errors = future.result()
            predictions.update(chunk_predictions)
            errors.update(chunk_errors)
            done += len(futures[future])
            if progress:
                progress(done, len(windows))
    return predictions, errors


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

    # Resolve each row's window first (without predicting) so identical windows
    # shared across rows are only predicted once, not once per row
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
        problem = check_kinase_helper()
        if problem:
            failed_windows = set(to_predict)
            first_error = problem
        else:
            predictions, errors = predict_kinase_windows(
                to_predict,
                progress=lambda done, total: _emit_progress(2, done / total * 100,
                                                            f"Kinase predictions: {done}/{total}"),
            )
            cache.update(predictions)
            failed_windows = set(errors)
            first_error = next(iter(errors.values()), "")
        if failed_windows:
            print(f"  Warning: kinase prediction failed for {len(failed_windows)} of {len(to_predict)} "
                  f"new windows ({first_error}) -- those sites get no prediction this run; they "
                  f"aren't cached, so they're retried next run")
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


_WITHIN_COL, _BEYOND_COL = _MUTATION_COLS


def _filter_split_mut_cols(df: pd.DataFrame, exclude_codes: set[str], ref_pos_col: str,
                           extra_cols: tuple[str, ...] = ()) -> None:
    """Filter the within-5/more-than-5 mutation strings in place and recompute
    their count, unique-position, and morethan5_linear_distance columns.

    ref_pos_col names the column whose first integer is the reference position
    linear distances are measured from (ptm_site or anchor_position).
    extra_cols are filtered too, but have no derived columns to recompute.
    """
    for col in (_WITHIN_COL, _BEYOND_COL, *extra_cols):
        if col in df.columns:
            df[col] = df[col].fillna("").apply(
                lambda s: _filter_mut_str(s, exclude_codes)
            )

    def _count(s: str) -> int:
        return len([e for e in s.split(", ") if e.strip()]) if s else 0

    def _uniq_pos(s: str) -> int:
        return len(set(_positions_from_str(s)))

    for prefix, col in [("within_5", _WITHIN_COL), ("more_than_5", _BEYOND_COL)]:
        df[f"mutation_count_{prefix}_positions"] = df[col].apply(_count)
        df[f"unique_mutation_position_count_{prefix}_positions"] = df[col].apply(_uniq_pos)

    # Recompute morethan5_linear_distance from filtered beyond-5 string
    if "morethan5_linear_distance" in df.columns and ref_pos_col in df.columns:
        def _linear_dists(row) -> str:
            ref_m = re.search(r"(\d+)", str(row.get(ref_pos_col, "")))
            if not ref_m:
                return ""
            ref_pos = int(ref_m.group(1))
            return ",".join(str(abs(p - ref_pos)) for p in _positions_from_str(row[_BEYOND_COL]))
        df["morethan5_linear_distance"] = df.apply(_linear_dists, axis=1)


def _has_remaining_mutations(df: pd.DataFrame) -> pd.Series:
    """Rows with at least one mutation left in either the within-5 or more-than-5 string."""
    return (
        (df[_WITHIN_COL].fillna("").str.len() > 0) |
        (df[_BEYOND_COL].fillna("").str.len() > 0)
    )


def apply_polyphen_filter(df: pd.DataFrame, exclude_classes: list[str]) -> pd.DataFrame:
    """Remove mutations of excluded PP classes from the wide-format proximity DB.

    Filters mutation strings, recomputes count/distance columns, and drops PTM
    rows with no qualifying mutations left. *_total_patient_count columns are
    left as pre-filter totals -- they can't be recomputed here.
    """
    exclude_codes = {_PP_CODE_MAP[c] for c in exclude_classes if c in _PP_CODE_MAP}
    if not exclude_codes:
        return df

    print(f"\nApplying PolyPhen filter — excluding: {', '.join(exclude_classes)}")
    print("  Note: *_total_patient_count columns retain pre-filter totals")

    _filter_split_mut_cols(df, exclude_codes, "ptm_site",
                           extra_cols=("confirmed_disrupting_mutations",))

    # Recompute mutation_at_ptm_site from filtered within-5 string
    if "mutation_at_ptm_site" in df.columns and "ptm_site" in df.columns:
        def _at_ptm(row) -> str:
            ptm_m = re.search(r"(\d+)", str(row.get("ptm_site", "")))
            if not ptm_m:
                return "no"
            return "yes" if int(ptm_m.group(1)) in _positions_from_str(row[_WITHIN_COL]) else "no"
        df["mutation_at_ptm_site"] = df.apply(_at_ptm, axis=1)

    before = len(df)
    df = df[_has_remaining_mutations(df)].reset_index(drop=True)
    removed = before - len(df)
    print(f"  Removed {removed} PTM rows with no qualifying mutations; "
          f"{len(df)} rows remaining")
    return df


def apply_polyphen_filter_cluster(df: pd.DataFrame, exclude_classes: list[str]) -> pd.DataFrame:
    """Remove mutations of excluded PP classes from the wide-format cluster DB.

    Filters the within-5/more-than-5 neighbor strings and recomputes their
    count/position/linear-distance columns, then drops the whole anchor row if
    the anchor's own class is excluded or no neighbors remain in either group.
    *_total_patient_count columns are left as pre-filter totals -- same caveat
    as apply_polyphen_filter.
    """
    exclude_codes = {_PP_CODE_MAP[c] for c in exclude_classes if c in _PP_CODE_MAP}
    if not exclude_codes:
        return df

    print(f"\nApplying PolyPhen filter — excluding: {', '.join(exclude_classes)}")
    print("  Note: *_total_patient_count columns retain pre-filter totals")

    _filter_split_mut_cols(df, exclude_codes, "anchor_position")

    before = len(df)
    anchor_excluded = df["anchor_polyphen_class"].isin(exclude_classes)
    df = df[~anchor_excluded & _has_remaining_mutations(df)].reset_index(drop=True)
    removed = before - len(df)
    print(f"  Removed {removed} anchor rows (excluded anchor class or no qualifying neighbors); "
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


def run_aiupred_phase(
    df: pd.DataFrame, phase_idx: int = 3, num_phases: int = _NUM_PHASES,
) -> dict[str, dict[str, dict[int, float]]]:
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
                _emit_progress(phase_idx, done / max(total_fetches, 1) * 100,
                               f"AIUPred {uid}: {done}/{total_fetches}", num_phases)
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

    Queries /entry/interpro/ (curated entries only) rather than /entry/all/,
    which returns many overlapping near-duplicates from individual member
    databases (Pfam, CDD, PROSITE, etc. each flagging the same domain).

    Returns a list of {"name", "type", "start", "end"} dicts, one per
    (entry, fragment) pair -- an entry can have multiple discontinuous fragments.
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


def run_interpro_phase(
    df: pd.DataFrame, phase_idx: int = 4, num_phases: int = _NUM_PHASES,
) -> dict[str, list[dict]]:
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
                _emit_progress(phase_idx, i / max(len(need), 1) * 100, f"InterPro {uid}: {i}/{len(need)}", num_phases)
        _interpro_save_cache(cache)

    return {uid: cache.get(uid, []) for uid in uniprots}


def find_domain_at_position(entries: list[dict], position: int | None) -> str:
    """Return "name (type, start-end)" for every InterPro entry containing
    *position*, semicolon-joined (a residue can fall inside more than one
    entry, e.g. a domain nested in a broader superfamily call).
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


def annotate_cluster_long_format(
    df: pd.DataFrame,
    pp_cache: dict,
    disorder_maps: dict | None = None,
    domain_maps: dict | None = None,
) -> None:
    """Fill mutation-level annotation columns in the long-format cluster table in-place.

    Cluster-mode counterpart of annotate_long_format's PolyPhen/AIUPred/domain
    blocks (no 14-3-3/kinase -- those are PTM-site-specific). Simpler than the
    PTM version since mutation_position is already a column here, not something
    to regex out of the mutation label.
    """
    pp_scores, pp_classes = [], []
    _ATYPES = ("general", "binding")
    mut_disorder: dict[str, list] = {t: [] for t in _ATYPES}
    mut_domains = []

    for _, row in df.iterrows():
        uid = str(row.get("UniProt", "") or "")
        gene = str(row.get("gene", "") or "")
        mutation = str(row.get("mutation", "") or "")
        clean_mut = mutation.replace("(isoform?)", "")
        try:
            raw_pos = str(row.get("mutation_position", "")).strip()
            mut_pos = int(float(raw_pos)) if raw_pos else None
        except ValueError:
            mut_pos = None

        pred, pp_score = pp_cache.get((gene, clean_mut), ("", ""))
        pp_scores.append(pp_score)
        pp_classes.append(_PP_CLASS.get(pred, ""))

        for t in _ATYPES:
            if disorder_maps is not None:
                ps = disorder_maps[t].get(uid, {})
                mut_disorder[t].append(f"{ps[mut_pos]:.3f}" if mut_pos in ps else "")
            else:
                mut_disorder[t].append("")

        if domain_maps is not None:
            mut_domains.append(find_domain_at_position(domain_maps.get(uid, []), mut_pos))
        else:
            mut_domains.append("")

    df["polyphen_score"] = pp_scores
    df["polyphen_class"] = pp_classes
    for t in _ATYPES:
        df[f"mut_aiupred_{t}"] = mut_disorder[t]
    df["mut_is_disordered"] = [
        "yes" if v and float(v) > 0.5 else "no" for v in mut_disorder["general"]
    ]
    df["mut_is_binding"] = [
        "yes" if v and float(v) > 0.5 else "no" for v in mut_disorder["binding"]
    ]
    df["mutation_domain"] = mut_domains


# ═════════════════════════════════════════════════════════════════════════════
# Main
# ═════════════════════════════════════════════════════════════════════════════

def _annotate_ptm_proximity(output_dir: Path, pp_exclude: list[str],
                            ptm_source: str = DEFAULT_PTM_SOURCE) -> None:
    """Run all 5 annotation phases on the PTM Proximity proximity database."""
    ptm_paths = ptm_output_paths(output_dir, ptm_source)
    proximity_db = ptm_paths["db"]
    print(f"Reading proximity DB: {proximity_db}")
    df = pd.read_csv(proximity_db, sep="\t", encoding="utf-16", dtype=str,
                     keep_default_na=False)
    print(f"{len(df)} rows, {df['UniProt'].nunique()} unique proteins\n")

    # The five phases run concurrently, each on its own copy of the columns it
    # reads; their new/updated columns are merged back below in the original
    # phase order, so the output matches running them one after another.
    df_1433 = df[["UniProt", "ptm_site"]].copy()
    df_pp = df[["gene", *[c for c in _MUTATION_COLS if c in df.columns]]].copy()
    df_kin = df[["UniProt", "ptm_site", "ptm_type"]].copy()
    results = _run_phases({
        "Phase 1 (14-3-3)": (0, _NUM_PHASES, lambda: run_1433_phase(df_1433)),
        "Phase 2 (PolyPhen-2)": (1, _NUM_PHASES, lambda: run_polyphen_phase(df_pp)),
        "Phase 3 (Kinase)": (2, _NUM_PHASES, lambda: run_kinase_phase(df_kin)),
        "Phase 4 (AIUPred)": (3, _NUM_PHASES, lambda: run_aiupred_phase(df[["UniProt"]].copy())),
        "Phase 5 (InterPro)": (4, _NUM_PHASES, lambda: run_interpro_phase(df[["UniProt"]].copy())),
    })
    score_maps, confirmed_sites = results["Phase 1 (14-3-3)"]
    for col in ("1433pred_binding_site", "1433pred_consensus", "1433_confirmed_site", "1433_confirmed_pmid"):
        df[col] = df_1433[col]
    pp_cache = results["Phase 2 (PolyPhen-2)"]
    for col in _MUTATION_COLS:
        if col in df_pp.columns:
            df[col] = df_pp[col]
    seq_maps, kin_cache = results["Phase 3 (Kinase)"]
    df["kinase_predictions"] = df_kin["kinase_predictions"]

    disorder_maps = results["Phase 4 (AIUPred)"]
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

    domain_maps = results["Phase 5 (InterPro)"]
    col_vals = []
    for _, row in df.iterrows():
        uid = str(row.get("UniProt", "") or "")
        ptm_site = str(row.get("ptm_site", "") or "")
        m = SITE_RE.match(ptm_site.strip()) if ptm_site else None
        pos = int(m.group(2)) if m else None
        col_vals.append(find_domain_at_position(domain_maps.get(uid, []), pos))
    df["ptm_domain"] = col_vals

    if pp_exclude:
        df = apply_polyphen_filter(df, pp_exclude)

    df.to_csv(proximity_db, sep="\t", index=False, encoding="utf-16")
    print(f"\nUpdated proximity DB written to: {proximity_db}")

    long_db = ptm_paths["long"]
    if long_db.exists():
        print(f"\nAnnotating long-format DB: {long_db}")
        df_long = pd.read_csv(long_db, sep="\t", encoding="utf-16", dtype=str,
                              keep_default_na=False)
        print(f"{len(df_long)} rows")
        annotate_long_format(
            df_long, score_maps, confirmed_sites, pp_cache, seq_maps, kin_cache,
            disorder_maps, domain_maps,
        )
        if pp_exclude:
            before = len(df_long)
            df_long = df_long[
                ~df_long["polyphen_class"].isin(pp_exclude)
            ].reset_index(drop=True)
            print(f"  Long format: removed {before - len(df_long)} rows by PolyPhen filter")
        df_long.to_csv(long_db, sep="\t", index=False, encoding="utf-16")
        print(f"Updated long-format DB written to: {long_db}")
    else:
        print(f"\nNo long-format DB found at {long_db} — skipping")


def _annotate_mutation_clustering(output_dir: Path, pp_exclude: list[str]) -> None:
    """Run PolyPhen-2/AIUPred/InterPro (no 14-3-3/kinase) on the Mutation Clustering database."""
    cluster_db = output_dir / "mutation_cluster_db.tsv"
    print(f"Reading mutation cluster DB: {cluster_db}")
    df = pd.read_csv(cluster_db, sep="\t", encoding="utf-16", dtype=str,
                     keep_default_na=False)
    print(f"{len(df)} rows, {df['UniProt'].nunique()} unique proteins\n")

    # PolyPhen, AIUPred and InterPro run concurrently -- see _annotate_ptm_proximity
    df_pp = df[["gene", "anchor_mutation", *[c for c in _MUTATION_COLS if c in df.columns]]].copy()
    results = _run_phases({
        "Phase 2 (PolyPhen-2)": (0, 3, lambda: run_polyphen_phase(
            df_pp, bare_mutation_cols=("anchor_mutation",), phase_idx=0, num_phases=3)),
        "Phase 4 (AIUPred)": (1, 3, lambda: run_aiupred_phase(df[["UniProt"]].copy(), phase_idx=1, num_phases=3)),
        "Phase 5 (InterPro)": (2, 3, lambda: run_interpro_phase(df[["UniProt"]].copy(), phase_idx=2, num_phases=3)),
    })
    pp_cache = results["Phase 2 (PolyPhen-2)"]
    for col in _MUTATION_COLS:
        if col in df_pp.columns:
            df[col] = df_pp[col]
    anchor_classes, anchor_scores = [], []
    for _, row in df.iterrows():
        pred, score = _pp_lookup_single(row.get("anchor_mutation", ""), row.get("gene", ""), pp_cache)
        anchor_classes.append(_PP_CLASS.get(pred, ""))
        anchor_scores.append(score)
    df["anchor_polyphen_class"] = anchor_classes
    df["anchor_polyphen_score"] = anchor_scores

    disorder_maps = results["Phase 4 (AIUPred)"]
    for atype in ("general", "binding"):
        col_vals = []
        for _, row in df.iterrows():
            uid = str(row.get("UniProt", "") or "")
            try:
                raw_pos = str(row.get("anchor_position", "")).strip()
                pos = int(float(raw_pos)) if raw_pos else None
            except ValueError:
                pos = None
            ps = disorder_maps[atype].get(uid, {})
            col_vals.append(f"{ps[pos]:.3f}" if pos in ps else "")
        df[f"anchor_aiupred_{atype}"] = col_vals
    df["anchor_is_disordered"] = df["anchor_aiupred_general"].apply(
        lambda v: "yes" if v and float(v) > 0.5 else "no"
    )
    df["anchor_is_binding"] = df["anchor_aiupred_binding"].apply(
        lambda v: "yes" if v and float(v) > 0.5 else "no"
    )

    domain_maps = results["Phase 5 (InterPro)"]
    col_vals = []
    for _, row in df.iterrows():
        uid = str(row.get("UniProt", "") or "")
        try:
            raw_pos = str(row.get("anchor_position", "")).strip()
            pos = int(float(raw_pos)) if raw_pos else None
        except ValueError:
            pos = None
        col_vals.append(find_domain_at_position(domain_maps.get(uid, []), pos))
    df["anchor_domain"] = col_vals

    if pp_exclude:
        df = apply_polyphen_filter_cluster(df, pp_exclude)

    df.to_csv(cluster_db, sep="\t", index=False, encoding="utf-16")
    print(f"\nUpdated mutation cluster DB written to: {cluster_db}")

    long_db = output_dir / "mutation_cluster_long.tsv"
    if long_db.exists():
        print(f"\nAnnotating long-format cluster DB: {long_db}")
        df_long = pd.read_csv(long_db, sep="\t", encoding="utf-16", dtype=str,
                              keep_default_na=False)
        print(f"{len(df_long)} rows")
        annotate_cluster_long_format(df_long, pp_cache, disorder_maps, domain_maps)
        if pp_exclude:
            before = len(df_long)
            kept_anchors = set(zip(df["UniProt"], df["anchor_mutation"]))
            anchor_kept = df_long.apply(
                lambda r: (r["UniProt"], r["anchor_mutation"]) in kept_anchors, axis=1
            )
            df_long = df_long[
                anchor_kept & ~df_long["polyphen_class"].isin(pp_exclude)
            ].reset_index(drop=True)
            print(f"  Long format: removed {before - len(df_long)} rows by PolyPhen filter")
        df_long.to_csv(long_db, sep="\t", index=False, encoding="utf-16")
        print(f"Updated long-format cluster DB written to: {long_db}")
    else:
        print(f"\nNo long-format cluster DB found at {long_db} — skipping")


def main() -> None:
    """Parse CLI args and annotate the selected pipeline mode's output database."""
    import argparse
    parser = argparse.ArgumentParser(description="Annotate the pipeline's output database.")
    parser.add_argument(
        "--mode",
        choices=["ptm-proximity", "mutation-clustering"],
        default="ptm-proximity",
        help="Which pipeline mode's output to annotate (default: ptm-proximity)",
    )
    parser.add_argument(
        "--ptm-source",
        choices=PTM_SOURCES,
        default=DEFAULT_PTM_SOURCE,
        help="Which PTM source's ptm-proximity output to annotate: 'ptmd' (default) or 'psp'",
    )
    parser.add_argument(
        "--output-dir",
        default=str(DEFAULT_OUTPUT_DIR),
        help="Directory containing the output database (default: Output/)",
    )
    parser.add_argument(
        "--pp-exclude",
        nargs="+",
        choices=["benign", "possibly_damaging", "probably_damaging"],
        default=[],
        metavar="CLASS",
        help="Exclude mutations with these PolyPhen-2 classes from the output",
    )
    parser.add_argument(
        "--prefetch-1433",
        action="store_true",
        help="Only fetch 14-3-3-Pred results into the cache for the ptm-proximity step-1 table's proteins "
             "(run in the background during steps 2-3), then exit",
    )
    args = parser.parse_args()

    if args.prefetch_1433:
        prefetch_1433(hotspots_tsv_path(PROJECT_ROOT, "ptm-proximity", args.ptm_source))
        return

    output_dir = Path(args.output_dir)
    if args.mode == "mutation-clustering":
        _annotate_mutation_clustering(output_dir, args.pp_exclude)
    else:
        _annotate_ptm_proximity(output_dir, args.pp_exclude, args.ptm_source)


if __name__ == "__main__":
    main()
