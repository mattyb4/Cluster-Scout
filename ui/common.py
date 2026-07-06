"""Shared constants and small helper functions used across the app's UI mixins.

Mirrors the role scripts/pipeline_utils.py plays for the backend pipeline
scripts: one shared-utilities module, imported by every ui/*.py file.
"""
from __future__ import annotations

import json
import re
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
SCRIPTS_DIR = PROJECT_ROOT / "scripts"
OUTPUT_DIR = PROJECT_ROOT / "Output"

sys.path.insert(0, str(SCRIPTS_DIR))
from pipeline_utils import (  # noqa: E402
    PTM_PROXIMITY_STEPS, MUTATION_CLUSTERING_STEPS,
    input_dir, resolve_input_file, extract_uniprot_from_cif,
    COSMIC_INPUT_DIR, PTMD_INPUT_DIR, INTERACTORS_1433_INPUT_DIR,
)

_INPUT_FOLDERS: dict[str, tuple[Path, tuple[str, ...], str]] = {
    "COSMIC": (
        input_dir(PROJECT_ROOT, COSMIC_INPUT_DIR),
        (".tsv",),
        "COSMIC Mutant Census TSV",
    ),
    "PTMD": (
        input_dir(PROJECT_ROOT, PTMD_INPUT_DIR),
        (".tsv",),
        "PTMD disease-associated PTMs TSV",
    ),
    "14-3-3": (
        input_dir(PROJECT_ROOT, INTERACTORS_1433_INPUT_DIR),
        (".xlsx", ".xls"),
        "14-3-3 confirmed interactors Excel",
    ),
}

_GRAY = "gray"
_BLUE = "#3a86ff"
_GREEN = "#2ecc71"
_RED = "#e74c3c"
_YELLOW = "#f1c40f"

_CACHE_DIR = PROJECT_ROOT / "data" / "cache"
_CACHE_ITEMS = [
    # (step_label, display_name, path, is_dir)
    ("Step 1", "UniProt gene mapping",   _CACHE_DIR / "uniprot_gene_mapping.tsv",    False),
    ("Step 1", "Gene → UniProt mapping", _CACHE_DIR / "gene_to_uniprot_mapping.tsv", False),
    ("Step 1", "Isoform safe lengths",   _CACHE_DIR / "isoform_safe_lengths.tsv",    False),
    ("Step 4", "14-3-3 predictions",     _CACHE_DIR / "1433pred",                    True),
    ("Step 4", "PolyPhen-2 scores",      _CACHE_DIR / "polyphen.tsv",                False),
    ("Step 4", "Kinase predictions",     _CACHE_DIR / "kinase_predictions.tsv",      False),
    ("Step 4", "AIUPred disorder",       _CACHE_DIR / "aiupred_disorder.tsv",        False),
]


def _cache_entry_count(path: Path, is_dir: bool) -> str:
    """Return a human-readable entry count string for a cache path."""
    if is_dir:
        if not path.is_dir():
            return "empty"
        n = sum(1 for f in path.iterdir() if f.is_file())
        return f"{n:,} entries" if n else "empty"
    if not path.exists():
        return "empty"
    try:
        with path.open(encoding="utf-8", errors="replace") as fh:
            n = sum(1 for _ in fh) - 1  # subtract header row
        return f"{max(0, n):,} entries" if n > 0 else "empty"
    except Exception:
        return "?"


# Results-tab helpers
_MUT_ENTRY_RE = re.compile(
    r"([A-Z]\d+[A-Z*](?:\(isoform\?\))?)"
    r"(?:\(PP:([DPB]),([0-9.]*)\))?"
    r"-([0-9.]+)Å"
    r"(?:\(PAE:([0-9.]+)\))?"
)
_PP_LABEL = {"D": "probably_damaging", "P": "possibly_damaging", "B": "benign"}
_PP_COLORS = {
    "probably_damaging": _RED,
    "possibly_damaging": _YELLOW,
    "benign": _GREEN,
}
_PTM_MARKER_COLOR = _BLUE
_NEEDLE_DEFAULT_COLOR = "#888888"

_PTM_TV_COLS = [
    ("#",              "#col",   32),
    ("UniProt",        "uniprot",70),
    ("Gene",           "gene",   58),
    ("PTM Site",       "site",   65),
    ("Type",           "type",  110),
    ("≤5 pos",         "near",   52),
    (">5 pos",         "far",    52),
    ("Unique pos",     "total",  68),
    ("Patients",       "pts",    65),
    ("COSMIC",         "cosmic", 65),
    ("At PTM",         "atptm",  52),
    ("14-3-3",         "pred14", 58),
    ("Disordered?",    "disord", 78),
    ("Binding?",       "bind",   58),
    ("Max lin. dist.", "maxlin", 90),
]

_MUT_TV_COLS = [
    ("#",          "#col",   32),
    ("Mutation",   "mut",    80),
    ("Seq dist",   "seqd",   62),
    ("Dist (Å)",   "dist",   62),
    ("Binding?",   "isbnd",  58),
    ("Disordered?","isdis",  78),
    ("PP Class",   "ppc",   115),
    ("Mut pLDDT",  "mpld",   72),
    ("PAE",        "pae",    48),
    ("Patients",   "pts",    62),
]


def _fmt_time(seconds: float) -> str:
    s = int(seconds)
    if s < 60:
        return f"{s}s"
    return f"{s // 60}m {s % 60:02d}s"


RUNTIMES_FILE = OUTPUT_DIR / "logs" / "pipeline_runtimes.json"
_CIF_DIR = PROJECT_ROOT / "cif_models"


def _detect_run_type() -> str:
    """Return 'cold' if key resources are missing, 'warm' if they are cached."""
    has_cifs = _CIF_DIR.exists() and any(_CIF_DIR.glob("*/*.cif"))
    has_cache = (_CACHE_DIR / "uniprot_gene_mapping.tsv").exists()
    return "warm" if (has_cifs and has_cache) else "cold"


def _load_runtimes(mode: str, run_type: str) -> list[float] | None:
    try:
        data = json.loads(RUNTIMES_FILE.read_text())
        return data.get(mode, {}).get(run_type)
    except Exception:
        return None


def _save_runtimes(mode: str, run_type: str, times: list[float]) -> None:
    try:
        data: dict = {}
        if RUNTIMES_FILE.exists():
            data = json.loads(RUNTIMES_FILE.read_text())
        data.setdefault(mode, {})[run_type] = times
        RUNTIMES_FILE.parent.mkdir(parents=True, exist_ok=True)
        RUNTIMES_FILE.write_text(json.dumps(data, indent=2))
    except Exception:
        pass
