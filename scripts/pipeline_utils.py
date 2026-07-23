"""Shared constants, regexes, and helper functions used across pipeline scripts.

Centralises definitions that were previously duplicated in multiple scripts
so that changes only need to be made in one place.
"""
from __future__ import annotations

import json
import re
from pathlib import Path

import numpy as np
import pandas as pd
from biotite.structure.io.pdbx import CIFFile, get_structure  # type: ignore[import-untyped]


def project_root(script_file: str) -> Path:
    """Derive the project root from any script's ``__file__``."""
    return Path(script_file).resolve().parent.parent


def fmt_time(seconds: float) -> str:
    """Format a duration as '12s' or '4m 05s', matching the app's runtime display."""
    s = int(seconds)
    if s < 60:
        return f"{s}s"
    return f"{s // 60}m {s % 60:02d}s"


# ── Input folder resolution ───────────────────────────────────────────────────

COSMIC_INPUT_DIR = "cosmic"
PTMD_INPUT_DIR = "ptmd"
INTERACTORS_1433_INPUT_DIR = "1433_interactors"


def input_dir(root: Path, subfolder: str) -> Path:
    """Return the input folder path for a given input type, creating it if needed."""
    d = root / "data" / "input" / subfolder
    d.mkdir(parents=True, exist_ok=True)
    return d


def resolve_input_file(
    folder: Path,
    extensions: tuple[str, ...] = (".tsv", ".csv", ".xlsx", ".xls"),
) -> Path:
    """Find the single input file in *folder*, erroring if zero or more than one are found."""
    if not folder.exists():
        raise FileNotFoundError(f"Input folder does not exist: {folder}")
    matches = [f for f in sorted(folder.iterdir())
               if f.is_file() and f.suffix.lower() in extensions]
    if len(matches) == 0:
        raise FileNotFoundError(
            f"No input file found in {folder}.\n"
            f"Place a file with one of these extensions in the folder: {', '.join(extensions)}"
        )
    if len(matches) > 1:
        names = [f.name for f in matches]
        raise RuntimeError(
            f"Multiple files found in {folder}: {names}\n"
            "Remove extras so only one input file remains."
        )
    return matches[0]


# ── Input file content validation ─────────────────────────────────────────────

# Only columns the pipeline scripts index directly (df["col"]); .get(...)-with-default columns are left out
COSMIC_REQUIRED_COLUMNS = (
    "GENE_SYMBOL", "MUTATION_AA", "COSMIC_SAMPLE_ID",
    "MUTATION_SOMATIC_STATUS", "TRANSCRIPT_ACCESSION",
)
PTMD_REQUIRED_COLUMNS = (
    "State", "UniProt", "Disease", "MutationSite", "Residue", "Position", "Type",
)
INTERACTORS_1433_REQUIRED_COLUMNS = ("Residue", "PMID")


def _peek_columns(path: Path, is_excel: bool) -> list[str]:
    """Read only the header row of *path*, so this stays fast even on multi-GB TSVs."""
    if is_excel:
        df = pd.read_excel(path, nrows=0)
    else:
        df = pd.read_csv(path, sep="\t", nrows=0)
    return [c.strip() for c in df.columns]


def _validate_columns(path: Path, required: tuple[str, ...], is_excel: bool = False) -> list[str]:
    """Check that *path* contains all of *required* as columns.

    Returns a list of human-readable problem descriptions; empty means valid.
    """
    try:
        columns = _peek_columns(path, is_excel)
    except Exception as exc:
        return [f"{path.name}: could not be read as a {'spreadsheet' if is_excel else 'TSV'} file ({exc})"]
    missing = [c for c in required if c not in columns]
    return [f"{path.name}: missing expected column '{c}'" for c in missing]


def validate_cosmic_file(path: Path) -> list[str]:
    """Validate that *path* looks like a COSMIC Mutant Census TSV."""
    return _validate_columns(path, COSMIC_REQUIRED_COLUMNS)


def validate_ptmd_file(path: Path) -> list[str]:
    """Validate that *path* looks like a PTMD disease-associated PTMs TSV."""
    return _validate_columns(path, PTMD_REQUIRED_COLUMNS)


def validate_1433_file(path: Path) -> list[str]:
    """Validate that *path* looks like a 14-3-3 confirmed interactors spreadsheet."""
    return _validate_columns(path, INTERACTORS_1433_REQUIRED_COLUMNS, is_excel=True)


# ── CIF metadata extraction ───────────────────────────────────────────────────

def extract_uniprot_from_cif(cif_path: Path) -> str | None:
    """Read the UniProt accession embedded in an AlphaFold CIF file, or None if not found."""
    try:
        with open(cif_path, encoding="utf-8", errors="replace") as f:
            for line in f:
                if line.startswith("_ma_target_ref_db_details.db_accession"):
                    return line.split()[-1].strip()
    except OSError:
        pass
    return None


# ── Amino-acid codes ──────────────────────────────────────────────────────────

AA3TO1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    "SEC": "U", "PYL": "O",
}


# ── Regex patterns ────────────────────────────────────────────────────────────

MUT_RE = re.compile(r"([A-Z])(\d+)([A-Z*])")
SITE_RE = re.compile(r"^([A-Z])(\d+)$")

COSMIC_SOMATIC_STATUSES = {
    "Confirmed somatic variant",
    "Reported in another cancer sample as somatic",
}


# ── AlphaFold CIF file helpers ────────────────────────────────────────────────

def _canonical_cif_re(uid: str) -> re.Pattern:
    return re.compile(rf"^AF-{re.escape(uid)}-F(\d+)-model_v\d+\.", re.IGNORECASE)


def find_canonical_cif(uniprot_dir: Path) -> Path | None:
    """Return the first canonical AlphaFold CIF for *uniprot_dir*, or None."""
    uid = uniprot_dir.name
    pat = _canonical_cif_re(uid)
    candidates = [p for p in sorted(uniprot_dir.glob("*.cif")) if pat.match(p.name)]
    return candidates[0] if candidates else None


def find_canonical_cifs(uniprot_dir: Path) -> list[Path]:
    """Return all canonical AlphaFold CIF fragments, sorted by fragment number."""
    uid = uniprot_dir.name
    pat = _canonical_cif_re(uid)
    hits = [(int(pat.match(p.name).group(1)), p)
            for p in uniprot_dir.glob("*.cif") if pat.match(p.name)]
    return [p for _, p in sorted(hits)]


def load_first_chain(model_file: Path):
    """Parse a CIF file and return the first chain as a biotite AtomArray, or None.

    Includes per-atom pLDDT in the ``b_factor`` attribute.
    """
    try:
        cif = CIFFile.read(str(model_file))
        structure = get_structure(cif, model=1, extra_fields=["b_factor"])
    except Exception:
        return None

    if structure is None or len(structure) == 0:
        return None

    chain_ids = list(dict.fromkeys(structure.chain_id))
    if not chain_ids:
        return None

    return structure[structure.chain_id == chain_ids[0]]


def get_plddt_map(chain) -> dict[int, float]:
    """Build a {residue_position: pLDDT} dict from a chain's CA atoms."""
    ca_mask = chain.atom_name == "CA"
    ca = chain[ca_mask]
    return {int(ca.res_id[i]): float(ca.b_factor[i]) for i in range(len(ca))}


def get_protein_length(uniprot_dir: Path) -> int | None:
    """Return the highest modeled residue position across all AlphaFold
    fragments for this protein, or None if no CIFs are found.

    AlphaFold fragment numbering is continuous in canonical UniProt coordinates
    (fragment 2 continues from where fragment 1 left off), so this is a
    reliable proxy for protein length without a separate UniProt lookup.
    """
    cif_files = find_canonical_cifs(uniprot_dir)
    if not cif_files:
        return None
    max_pos = 0
    for cif_file in cif_files:
        chain = load_first_chain(cif_file)
        if chain is None:
            continue
        ca_mask = chain.atom_name == "CA"
        ca = chain[ca_mask]
        if len(ca):
            max_pos = max(max_pos, int(ca.res_id.max()))
    return max_pos or None


def load_pae_matrix(uniprot_dir: Path):
    """Load the PAE matrix JSON for the canonical model, or return None."""
    uid = uniprot_dir.name
    pat = re.compile(rf"^AF-{re.escape(uid)}-F\d+-predicted_aligned_error_v\d+\.",
                     re.IGNORECASE)
    candidates = [p for p in sorted(uniprot_dir.glob("*.json")) if pat.match(p.name)]
    if not candidates:
        return None
    with candidates[0].open() as f:
        data = json.load(f)
    if isinstance(data, list):
        data = data[0]
    matrix = data.get("predicted_aligned_error")
    return np.array(matrix) if matrix else None


# ── Structural-hotspot significance testing (prototype) ──────────────────────
#
# Permutation-test + FDR significance for PTM-mutation 3D proximity, adapted
# from HotMAPS's method for cancer mutation hotspots (Tokheim et al. 2016,
# https://doi.org/10.1158/0008-5472.CAN-15-3190) and Benjamini-Hochberg FDR
# (Benjamini & Hochberg 1995).
#
# Null hypothesis per PTM site: if the same number of distinct mutated
# positions had instead been placed uniformly at random among the protein's
# structurally-resolved residues, how often would at least as many land within
# `cutoff` Angstroms of the site as were actually observed? Normalizes for
# protein size/shape, unlike a fixed Angstrom cutoff alone.
#
# Prototype status: standalone primitives, not wired into the production
# pipeline (3_find_nearby_mutations.py still uses the fixed-cutoff filter).
# See scripts/prototype_hotspot_significance.py for a runnable demo.

def sample_permutation_indices(
    n_residues: int, n_mutations: int, n_permutations: int, rng: np.random.Generator,
) -> np.ndarray:
    """Return an (n_permutations, n_mutations) array of residue indices, each row
    an independent random sample without replacement from range(n_residues) --
    i.e. n_permutations simulated "random mutation sets".

    Vectorized via argsort-of-random-keys instead of a per-trial
    ``rng.choice(..., replace=False)`` loop, since the latter is the dominant
    cost when generating thousands of samples. Safe to reuse across every PTM
    site in the same protein: n_mutations doesn't vary by site.
    """
    n_mutations = min(n_mutations, n_residues)
    random_keys = rng.random((n_permutations, n_residues))
    return np.argsort(random_keys, axis=1)[:, :n_mutations]


def permutation_pvalue(
    site_coord: np.ndarray,
    residue_coords: np.ndarray,
    observed_count: int,
    cutoff: float,
    sampled_idx: np.ndarray,
) -> tuple[float, np.ndarray]:
    """Empirical one-sided permutation p-value for one PTM site: how often
    does a random placement of the same number of mutations produce at
    least as many within `cutoff` of `site_coord` as were actually observed?

    Returns (p_value, null_counts); null_counts is the simulated null
    distribution, kept for inspection/plotting.
    """
    dists = np.linalg.norm(residue_coords - site_coord, axis=1)
    sampled_dists = dists[sampled_idx]
    null_counts = np.sum(sampled_dists <= cutoff, axis=1)
    n_permutations = len(null_counts)
    # +1 smoothing avoids a meaningless p=0 from finite resampling.
    # (Davison AC, Hinkley DV. "Bootstrap Methods and Their Application."
    # Cambridge University Press, 1997, section 4.2.)
    p = (np.sum(null_counts >= observed_count) + 1) / (n_permutations + 1)
    return float(p), null_counts


def benjamini_hochberg(pvalues) -> np.ndarray:
    """Benjamini-Hochberg FDR-adjusted q-values (Benjamini & Hochberg, 1995).

    Matches R's ``p.adjust(method="BH")`` / statsmodels' ``fdr_bh``.
    """
    p = np.asarray(pvalues, dtype=float)
    m = len(p)
    order = np.argsort(p)
    ranked = p[order] * m / (np.arange(m) + 1)
    q_sorted = np.minimum.accumulate(ranked[::-1])[::-1]  # enforce monotonicity
    q_sorted = np.clip(q_sorted, 0, 1)
    q = np.empty(m, dtype=float)
    q[order] = q_sorted
    return q


# ── Pipeline step labels (single source of truth) ────────────────────────────

# Each step is (panel_label, log_label):
#   panel_label — declarative, shown in the app's steps panel
#   log_label   — present tense, shown in the log while running
PTM_PROXIMITY_STEPS = [
    ("Filter and merge PTMD + COSMIC data",
     "Filtering and merging PTMD + COSMIC data - this may take a moment"),
    ("Download AlphaFold CIF models and PAE files",
     "Downloading AlphaFold CIF models and PAE files"),
    ("Find nearby mutations and compute distances",
     "Finding nearby mutations and computing distances"),
    ("Annotate results (14-3-3, PolyPhen-2, kinase, AIUPred, InterPro predictions)",
     "Annotating results (14-3-3, PolyPhen-2, kinase, AIUPred, InterPro predictions)"),
]

MUTATION_CLUSTERING_STEPS = [
    ("Filter COSMIC hotspot mutations",
     "Filtering COSMIC hotspot mutations"),
    ("Download AlphaFold CIF models and PAE files",
     "Downloading AlphaFold CIF models and PAE files"),
    ("Find mutation clusters in 3D space",
     "Finding mutation clusters in 3D space"),
]
