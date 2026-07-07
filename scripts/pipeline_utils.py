"""Shared constants, regexes, and helper functions used across pipeline scripts.

Centralises definitions that were previously duplicated in multiple scripts
so that changes only need to be made in one place.
"""
from __future__ import annotations

import json
import re
from pathlib import Path

import numpy as np
from biotite.structure.io.pdbx import CIFFile, get_structure  # type: ignore[import-untyped]


def project_root(script_file: str) -> Path:
    """Derive the project root from any script's ``__file__``."""
    return Path(script_file).resolve().parent.parent


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
    ("Annotate results (14-3-3, PolyPhen-2, kinase, AIUPred predictions)",
     "Annotating results (14-3-3, PolyPhen-2, kinase, AIUPred predictions)"),
]

MUTATION_CLUSTERING_STEPS = [
    ("Filter COSMIC hotspot mutations",
     "Filtering COSMIC hotspot mutations"),
    ("Download AlphaFold CIF models and PAE files",
     "Downloading AlphaFold CIF models and PAE files"),
    ("Find mutation clusters in 3D space",
     "Finding mutation clusters in 3D space"),
]
