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
    """Parse a CIF file and return the first chain as a biotite AtomArray, or None."""
    try:
        cif = CIFFile.read(str(model_file))
        structure = get_structure(cif, model=1)
    except Exception:
        return None

    if structure is None or len(structure) == 0:
        return None

    chain_ids = list(dict.fromkeys(structure.chain_id))
    if not chain_ids:
        return None

    return structure[structure.chain_id == chain_ids[0]]


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

PTM_PROXIMITY_STEPS = [
    "Filter and merge PTMD + COSMIC data",
    "Download AlphaFold CIF models and PAE files",
    "Find nearby mutations and compute distances",
    "Annotate 14-3-3-Pred binding-site predictions",
    "Annotate mutations with PolyPhen-2 scores",
]

MUTATION_CLUSTERING_STEPS = [
    "Filter COSMIC hotspot mutations",
    "Download AlphaFold CIF models and PAE files",
    "Find mutation clusters in 3D space",
]
