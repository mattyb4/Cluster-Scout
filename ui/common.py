"""Shared constants and small helper functions used across the app's UI mixins.

Mirrors the role scripts/pipeline_utils.py plays for the backend pipeline
scripts: one shared-utilities module, imported by every ui/*.py file.
"""
from __future__ import annotations

import json
import re
import sys
import tkinter as tk
from pathlib import Path

import customtkinter as ctk

PROJECT_ROOT = Path(__file__).resolve().parent.parent
SCRIPTS_DIR = PROJECT_ROOT / "scripts"
OUTPUT_DIR = PROJECT_ROOT / "Output"

sys.path.insert(0, str(SCRIPTS_DIR))
from pipeline_utils import (  # noqa: E402
    PTM_PROXIMITY_STEPS, MUTATION_CLUSTERING_STEPS,
    input_dir, resolve_input_file, extract_uniprot_from_cif,
    COSMIC_INPUT_DIR, PTMD_INPUT_DIR, INTERACTORS_1433_INPUT_DIR,
    COSMIC_SOMATIC_STATUSES, fmt_time as _fmt_time,
    validate_cosmic_file, validate_ptmd_file, validate_1433_file,
)

# The 14-3-3 confirmed-interactors file isn't listed here: unlike COSMIC/PTMD,
# it's small, rarely updated, and bundled with the app (see
# data/input/1433_interactors/) rather than something the user is expected to
# provide — scripts/4_annotate.py still reads it from that same folder via
# INTERACTORS_1433_INPUT_DIR, this just keeps it out of the Pipeline tab's
# input-file browse/status UI.
_INPUT_FOLDERS: dict[str, tuple[Path, tuple[str, ...], str, object]] = {
    "COSMIC": (
        input_dir(PROJECT_ROOT, COSMIC_INPUT_DIR),
        (".tsv",),
        "COSMIC Mutant Census TSV",
        validate_cosmic_file,
    ),
    "PTMD": (
        input_dir(PROJECT_ROOT, PTMD_INPUT_DIR),
        (".tsv",),
        "PTMD disease-associated PTMs TSV",
        validate_ptmd_file,
    ),
}


# ── Hover-tooltip help icons ──────────────────────────────────────────────────

class _Tooltip:
    """Hover pop-up bubble anchored to a single widget.

    A plain tk.Toplevel/tk.Label rather than CTk widgets: it's a short-lived,
    unmanaged popup outside the normal widget tree, so there's nothing to gain
    from CTk's theming machinery here — just fixed colors that match the
    app's dark theme (see app.py's `ctk.set_appearance_mode("dark")`).
    """
    _DELAY_MS = 400

    def __init__(self, widget, text: str):
        self._widget = widget
        self._text = text
        self._after_id: str | None = None
        self._popup: tk.Toplevel | None = None
        widget.bind("<Enter>", self._schedule)
        widget.bind("<Leave>", self._hide)

    def _schedule(self, _event=None) -> None:
        self._cancel()
        self._after_id = self._widget.after(self._DELAY_MS, self._show)

    def _cancel(self) -> None:
        if self._after_id is not None:
            self._widget.after_cancel(self._after_id)
            self._after_id = None

    def _show(self) -> None:
        if self._popup is not None:
            return
        x = self._widget.winfo_rootx() + self._widget.winfo_width() // 2
        y = self._widget.winfo_rooty() + self._widget.winfo_height() + 6

        self._popup = tk.Toplevel(self._widget)
        self._popup.wm_overrideredirect(True)
        self._popup.wm_geometry(f"+{x}+{y}")
        self._popup.attributes("-topmost", True)

        tk.Label(
            self._popup, text=self._text, justify="left",
            background="#2b2b2b", foreground="#dce4ee",
            font=("Segoe UI", 11), padx=10, pady=6,
            wraplength=280, borderwidth=1, relief="solid",
        ).pack()

    def _hide(self, _event=None) -> None:
        self._cancel()
        if self._popup is not None:
            self._popup.destroy()
            self._popup = None


def help_icon(parent, text: str) -> ctk.CTkLabel:
    """A small "?" badge that shows *text* in a hover tooltip.

    Pack/grid the returned label right next to whatever it explains.
    """
    badge = ctk.CTkLabel(
        parent, text="?", width=16, height=16, corner_radius=8,
        fg_color="gray30", text_color="#cfcfcf",
        font=ctk.CTkFont(size=10, weight="bold"),
    )
    _Tooltip(badge, text)
    return badge


_GRAY = "gray"
_BLUE = "#3a86ff"
_GREEN = "#2ecc71"
_RED = "#e74c3c"
_YELLOW = "#f1c40f"

# Ctrl+scroll UI zoom (see App._on_ctrl_scroll_zoom)
MIN_UI_SCALE = 0.6
MAX_UI_SCALE = 2.0
UI_SCALE_STEP = 0.1

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

# Full column registries for the Results-tab treeviews: (label, col_id, width,
# numeric, default_visible). default_visible=True columns are shown out of the
# box; the rest are available via the Columns picker (ResultsTabMixin
# ._open_column_picker). Order here is the fixed display order when visible.
_PTM_TV_COLS = [
    ("#",                    "#col",              32, True,  True),
    ("UniProt",              "uniprot",           70, False, True),
    ("Gene",                 "gene",              58, False, True),
    ("PTM Site",             "site",              65, False, True),
    ("Type",                 "type",             110, False, True),
    ("≤5 pos",               "near",              52, True,  True),
    (">5 pos",               "far",               52, True,  True),
    ("≤5 pos patients",      "near_pts",          90, True,  False),
    (">5 pos patients",      "far_pts",           90, True,  False),
    ("≤5 unique pos",        "near_unique",       90, True,  False),
    (">5 unique pos",        "far_unique",        90, True,  False),
    ("Unique pos",           "total",             68, True,  True),
    ("Patients",             "pts",               65, True,  True),
    ("COSMIC",               "cosmic",            65, True,  True),
    ("At PTM",               "atptm",             52, False, True),
    ("Confirmed disrupting", "confirmed_disrupt",110, False, False),
    ("PTM diseases",         "diseases",         140, False, False),
    ("14-3-3",               "pred14",            58, False, True),
    ("14-3-3 consensus",     "pred14_consensus",  90, False, False),
    ("14-3-3 confirmed",     "conf14",            90, False, False),
    ("14-3-3 PMID",          "conf14_pmid",       90, False, False),
    ("Kinases",              "kinases",          160, False, False),
    ("PTM AIUPred gen.",     "aiupred_gen",       90, True,  False),
    ("PTM AIUPred bind.",    "aiupred_bind",      90, True,  False),
    ("Disordered?",          "disord",            78, False, True),
    ("Binding?",             "bind",              58, False, True),
    ("Max lin. dist.",       "maxlin",            90, True,  True),
    ("≤5 pos mutations",     "near_muts_raw",    200, False, False),
    (">5 pos mutations",     "far_muts_raw",     200, False, False),
    ("Linear distances",     "lin_dist_raw",     150, False, False),
]

_MUT_TV_COLS = [
    ("#",                    "#col",              32, True,  True),
    ("Mutation",             "mut",               80, False, True),
    ("Seq dist",             "seqd",              62, True,  True),
    ("Dist (Å)",             "dist",              62, True,  True),
    ("Binding?",             "isbnd",             58, False, True),
    ("Disordered?",          "isdis",             78, False, True),
    ("PP Class",             "ppc",              115, False, True),
    ("PP Score",             "pps",               62, True,  False),
    ("Mut pLDDT",            "mpld",              72, True,  True),
    ("PAE",                  "pae",               48, True,  True),
    ("Patients",             "pts",               62, True,  True),
    ("Total nearby pts",     "total_near_pts",   100, True,  False),
    ("COSMIC",               "cosmic",            65, True,  False),
    ("Nearby mut count",     "near_mut_count",   100, True,  False),
    ("Confirmed disrupting", "confirmed_disrupt",110, False, False),
    ("PTM diseases",         "diseases",         140, False, False),
    ("14-3-3",               "pred14",            58, False, False),
    ("14-3-3 consensus",     "pred14_consensus",  90, False, False),
    ("14-3-3 confirmed",     "conf14",            90, False, False),
    ("Kinases",              "kinases",          160, False, False),
    ("PTM AIUPred gen.",     "ptm_aiupred_gen",   90, True,  False),
    ("PTM AIUPred bind.",    "ptm_aiupred_bind",  90, True,  False),
    ("PTM disordered?",      "ptm_disord",        90, False, False),
    ("PTM binding?",         "ptm_bind",          80, False, False),
    ("Mut AIUPred gen.",     "mut_aiupred_gen",   90, True,  False),
    ("Mut AIUPred bind.",    "mut_aiupred_bind",  90, True,  False),
    ("Gene",                 "gene",              58, False, False),
    ("UniProt",              "uniprot",           70, False, False),
    ("PTM Site",             "ptm_position",      65, False, False),
    ("PTM Type",             "ptm_type_l",       110, False, False),
    ("PTM pLDDT",            "ptm_plddt",         72, True,  False),
]

# Column-picker hover help, keyed by col_id. Two separate dicts (not one
# shared by col_id) because a few ids mean different things in each table —
# e.g. "pts" is the PTM table's total patient count across ALL nearby
# mutations, but the per-mutation patient count in the mutation table.
_PTM_COL_HELP: dict[str, str] = {
    "uniprot": "UniProt accession for this protein.",
    "gene": "Gene symbol for this protein.",
    "site": "The modified residue and its position (e.g. S557), from PTMD.",
    "type": "The type of post-translational modification at this site "
            "(e.g. Phosphorylation, Ubiquitination).",
    "near": "Number of nearby mutations (within the distance cutoff) that are "
            "also within 5 residues of this site in the linear sequence - "
            "likely to directly disrupt the modified residue itself.",
    "far": "Number of nearby mutations (within the distance cutoff) that are "
           "more than 5 residues away in the linear sequence - close in 3D "
           "space but not sequence-adjacent, suggesting fold-mediated or "
           "allosteric proximity rather than direct disruption.",
    "near_pts": "Total COSMIC patient count summed across the ≤ 5 pos "
                "(sequence-adjacent) nearby mutations.",
    "far_pts": "Total COSMIC patient count summed across the > 5 pos "
               "(sequence-distant) nearby mutations.",
    "near_unique": "Number of distinct mutated positions among the ≤ 5 pos "
                   "group (vs. ≤ 5 pos itself, which counts every mutation, "
                   "including multiple substitutions at the same position).",
    "far_unique": "Number of distinct mutated positions among the > 5 pos group.",
    "total": "Total number of distinct nearby mutation positions "
             "(≤ 5 pos + > 5 pos, unique positions only).",
    "pts": "Total COSMIC patient count across every nearby mutation for this "
           "PTM site (≤ 5 pos + > 5 pos combined).",
    "cosmic": "Total number of COSMIC patients with any missense mutation in "
              "this gene, regardless of distance to this PTM site - for "
              "context on how mutated the gene is overall.",
    "atptm": "Whether any nearby mutation occurs exactly at the modified "
             "residue itself, not just nearby - the most direct possible "
             "disruption.",
    "confirmed_disrupt": "Nearby mutations experimentally confirmed, in "
                          "PTMD's literature, to disrupt this specific PTM site.",
    "diseases": "Cancer-related diseases associated with this PTM site in "
                "PTMD's literature-curated data.",
    "pred14": "Predicted 14-3-3 binding at this site (14-3-3 proteins often "
              "bind phosphorylated motifs).",
    "pred14_consensus": "Agreement across multiple 14-3-3 binding predictors "
                         "for this site.",
    "conf14": "Whether this site is a literature-confirmed 14-3-3 binding "
              "site (from the bundled confirmed-interactors reference).",
    "conf14_pmid": "PubMed ID citing the literature confirmation for this "
                   "14-3-3 site.",
    "kinases": "Kinases predicted to phosphorylate this site.",
    "aiupred_gen": "AIUPred intrinsic disorder score (0-1) at this PTM "
                   "residue; above 0.5 is classified \"Disordered.\"",
    "aiupred_bind": "AIUPred binding-region disorder score (0-1) at this PTM "
                    "residue - disorder specifically linked to protein-binding "
                    "regions; above 0.5 is classified \"Binding.\"",
    "disord": "Yes/no: is this PTM residue predicted to be intrinsically "
              "disordered (AIUPred general score > 0.5)?",
    "bind": "Yes/no: is this PTM residue predicted to be a disordered "
            "binding region (AIUPred binding score > 0.5)?",
    "maxlin": "The largest linear (sequence) distance among the > 5 pos "
              "nearby mutations - how far the most sequence-distant-but-"
              "3D-close mutation actually is.",
    "near_muts_raw": "The individual ≤ 5 pos mutations, each with its 3D "
                     "distance (and PAE, if available) to this PTM site.",
    "far_muts_raw": "The individual > 5 pos mutations, each with its 3D "
                    "distance (and PAE, if available) to this PTM site.",
    "lin_dist_raw": "Linear (sequence) distance from this PTM site to each "
                    "individual > 5 pos mutation.",
}

_MUT_COL_HELP: dict[str, str] = {
    "mut": "The specific mutation (e.g. R175H) shown in this row.",
    "seqd": "Linear (sequence) distance, in residues, between this mutation "
            "and the PTM site.",
    "dist": "3D spatial distance, in Ångströms, between this mutation and "
            "the PTM site in the AlphaFold structure - the core proximity "
            "measurement.",
    "isbnd": "Yes/no: is this mutation's residue predicted to be a "
             "disordered binding region (AIUPred binding score > 0.5)?",
    "isdis": "Yes/no: is this mutation's residue predicted to be "
             "intrinsically disordered (AIUPred general score > 0.5)?",
    "ppc": "PolyPhen-2's classification of this mutation's predicted effect "
           "on protein function: benign, possibly damaging, or probably "
           "damaging.",
    "pps": "PolyPhen-2's raw score (0-1) for this mutation; higher means "
           "more likely to be damaging.",
    "mpld": "AlphaFold's per-residue confidence (pLDDT, 0-100) at this "
            "mutation's position.",
    "pae": "Predicted Aligned Error (Å), between the PTM site and this "
           "mutation - AlphaFold's confidence in their relative 3D position, "
           "independent of each residue's individual confidence. High PAE "
           "means the distance shown may not be reliable even if both "
           "residues have good pLDDT.",
    "pts": "Number of distinct COSMIC patient samples carrying this "
           "specific mutation.",
    "total_near_pts": "Total COSMIC patient count summed across every "
                       "mutation nearby this PTM site, not just this row - "
                       "for context.",
    "cosmic": "Total number of COSMIC patients with any missense mutation in "
              "this gene, regardless of distance to the PTM site - for "
              "context on how mutated the gene is overall.",
    "near_mut_count": "Total number of distinct mutations found nearby this "
                      "PTM site, not just this row - for context.",
    "confirmed_disrupt": "Nearby mutations experimentally confirmed, in "
                          "PTMD's literature, to disrupt this PTM site.",
    "diseases": "Cancer-related diseases associated with this PTM site in "
                "PTMD's literature-curated data.",
    "pred14": "Predicted 14-3-3 binding at this PTM site (14-3-3 proteins "
              "often bind phosphorylated motifs).",
    "pred14_consensus": "Agreement across multiple 14-3-3 binding predictors "
                         "for this PTM site.",
    "conf14": "Whether this PTM site is a literature-confirmed 14-3-3 "
              "binding site (from the bundled confirmed-interactors reference).",
    "kinases": "Kinases predicted to phosphorylate this PTM site.",
    "ptm_aiupred_gen": "AIUPred intrinsic disorder score (0-1) at the PTM "
                       "residue.",
    "ptm_aiupred_bind": "AIUPred binding-region disorder score (0-1) at the "
                        "PTM residue.",
    "ptm_disord": "Yes/no: is the PTM residue predicted to be intrinsically "
                  "disordered (AIUPred general score > 0.5)?",
    "ptm_bind": "Yes/no: is the PTM residue predicted to be a disordered "
                "binding region (AIUPred binding score > 0.5)?",
    "mut_aiupred_gen": "AIUPred intrinsic disorder score (0-1) at this "
                       "mutation's residue.",
    "mut_aiupred_bind": "AIUPred binding-region disorder score (0-1) at "
                        "this mutation's residue.",
    "gene": "Gene symbol for this protein.",
    "uniprot": "UniProt accession for this protein.",
    "ptm_position": "The modified residue and its position (e.g. S557) that "
                    "this mutation is being compared against.",
    "ptm_type_l": "The type of post-translational modification at the PTM "
                  "site (e.g. Phosphorylation, Ubiquitination).",
    "ptm_plddt": "AlphaFold's per-residue confidence (pLDDT, 0-100) at the "
                 "PTM site's position.",
}

# df_long column names for every _MUT_TV_COLS entry that's a direct pass-through
# (i.e. everything except "#col", the synthetic row index).
_MUT_LONG_SRC_MAP = {
    "mut": "mutation",
    "seqd": "sequence_distance",
    "dist": "distance_angstrom",
    "isbnd": "mut_is_binding",
    "isdis": "mut_is_disordered",
    "ppc": "polyphen_class",
    "pps": "polyphen_score",
    "mpld": "mutation_plddt",
    "pae": "pair_pae",
    "pts": "patient_count",
    "total_near_pts": "total_nearby_patient_count",
    "cosmic": "total_cosmic_missense_patients",
    "near_mut_count": "nearby_mutation_count",
    "confirmed_disrupt": "confirmed_disrupting_mutation",
    "diseases": "ptm_diseases",
    "pred14": "1433_predicted",
    "pred14_consensus": "1433_predicted_consensus",
    "conf14": "1433_confirmed",
    "kinases": "kinase_predictions",
    "ptm_aiupred_gen": "ptm_aiupred_general",
    "ptm_aiupred_bind": "ptm_aiupred_binding",
    "ptm_disord": "ptm_is_disordered",
    "ptm_bind": "ptm_is_binding",
    "mut_aiupred_gen": "mut_aiupred_general",
    "mut_aiupred_bind": "mut_aiupred_binding",
    "gene": "gene",
    "uniprot": "uniprot_id",
    "ptm_position": "ptm_position",
    "ptm_type_l": "ptm_type",
    "ptm_plddt": "ptm_plddt",
}


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


# Results-tab column visibility preferences (persisted across app launches).
COLUMN_PREFS_FILE = OUTPUT_DIR / "logs" / "results_column_prefs.json"


def _load_column_prefs(which: str) -> list[str] | None:
    try:
        data = json.loads(COLUMN_PREFS_FILE.read_text())
        cols = data.get(which)
        return cols if isinstance(cols, list) and cols else None
    except Exception:
        return None


def _save_column_prefs(which: str, col_ids: list[str]) -> None:
    try:
        data: dict = {}
        if COLUMN_PREFS_FILE.exists():
            data = json.loads(COLUMN_PREFS_FILE.read_text())
        data[which] = col_ids
        COLUMN_PREFS_FILE.parent.mkdir(parents=True, exist_ok=True)
        COLUMN_PREFS_FILE.write_text(json.dumps(data, indent=2))
    except Exception:
        pass
