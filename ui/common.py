"""Shared constants and small helper functions used across the app's UI mixins.

Mirrors the role scripts/pipeline_utils.py plays for the backend pipeline
scripts: one shared-utilities module, imported by every ui/*.py file.
"""
from __future__ import annotations

import csv
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
    get_protein_length,
)

# 14-3-3 confirmed-interactors file isn't listed here: it's bundled with the app
# (data/input/1433_interactors/), not something the user provides, so it's kept
# out of the Pipeline tab's input-file browse/status UI even though
# scripts/4_annotate.py still reads it via INTERACTORS_1433_INPUT_DIR.
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

    Plain tk.Toplevel/tk.Label rather than CTk widgets, since it's a short-lived
    unmanaged popup with no benefit from CTk's theming -- just fixed colors
    matching the app's dark theme.
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


def add_resize_grip(widget, min_height: int = 60, max_height: int = 800) -> ctk.CTkFrame:
    """A thin draggable handle that resizes *widget* vertically, the same way
    an OS window's edge can be dragged.

    *widget* must support ``.configure(height=...)`` / ``.cget("height")``
    (CTkTextbox does). Pack or grid the returned frame immediately after
    *widget*, in the same parent, so it reads as an edge of *widget*.
    """
    grip = ctk.CTkFrame(
        widget.master, height=7, fg_color="gray25", cursor="sb_v_double_arrow",
    )
    grip.grid_propagate(False)
    grip.pack_propagate(False)

    drag_state: dict[str, int] = {}

    def _start_drag(event) -> None:
        drag_state["y"] = event.y_root
        drag_state["height"] = int(widget.cget("height"))

    def _do_drag(event) -> None:
        delta = event.y_root - drag_state["y"]
        new_height = max(min_height, min(max_height, drag_state["height"] + delta))
        widget.configure(height=new_height)

    grip.bind("<Button-1>", _start_drag)
    grip.bind("<B1-Motion>", _do_drag)
    return grip


def isolate_textbox_scroll(textbox: ctk.CTkTextbox) -> None:
    """Keep mouse-wheel scrolling over *textbox* from also scrolling the
    outer page it sits on.

    CTkScrollableFrame binds <MouseWheel> at bind_all and scrolls itself
    whenever the event's target has its canvas in its ancestor chain -- true
    for any textbox on a scrollable tab, so a normal wheel-scroll runs both
    the textbox's own scroll and the page's. Binding directly on the real
    tkinter.Text and returning "break" stops the event before bind_all sees it.
    """
    real_text = textbox._textbox

    def _on_wheel(event):
        if getattr(event, "num", None) == 4:
            real_text.yview_scroll(-1, "units")
        elif getattr(event, "num", None) == 5:
            real_text.yview_scroll(1, "units")
        elif sys.platform == "darwin":
            # macOS reports small per-tick deltas directly (no /120 scaling,
            # unlike Windows, where a notch is always a multiple of 120).
            real_text.yview_scroll(int(-1 * event.delta), "units")
        else:
            real_text.yview_scroll(int(-1 * (event.delta / 120)), "units")
        return "break"

    real_text.bind("<MouseWheel>", _on_wheel)
    real_text.bind("<Button-4>", _on_wheel)
    real_text.bind("<Button-5>", _on_wheel)


# Pipeline tab mode-selector hover help, keyed by mode value.
_MODE_HELP: dict[str, str] = {
    "ptm-proximity": "Find cancer mutations that cluster near, or directly "
                      "disrupt, known PTM (post-translational modification) "
                      "sites in 3D protein structure. Merges PTMD's PTM-site "
                      "data with recurrent COSMIC hotspot mutations, then "
                      "annotates results with PolyPhen-2, kinase, 14-3-3, "
                      "and AIUPred predictions.",
    "mutation-clustering": "Find recurrent COSMIC hotspot mutations that "
                           "cluster together in 3D space, independent of PTM "
                           "sites - reveals spatial mutation hotspots a "
                           "linear sequence view wouldn't show.",
    "single-protein": "Run proximity analysis on one CIF structure file you "
                      "provide directly, without running the full pipeline "
                      "- useful for a quick, one-off look at a single protein.",
    "ca-coordinates": "Export alpha-carbon coordinates for every residue of "
                      "one protein (by UniProt ID or gene), along with its "
                      "COSMIC missense mutation positions, for use in "
                      "external visualization tools.",
}

# Analysis Tools tab hover help, keyed by field name.
_RADIUS_SWEEP_HELP: dict[str, str] = {
    "genes": "Add one or more proteins, by gene symbol or UniProt accession, "
             "to test. Each one is validated upfront - checked for hotspot "
             "mutation data and a downloaded AlphaFold structure - before "
             "being added.",
    "radius_range": "The range of distance cutoffs, in Ångströms, to test - "
                    "start, stop, and step size. The sweep re-runs the "
                    "proximity search at each radius in this range, so you "
                    "can see how results change as the cutoff changes.",
    "unfiltered": "Also run the sweep against every COSMIC missense mutation "
                  "for these genes, not just the recurrent hotspot ones - "
                  "lets you compare hotspot-filtered results against the "
                  "full, unfiltered mutation set at each radius.",
}

_CIF_VARIANCE_HELP: dict[str, str] = {
    "input_dir": "Folder of multiple CIF files for the SAME protein to "
                 "compare - e.g. different AlphaFold seeds or model "
                 "versions. Structures are aligned and compared to compute "
                 "per-residue positional variance.",
    "top_n": "Number of most structurally-variable residues to report.",
    "report_range": "Restrict the reported output to this residue range "
                    "(e.g. 0-630). Leave blank to report all residues.",
    "align_range": "Use only these residues for structural alignment, "
                   "rather than the whole protein. Defaults to the Report "
                   "range if left blank. Useful for excluding disordered "
                   "regions from alignment while still reporting their "
                   "variance.",
    "uniprot_override": "UniProt accession to use for PTM/mutation "
                        "cross-referencing. Auto-detected from the CIF file "
                        "if left blank.",
    "gene_override": "Gene symbol used to look up the UniProt ID from the "
                     "pipeline's intermediate data, if the UniProt override "
                     "above is left blank.",
}


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
    ("Step 4", "InterPro domains",       _CACHE_DIR / "interpro_domains.tsv",        False),
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

# Domain-map diagram (Visualization tab): color and lane per InterPro entry
# type. Lane groups by specificity (0 = broadest, rendered lowest) so a domain
# nested inside a broader family/superfamily call stays visually distinct
# instead of overdrawing the same strip. Unlisted types fall back to lane 1.
_DOMAIN_TYPE_COLORS: dict[str, str] = {
    "homologous_superfamily": _GRAY,
    "family": _GREEN,
    "domain": _BLUE,
    "repeat": "#e67e22",
    "conserved_site": _YELLOW,
    "active_site": _RED,
    "binding_site": "#9b59b6",
    "ptm": _PTM_MARKER_COLOR,
}
_DOMAIN_TYPE_LANES: dict[str, int] = {
    "homologous_superfamily": 0,
    "family": 0,
    "domain": 1,
    "repeat": 1,
    "conserved_site": 2,
    "active_site": 2,
    "binding_site": 2,
    "ptm": 2,
}
_DOMAIN_TYPE_FALLBACK_COLOR = "#5a5a5a"
_DOMAIN_TYPE_FALLBACK_LANE = 1


def _load_interpro_entries(uid: str) -> list[dict]:
    """Read cached InterPro domain/family/site entries for one protein.

    Returns [] if there's no cached row for this UniProt ID -- callers should
    treat that as "no domain data available", not an error.
    """
    cache_file = _CACHE_DIR / "interpro_domains.tsv"
    if not cache_file.exists():
        return []
    with cache_file.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            if row.get("uniprot_id") == uid:
                try:
                    return json.loads(row.get("entries_json", "[]"))
                except Exception:
                    return []
    return []

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
    ("≤5 pos patients",      "near_pts",          100, True,  False),
    (">5 pos patients",      "far_pts",           100, True,  False),
    ("≤5 unique pos",        "near_unique",       90, True,  False),
    (">5 unique pos",        "far_unique",        90, True,  False),
    ("Unique pos",           "total",             75, True,  True),
    ("Patients",             "pts",               65, True,  True),
    ("Total gene pts",       "cosmic",            95, True,  True),
    ("At PTM",               "atptm",             52, False, True),
    ("Confirmed disrupted", "confirmed_disrupt",150, False, False),
    ("PTM diseases",         "diseases",         140, False, False),
    ("14-3-3",               "pred14",            58, False, True),
    ("14-3-3 consensus score", "pred14_consensus",  150, False, False),
    ("14-3-3 confirmed",     "conf14",            120, False, False),
    ("14-3-3 PMID",          "conf14_pmid",       90, False, False),
    ("Predicted kinases",    "kinases",          160, False, False),
    ("PTM AIUPred gen.",     "aiupred_gen",       120, True,  False),
    ("PTM AIUPred bind.",    "aiupred_bind",      120, True,  False),
    ("Disordered?",          "disord",            78, False, True),
    ("Binding?",             "bind",              58, False, True),
    ("Max lin. dist.",       "maxlin",            90, True,  True),
    ("PTM pLDDT",            "ptm_plddt",         90, True,  False),
    ("Linear distances",     "lin_dist_raw",     150, False, False),
    ("PTM domain",           "ptm_domain",       180, False, False),
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
    ("Mut pLDDT",            "mpld",              75, True,  True),
    ("PAE",                  "pae",               48, True,  True),
    ("Patients",             "pts",               62, True,  True),
    ("Confirmed disrupting", "confirmed_disrupt",150, False, False),
    ("Mut AIUPred gen.",     "mut_aiupred_gen",   120, True,  False),
    ("Mut AIUPred bind.",    "mut_aiupred_bind",  120, True,  False),
    ("Mutation domain",      "mut_domain",       180, False, False),
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
    "pred14_consensus": "14-3-3-Pred's combined \"Consensus\" score for this "
                         "site, combining its individual prediction methods "
                         "(ANN, PSSM) into one number, roughly -1 to +1.5 in "
                         "practice. Positive is what the separate \"14-3-3\" "
                         "Yes/No column is based on - the more positive, the "
                         "stronger the predicted binding signal; negative "
                         "means predicted non-binding.",
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
    "lin_dist_raw": "Linear (sequence) distance from this PTM site to each "
                    "individual > 5 pos mutation.",
    "ptm_domain": "InterPro functional domain(s) (name, type, and residue "
                  "range) containing this PTM site's position, if any. A "
                  "residue can fall inside more than one entry, e.g. a "
                  "specific domain nested inside a broader superfamily call "
                  "- all are shown, semicolon-separated.",
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
    "confirmed_disrupt": "Yes/no: is this specific mutation experimentally "
                          "confirmed, in PTMD's literature, to disrupt this "
                          "PTM site?",
    "mut_aiupred_gen": "AIUPred intrinsic disorder score (0-1) at this "
                       "mutation's residue.",
    "mut_aiupred_bind": "AIUPred binding-region disorder score (0-1) at "
                        "this mutation's residue.",
    "ptm_plddt": "AlphaFold's per-residue confidence (pLDDT, 0-100) at the "
                 "PTM site's position - not shown in the PTM Sites table, "
                 "so it's kept here.",
    "mut_domain": "InterPro functional domain(s) (name, type, and residue "
                  "range) containing this mutation's position, if any. A "
                  "residue can fall inside more than one entry, e.g. a "
                  "specific domain nested inside a broader superfamily call "
                  "- all are shown, semicolon-separated.",
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
    "confirmed_disrupt": "confirmed_disrupting_mutation",
    "mut_aiupred_gen": "mut_aiupred_general",
    "mut_aiupred_bind": "mut_aiupred_binding",
    "ptm_plddt": "ptm_plddt",
    "mut_domain": "mutation_domain",
}

# Mutation-Clustering mode's equivalent of _PTM_TV_COLS/_MUT_TV_COLS. No
# 14-3-3/kinase columns -- those are PTM-site-specific and this mode has no
# concept of a PTM site (the anchor is itself a mutation).
_ANCHOR_TV_COLS = [
    ("#",                "#col",        32,  True,  True),
    ("UniProt",          "uniprot",     70,  False, True),
    ("Gene",             "gene",        58,  False, True),
    ("Anchor mutation",  "anchor",     100,  False, True),
    ("Binding?",         "isbnd",       58,  False, True),
    ("Disordered?",      "isdis",       78,  False, True),
    ("PP Class",         "ppc",        115,  False, True),
    ("PP Score",         "pps",         62,  True,  False),
    ("Anchor pLDDT",     "anchor_plddt", 90, True,  True),
    ("Nearby count",     "near_count",  95,  True,  True),
    ("Unique positions", "uniq_pos",   100,  True,  True),
    ("Nearby patients",  "near_pts",   100,  True,  True),
    ("Anchor AIUPred gen.",  "anchor_aiupred_gen",  120, True,  False),
    ("Anchor AIUPred bind.", "anchor_aiupred_bind", 120, True,  False),
    ("Anchor domain",    "anchor_domain", 180, False, False),
]

_NEARBY_TV_COLS = [
    ("#",                    "#col", 32, True,  True),
    ("Mutation",             "mut",  80, False, True),
    ("Seq dist",             "seqd", 62, True,  True),
    ("Dist (Å)",             "dist", 62, True,  True),
    ("Binding?",             "isbnd",             58, False, True),
    ("Disordered?",          "isdis",             78, False, True),
    ("PP Class",             "ppc",              115, False, True),
    ("PP Score",             "pps",               62, True,  False),
    ("PAE",                  "pae",  48, True,  True),
    ("Mut pLDDT",            "mpld", 75, True,  True),
    ("Patients",             "pts",  62, True,  True),
    ("Mut AIUPred gen.",     "mut_aiupred_gen",   120, True,  False),
    ("Mut AIUPred bind.",    "mut_aiupred_bind",  120, True,  False),
    ("Mutation domain",      "mut_domain",       180, False, False),
]

_ANCHOR_COL_HELP: dict[str, str] = {
    "uniprot": "UniProt accession for this protein.",
    "gene": "Gene symbol for this protein.",
    "anchor": "The mutation this cluster is centered on. Clustering is "
              "symmetric - every mutation with at least one 3D neighbor gets "
              "its own anchor row, so a pair (A, B) appears twice: once "
              "anchored on A, once anchored on B.",
    "anchor_plddt": "AlphaFold's per-residue confidence (pLDDT, 0-100) at "
                    "the anchor mutation's position.",
    "near_count": "Number of other recurrent mutations within the distance "
                  "cutoff of this anchor mutation.",
    "uniq_pos": "Number of distinct mutated positions among the nearby "
                "mutations (vs. Nearby count, which counts every mutation, "
                "including multiple substitutions at the same position).",
    "near_pts": "Total COSMIC patient count summed across every nearby "
                "mutation for this anchor.",
    "isbnd": "Yes/no: is the anchor mutation's residue predicted to be a "
             "disordered binding region (AIUPred binding score > 0.5)?",
    "isdis": "Yes/no: is the anchor mutation's residue predicted to be "
             "intrinsically disordered (AIUPred general score > 0.5)?",
    "ppc": "PolyPhen-2's classification of the anchor mutation's predicted "
           "effect on protein function: benign, possibly damaging, or "
           "probably damaging.",
    "pps": "PolyPhen-2's raw score (0-1) for the anchor mutation; higher "
           "means more likely to be damaging.",
    "anchor_aiupred_gen": "AIUPred intrinsic disorder score (0-1) at the "
                          "anchor mutation's residue.",
    "anchor_aiupred_bind": "AIUPred binding-region disorder score (0-1) at "
                           "the anchor mutation's residue.",
    "anchor_domain": "InterPro functional domain(s) (name, type, and residue "
                     "range) containing the anchor mutation's position, if "
                     "any. A residue can fall inside more than one entry - "
                     "all are shown, semicolon-separated.",
}

_NEARBY_COL_HELP: dict[str, str] = {
    "mut": "The specific nearby mutation (e.g. R175H) shown in this row.",
    "seqd": "Linear (sequence) distance, in residues, between this mutation "
            "and the anchor mutation.",
    "dist": "3D spatial distance, in Ångströms, between this mutation and "
            "the anchor mutation in the AlphaFold structure.",
    "isbnd": "Yes/no: is this mutation's residue predicted to be a "
             "disordered binding region (AIUPred binding score > 0.5)?",
    "isdis": "Yes/no: is this mutation's residue predicted to be "
             "intrinsically disordered (AIUPred general score > 0.5)?",
    "ppc": "PolyPhen-2's classification of this mutation's predicted effect "
           "on protein function: benign, possibly damaging, or probably "
           "damaging.",
    "pps": "PolyPhen-2's raw score (0-1) for this mutation; higher means "
           "more likely to be damaging.",
    "pae": "Predicted Aligned Error (Å) between the anchor and this "
           "mutation - AlphaFold's confidence in their relative 3D position.",
    "mpld": "AlphaFold's per-residue confidence (pLDDT, 0-100) at this "
            "mutation's position.",
    "pts": "Number of distinct COSMIC patient samples carrying this "
           "specific mutation.",
    "mut_aiupred_gen": "AIUPred intrinsic disorder score (0-1) at this "
                       "mutation's residue.",
    "mut_aiupred_bind": "AIUPred binding-region disorder score (0-1) at "
                        "this mutation's residue.",
    "mut_domain": "InterPro functional domain(s) (name, type, and residue "
                  "range) containing this mutation's position, if any. A "
                  "residue can fall inside more than one entry - all are "
                  "shown, semicolon-separated.",
}

# mutation_cluster_long.tsv column names for every _NEARBY_TV_COLS entry
# that's a direct pass-through (i.e. everything except "#col").
_CLUSTER_LONG_SRC_MAP = {
    "mut": "mutation",
    "seqd": "sequence_distance",
    "dist": "distance_angstrom",
    "isbnd": "mut_is_binding",
    "isdis": "mut_is_disordered",
    "ppc": "polyphen_class",
    "pps": "polyphen_score",
    "pae": "pair_pae",
    "mpld": "mutation_plddt",
    "pts": "patient_count",
    "mut_aiupred_gen": "mut_aiupred_general",
    "mut_aiupred_bind": "mut_aiupred_binding",
    "mut_domain": "mutation_domain",
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
