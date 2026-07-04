"""Desktop GUI for the Mutation Cluster Proximity pipeline.

Run with:  uv run app.py
"""
from __future__ import annotations

import json
import platform
import queue
import re
import shutil
import subprocess
import sys
import threading
import time
from pathlib import Path
from tkinter import filedialog

import psutil

import customtkinter as ctk
from PIL import Image

PROJECT_ROOT = Path(__file__).resolve().parent
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

ctk.set_appearance_mode("dark")
ctk.set_default_color_theme("blue")

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
    ("#",          "#col",   32),
    ("Gene",       "gene",   58),
    ("PTM Site",   "site",   65),
    ("Type",       "type",  110),
    ("≤5 pos",     "near",   52),
    (">5 pos",     "far",    52),
    ("Unique pos", "total",  68),
    ("Patients",   "pts",    65),
    ("COSMIC",     "cosmic", 65),
    ("At PTM",     "atptm",  52),
    ("Disrupting", "disrupt",68),
    ("14-3-3",     "pred14", 58),
    ("Conf.",      "conf14", 52),
    ("Diseases",   "dis",   200),
]

_MUT_TV_COLS = [
    ("#",          "#col",   32),
    ("Mutation",   "mut",    80),
    ("Seq dist",   "seqd",   62),
    ("Dist (Å)",   "dist",   62),
    ("Binding",    "bind",   62),
    ("Binding?",   "isbnd",  58),
    ("Disorder",   "dsord",  62),
    ("Disordr?",   "isdis",  58),
    ("PP Class",   "ppc",   115),
    ("PP Score",   "pps",    62),
    ("Mut pLDDT",  "mpld",   72),
    ("PAE",        "pae",    48),
    ("Patients",   "pts",    62),
    ("Disrupting", "dis",    68),
]


def _fmt_time(seconds: float) -> str:
    s = int(seconds)
    if s < 60:
        return f"{s}s"
    return f"{s // 60}m {s % 60:02d}s"


RUNTIMES_FILE = OUTPUT_DIR / "logs" / "pipeline_runtimes.json"
_CIF_DIR = PROJECT_ROOT / "cif_models"
_CACHE_DIR = PROJECT_ROOT / "data" / "cache"


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


class App(ctk.CTk):
    def __init__(self):
        super().__init__()
        self.title("Cluster-Scout")
        self.geometry("1100x820")
        self.minsize(900, 560)
        ico = PROJECT_ROOT / "cluster_scout.ico"
        if ico.exists():
            self.iconbitmap(str(ico))

        self._queue: queue.Queue[tuple] = queue.Queue()
        self._radius_sweep_result = None
        self._cif_variance_result = None
        self._running = False
        self._stop_requested = False
        self._suspended = False
        self._current_proc: subprocess.Popen | None = None
        self._step_status_labels: list[ctk.CTkLabel] = []
        self._pipeline_start: float | None = None
        self._step_start: float | None = None
        self._steps_done = 0
        self._total_steps = 0
        self._step_times: list[float] = []
        self._historical_times: list[float] | None = None

        self._build_ui()
        self._refresh_file_status()
        self._poll_queue()

    # ── UI construction ──────────────────────────────────────────────────────

    def _build_ui(self):
        self.grid_columnconfigure(0, weight=1)
        self.grid_rowconfigure(0, weight=1)

        # Tab view
        self._tabview = ctk.CTkTabview(self, command=self._on_tab_change)
        self._tabview.grid(row=0, column=0, padx=12, pady=12, sticky="nsew")
        pipeline_tab = self._tabview.add("Pipeline")
        results_tab = self._tabview.add("Results")
        viz_tab = self._tabview.add("Visualization")
        radius_sweep_tab = self._tabview.add("Radius Sweep")
        cif_variance_tab = self._tabview.add("CIF Variance")
        help_tab = self._tabview.add("Help / Documentation")

        # ── Pipeline tab ──
        pipeline_tab.grid_columnconfigure(0, weight=1)
        pipeline_tab.grid_rowconfigure(0, weight=1)

        scroll = ctk.CTkScrollableFrame(pipeline_tab, fg_color="transparent")
        scroll.grid(row=0, column=0, sticky="nsew")
        scroll.grid_columnconfigure(0, weight=1)

        p = scroll  # all pipeline widgets go in the scrollable frame

        # Title
        _logo_path_dark = PROJECT_ROOT / "cluster_scout_logo_dark.png"
        if _logo_path_dark.exists():
            _pil_dark = Image.open(_logo_path_dark)
            _h = 200
            _w = int(_pil_dark.width * _h / _pil_dark.height)
            _logo_img = ctk.CTkImage(
                light_image=_pil_dark,
                dark_image=_pil_dark,
                size=(_w, _h),
            )
            ctk.CTkLabel(p, image=_logo_img, text="").grid(
                row=0, column=0, padx=24, pady=(12, 4), sticky="w"
            )
        else:
            ctk.CTkLabel(
                p,
                text="Cluster-Scout",
                font=ctk.CTkFont(size=22, weight="bold"),
            ).grid(row=0, column=0, padx=24, pady=(12, 4), sticky="w")

        # Data-file status bar with Browse buttons
        self._file_frame = ctk.CTkFrame(p)
        self._file_frame.grid(row=1, column=0, padx=24, pady=4, sticky="ew")
        self._file_frame.grid_columnconfigure(1, weight=1)
        self._file_indicators: dict[str, ctk.CTkLabel] = {}
        self._file_buttons: dict[str, ctk.CTkButton] = {}

        ctk.CTkLabel(
            self._file_frame,
            text="Input files:",
            font=ctk.CTkFont(weight="bold"),
        ).grid(row=0, column=0, columnspan=3, padx=12, pady=(8, 2), sticky="w")

        for i, (name, (folder, exts, desc)) in enumerate(_INPUT_FOLDERS.items(), 1):
            lbl = ctk.CTkLabel(self._file_frame, text=f"{name} …", anchor="w")
            lbl.grid(row=i, column=0, columnspan=2, padx=(12, 6), pady=3, sticky="ew")
            self._file_indicators[name] = lbl

            filetypes = [(desc, " ".join(f"*{e}" for e in exts))]
            btn = ctk.CTkButton(
                self._file_frame,
                text="Browse",
                width=70,
                height=26,
                font=ctk.CTkFont(size=12),
                command=lambda n=name, f=folder, ft=filetypes: self._browse_file(n, f, ft),
            )
            btn.grid(row=i, column=2, padx=12, pady=3, sticky="e")
            self._file_buttons[name] = btn

        # Mode selection
        mode_frame = ctk.CTkFrame(p)
        mode_frame.grid(row=2, column=0, padx=24, pady=4, sticky="ew")

        ctk.CTkLabel(
            mode_frame, text="Mode:", font=ctk.CTkFont(weight="bold")
        ).pack(side="left", padx=(12, 8), pady=10)

        self._mode = ctk.StringVar(value="ptm-proximity")
        for label, value in [
            ("PTM Proximity", "ptm-proximity"),
            ("Mutation Clustering", "mutation-clustering"),
            ("Single Protein", "single-protein"),
            ("Radius Sweep", "radius-sweep"),
            ("CIF Variance", "cif-variance"),
            ("CA Coordinates", "ca-coordinates"),
        ]:
            ctk.CTkRadioButton(
                mode_frame,
                text=label,
                variable=self._mode,
                value=value,
                command=self._rebuild_step_rows,
            ).pack(side="left", padx=8, pady=10)

        # Output folder selector
        out_frame = ctk.CTkFrame(p)
        out_frame.grid(row=3, column=0, padx=24, pady=4, sticky="ew")
        out_frame.grid_columnconfigure(1, weight=1)

        ctk.CTkLabel(
            out_frame, text="Output folder:", font=ctk.CTkFont(weight="bold"),
        ).grid(row=0, column=0, padx=(12, 6), pady=8, sticky="w")

        self._output_dir_var = ctk.StringVar(value=str(OUTPUT_DIR))
        self._output_dir_entry = ctk.CTkEntry(
            out_frame, textvariable=self._output_dir_var, state="readonly",
        )
        self._output_dir_entry.grid(row=0, column=1, padx=6, pady=8, sticky="ew")

        ctk.CTkButton(
            out_frame, text="Change", width=70, height=26,
            font=ctk.CTkFont(size=12),
            command=self._browse_output_dir,
        ).grid(row=0, column=2, padx=(6, 4), pady=8, sticky="e")

        ctk.CTkButton(
            out_frame, text="Reset", width=60, height=26,
            font=ctk.CTkFont(size=12), fg_color="gray30", hover_color="gray40",
            command=lambda: self._output_dir_var.set(str(OUTPUT_DIR)),
        ).grid(row=0, column=3, padx=(0, 12), pady=8, sticky="e")

        # Pipeline settings
        settings_frame = ctk.CTkFrame(p)
        settings_frame.grid(row=4, column=0, padx=24, pady=4, sticky="ew")

        ctk.CTkLabel(
            settings_frame, text="Settings:", font=ctk.CTkFont(weight="bold"),
        ).pack(side="left", padx=(12, 8), pady=8)

        ctk.CTkLabel(settings_frame, text="Cutoff (Å):").pack(side="left", padx=(8, 4), pady=8)
        self._cutoff_var = ctk.StringVar(value="10.0")
        ctk.CTkEntry(
            settings_frame, textvariable=self._cutoff_var, width=60,
        ).pack(side="left", padx=(0, 16), pady=8)

        ctk.CTkLabel(settings_frame, text="Min samples:").pack(side="left", padx=(8, 4), pady=8)
        self._min_samples_var = ctk.StringVar(value="3")
        ctk.CTkEntry(
            settings_frame, textvariable=self._min_samples_var, width=60,
        ).pack(side="left", padx=(0, 16), pady=8)

        ctk.CTkLabel(settings_frame, text="Min pLDDT:").pack(side="left", padx=(8, 4), pady=8)
        self._min_plddt_var = ctk.StringVar(value="")
        ctk.CTkEntry(
            settings_frame, textvariable=self._min_plddt_var, width=60,
            placeholder_text="off",
        ).pack(side="left", padx=(0, 16), pady=8)

        ctk.CTkLabel(settings_frame, text="Max PAE:").pack(side="left", padx=(8, 4), pady=8)
        self._max_pae_var = ctk.StringVar(value="")
        ctk.CTkEntry(
            settings_frame, textvariable=self._max_pae_var, width=60,
            placeholder_text="off",
        ).pack(side="left", padx=(0, 12), pady=8)

        # PolyPhen filter
        pp_frame = ctk.CTkFrame(p)
        pp_frame.grid(row=5, column=0, padx=24, pady=4, sticky="ew")

        ctk.CTkLabel(
            pp_frame, text="PolyPhen filter:", font=ctk.CTkFont(weight="bold"),
        ).pack(side="left", padx=(12, 4), pady=8)

        self._pp_benign_var = ctk.BooleanVar(value=True)
        ctk.CTkCheckBox(
            pp_frame, text="Benign",
            variable=self._pp_benign_var,
            checkbox_width=18, checkbox_height=18,
        ).pack(side="left", padx=8, pady=8)

        self._pp_possibly_var = ctk.BooleanVar(value=True)
        ctk.CTkCheckBox(
            pp_frame, text="Possibly damaging",
            variable=self._pp_possibly_var,
            checkbox_width=18, checkbox_height=18,
        ).pack(side="left", padx=8, pady=8)

        self._pp_probably_var = ctk.BooleanVar(value=True)
        ctk.CTkCheckBox(
            pp_frame, text="Probably damaging",
            variable=self._pp_probably_var,
            checkbox_width=18, checkbox_height=18,
        ).pack(side="left", padx=(8, 12), pady=8)

        # Steps panel
        self._steps_outer = ctk.CTkFrame(p)
        self._steps_outer.grid(row=7, column=0, padx=24, pady=4, sticky="ew")
        self._steps_outer.grid_columnconfigure(1, weight=1)
        self._rebuild_step_rows()

        # Buttons
        btn_frame = ctk.CTkFrame(p, fg_color="transparent")
        btn_frame.grid(row=8, column=0, padx=24, pady=8, sticky="ew")

        self._run_btn = ctk.CTkButton(
            btn_frame,
            text="▶  Run Pipeline",
            command=self._start_pipeline,
            width=160,
            height=44,
            font=ctk.CTkFont(size=15, weight="bold"),
        )
        self._run_btn.pack(side="left", padx=(0, 8))

        self._stop_btn = ctk.CTkButton(
            btn_frame,
            text="■  Stop",
            command=self._stop_pipeline,
            width=110,
            height=44,
            font=ctk.CTkFont(size=15, weight="bold"),
            state="disabled",
            fg_color="gray30",
            hover_color=_RED,
        )
        self._stop_btn.pack(side="left", padx=(0, 12))

        self._open_btn = ctk.CTkButton(
            btn_frame,
            text="Open Output Folder",
            command=self._open_output_folder,
            width=180,
            height=44,
        )
        self._open_btn.pack(side="left")

        ctk.CTkButton(
            btn_frame,
            text="Manage Cache",
            command=self._manage_cache,
            width=140,
            height=44,
            fg_color="gray30",
            hover_color="gray40",
        ).pack(side="left", padx=(8, 0))

        self._timer_label = ctk.CTkLabel(
            btn_frame,
            text="",
            text_color=_GRAY,
            font=ctk.CTkFont(size=13),
        )
        self._timer_label.pack(side="right", padx=12)

        import tkinter as tk
        self._activity_width = 40
        self._activity_height = 6
        self._activity_chunk = 12
        self._activity_canvas = tk.Canvas(
            btn_frame, width=self._activity_width, height=self._activity_height,
            highlightthickness=0, bg="#333333", bd=0,
        )
        self._activity_canvas.pack(side="right", padx=(0, 4))
        self._activity_canvas.pack_forget()
        self._activity_animating = False
        self._activity_pos = -self._activity_chunk

        # Log (collapsible)
        self._log_visible = False
        self._log_toggle = ctk.CTkButton(
            p,
            text="Show Details",
            width=120,
            height=28,
            font=ctk.CTkFont(size=12),
            fg_color="gray30",
            hover_color="gray40",
            command=self._toggle_log,
        )
        self._log_toggle.grid(row=9, column=0, padx=24, pady=(8, 0), sticky="w")

        self._log = ctk.CTkTextbox(
            p,
            height=140,
            font=ctk.CTkFont(family="Courier New", size=12),
            wrap="word",
            state="disabled",
        )
        self._toggle_log()  # visible by default

        # Hide the scrollbar when content fits; show it only when scrolling is needed
        def _update_scrollbar(*_):
            try:
                canvas = scroll._parent_canvas
                sr = canvas.cget("scrollregion")
                if not sr:
                    scroll._scrollbar.grid_remove()
                    return
                content_h = int(float(sr.split()[3]))
                if content_h > canvas.winfo_height():
                    scroll._scrollbar.grid()
                else:
                    scroll._scrollbar.grid_remove()
            except Exception:
                pass

        # add="+" is required here: CTkScrollableFrame already binds its own
        # <Configure> handler (to keep the canvas scrollregion in sync with
        # content size) — a plain .bind() would replace it instead of adding
        # to it, freezing the scrollregion at whatever it was when this ran.
        scroll.bind("<Configure>", _update_scrollbar, add="+")
        pipeline_tab.bind("<Configure>", _update_scrollbar, add="+")
        self.after(200, _update_scrollbar)

        # ── Results tab ──
        self._build_results_tab(results_tab)

        # ── Visualization tab ──
        self._build_viz_tab(viz_tab)

        # ── Radius Sweep tab ──
        self._build_radius_sweep_tab(radius_sweep_tab)

        # ── CIF Variance tab ──
        self._build_cif_variance_tab(cif_variance_tab)

        # ── Help / Documentation tab ──
        self._build_help_tab(help_tab)

    def _build_help_tab(self, tab):
        """Load docs/help.md, convert to HTML, and display in the Help tab."""
        import markdown

        tab.grid_columnconfigure(0, weight=1)
        tab.grid_rowconfigure(0, weight=1)

        help_md = PROJECT_ROOT / "docs" / "help.md"
        try:
            md_text = help_md.read_text(encoding="utf-8")
        except FileNotFoundError:
            md_text = "# Help\n\nDocumentation file not found at `docs/help.md`."

        html_body = markdown.markdown(md_text, extensions=["tables", "fenced_code"])
        css = """
        <style>
            body { font-family: Segoe UI, Arial, sans-serif; padding: 24px;
                   background: #2b2b2b; color: #dcdcdc; line-height: 1.6; }
            h1 { color: #3a86ff; border-bottom: 2px solid #3a86ff; padding-bottom: 6px; }
            h2 { color: #6cb4ee; margin-top: 28px; }
            h3 { color: #a0cfff; margin-top: 20px; }
            table { border-collapse: collapse; margin: 12px 0; width: 100%; }
            th, td { border: 1px solid #555; padding: 6px 10px; text-align: left; }
            th { background: #3a3a3a; color: #dcdcdc; }
            tr:nth-child(even) { background: #333; }
            code { background: #3a3a3a; padding: 2px 5px; border-radius: 3px; }
            hr { border: 1px solid #555; margin: 20px 0; }
            strong { color: #f0f0f0; }
            ul, ol { padding-left: 24px; }
            li { margin: 4px 0; }
        </style>
        """
        full_html = f"<html><head>{css}</head><body>{html_body}</body></html>"

        try:
            from tkinterweb import HtmlFrame
            frame = HtmlFrame(tab, messages_enabled=False)
            frame.load_html(full_html)
            frame.grid(row=0, column=0, sticky="nsew")
        except ImportError:
            fallback = ctk.CTkTextbox(tab, wrap="word", font=ctk.CTkFont(size=13))
            fallback.insert("1.0", md_text)
            fallback.configure(state="disabled")
            fallback.grid(row=0, column=0, sticky="nsew")

    def _rebuild_step_rows(self):
        for w in self._steps_outer.winfo_children():
            w.destroy()
        self._step_status_labels = []
        self._step_progress_bars: list[ctk.CTkProgressBar] = []
        self._steps_outer.grid_columnconfigure(1, weight=1)

        mode = self._mode.get()

        if mode == "single-protein":
            self._build_single_protein_panel()
            return
        if mode == "radius-sweep":
            self._build_radius_sweep_panel()
            return
        if mode == "cif-variance":
            self._build_cif_variance_panel()
            return
        if mode == "ca-coordinates":
            self._build_ca_coordinates_panel()
            return

        steps = PTM_PROXIMITY_STEPS if mode == "ptm-proximity" else MUTATION_CLUSTERING_STEPS

        ctk.CTkLabel(
            self._steps_outer,
            text="Steps",
            font=ctk.CTkFont(weight="bold"),
        ).grid(row=0, column=0, columnspan=4, padx=12, pady=(8, 2), sticky="w")

        for i, (panel_label, _log_label) in enumerate(steps, 1):
            ctk.CTkLabel(self._steps_outer, text=f"  {i}.", width=28).grid(
                row=i, column=0, padx=(12, 0), pady=5, sticky="w"
            )
            ctk.CTkLabel(self._steps_outer, text=panel_label, anchor="w").grid(
                row=i, column=1, padx=6, pady=5, sticky="ew"
            )

            bar = ctk.CTkProgressBar(self._steps_outer, width=120, height=14)
            bar.set(0)
            bar.grid(row=i, column=2, padx=6, pady=5, sticky="e")
            bar.grid_remove()
            self._step_progress_bars.append(bar)

            status = ctk.CTkLabel(
                self._steps_outer,
                text="●  Waiting",
                width=100,
                anchor="e",
                text_color=_GRAY,
            )
            status.grid(row=i, column=3, padx=12, pady=5, sticky="e")
            self._step_status_labels.append(status)

    def _build_single_protein_panel(self):
        """Build the input fields for single-protein analysis mode."""
        ctk.CTkLabel(
            self._steps_outer,
            text="Single Protein Analysis",
            font=ctk.CTkFont(weight="bold"),
        ).grid(row=0, column=0, columnspan=3, padx=12, pady=(8, 2), sticky="w")

        # CIF file picker
        ctk.CTkLabel(self._steps_outer, text="CIF file:", anchor="w").grid(
            row=1, column=0, padx=(12, 6), pady=6, sticky="w"
        )
        self._single_cif_var = ctk.StringVar(value="")
        self._single_cif_entry = ctk.CTkEntry(
            self._steps_outer, textvariable=self._single_cif_var, width=400,
        )
        self._single_cif_entry.grid(row=1, column=1, padx=6, pady=6, sticky="ew")
        ctk.CTkButton(
            self._steps_outer, text="Browse", width=70, height=26,
            font=ctk.CTkFont(size=12),
            command=self._browse_cif,
        ).grid(row=1, column=2, padx=12, pady=6, sticky="e")

        # UniProt ID
        ctk.CTkLabel(self._steps_outer, text="UniProt ID:", anchor="w").grid(
            row=2, column=0, padx=(12, 6), pady=6, sticky="w"
        )
        self._single_uniprot_var = ctk.StringVar(value="")
        ctk.CTkEntry(
            self._steps_outer, textvariable=self._single_uniprot_var, width=200,
            placeholder_text="Edit if CIF is not in a UniProt-named folder",
        ).grid(row=2, column=1, padx=6, pady=6, sticky="w")

        # Status label (reuse the step status pattern so the pipeline runner can update it)
        status = ctk.CTkLabel(
            self._steps_outer, text="●  Ready", width=100,
            anchor="e", text_color=_GRAY,
        )
        status.grid(row=3, column=1, columnspan=2, padx=12, pady=6, sticky="e")
        self._step_status_labels.append(status)

        bar = ctk.CTkProgressBar(self._steps_outer, width=120, height=14)
        bar.set(0)
        bar.grid(row=3, column=0, padx=12, pady=6, sticky="w")
        bar.grid_remove()
        self._step_progress_bars.append(bar)

    def _browse_cif(self):
        """Open a file dialog for selecting a CIF file and auto-fill the UniProt ID."""
        path = filedialog.askopenfilename(
            title="Select CIF structure file",
            filetypes=[("CIF files", "*.cif"), ("All files", "*.*")],
        )
        if not path:
            return
        self._single_cif_var.set(path)
        uid = extract_uniprot_from_cif(Path(path))
        if uid:
            self._single_uniprot_var.set(uid)
        else:
            parent_name = Path(path).parent.name
            self._single_uniprot_var.set(parent_name)

    def _build_radius_sweep_panel(self):
        """Build the input fields for radius-sweep analysis mode."""
        ctk.CTkLabel(
            self._steps_outer,
            text="Radius Sweep Analysis",
            font=ctk.CTkFont(weight="bold"),
        ).grid(row=0, column=0, columnspan=4, padx=12, pady=(8, 2), sticky="w")

        # Genes
        ctk.CTkLabel(self._steps_outer, text="Genes:", anchor="w").grid(
            row=1, column=0, padx=(12, 6), pady=6, sticky="w"
        )
        if not hasattr(self, "_radius_genes_var"):
            from radius_sweep import DEFAULT_GENES
            self._radius_genes_var = ctk.StringVar(value=" ".join(DEFAULT_GENES))
        ctk.CTkEntry(
            self._steps_outer, textvariable=self._radius_genes_var, width=400,
            placeholder_text="Space-separated gene symbols",
        ).grid(row=1, column=1, columnspan=3, padx=6, pady=6, sticky="ew")

        # Radius range
        ctk.CTkLabel(self._steps_outer, text="Radius range (Å):", anchor="w").grid(
            row=2, column=0, padx=(12, 6), pady=6, sticky="w"
        )
        range_frame = ctk.CTkFrame(self._steps_outer, fg_color="transparent")
        range_frame.grid(row=2, column=1, columnspan=3, padx=6, pady=6, sticky="w")
        if not hasattr(self, "_radius_start_var"):
            self._radius_start_var = ctk.StringVar(value="4")
            self._radius_stop_var = ctk.StringVar(value="20")
            self._radius_step_var = ctk.StringVar(value="1")
        for label, var in [
            ("start", self._radius_start_var),
            ("stop", self._radius_stop_var),
            ("step", self._radius_step_var),
        ]:
            ctk.CTkLabel(range_frame, text=label).pack(side="left", padx=(0, 4))
            ctk.CTkEntry(range_frame, textvariable=var, width=50).pack(side="left", padx=(0, 12))

        # Unfiltered comparison
        if not hasattr(self, "_radius_unfiltered_var"):
            self._radius_unfiltered_var = ctk.BooleanVar(value=False)
        ctk.CTkCheckBox(
            self._steps_outer, text="Include unfiltered COSMIC comparison",
            variable=self._radius_unfiltered_var,
            checkbox_width=18, checkbox_height=18,
        ).grid(row=3, column=0, columnspan=4, padx=12, pady=6, sticky="w")

        # Status label + progress bar (reuse the step status pattern)
        status = ctk.CTkLabel(
            self._steps_outer, text="●  Ready", width=100,
            anchor="e", text_color=_GRAY,
        )
        status.grid(row=4, column=1, columnspan=3, padx=12, pady=6, sticky="e")
        self._step_status_labels.append(status)

        bar = ctk.CTkProgressBar(self._steps_outer, width=120, height=14)
        bar.set(0)
        bar.grid(row=4, column=0, padx=12, pady=6, sticky="w")
        bar.grid_remove()
        self._step_progress_bars.append(bar)

    def _build_cif_variance_panel(self):
        """Build the input fields for CIF-variance analysis mode."""
        ctk.CTkLabel(
            self._steps_outer,
            text="CIF Variance Analysis",
            font=ctk.CTkFont(weight="bold"),
        ).grid(row=0, column=0, columnspan=3, padx=12, pady=(8, 2), sticky="w")

        # Input folder
        ctk.CTkLabel(self._steps_outer, text="Input folder:", anchor="w").grid(
            row=1, column=0, padx=(12, 6), pady=6, sticky="w"
        )
        if not hasattr(self, "_variance_input_dir_var"):
            from cif_variance import DEFAULT_INPUT_DIR
            self._variance_input_dir_var = ctk.StringVar(value=str(DEFAULT_INPUT_DIR))
            self._variance_input_dir_var.trace_add("write", self._update_variance_cif_count)
        ctk.CTkEntry(
            self._steps_outer, textvariable=self._variance_input_dir_var, width=340,
        ).grid(row=1, column=1, padx=6, pady=6, sticky="ew")
        ctk.CTkButton(
            self._steps_outer, text="Browse", width=70, height=26,
            font=ctk.CTkFont(size=12),
            command=self._browse_variance_input_dir,
        ).grid(row=1, column=2, padx=12, pady=6, sticky="e")

        self._variance_cif_count_label = ctk.CTkLabel(
            self._steps_outer, text="", anchor="w", font=ctk.CTkFont(size=11),
        )
        self._variance_cif_count_label.grid(row=2, column=1, columnspan=2, padx=6, pady=(0, 6), sticky="w")
        self._update_variance_cif_count()

        # Top N
        ctk.CTkLabel(self._steps_outer, text="Top N residues:", anchor="w").grid(
            row=3, column=0, padx=(12, 6), pady=6, sticky="w"
        )
        if not hasattr(self, "_variance_top_var"):
            self._variance_top_var = ctk.StringVar(value="10")
        ctk.CTkEntry(
            self._steps_outer, textvariable=self._variance_top_var, width=60,
        ).grid(row=3, column=1, padx=6, pady=6, sticky="w")

        # Report range
        ctk.CTkLabel(self._steps_outer, text="Report range:", anchor="w").grid(
            row=4, column=0, padx=(12, 6), pady=6, sticky="w"
        )
        report_frame = ctk.CTkFrame(self._steps_outer, fg_color="transparent")
        report_frame.grid(row=4, column=1, columnspan=2, padx=6, pady=6, sticky="w")
        if not hasattr(self, "_variance_range_start_var"):
            self._variance_range_start_var = ctk.StringVar(value="")
            self._variance_range_end_var = ctk.StringVar(value="")
        ctk.CTkEntry(report_frame, textvariable=self._variance_range_start_var, width=70,
                     placeholder_text="start").pack(side="left", padx=(0, 6))
        ctk.CTkEntry(report_frame, textvariable=self._variance_range_end_var, width=70,
                     placeholder_text="end (blank = all)").pack(side="left")

        # Align range
        ctk.CTkLabel(self._steps_outer, text="Align range:", anchor="w").grid(
            row=5, column=0, padx=(12, 6), pady=6, sticky="w"
        )
        align_frame = ctk.CTkFrame(self._steps_outer, fg_color="transparent")
        align_frame.grid(row=5, column=1, columnspan=2, padx=6, pady=6, sticky="w")
        if not hasattr(self, "_variance_align_start_var"):
            self._variance_align_start_var = ctk.StringVar(value="")
            self._variance_align_end_var = ctk.StringVar(value="")
        ctk.CTkEntry(align_frame, textvariable=self._variance_align_start_var, width=70,
                     placeholder_text="start").pack(side="left", padx=(0, 6))
        ctk.CTkEntry(align_frame, textvariable=self._variance_align_end_var, width=70,
                     placeholder_text="end (blank = same as report range)").pack(side="left")

        # UniProt / gene overrides
        ctk.CTkLabel(self._steps_outer, text="UniProt override:", anchor="w").grid(
            row=6, column=0, padx=(12, 6), pady=6, sticky="w"
        )
        if not hasattr(self, "_variance_uniprot_var"):
            self._variance_uniprot_var = ctk.StringVar(value="")
        ctk.CTkEntry(
            self._steps_outer, textvariable=self._variance_uniprot_var, width=150,
            placeholder_text="auto-detected from CIF",
        ).grid(row=6, column=1, padx=6, pady=6, sticky="w")

        ctk.CTkLabel(self._steps_outer, text="Gene override:", anchor="w").grid(
            row=7, column=0, padx=(12, 6), pady=6, sticky="w"
        )
        if not hasattr(self, "_variance_gene_var"):
            self._variance_gene_var = ctk.StringVar(value="")
        ctk.CTkEntry(
            self._steps_outer, textvariable=self._variance_gene_var, width=150,
            placeholder_text="optional, for UniProt lookup",
        ).grid(row=7, column=1, padx=6, pady=6, sticky="w")

        # Status label + progress bar (reuse the step status pattern)
        status = ctk.CTkLabel(
            self._steps_outer, text="●  Ready", width=100,
            anchor="e", text_color=_GRAY,
        )
        status.grid(row=8, column=1, columnspan=2, padx=12, pady=6, sticky="e")
        self._step_status_labels.append(status)

        bar = ctk.CTkProgressBar(self._steps_outer, width=120, height=14)
        bar.set(0)
        bar.grid(row=8, column=0, padx=12, pady=6, sticky="w")
        bar.grid_remove()
        self._step_progress_bars.append(bar)

    def _browse_variance_input_dir(self):
        """Open a folder dialog for selecting the CIF-comparison input directory."""
        path = filedialog.askdirectory(
            title="Select CIF comparison folder",
            initialdir=self._variance_input_dir_var.get(),
        )
        if path:
            self._variance_input_dir_var.set(path)

    def _update_variance_cif_count(self, *_args):
        """Update the live .cif file count label for the CIF Variance input folder."""
        if not hasattr(self, "_variance_cif_count_label"):
            return
        path = Path(self._variance_input_dir_var.get().strip())
        n = len(list(path.glob("*.cif"))) if path.is_dir() else 0
        if n >= 2:
            self._variance_cif_count_label.configure(text=f"✓  {n} .cif files found", text_color=_GREEN)
        else:
            self._variance_cif_count_label.configure(
                text=f"⚠  {n} .cif file(s) found — need at least 2", text_color=_YELLOW,
            )

    def _build_ca_coordinates_panel(self):
        """Build the input fields for CA-coordinate export mode."""
        ctk.CTkLabel(
            self._steps_outer,
            text="Export Alpha-Carbon Coordinates",
            font=ctk.CTkFont(weight="bold"),
        ).grid(row=0, column=0, columnspan=3, padx=12, pady=(8, 2), sticky="w")

        # UniProt ID
        ctk.CTkLabel(self._steps_outer, text="UniProt ID:", anchor="w").grid(
            row=1, column=0, padx=(12, 6), pady=6, sticky="w"
        )
        if not hasattr(self, "_ca_uniprot_var"):
            self._ca_uniprot_var = ctk.StringVar(value="")
        ctk.CTkEntry(
            self._steps_outer, textvariable=self._ca_uniprot_var, width=200,
            placeholder_text="e.g. P04637",
        ).grid(row=1, column=1, padx=6, pady=6, sticky="w")

        # Gene override
        ctk.CTkLabel(self._steps_outer, text="Gene override:", anchor="w").grid(
            row=2, column=0, padx=(12, 6), pady=6, sticky="w"
        )
        if not hasattr(self, "_ca_gene_var"):
            self._ca_gene_var = ctk.StringVar(value="")
        ctk.CTkEntry(
            self._steps_outer, textvariable=self._ca_gene_var, width=200,
            placeholder_text="optional, skips UniProt API lookup",
        ).grid(row=2, column=1, padx=6, pady=6, sticky="w")

        # Status label + progress bar (reuse the step status pattern)
        status = ctk.CTkLabel(
            self._steps_outer, text="●  Ready", width=100,
            anchor="e", text_color=_GRAY,
        )
        status.grid(row=3, column=1, columnspan=2, padx=12, pady=6, sticky="e")
        self._step_status_labels.append(status)

        bar = ctk.CTkProgressBar(self._steps_outer, width=120, height=14)
        bar.set(0)
        bar.grid(row=3, column=0, padx=12, pady=6, sticky="w")
        bar.grid_remove()
        self._step_progress_bars.append(bar)

    # ── Timer ────────────────────────────────────────────────────────────────

    def _tick(self) -> None:
        if not self._running or self._pipeline_start is None:
            return
        if self._suspended:
            return
        elapsed = time.time() - self._pipeline_start

        text = f"Elapsed: {_fmt_time(elapsed)}"

        est = getattr(self, "_precheck_estimate", None)
        hist = self._historical_times
        if est is not None:
            text += f"  |  Est. total: ~{_fmt_time(est)}"
        elif hist and len(hist) == self._total_steps:
            text += f"  |  Est. total: ~{_fmt_time(sum(hist))}"
        elif self._step_times:
            avg = sum(self._step_times) / len(self._step_times)
            text += f"  |  Est. total: ~{_fmt_time(avg * self._total_steps)}"

        self._timer_label.configure(text=text)
        self.after(1000, self._tick)

    # ── Activity animation ─────────────────────────────────────────────────

    def _start_activity_animation(self):
        self._activity_animating = True
        self._activity_pos = -self._activity_chunk
        self._animate_activity()

    def _stop_activity_animation(self):
        self._activity_animating = False
        self._activity_canvas.delete("chunk")

    def _animate_activity(self):
        if not self._activity_animating:
            return
        self._activity_canvas.delete("chunk")
        x1 = max(0, self._activity_pos)
        x2 = min(self._activity_width, self._activity_pos + self._activity_chunk)
        if x2 > x1:
            self._activity_canvas.create_rectangle(
                x1, 0, x2, self._activity_height,
                fill=_BLUE, outline="", tags="chunk",
            )
        self._activity_pos += 1
        if self._activity_pos > self._activity_width:
            self._activity_pos = -self._activity_chunk
        self.after(25, self._animate_activity)

    # ── File-status bar ──────────────────────────────────────────────────────

    def _refresh_file_status(self):
        """Update the status indicator for each input folder."""
        for name, (folder, exts, _desc) in _INPUT_FOLDERS.items():
            lbl = self._file_indicators[name]
            try:
                f = resolve_input_file(folder, exts)
                lbl.configure(text=f"✓  {name}: {f.name}", text_color=_GREEN)
            except FileNotFoundError:
                lbl.configure(text=f"✗  {name}: no file", text_color=_RED)
            except RuntimeError:
                lbl.configure(text=f"⚠  {name}: multiple files", text_color=_YELLOW)

    def _browse_file(self, name: str, folder: Path, filetypes: list) -> None:
        """Open a file dialog, copy the selected file into the input folder, and refresh status."""
        path = filedialog.askopenfilename(
            title=f"Select {name} input file",
            filetypes=filetypes + [("All files", "*.*")],
        )
        if not path:
            return

        src = Path(path)
        folder.mkdir(parents=True, exist_ok=True)

        for existing in folder.iterdir():
            if existing.is_file():
                existing.unlink()

        shutil.copy2(src, folder / src.name)
        self._refresh_file_status()

    def _toggle_log(self):
        """Show or hide the raw log output panel (a normal part of the scrollable content)."""
        if self._log_visible:
            self._log.grid_remove()
            self._log_toggle.configure(text="Show Details")
            self._log_visible = False
        else:
            self._log.grid(row=10, column=0, padx=24, pady=(6, 12), sticky="ew")
            self._log_toggle.configure(text="Hide Details")
            self._log_visible = True

    def _browse_output_dir(self):
        """Let the user pick a custom output folder."""
        path = filedialog.askdirectory(
            title="Select output folder",
            initialdir=self._output_dir_var.get(),
        )
        if path:
            self._output_dir_var.set(path)

    @property
    def _output_dir(self) -> Path:
        return Path(self._output_dir_var.get())

    # ── Pipeline execution ───────────────────────────────────────────────────

    _DEFAULT_OUTPUT_FILES = [
        "ptm_mutation_proximity_db.tsv",
        "ptm_mutation_proximity_long.tsv",
        "mutation_cluster_db.tsv",
    ]

    @property
    def _output_files(self) -> list[Path]:
        return [self._output_dir / name for name in self._DEFAULT_OUTPUT_FILES]

    def _find_locked_output_file(self) -> Path | None:
        """Return the first existing output file that cannot be opened for writing, or None."""
        for path in self._output_files:
            if not path.exists():
                continue
            try:
                with path.open("r+b"):
                    pass
            except PermissionError:
                return path
        return None

    def _backup_outputs(self) -> list[Path]:
        """Copy existing output files to .bak before the pipeline overwrites them."""
        backups = []
        for path in self._output_files:
            if path.exists():
                bak = path.with_suffix(path.suffix + ".bak")
                shutil.copy2(path, bak)
                backups.append(bak)
        return backups

    @staticmethod
    def _restore_backups(backups: list[Path]) -> None:
        """Restore .bak files over their originals (undo a partial pipeline run)."""
        for bak in backups:
            if bak.exists():
                original = bak.with_suffix("")  # strip .bak
                shutil.copy2(bak, original)

    def _remove_backups(self, backups: list[Path]) -> None:
        for bak in backups:
            try:
                if bak.exists():
                    bak.unlink()
            except Exception as exc:
                self._q("log", f"Warning: could not delete {bak.name}: {exc}")

    def _start_pipeline(self):
        if self._running:
            return

        locked = self._find_locked_output_file()
        if locked:
            from tkinter import messagebox
            messagebox.showerror(
                "Output File Not Writable",
                f"Cannot start the pipeline — the output file\n\n"
                f"    {locked.name}\n\n"
                f"cannot be written to. Common causes:\n"
                f"  • Open in another program (e.g. Excel)\n"
                f"  • File is marked read-only\n"
                f"  • A cloud sync service (OneDrive, Dropbox) has it locked\n"
                f"  • You do not have write permission on the output folder\n\n"
                f"Close or resolve the issue above and try again.",
            )
            return

        mode = self._mode.get()
        if mode == "single-protein":
            cif = getattr(self, "_single_cif_var", None)
            if not cif or not cif.get().strip():
                from tkinter import messagebox
                messagebox.showwarning("Missing input", "Please select a CIF file first.")
                return
        elif mode == "radius-sweep":
            from tkinter import messagebox
            genes = self._radius_genes_var.get().split()
            if not genes:
                messagebox.showwarning("Missing input", "Enter at least one gene symbol.")
                return
            try:
                start = float(self._radius_start_var.get())
                stop = float(self._radius_stop_var.get())
                step = float(self._radius_step_var.get())
                if step <= 0 or stop <= start:
                    raise ValueError
            except ValueError:
                messagebox.showwarning(
                    "Invalid input",
                    "Radius start/stop/step must be numeric, with stop > start and step > 0.",
                )
                return
            ptm_tsv = PROJECT_ROOT / "data" / "steps" / "PTMD_TCGA_hotspots_by_protein.tsv"
            if not ptm_tsv.exists():
                messagebox.showerror(
                    "Missing data",
                    f"Required file not found:\n{ptm_tsv}\n\n"
                    f"Run the PTM Proximity or Mutation Clustering pipeline (step 1) first.",
                )
                return
        elif mode == "cif-variance":
            from tkinter import messagebox
            input_dir = Path(self._variance_input_dir_var.get().strip())
            n_cifs = len(list(input_dir.glob("*.cif"))) if input_dir.is_dir() else 0
            if n_cifs < 2:
                messagebox.showwarning(
                    "Missing input", f"Need at least 2 .cif files in\n{input_dir}\n\nFound {n_cifs}.",
                )
                return
            if self._variance_top_var.get().strip() and not self._variance_top_var.get().strip().isdigit():
                messagebox.showwarning("Invalid input", "Top N must be a whole number.")
                return
            for start_var, end_var, label in [
                (self._variance_range_start_var, self._variance_range_end_var, "Report range"),
                (self._variance_align_start_var, self._variance_align_end_var, "Align range"),
            ]:
                s, e = start_var.get().strip(), end_var.get().strip()
                if s or e:
                    if not (s.isdigit() and e.isdigit()) or int(s) >= int(e):
                        messagebox.showwarning(
                            "Invalid input", f"{label} must be two whole numbers with start < end.",
                        )
                        return
        elif mode == "ca-coordinates":
            from tkinter import messagebox
            if not self._ca_uniprot_var.get().strip():
                messagebox.showwarning("Missing input", "Enter a UniProt accession.")
                return

        self._running = True
        self._stop_requested = False
        self._suspended = False
        self._run_btn.configure(state="disabled")
        self._stop_btn.configure(state="normal", fg_color=_RED)
        self._timer_label.configure(text="")

        self._log.configure(state="normal")
        self._log.delete("1.0", "end")
        self._log.configure(state="disabled")

        for lbl in self._step_status_labels:
            lbl.configure(text="●  Waiting", text_color=_GRAY)
        for bar in self._step_progress_bars:
            bar.set(0)
            bar.grid_remove()

        mode = self._mode.get()
        threading.Thread(target=self._run_pipeline, args=(mode,), daemon=True).start()

    def _stop_pipeline(self):
        """Suspend the running subprocess immediately and show Resume/Cancel options."""
        if not self._running:
            return
        proc = self._current_proc
        actually_suspended = False
        if proc and proc.poll() is None:
            try:
                psutil.Process(proc.pid).suspend()
                self._suspended = True
                actually_suspended = True
            except psutil.NoSuchProcess:
                pass
        self._stop_btn.configure(state="disabled", fg_color="gray30")
        if actually_suspended:
            self._paused_elapsed = time.time() - self._pipeline_start
            self._stop_activity_animation()
            self._activity_canvas.pack_forget()
            self._run_btn.configure(
                state="normal", text="▶  Resume", command=self._resume_pipeline,
            )
            self._cancel_btn_ref = ctk.CTkButton(
                self._run_btn.master,
                text="✗  Cancel",
                command=self._cancel_pipeline,
                width=110,
                height=44,
                font=ctk.CTkFont(size=15, weight="bold"),
                fg_color=_RED,
                hover_color="#b82020",
            )
            self._cancel_btn_ref.pack(side="left", padx=(8, 0))
            self._q("log", "Pipeline paused.")

    def _resume_pipeline(self):
        """Resume the suspended subprocess."""
        proc = self._current_proc
        if proc and self._suspended:
            # Adjust start time so elapsed doesn't include the pause duration
            pause_duration = time.time() - (self._pipeline_start + self._paused_elapsed)
            self._pipeline_start += pause_duration
            try:
                psutil.Process(proc.pid).resume()
            except psutil.NoSuchProcess:
                pass
            self._suspended = False
        if hasattr(self, "_cancel_btn_ref"):
            self._cancel_btn_ref.destroy()
            del self._cancel_btn_ref
        self._run_btn.configure(
            state="disabled", text="▶  Run Pipeline", command=self._start_pipeline,
        )
        self._stop_btn.configure(state="normal", fg_color=_RED)
        self._activity_canvas.pack(side="right", padx=(0, 4))
        self._start_activity_animation()
        self._tick()
        self._q("log", "Pipeline resumed.")

    def _cancel_pipeline(self):
        """Kill the suspended subprocess and restore previous output."""
        self._stop_requested = True
        proc = self._current_proc
        if proc and proc.poll() is None:
            if self._suspended:
                try:
                    psutil.Process(proc.pid).resume()
                except psutil.NoSuchProcess:
                    pass
                self._suspended = False
            proc.terminate()
        if hasattr(self, "_cancel_btn_ref"):
            self._cancel_btn_ref.destroy()
            del self._cancel_btn_ref
        self._run_btn.configure(
            state="normal", text="▶  Run Pipeline", command=self._start_pipeline,
        )
        self._stop_btn.configure(state="disabled", fg_color="gray30")
        self._open_btn.configure(state="normal", fg_color=_BLUE, hover_color="#2563eb")
        for lbl in self._step_status_labels:
            lbl.configure(text="●  Waiting", text_color=_GRAY)
        for bar in self._step_progress_bars:
            bar.set(0)
            bar.grid_remove()

    # Per-item time estimates (seconds) for runtime calculation
    _TIME_PER_CIF_DOWNLOAD = 2.5
    _TIME_PER_UNIPROT_BATCH = 1.5       # ~100 IDs per batch
    _TIME_PER_1433_FETCH = 0.2          # 5 concurrent workers
    _TIME_PER_PP_FETCH = 0.05           # 10 concurrent workers
    _TIME_PER_KINASE_PREDICT = 0.5
    _TIME_PER_AIUPRED_FETCH = 2.0       # per-protein REST call, sequential
    _TIME_PER_PROTEIN_STEP3 = 0.35
    _TIME_STEP1_BASE = 40               # base time for filtering/merging
    _TIME_STEP4_BASE = 40               # base time for reading/writing the proximity DB

    def _run_precheck(self, mode: str) -> float:
        """Analyze caches and data to estimate total pipeline runtime. Returns seconds."""
        import pandas as pd

        self._q("log", "Initializing pipeline...")
        self._q("log", "")

        input_tsv = PROJECT_ROOT / "data" / "steps" / "PTMD_TCGA_hotspots_by_protein.tsv"
        cache_dir = PROJECT_ROOT / "data" / "cache"
        models_dir = PROJECT_ROOT / "cif_models"

        # Count proteins from intermediate TSV (if exists) or estimate from input
        n_proteins = 0
        if input_tsv.exists():
            try:
                df = pd.read_csv(input_tsv, sep="\t", usecols=["uniprot_id"], dtype=str)
                n_proteins = df["uniprot_id"].nunique()
            except Exception:
                pass

        # Step 1: UniProt API cache
        gene_cache = cache_dir / "uniprot_gene_mapping.tsv"
        cached_genes = 0
        if gene_cache.exists():
            try:
                cached_genes = len(pd.read_csv(gene_cache, sep="\t", dtype=str))
            except Exception:
                pass
        uncached_batches = max(0, (n_proteins - cached_genes)) // 100 + 1 if n_proteins > cached_genes else 0
        step1_est = self._TIME_STEP1_BASE + uncached_batches * self._TIME_PER_UNIPROT_BATCH
        self._q("log", f"Step 1: {cached_genes} UniProt gene mappings cached")

        # Step 2: CIF downloads
        cifs_present = 0
        if models_dir.exists():
            cifs_present = sum(1 for d in models_dir.iterdir()
                               if d.is_dir() and any(d.glob("*.cif")))
        cifs_needed = max(0, n_proteins - cifs_present)
        step2_est = cifs_needed * self._TIME_PER_CIF_DOWNLOAD + 5
        self._q("log", f"Step 2: {cifs_present} CIF files present"
                + (f", ~{cifs_needed} to download" if cifs_needed else " (all cached)"))

        # Step 3: local computation
        step3_est = n_proteins * self._TIME_PER_PROTEIN_STEP3 if n_proteins else 60
        self._q("log", f"Step 3: {n_proteins} proteins to process")

        # Step 4: annotation caches (only for ptm-proximity)
        step4_est = 0
        if mode == "ptm-proximity":
            # 14-3-3
            cache_1433 = cache_dir / "1433pred"
            cached_1433 = len(list(cache_1433.glob("*.json"))) if cache_1433.exists() else 0
            uncached_1433 = max(0, n_proteins - cached_1433)

            # PolyPhen
            pp_cache = cache_dir / "polyphen.tsv"
            cached_pp = 0
            if pp_cache.exists():
                try:
                    cached_pp = len(pd.read_csv(pp_cache, sep="\t", dtype=str))
                except Exception:
                    pass

            # Kinase
            kin_cache = cache_dir / "kinase_predictions.tsv"
            cached_kin = 0
            if kin_cache.exists():
                try:
                    cached_kin = len(pd.read_csv(kin_cache, sep="\t", dtype=str))
                except Exception:
                    pass

            # AIUPred — one binding call per protein (yields general + binding scores)
            cached_aiupred = 0
            aiupred_cache = cache_dir / "aiupred_disorder.tsv"
            if aiupred_cache.exists():
                try:
                    df_ac = pd.read_csv(aiupred_cache, sep="\t", dtype=str,
                                        keep_default_na=False)
                    cached_aiupred = int((df_ac["analysis_type"] == "binding").sum())
                except Exception:
                    pass
            uncached_aiupred = max(0, n_proteins - cached_aiupred)

            self._q("log", f"Step 4: {cached_1433}/{n_proteins} 14-3-3 predictions cached, "
                    f"{cached_pp} PolyPhen pairs cached, {cached_kin} kinase windows cached, "
                    f"{cached_aiupred}/{n_proteins} AIUPred cached")

            step4_est = (uncached_1433 * self._TIME_PER_1433_FETCH
                         + max(0, n_proteins * 2 - cached_pp) * self._TIME_PER_PP_FETCH
                         + max(0, n_proteins - cached_kin) * self._TIME_PER_KINASE_PREDICT
                         + uncached_aiupred * self._TIME_PER_AIUPRED_FETCH
                         + self._TIME_STEP4_BASE)

        total_est = step1_est + step2_est + step3_est + step4_est
        self._q("log", "")
        self._q("log", f"Estimated runtime: ~{_fmt_time(total_est)}")
        self._q("log", "")

        return total_est

    def _run_pipeline(self, mode: str):
        if mode == "single-protein":
            self._run_single_protein()
            return
        if mode == "radius-sweep":
            self._run_radius_sweep()
            return
        if mode == "cif-variance":
            self._run_cif_variance()
            return
        if mode == "ca-coordinates":
            self._run_ca_coordinates()
            return

        python = [sys.executable, "-u"]
        input_tsv = PROJECT_ROOT / "data" / "steps" / "PTMD_TCGA_hotspots_by_protein.tsv"
        models_dir = PROJECT_ROOT / "cif_models"

        cutoff = self._cutoff_var.get().strip() or "10.0"
        min_samples = self._min_samples_var.get().strip() or "3"
        min_plddt = self._min_plddt_var.get().strip()
        max_pae = self._max_pae_var.get().strip()
        pp_exclude = []
        if not self._pp_benign_var.get():
            pp_exclude.append("benign")
        if not self._pp_possibly_var.get():
            pp_exclude.append("possibly_damaging")
        if not self._pp_probably_var.get():
            pp_exclude.append("probably_damaging")

        cmds = [
            [*python, str(SCRIPTS_DIR / "1_filter.py"), "--mode", mode,
             "--min-samples", min_samples],
            [
                *python, str(SCRIPTS_DIR / "2_download_structures.py"),
                str(input_tsv),
                "--id_column", "uniprot_id",
                "--out_dir", str(models_dir),
                "--prefer", "cif",
                "--delay", "0.1",
                "--also_pae",
                "--logs_dir", str(self._output_dir / "logs"),
            ],
            [*python, str(SCRIPTS_DIR / "3_find_nearby_mutations.py"), "--mode", mode,
             "--output-dir", str(self._output_dir), "--cutoff", cutoff,
             *(["--min-plddt", min_plddt] if min_plddt else []),
             *(["--max-pae", max_pae] if max_pae else [])],
        ]
        if mode == "ptm-proximity":
            cmds.append([*python, str(SCRIPTS_DIR / "4_annotate.py"),
                         "--output-dir", str(self._output_dir),
                         *(["--pp-exclude"] + pp_exclude if pp_exclude else [])])

        steps = PTM_PROXIMITY_STEPS if mode == "ptm-proximity" else MUTATION_CLUSTERING_STEPS
        run_type = _detect_run_type()
        self._q("pipeline_start", len(steps), mode, run_type)

        estimated_total = self._run_precheck(mode)
        self._q("set_estimate", estimated_total)

        backups = self._backup_outputs()

        try:
            all_ok = True
            cancelled = False
            step_times: list[float] = []
            for i, (step_def, cmd) in enumerate(zip(steps, cmds)):
                _panel_label, log_label = step_def
                if self._stop_requested:
                    self._q("status", i, "●  Cancelled", _YELLOW)
                    cancelled = True
                    break

                self._q("status", i, "▶  Running…", _BLUE)
                self._q("show_progress", i)
                self._q("log", f"Step {i + 1}/{len(steps)}: {log_label}")
                self._q("step_start")

                t0 = time.time()
                ok = self._stream_cmd(cmd, i)
                elapsed = time.time() - t0

                if self._stop_requested:
                    self._q("hide_progress", i)
                    self._q("status", i, "■  Stopped", _YELLOW)
                    self._q("log", f"Step {i + 1} stopped by user")
                    cancelled = True
                    break

                if ok:
                    step_times.append(elapsed)
                    self._q("progress", i, 1.0, "Done")
                    self._q("hide_progress", i)
                    self._q("status", i, f"✓  {_fmt_time(elapsed)}", _GREEN)
                    self._q("step_complete", elapsed)
                else:
                    self._q("hide_progress", i)
                    self._q("status", i, "✗  Failed", _RED)
                    self._q("log", f"Step {i + 1} FAILED after {_fmt_time(elapsed)}")
                    all_ok = False
                    break

            restore_ok = False
            if cancelled:
                try:
                    self._restore_backups(backups)
                    self._q("log", "Pipeline cancelled — previous output restored.")
                    restore_ok = True
                except Exception:
                    self._q("log", "Pipeline cancelled. Output file may be locked — "
                            "close it and rename the .bak file to restore your previous output.")
            elif not all_ok:
                try:
                    self._restore_backups(backups)
                    self._q("log", "Pipeline failed — previous output restored.")
                    restore_ok = True
                except Exception:
                    self._q("log", "Pipeline failed. Output file may be locked — "
                            "close it and rename the .bak file to restore your previous output.")
            else:
                restore_ok = True
                self._q("save_runtimes", mode, run_type, step_times)
                self._q("log", "Pipeline complete! Output saved to the Output/ folder.")
        except Exception as exc:
            restore_ok = False
            self._q("log", f"Pipeline error: {exc}")
        finally:
            if restore_ok:
                self._remove_backups(backups)
            self._q("finished")

    def _run_single_protein(self):
        """Run the single-protein CIF analysis in the background thread."""
        cif_path = self._single_cif_var.get().strip()
        uniprot = self._single_uniprot_var.get().strip()

        if not cif_path:
            self._q("log", "Error: no CIF file selected.")
            self._q("finished")
            return

        cmd = [sys.executable, "-u", str(SCRIPTS_DIR / "analyze_single_cif_nearby_mutations.py"), cif_path]
        if uniprot:
            cmd.extend(["--uniprot", uniprot])

        self._q("pipeline_start", 1, "single-protein", "warm")
        self._q("show_log")
        self._q("status", 0, "▶  Analyzing…", _BLUE)
        self._q("show_progress", 0)
        self._q("log", f"Analyzing CIF: {Path(cif_path).name}")
        if uniprot:
            self._q("log", f"UniProt ID: {uniprot}")
        self._q("log", "")

        t0 = time.time()
        ok = self._stream_cmd(cmd, 0)
        elapsed = time.time() - t0

        if self._stop_requested:
            self._q("hide_progress", 0)
            self._q("status", 0, "■  Stopped", _YELLOW)
            self._q("finished")
            return

        if not ok:
            self._q("hide_progress", 0)
            self._q("status", 0, "✗  Failed", _RED)
            self._q("log", f"Analysis failed after {_fmt_time(elapsed)}")
            self._q("finished")
            return

        self._q("progress", 0, 1.0, "Done")
        self._q("hide_progress", 0)
        self._q("status", 0, f"✓  {_fmt_time(elapsed)}", _GREEN)
        self._q("log", "")
        self._q("log", f"Analysis complete in {_fmt_time(elapsed)}.")

        output_db = self._output_dir / "ptm_mutation_proximity_db.tsv"
        if output_db.exists():
            self._q("log", "Review the results above, then choose whether to append to the existing database.")
            self._q("offer_append", cif_path, uniprot)
        else:
            self._q("log", "No existing proximity database found — results printed above.")

        self._q("enable_open")
        self._q("finished")

    def _run_radius_sweep(self):
        """Run the radius-sweep analysis in-process, in the background thread."""
        import numpy as np

        genes = self._radius_genes_var.get().split()
        try:
            start = float(self._radius_start_var.get())
            stop = float(self._radius_stop_var.get())
            step = float(self._radius_step_var.get())
        except ValueError:
            self._q("log", "Error: invalid radius range.")
            self._q("finished")
            return
        radii = list(np.arange(start, stop + step / 2, step))
        unfiltered = self._radius_unfiltered_var.get()

        self._q("pipeline_start", 1, "radius-sweep", "warm")
        self._q("show_log")
        self._q("status", 0, "▶  Running sweep…", _BLUE)
        self._q("show_progress", 0)
        self._q("log", f"Testing radii {radii[0]:.0f}-{radii[-1]:.0f} Å in {step:.0f} Å steps "
                        f"for {len(genes)} gene(s)")
        self._q("log", "")

        t0 = time.time()
        try:
            from radius_sweep import run_sweep
            result = run_sweep(
                genes, radii, unfiltered=unfiltered,
                output_tsv_path=self._output_dir / "radius_sweep.tsv",
                log_cb=lambda line: self._q("log", line),
            )
        except ImportError as exc:
            self._q("status", 0, "✗  Failed", _RED)
            self._q("log", f"Missing dependency: {exc}. Run: uv sync")
            self._q("hide_progress", 0)
            self._q("finished")
            return
        except (FileNotFoundError, ValueError) as exc:
            self._q("status", 0, "✗  Failed", _RED)
            self._q("log", f"Error: {exc}")
            self._q("hide_progress", 0)
            self._q("finished")
            return
        except Exception as exc:
            self._q("status", 0, "✗  Failed", _RED)
            self._q("log", f"Unexpected error: {exc}")
            self._q("hide_progress", 0)
            self._q("finished")
            return

        elapsed = time.time() - t0
        self._q("hide_progress", 0)
        self._q("status", 0, f"✓  {_fmt_time(elapsed)}", _GREEN)
        self._q("log", "")
        self._q("log", f"Sweep complete in {_fmt_time(elapsed)}.")
        self._q("viz_data", "radius-sweep", result)
        self._q("enable_open")
        self._q("finished")

    def _run_cif_variance(self):
        """Run the CIF variance analysis in-process, in the background thread."""
        input_dir = Path(self._variance_input_dir_var.get().strip())
        top = int(self._variance_top_var.get().strip() or 10)

        def _range_or_none(start_var, end_var):
            s, e = start_var.get().strip(), end_var.get().strip()
            return (int(s), int(e)) if s and e else None

        range_ = _range_or_none(self._variance_range_start_var, self._variance_range_end_var)
        align_range = _range_or_none(self._variance_align_start_var, self._variance_align_end_var)
        uniprot = self._variance_uniprot_var.get().strip() or None
        gene = self._variance_gene_var.get().strip() or None

        self._q("pipeline_start", 1, "cif-variance", "warm")
        self._q("show_log")
        self._q("status", 0, "▶  Running analysis…", _BLUE)
        self._q("show_progress", 0)
        self._q("log", f"Comparing CIF files in {input_dir}")
        self._q("log", "")

        t0 = time.time()
        try:
            from cif_variance import run_variance_analysis
            result = run_variance_analysis(
                input_dir=input_dir,
                output_dir=self._output_dir / "cif_variance",
                top=top,
                range_=range_,
                align_range=align_range,
                uniprot=uniprot,
                gene=gene,
                log_cb=lambda line: self._q("log", line),
            )
        except ImportError as exc:
            self._q("status", 0, "✗  Failed", _RED)
            self._q("log", f"Missing dependency: {exc}. Run: uv sync")
            self._q("hide_progress", 0)
            self._q("finished")
            return
        except (FileNotFoundError, ValueError) as exc:
            self._q("status", 0, "✗  Failed", _RED)
            self._q("log", f"Error: {exc}")
            self._q("hide_progress", 0)
            self._q("finished")
            return
        except Exception as exc:
            self._q("status", 0, "✗  Failed", _RED)
            self._q("log", f"Unexpected error: {exc}")
            self._q("hide_progress", 0)
            self._q("finished")
            return

        elapsed = time.time() - t0
        self._q("hide_progress", 0)
        self._q("status", 0, f"✓  {_fmt_time(elapsed)}", _GREEN)
        self._q("log", "")
        self._q("log", f"Analysis complete in {_fmt_time(elapsed)}.")
        self._q("viz_data", "cif-variance", result)
        self._q("enable_open")
        self._q("finished")

    def _run_ca_coordinates(self):
        """Run the CA-coordinate export in-process, in the background thread."""
        uniprot = self._ca_uniprot_var.get().strip()
        gene = self._ca_gene_var.get().strip() or None

        self._q("pipeline_start", 1, "ca-coordinates", "warm")
        self._q("show_log")
        self._q("status", 0, "▶  Exporting…", _BLUE)
        self._q("show_progress", 0)
        self._q("log", f"Exporting CA coordinates for {uniprot}")
        self._q("log", "")

        t0 = time.time()
        try:
            from export_ca_coordinates import run_export
            run_export(
                uniprot, gene=gene,
                output_dir=self._output_dir / "coordinates",
                log_cb=lambda line: self._q("log", line),
            )
        except ImportError as exc:
            self._q("status", 0, "✗  Failed", _RED)
            self._q("log", f"Missing dependency: {exc}. Run: uv sync")
            self._q("hide_progress", 0)
            self._q("finished")
            return
        except (FileNotFoundError, ValueError, RuntimeError) as exc:
            self._q("status", 0, "✗  Failed", _RED)
            self._q("log", f"Error: {exc}")
            self._q("hide_progress", 0)
            self._q("finished")
            return
        except Exception as exc:
            self._q("status", 0, "✗  Failed", _RED)
            self._q("log", f"Unexpected error: {exc}")
            self._q("hide_progress", 0)
            self._q("finished")
            return

        elapsed = time.time() - t0
        self._q("hide_progress", 0)
        self._q("status", 0, f"✓  {_fmt_time(elapsed)}", _GREEN)
        self._q("enable_open")
        self._q("finished")

    _TQDM_RE = re.compile(r"^(.*?):\s+(\d+)%\|[^|]*\|\s+(\d+)/(\d+)")
    _OVERALL_RE = re.compile(r"^##PROGRESS##\s+(\d+)\s+(.+)$")
    _using_overall_progress = False

    def _parse_tqdm(self, text: str, step_idx: int) -> bool:
        """Try to parse a progress line (custom ##PROGRESS## or tqdm). Returns True if matched."""
        m_overall = self._OVERALL_RE.search(text)
        if m_overall:
            self._using_overall_progress = True
            pct = int(m_overall.group(1))
            desc = m_overall.group(2)
            self._q("progress", step_idx, pct / 100, f"{pct}%")
            self._q("progress_log", desc)
            return True
        m = self._TQDM_RE.search(text)
        if not m:
            return False
        if self._using_overall_progress:
            return True
        desc, pct, current, total = m.group(1), int(m.group(2)), m.group(3), m.group(4)
        self._q("progress", step_idx, pct / 100, f"{pct}% ({current}/{total})")
        self._q("progress_log", f"{desc}: {current}/{total}")
        return True

    def _stream_cmd(self, cmd: list[str], step_idx: int) -> bool:
        """Run a subprocess, parsing tqdm output for progress bar updates."""
        env = {**__import__("os").environ, "PYTHONIOENCODING": "utf-8"}
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, env=env)
        self._current_proc = proc
        buf = b""
        while True:
            byte = proc.stdout.read(1)
            if not byte:
                break
            if self._suspended:
                buf = b""
                continue
            if byte == b"\r":
                text = buf.decode("utf-8", errors="replace")
                if self._parse_tqdm(text, step_idx):
                    buf = b""
            elif byte == b"\n":
                text = buf.decode("utf-8", errors="replace").strip()
                if text and not self._parse_tqdm(text, step_idx):
                    self._q("log", text)
                buf = b""
            else:
                buf += byte

        if buf:
            text = buf.decode("utf-8", errors="replace").strip()
            if text:
                self._q("log", text)

        proc.wait()
        return proc.returncode == 0

    # ── Queue helpers ────────────────────────────────────────────────────────

    def _q(self, *args):
        self._queue.put(args)

    def _poll_queue(self):
        try:
            while True:
                msg = self._queue.get_nowait()
                kind = msg[0]
                if kind == "log":
                    self._append_log(msg[1])
                elif kind == "show_log":
                    if not self._log_visible:
                        self._toggle_log()
                elif kind == "status":
                    _, idx, text, color = msg
                    if 0 <= idx < len(self._step_status_labels):
                        self._step_status_labels[idx].configure(text=text, text_color=color)
                elif kind == "progress":
                    _, idx, pct, status_text = msg
                    if 0 <= idx < len(self._step_progress_bars):
                        self._step_progress_bars[idx].set(pct)
                        self._step_status_labels[idx].configure(
                            text=status_text, text_color=_BLUE
                        )
                elif kind == "progress_log":
                    self._update_progress_line(msg[1])
                elif kind == "show_progress":
                    _, idx = msg
                    self._using_overall_progress = False
                    if 0 <= idx < len(self._step_progress_bars):
                        self._step_progress_bars[idx].set(0)
                        self._step_progress_bars[idx].grid()
                elif kind == "hide_progress":
                    _, idx = msg
                    if 0 <= idx < len(self._step_progress_bars):
                        self._step_progress_bars[idx].grid_remove()
                elif kind == "pipeline_start":
                    _, total, mode, run_type = msg
                    self._pipeline_start = time.time()
                    self._step_start = None
                    self._total_steps = total
                    self._steps_done = 0
                    self._step_times = []
                    self._historical_times = _load_runtimes(mode, run_type)
                    self._precheck_estimate = None
                    self._activity_canvas.pack(side="right", padx=(0, 8))
                    self._start_activity_animation()
                    self._tick()
                elif kind == "set_estimate":
                    self._precheck_estimate = msg[1]
                elif kind == "step_start":
                    self._step_start = time.time()
                elif kind == "step_complete":
                    _, elapsed = msg
                    self._steps_done += 1
                    self._step_times.append(elapsed)
                    self._step_start = None
                elif kind == "save_runtimes":
                    _, mode, run_type, times = msg
                    _save_runtimes(mode, run_type, times)
                elif kind == "offer_append":
                    _, cif_path, uniprot = msg
                    self._show_append_dialog(cif_path, uniprot)
                elif kind == "viz_data":
                    _, which, result = msg
                    if which == "radius-sweep":
                        self._radius_sweep_result = result
                        self._tabview.set("Radius Sweep")
                        self._draw_radius_sweep_plot()
                    elif which == "cif-variance":
                        self._cif_variance_result = result
                        self._tabview.set("CIF Variance")
                        self._draw_cif_variance_plot()
                elif kind == "enable_open":
                    self._open_btn.configure(state="normal", fg_color=_BLUE, hover_color="#2563eb")
                elif kind == "finished":
                    was_cancelled = self._stop_requested
                    self._running = False
                    self._stop_requested = False
                    self._suspended = False
                    self._current_proc = None
                    self._stop_activity_animation()
                    self._activity_canvas.pack_forget()
                    if hasattr(self, "_cancel_btn_ref"):
                        self._cancel_btn_ref.destroy()
                        del self._cancel_btn_ref
                    self._run_btn.configure(
                        state="normal", text="▶  Run Pipeline",
                        command=self._start_pipeline,
                    )
                    self._stop_btn.configure(state="disabled", fg_color="gray30")
                    self._open_btn.configure(state="normal", fg_color=_BLUE, hover_color="#2563eb")
                    if was_cancelled:
                        for lbl in self._step_status_labels:
                            lbl.configure(text="●  Waiting", text_color=_GRAY)
                        for bar in self._step_progress_bars:
                            bar.set(0)
                            bar.grid_remove()
                    if self._pipeline_start is not None:
                        total = time.time() - self._pipeline_start
                        self._timer_label.configure(
                            text=f"Total: {_fmt_time(total)}", text_color=_GREEN
                        )
        except queue.Empty:
            pass
        self.after(100, self._poll_queue)

    def _append_log(self, line: str):
        self._log.configure(state="normal")
        self._log.insert("end", line + "\n")
        self._log.see("end")
        self._log.configure(state="disabled")
        self._last_log_was_progress = False

    def _update_progress_line(self, text: str):
        """Replace the last line in the log if it was a progress update, otherwise append."""
        self._log.configure(state="normal")
        if getattr(self, "_last_log_was_progress", False):
            self._log.delete("end-2l", "end-1l")
        self._log.insert("end", text + "\n")
        self._log.see("end")
        self._log.configure(state="disabled")
        self._last_log_was_progress = True

    # ── Single-protein append dialog ────────────────────────────────────────

    def _check_existing_uniprot(self, uniprot: str) -> int:
        """Return the number of rows in the proximity DB that match this UniProt ID."""
        output_db = self._output_dir / "ptm_mutation_proximity_db.tsv"
        if not output_db.exists() or not uniprot:
            return 0
        try:
            import pandas as pd
            df = pd.read_csv(output_db, sep="\t", encoding="utf-16", usecols=["UniProt"],
                             dtype=str, keep_default_na=False)
            return int((df["UniProt"] == uniprot).sum())
        except Exception:
            return 0

    def _remove_uniprot_rows(self, uniprot: str) -> None:
        """Remove all rows for a UniProt ID from the proximity DB."""
        output_db = self._output_dir / "ptm_mutation_proximity_db.tsv"
        if not output_db.exists():
            return
        try:
            import pandas as pd
            df = pd.read_csv(output_db, sep="\t", encoding="utf-16", dtype=str,
                             keep_default_na=False)
            before = len(df)
            df = df[df["UniProt"] != uniprot]
            df.to_csv(output_db, sep="\t", index=False, encoding="utf-16")
            self._q("log", f"Removed {before - len(df)} existing row(s) for {uniprot}.")
        except Exception as exc:
            self._q("log", f"Warning: could not remove existing rows: {exc}")

    def _show_append_dialog(self, cif_path: str, uniprot: str):
        """Ask the user whether to append, replace, or skip based on existing data."""
        existing_count = self._check_existing_uniprot(uniprot)

        dialog = ctk.CTkToplevel(self)
        dialog.title("Save to Database?")
        dialog.resizable(False, False)
        dialog.grab_set()

        dialog.geometry("500x190")

        if existing_count > 0:
            msg = (f"The output database already contains {existing_count} row(s) "
                   f"for {uniprot}.\nHow would you like to proceed?")
        else:
            msg = (f"The output database has no existing data for {uniprot}.\n"
                   f"How would you like to proceed?")

        ctk.CTkLabel(
            dialog, text=msg, font=ctk.CTkFont(size=14), justify="center",
        ).pack(pady=(20, 16))

        btn_frame = ctk.CTkFrame(dialog, fg_color="transparent")
        btn_frame.pack()

        def do_replace():
            dialog.destroy()
            threading.Thread(
                target=self._replace_and_append, args=(cif_path, uniprot),
                daemon=True,
            ).start()

        def do_append():
            dialog.destroy()
            threading.Thread(
                target=self._append_single_protein, args=(cif_path, uniprot),
                daemon=True,
            ).start()

        if existing_count > 0:
            ctk.CTkButton(
                btn_frame, text="Replace existing", width=140,
                command=do_replace,
            ).pack(side="left", padx=6)

        ctk.CTkButton(
            btn_frame, text="Append to database", width=140,
            command=do_append,
        ).pack(side="left", padx=6)

        ctk.CTkButton(
            btn_frame, text="Skip", width=100,
            fg_color="gray30", hover_color="gray40",
            command=dialog.destroy,
        ).pack(side="left", padx=6)

    def _replace_and_append(self, cif_path: str, uniprot: str):
        """Remove existing rows for this UniProt ID, then append new results."""
        self._remove_uniprot_rows(uniprot)
        self._append_single_protein(cif_path, uniprot)

    def _append_single_protein(self, cif_path: str, uniprot: str):
        """Re-run the single-protein script with --append-to-db."""
        cmd = [
            sys.executable, "-u",
            str(SCRIPTS_DIR / "analyze_single_cif_nearby_mutations.py"),
            cif_path, "--append-to-db",
        ]
        if uniprot:
            cmd.extend(["--uniprot", uniprot])
        self._q("log", "Appending results to proximity database...")
        proc = subprocess.Popen(
            cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
            text=True, encoding="utf-8", errors="replace",
        )
        for line in proc.stdout:
            self._q("log", line.rstrip())
        proc.wait()
        if proc.returncode == 0:
            self._q("log", "Results appended successfully.")
        else:
            self._q("log", "Failed to append results.")

    # ── Tab change handler ───────────────────────────────────────────────────

    def _on_tab_change(self) -> None:
        tab = self._tabview.get()
        if tab == "Results":
            self._load_results()
        elif tab == "Visualization" and self._results_df_wide is None:
            self._load_results()

    # ── Results tab ──────────────────────────────────────────────────────────

    def _setup_treeview_style(self) -> None:
        import tkinter.ttk as ttk
        style = ttk.Style()
        style.theme_use("default")
        bg, fg = "#2b2b2b", "#dcdcdc"
        heading_bg = "#3a3a3a"
        for name in ("Results.Treeview",):
            style.configure(name,
                background=bg, foreground=fg, rowheight=24,
                fieldbackground=bg, borderwidth=0, relief="flat",
                font=("Segoe UI", 10),
            )
            style.configure(f"{name}.Heading",
                background=heading_bg, foreground=fg,
                relief="flat", borderwidth=1,
                font=("Segoe UI", 10, "bold"),
            )
            style.map(name,
                background=[("selected", _BLUE)],
                foreground=[("selected", "white")],
            )
            style.map(f"{name}.Heading",
                background=[("active", "#4a4a4a")],
            )

    def _make_treeview(self, parent, col_defs: list):
        import tkinter as tk
        import tkinter.ttk as ttk

        frame = tk.Frame(parent, bg="#2b2b2b")
        frame.pack(fill=tk.BOTH, expand=True, padx=8, pady=(0, 8))
        frame.grid_rowconfigure(0, weight=1)
        frame.grid_columnconfigure(0, weight=1)

        tv = ttk.Treeview(frame, style="Results.Treeview", show="headings",
                           selectmode="browse")
        col_ids = [c[1] for c in col_defs]
        tv["columns"] = col_ids

        numeric_cols = {"#col", "near", "far", "total", "pts", "cosmic", "seqd", "dist", "pps", "pae", "mpld"}
        for display, col_id, width in col_defs:
            anchor = "e" if col_id in numeric_cols else "w"
            tv.heading(col_id, text=display,
                       command=lambda c=col_id: self._sort_tv(tv, c, False))
            stretch = col_id == col_ids[-1]
            tv.column(col_id, width=width, minwidth=30, stretch=stretch, anchor=anchor)

        vsb = ttk.Scrollbar(frame, orient="vertical", command=tv.yview)
        hsb = ttk.Scrollbar(frame, orient="horizontal", command=tv.xview)
        tv.configure(yscrollcommand=vsb.set, xscrollcommand=hsb.set)

        tv.tag_configure("odd",  background="#2b2b2b")
        tv.tag_configure("even", background="#313131")

        tv.grid(row=0, column=0, sticky="nsew")
        vsb.grid(row=0, column=1, sticky="ns")
        hsb.grid(row=1, column=0, sticky="ew")

        return tv

    def _sort_tv(self, tv, col: str, reverse: bool) -> None:
        def _key(v: str):
            try:
                return (0, float(v))
            except (ValueError, TypeError):
                return (1, str(v).lower())

        rows = [(tv.set(k, col), k) for k in tv.get_children("")]
        rows.sort(key=lambda x: _key(x[0]), reverse=reverse)
        for idx, (_, k) in enumerate(rows):
            tv.move(k, "", idx)
        if "#col" in tv["columns"]:
            for idx, k in enumerate(tv.get_children(""), 1):
                tv.set(k, "#col", idx)
        tv.heading(col, command=lambda: self._sort_tv(tv, col, not reverse))

    def _clear_treeview_fully(self, tv, all_rows: list) -> None:
        """Delete every row a treeview has ever held, including ones currently
        hidden by a search filter (plain tv.delete(*tv.get_children()) would
        miss detached rows, leaking them and causing iid collisions on reinsert).
        """
        all_iids = set(tv.get_children("")) | {iid for iid, _, _ in all_rows}
        if all_iids:
            tv.delete(*all_iids)

    def _capture_tv_rows(self, tv) -> list:
        """Snapshot (iid, values, tags) for every row, to filter against later.

        Must be called immediately after a full (unfiltered) populate, while
        every row is still attached.
        """
        return [(iid, tv.item(iid, "values"), tv.item(iid, "tags")) for iid in tv.get_children("")]

    def _filter_treeview(self, tv, all_rows: list, query: str) -> None:
        """Show only rows whose displayed values contain *query* (case-insensitive).

        Uses detach()/move() rather than delete(), so hidden rows keep their
        iid and can reappear — this preserves the iid-as-dataframe-position
        scheme that selection handlers rely on.
        """
        query = query.strip().lower()
        attached = set(tv.get_children(""))
        shown = 0
        for iid, values, _tags in all_rows:
            haystack = " ".join(str(v) for v in values).lower()
            if not query or query in haystack:
                tv.move(iid, "", shown)
                shown += 1
            elif iid in attached:
                tv.detach(iid)

    def _build_results_tab(self, tab) -> None:
        import tkinter as tk

        self._results_df_wide = None
        self._results_df_long = None
        self._ptm_tv_all_rows: list = []
        self._mut_tv_all_rows: list = []
        self._setup_treeview_style()

        tab.grid_columnconfigure(0, weight=1)
        tab.grid_rowconfigure(0, weight=1)

        outer = ctk.CTkFrame(tab, fg_color="transparent")
        outer.grid(row=0, column=0, sticky="nsew", padx=4, pady=4)
        outer.grid_columnconfigure(0, weight=1)
        outer.grid_rowconfigure(1, weight=3)
        outer.grid_rowconfigure(3, weight=2)

        # ── Top panel header ──
        top_header = ctk.CTkFrame(outer, fg_color="transparent")
        top_header.grid(row=0, column=0, sticky="ew", padx=8, pady=(8, 2))
        ctk.CTkLabel(top_header, text="PTM Sites",
                     font=ctk.CTkFont(size=14, weight="bold")).pack(side=tk.LEFT)
        self._results_status = ctk.CTkLabel(top_header, text="",
                                             text_color="gray60",
                                             font=ctk.CTkFont(size=11))
        self._results_status.pack(side=tk.LEFT, padx=(12, 0))
        ctk.CTkButton(top_header, text="↺  Refresh", width=90, height=28,
                       font=ctk.CTkFont(size=12),
                       command=self._load_results).pack(side=tk.RIGHT)
        ctk.CTkButton(top_header, text="📈  Visualize", width=100, height=28,
                       font=ctk.CTkFont(size=12),
                       command=self._visualize_selected_ptm).pack(side=tk.RIGHT, padx=(0, 6))
        self._results_ptm_search_var = ctk.StringVar(value="")
        self._ptm_search_entry = ctk.CTkEntry(
            top_header, textvariable=self._results_ptm_search_var,
            width=220, placeholder_text="🔍 Search PTM sites… (Ctrl+F)",
        )
        self._ptm_search_entry.pack(side=tk.RIGHT, padx=(0, 12))
        self._results_ptm_search_var.trace_add("write", self._filter_ptm_tv)

        # ── PTM site treeview ──
        top_frame = ctk.CTkFrame(outer)
        top_frame.grid(row=1, column=0, sticky="nsew", padx=8, pady=(0, 4))
        top_frame.grid_rowconfigure(0, weight=1)
        top_frame.grid_columnconfigure(0, weight=1)
        self._ptm_tv = self._make_treeview(top_frame, _PTM_TV_COLS)
        self._ptm_tv.bind("<<TreeviewSelect>>", self._on_ptm_select)
        self._ptm_tv.bind("<Double-Button-1>", lambda _e: self._visualize_selected_ptm())

        # ── Detail panel header ──
        bot_header = ctk.CTkFrame(outer, fg_color="transparent")
        bot_header.grid(row=2, column=0, sticky="ew", padx=8, pady=(4, 2))
        ctk.CTkLabel(bot_header, text="Mutation Details",
                     font=ctk.CTkFont(size=14, weight="bold")).pack(side=tk.LEFT, padx=(8, 0))
        self._results_mut_search_var = ctk.StringVar(value="")
        self._mut_search_entry = ctk.CTkEntry(
            bot_header, textvariable=self._results_mut_search_var,
            width=220, placeholder_text="🔍 Search mutations…",
        )
        self._mut_search_entry.pack(side=tk.RIGHT, padx=(0, 6))
        self._results_mut_search_var.trace_add("write", self._filter_mut_tv)

        # ── Mutation detail treeview ──
        bot_frame = ctk.CTkFrame(outer)
        bot_frame.grid(row=3, column=0, sticky="nsew", padx=8, pady=(0, 8))
        bot_frame.grid_rowconfigure(0, weight=1)
        bot_frame.grid_columnconfigure(0, weight=1)
        self._mut_tv = self._make_treeview(bot_frame, _MUT_TV_COLS)

        # Ctrl+F focuses the PTM search box whenever the Results tab is active
        self.bind_all("<Control-f>", self._focus_results_search)

    def _build_viz_tab(self, tab) -> None:
        import matplotlib
        matplotlib.use("TkAgg")
        from matplotlib.figure import Figure
        from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

        tab.grid_columnconfigure(0, weight=1)
        tab.grid_rowconfigure(2, weight=1)

        self._viz_ptm_rows: dict[str, int] = {}
        self._viz_all_labels: list[str] = []

        # ── Controls row 1: PTM selection ──
        controls = ctk.CTkFrame(tab)
        controls.grid(row=0, column=0, padx=12, pady=(12, 6), sticky="ew")

        ctk.CTkLabel(
            controls, text="PTM site:", font=ctk.CTkFont(weight="bold"),
        ).pack(side="left", padx=(12, 6), pady=10)

        self._viz_search_var = ctk.StringVar(value="")
        search_entry = ctk.CTkEntry(
            controls, textvariable=self._viz_search_var, width=160,
            placeholder_text="Search gene or site…",
        )
        search_entry.pack(side="left", padx=(0, 6), pady=10)
        self._viz_search_var.trace_add("write", self._on_viz_search)

        self._viz_combo = ctk.CTkComboBox(controls, width=220, values=[])
        self._viz_combo.pack(side="left", padx=(0, 12), pady=10)

        ctk.CTkLabel(controls, text="±aa window:").pack(side="left", padx=(8, 4), pady=10)
        self._viz_window_var = ctk.StringVar(value="15")
        ctk.CTkEntry(
            controls, textvariable=self._viz_window_var, width=50,
        ).pack(side="left", padx=(0, 12), pady=10)

        # ── Controls row 2: display mode + actions ──
        controls2 = ctk.CTkFrame(tab)
        controls2.grid(row=1, column=0, padx=12, pady=(0, 6), sticky="ew")

        ctk.CTkLabel(controls2, text="Show:", font=ctk.CTkFont(weight="bold")).pack(
            side="left", padx=(12, 6), pady=10,
        )
        self._viz_mode_var = ctk.StringVar(value="All mutations")
        ctk.CTkSegmentedButton(
            controls2, values=["All mutations", "Unique per position"],
            variable=self._viz_mode_var,
            command=lambda _v: self._generate_lollipop_plot(),
        ).pack(side="left", padx=(0, 12), pady=10)

        ctk.CTkButton(
            controls2, text="Generate", width=100,
            command=self._generate_lollipop_plot,
        ).pack(side="left", padx=(0, 8), pady=10)

        ctk.CTkButton(
            controls2, text="Save PNG", width=100,
            fg_color="gray30", hover_color="gray40",
            command=self._save_viz_plot,
        ).pack(side="left", padx=(0, 12), pady=10)

        self._viz_status = ctk.CTkLabel(
            controls2, text="Load results to select a PTM site.",
            text_color="gray60", font=ctk.CTkFont(size=11),
        )
        self._viz_status.pack(side="left", padx=(4, 12), pady=10)

        # ── Plot canvas ──
        canvas_frame = ctk.CTkFrame(tab, fg_color="#2b2b2b")
        canvas_frame.grid(row=2, column=0, padx=12, pady=(0, 12), sticky="nsew")
        canvas_frame.grid_columnconfigure(0, weight=1)
        canvas_frame.grid_rowconfigure(0, weight=1)

        self._viz_fig = Figure(figsize=(10, 6), dpi=100, facecolor="#2b2b2b")
        self._viz_canvas = FigureCanvasTkAgg(self._viz_fig, master=canvas_frame)
        self._viz_canvas.get_tk_widget().grid(row=0, column=0, sticky="nsew")

    # ── Visualization: PTM selection ────────────────────────────────────────

    def _refresh_viz_selector(self, df) -> None:
        """(Re)populate the PTM-site combobox from the currently loaded results."""
        self._viz_ptm_rows = {}
        labels: list[str] = []
        for idx, row in df.iterrows():
            gene = row.get("gene", "?")
            site = row.get("ptm_site", "?")
            uid = row.get("UniProt", "?")
            label = f"{gene}  {site}  ({uid})"
            labels.append(label)
            self._viz_ptm_rows[label] = idx
        self._viz_all_labels = labels
        current = self._viz_combo.get()
        self._viz_combo.configure(values=labels)
        if labels and current not in labels:
            self._viz_combo.set(labels[0])
        elif not labels:
            self._viz_combo.set("")

    def _on_viz_search(self, *_args) -> None:
        query = self._viz_search_var.get().strip().lower()
        filtered = (
            [label for label in self._viz_all_labels if query in label.lower()]
            if query else self._viz_all_labels
        )
        self._viz_combo.configure(values=filtered)
        if filtered and self._viz_combo.get() not in filtered:
            self._viz_combo.set(filtered[0])

    # ── Visualization: data prep ─────────────────────────────────────────────

    def _get_viz_mutation_df(self, uid: str, ptm_site: str, wide_row):
        """Return (DataFrame, note) of nearby mutations for a PTM site.

        Prefers the long-format table (has per-mutation patient counts);
        falls back to parsing the wide-format summary columns otherwise.
        """
        import pandas as pd

        if self._results_df_long is not None:
            df = self._results_df_long
            mask = (
                (df.get("uniprot_id", "") == uid) &
                (df.get("ptm_position", "") == ptm_site)
            )
            sub = df[mask].copy()
            if not sub.empty:
                def _pos(m):
                    match = re.search(r"\d+", str(m))
                    return int(match.group()) if match else None

                sub["mutation_position"] = sub["mutation"].apply(_pos)
                sub = sub.dropna(subset=["mutation_position"])
                sub["patient_count"] = pd.to_numeric(
                    sub.get("patient_count", 0), errors="coerce"
                ).fillna(0)
                sub["distance_angstrom"] = pd.to_numeric(
                    sub.get("distance_angstrom", 0), errors="coerce"
                ).fillna(0)
                if "polyphen_class" not in sub.columns:
                    sub["polyphen_class"] = ""
                return sub[[
                    "mutation", "mutation_position", "patient_count",
                    "distance_angstrom", "polyphen_class",
                ]], None

        # Fallback: wide-format columns only (no per-mutation patient counts)
        rows = []
        for col_key in ("mutations_within_5_positions", "mutations_more_than_5_positions"):
            for entry in (wide_row.get(col_key, "") or "").split(", "):
                entry = entry.strip()
                if not entry:
                    continue
                m = _MUT_ENTRY_RE.match(entry)
                if not m:
                    continue
                pos_m = re.search(r"\d+", m.group(1))
                if not pos_m:
                    continue
                rows.append({
                    "mutation": m.group(1),
                    "mutation_position": int(pos_m.group()),
                    "patient_count": 1,
                    "distance_angstrom": float(m.group(4)),
                    "polyphen_class": _PP_LABEL.get(m.group(2) or "", ""),
                })
        df_fallback = pd.DataFrame(rows)
        note = None
        if not df_fallback.empty:
            note = "patient counts unavailable — re-run the pipeline to generate the long-format table"
        return df_fallback, note

    # ── Visualization: plotting ──────────────────────────────────────────────

    def _style_lollipop_axis(self, ax) -> None:
        ax.set_facecolor("#2b2b2b")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["left"].set_color("#555555")
        ax.spines["bottom"].set_color("#555555")
        ax.tick_params(colors="#dcdcdc")
        ax.yaxis.label.set_color("#dcdcdc")
        ax.xaxis.label.set_color("#dcdcdc")

    def _draw_mutation_needles(self, ax, df, jitter: bool, x_col: str = "mutation_position") -> None:
        from collections import defaultdict

        positions = df[x_col].tolist()
        groups: dict = defaultdict(list)
        for i, pos in enumerate(positions):
            groups[round(pos)].append(i)

        x_plot = list(positions)
        if jitter:
            for pos, idxs in groups.items():
                n = len(idxs)
                if n > 1:
                    for k, i in enumerate(idxs):
                        x_plot[i] = pos + (k - (n - 1) / 2) * 0.6

        for i, (_, r) in enumerate(df.iterrows()):
            x = x_plot[i]
            count = r["patient_count"]
            color = _PP_COLORS.get(str(r.get("polyphen_class", "")).strip(), _NEEDLE_DEFAULT_COLOR)
            ax.vlines(x, 0, count, color=color, linewidth=2, zorder=2)
            ax.scatter([x], [count], s=90, color=color, edgecolor="black", zorder=3)
            # The far (categorical) panel already names each mutation via its x-tick
            # label, so only annotate above the marker for the true-scale local panel.
            if jitter:
                ax.annotate(
                    str(r["mutation"]), (x, count), xytext=(4, 10),
                    textcoords="offset points", ha="left", va="bottom",
                    fontsize=8, color="#dcdcdc", rotation=45,
                )

    def _draw_lollipop(self, gene: str, ptm_site: str, ptm_pos: int, mut_df, local_window: float) -> None:
        from matplotlib.lines import Line2D
        from matplotlib.patches import Patch

        fig = self._viz_fig
        fig.clf()
        fig.patch.set_facecolor("#2b2b2b")

        local_df = mut_df[(mut_df["mutation_position"] - ptm_pos).abs() <= local_window] \
            .sort_values("mutation_position")
        far_df = mut_df[(mut_df["mutation_position"] - ptm_pos).abs() > local_window] \
            .sort_values("mutation_position").reset_index(drop=True)
        has_far = not far_df.empty

        max_count = max(float(mut_df["patient_count"].max()), 1.0)
        ptm_height = max_count * 1.15
        y_top = ptm_height * 1.45

        if has_far:
            gs = fig.add_gridspec(
                1, 3, width_ratios=[0.68, 0.05, 0.27], wspace=0.05,
                left=0.08, right=0.97, top=0.86, bottom=0.28,
            )
            ax_local = fig.add_subplot(gs[0])
            ax_gap = fig.add_subplot(gs[1])
            ax_far = fig.add_subplot(gs[2])
        else:
            gs = fig.add_gridspec(1, 1, left=0.08, right=0.97, top=0.86, bottom=0.28)
            ax_local = fig.add_subplot(gs[0])
            ax_gap = None
            ax_far = None

        # ── Local panel (true sequence scale) ──
        self._style_lollipop_axis(ax_local)
        x_min, x_max = ptm_pos - local_window, ptm_pos + local_window
        ax_local.hlines(0, x_min, x_max, color=_NEEDLE_DEFAULT_COLOR, linewidth=2, zorder=1)
        ax_local.set_xlim(x_min - 0.5, x_max + 0.5)
        ax_local.set_ylim(0, y_top)
        ax_local.set_ylabel("Patient count")
        ax_local.set_xlabel("Residue position")
        if has_far:
            ax_local.spines["right"].set_visible(False)

        ax_local.vlines(ptm_pos, 0, ptm_height, color=_PTM_MARKER_COLOR, linewidth=2, zorder=2)
        ax_local.scatter(
            [ptm_pos], [ptm_height], marker="*", s=420,
            color=_PTM_MARKER_COLOR, edgecolor="white", linewidth=0.8, zorder=4,
        )
        ax_local.annotate(
            str(ptm_site), (ptm_pos, ptm_height), xytext=(0, 10),
            textcoords="offset points", ha="center", color="white",
            fontsize=11, fontweight="bold",
        )

        self._draw_mutation_needles(ax_local, local_df, jitter=True)

        # ── Broken-axis gap + far panel (categorical scale) ──
        if has_far:
            ax_gap.axis("off")
            ax_gap.set_xlim(0, 1)
            ax_gap.set_ylim(0, 1)
            for xc in (0.3, 0.7):
                ax_gap.plot([xc - 0.12, xc + 0.12], [0.42, 0.52],
                            color="#888888", lw=1.2, transform=ax_gap.transAxes, clip_on=False)
                ax_gap.plot([xc - 0.12, xc + 0.12], [0.48, 0.58],
                            color="#888888", lw=1.2, transform=ax_gap.transAxes, clip_on=False)

            self._style_lollipop_axis(ax_far)
            far_df["_x"] = range(len(far_df))
            n_far = len(far_df)
            ax_far.hlines(0, -0.5, n_far - 0.5, color=_NEEDLE_DEFAULT_COLOR, linewidth=2, zorder=1)
            ax_far.set_xlim(-0.5, n_far - 0.5)
            ax_far.set_ylim(0, y_top)
            ax_far.set_yticklabels([])
            ax_far.spines["left"].set_visible(False)

            self._draw_mutation_needles(ax_far, far_df, jitter=False, x_col="_x")
            ax_far.set_xticks(list(range(n_far)))
            ax_far.set_xticklabels(far_df["mutation"], rotation=45, ha="right", fontsize=8)
            ax_far.set_title(f"> {local_window:g} aa away", fontsize=9, color="gray")

        legend_handles = [
            Line2D([0], [0], marker="*", color="none", markerfacecolor=_PTM_MARKER_COLOR,
                   markeredgecolor="white", markersize=14, label="PTM site"),
            Patch(facecolor=_PP_COLORS["benign"], label="Benign"),
            Patch(facecolor=_PP_COLORS["possibly_damaging"], label="Possibly damaging"),
            Patch(facecolor=_PP_COLORS["probably_damaging"], label="Probably damaging"),
            Patch(facecolor=_NEEDLE_DEFAULT_COLOR, label="Unknown / not scored"),
        ]
        ax_local.legend(
            handles=legend_handles, loc="upper left", fontsize=8,
            facecolor="#3a3a3a", edgecolor="#555555", labelcolor="#dcdcdc",
        )

        fig.suptitle(f"{gene} — {ptm_site} and nearby mutations", color="white",
                     fontsize=14, fontweight="bold")

        self._viz_canvas.draw()

    def _generate_lollipop_plot(self) -> None:
        if self._results_df_wide is None:
            self._viz_status.configure(
                text="No results loaded — run the pipeline or open the Results tab first.",
                text_color=_RED,
            )
            return

        label = self._viz_combo.get()
        idx = self._viz_ptm_rows.get(label)
        if idx is None:
            self._viz_status.configure(text="Select a PTM site to plot.", text_color=_RED)
            return

        row = self._results_df_wide.loc[idx]
        gene = row.get("gene", "?")
        ptm_site = row.get("ptm_site", "?")
        uid = row.get("UniProt", "?")

        ptm_m = re.search(r"\d+", str(ptm_site))
        if not ptm_m:
            self._viz_status.configure(
                text=f"Could not parse a residue position from '{ptm_site}'.", text_color=_RED,
            )
            return
        ptm_pos = int(ptm_m.group())

        mut_df, note = self._get_viz_mutation_df(uid, ptm_site, row)
        if mut_df.empty:
            self._viz_fig.clf()
            self._viz_canvas.draw()
            self._viz_status.configure(
                text=f"No nearby mutations found for {gene} {ptm_site}.", text_color=_YELLOW,
            )
            return

        unique_only = self._viz_mode_var.get() == "Unique per position"
        if unique_only:
            mut_df = self._collapse_to_unique_positions(mut_df)

        try:
            local_window = float(self._viz_window_var.get().strip() or "15")
        except ValueError:
            local_window = 15.0

        self._draw_lollipop(gene, ptm_site, ptm_pos, mut_df, local_window)

        kind = "unique position(s)" if unique_only else "nearby mutation(s)"
        status = f"{len(mut_df)} {kind} for {gene} {ptm_site}"
        if note:
            status += f" — {note}"
        self._viz_status.configure(text=status, text_color="gray60")

    def _collapse_to_unique_positions(self, df):
        """Merge mutations that share a residue position into one entry per position.

        Patient counts are summed across substitutions at that residue; the
        displayed color reflects the most severe PolyPhen class present.
        """
        import pandas as pd

        _SEVERITY = {"probably_damaging": 3, "possibly_damaging": 2, "benign": 1}

        def _worst_class(classes) -> str:
            best, best_rank = "", -1
            for c in classes:
                c = str(c).strip()
                rank = _SEVERITY.get(c, 0)
                if rank > best_rank:
                    best, best_rank = c, rank
            return best

        def _position_label(mutation: str) -> str:
            m = re.match(r"^([A-Za-z]+\d+)", str(mutation))
            return m.group(1) if m else str(mutation)

        rows = []
        for pos, grp in df.groupby("mutation_position"):
            label = _position_label(grp["mutation"].iloc[0])
            if len(grp) > 1:
                label = f"{label} (×{len(grp)})"
            rows.append({
                "mutation": label,
                "mutation_position": pos,
                "patient_count": grp["patient_count"].sum(),
                "distance_angstrom": grp["distance_angstrom"].min(),
                "polyphen_class": _worst_class(grp["polyphen_class"]),
            })
        return pd.DataFrame(rows)

    def _save_viz_plot(self) -> None:
        if not self._viz_fig.axes:
            self._viz_status.configure(text="Generate a plot before saving.", text_color=_RED)
            return
        label = self._viz_combo.get()
        safe = re.sub(r"[^\w-]+", "_", label).strip("_") or "lollipop"
        out_dir = self._output_dir
        out_dir.mkdir(parents=True, exist_ok=True)
        out_path = out_dir / f"lollipop_{safe}.png"
        self._viz_fig.savefig(out_path, dpi=200, facecolor=self._viz_fig.get_facecolor(), bbox_inches="tight")
        self._viz_status.configure(text=f"Saved to {out_path}", text_color=_GREEN)

    # ── Radius Sweep / CIF Variance: shared dark-theme styling ────────────────

    def _style_dark_figure(self, fig) -> None:
        """Apply the app's dark theme to every axes of a Figure built by an external script."""
        fig.patch.set_facecolor("#2b2b2b")
        for ax in fig.axes:
            self._style_lollipop_axis(ax)
            if ax.get_title():
                ax.title.set_color("#dcdcdc")
            legend = ax.get_legend()
            if legend is not None:
                legend.get_frame().set_facecolor("#3a3a3a")
                legend.get_frame().set_edgecolor("#555555")
                for text in legend.get_texts():
                    text.set_color("#dcdcdc")

    # ── Radius Sweep tab ────────────────────────────────────────────────────

    def _build_radius_sweep_tab(self, tab) -> None:
        import matplotlib
        matplotlib.use("TkAgg")
        from matplotlib.figure import Figure
        from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

        tab.grid_columnconfigure(0, weight=1)
        tab.grid_rowconfigure(1, weight=1)

        controls = ctk.CTkFrame(tab)
        controls.grid(row=0, column=0, padx=12, pady=(12, 6), sticky="ew")

        ctk.CTkButton(
            controls, text="Save PNG", width=100,
            fg_color="gray30", hover_color="gray40",
            command=self._save_radius_sweep_plot,
        ).pack(side="left", padx=(12, 12), pady=10)

        self._radius_sweep_status = ctk.CTkLabel(
            controls, text="Run Radius Sweep mode from the Pipeline tab to see results here.",
            text_color="gray60", font=ctk.CTkFont(size=11),
        )
        self._radius_sweep_status.pack(side="left", padx=(0, 12), pady=10)

        canvas_frame = ctk.CTkFrame(tab, fg_color="#2b2b2b")
        canvas_frame.grid(row=1, column=0, padx=12, pady=(0, 12), sticky="nsew")
        canvas_frame.grid_columnconfigure(0, weight=1)
        canvas_frame.grid_rowconfigure(0, weight=1)

        self._radius_sweep_fig = Figure(figsize=(14, 9), dpi=100, facecolor="#2b2b2b")
        self._radius_sweep_canvas = FigureCanvasTkAgg(self._radius_sweep_fig, master=canvas_frame)
        self._radius_sweep_canvas.get_tk_widget().grid(row=0, column=0, sticky="nsew")

    def _draw_radius_sweep_plot(self) -> None:
        from radius_sweep import build_sweep_figure

        self._radius_sweep_fig.clf()
        build_sweep_figure(self._radius_sweep_result, fig=self._radius_sweep_fig)
        self._style_dark_figure(self._radius_sweep_fig)
        self._radius_sweep_canvas.draw()
        self._radius_sweep_status.configure(text="Sweep complete.", text_color=_GREEN)

    def _save_radius_sweep_plot(self) -> None:
        if not self._radius_sweep_fig.axes:
            self._radius_sweep_status.configure(text="Run a sweep before saving.", text_color=_RED)
            return
        out_dir = self._output_dir
        out_dir.mkdir(parents=True, exist_ok=True)
        out_path = out_dir / "radius_sweep_plot.png"
        self._radius_sweep_fig.savefig(
            out_path, dpi=200, facecolor=self._radius_sweep_fig.get_facecolor(), bbox_inches="tight",
        )
        self._radius_sweep_status.configure(text=f"Saved to {out_path}", text_color=_GREEN)

    # ── CIF Variance tab ────────────────────────────────────────────────────

    def _build_cif_variance_tab(self, tab) -> None:
        import matplotlib
        matplotlib.use("TkAgg")
        from matplotlib.figure import Figure
        from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

        tab.grid_columnconfigure(0, weight=1)
        tab.grid_rowconfigure(1, weight=1)

        controls = ctk.CTkFrame(tab)
        controls.grid(row=0, column=0, padx=12, pady=(12, 6), sticky="ew")

        ctk.CTkButton(
            controls, text="Save PNG", width=100,
            fg_color="gray30", hover_color="gray40",
            command=self._save_cif_variance_plot,
        ).pack(side="left", padx=(12, 12), pady=10)

        self._cif_variance_status = ctk.CTkLabel(
            controls, text="Run CIF Variance mode from the Pipeline tab to see results here.",
            text_color="gray60", font=ctk.CTkFont(size=11),
        )
        self._cif_variance_status.pack(side="left", padx=(0, 12), pady=10)

        canvas_frame = ctk.CTkFrame(tab, fg_color="#2b2b2b")
        canvas_frame.grid(row=1, column=0, padx=12, pady=(0, 12), sticky="nsew")
        canvas_frame.grid_columnconfigure(0, weight=1)
        canvas_frame.grid_rowconfigure(0, weight=1)

        self._cif_variance_fig = Figure(figsize=(14, 8), dpi=100, facecolor="#2b2b2b")
        self._cif_variance_canvas = FigureCanvasTkAgg(self._cif_variance_fig, master=canvas_frame)
        self._cif_variance_canvas.get_tk_widget().grid(row=0, column=0, sticky="nsew")

    def _draw_cif_variance_plot(self) -> None:
        from cif_variance import build_variance_figure

        self._cif_variance_fig.clf()
        build_variance_figure(self._cif_variance_result, fig=self._cif_variance_fig)
        self._style_dark_figure(self._cif_variance_fig)
        self._cif_variance_canvas.draw()
        self._cif_variance_status.configure(text="Analysis complete.", text_color=_GREEN)

    def _save_cif_variance_plot(self) -> None:
        if not self._cif_variance_fig.axes:
            self._cif_variance_status.configure(text="Run an analysis before saving.", text_color=_RED)
            return
        out_dir = self._output_dir
        out_dir.mkdir(parents=True, exist_ok=True)
        out_path = out_dir / "cif_variance_plot.png"
        self._cif_variance_fig.savefig(
            out_path, dpi=200, facecolor=self._cif_variance_fig.get_facecolor(), bbox_inches="tight",
        )
        self._cif_variance_status.configure(text=f"Saved to {out_path}", text_color=_GREEN)

    # ── Results tab: search/filter ───────────────────────────────────────────

    def _filter_ptm_tv(self, *_args) -> None:
        self._filter_treeview(self._ptm_tv, self._ptm_tv_all_rows, self._results_ptm_search_var.get())

    def _filter_mut_tv(self, *_args) -> None:
        self._filter_treeview(self._mut_tv, self._mut_tv_all_rows, self._results_mut_search_var.get())

    def _focus_results_search(self, event=None):
        """Ctrl+F: focus the PTM search box, but only while the Results tab is showing."""
        if self._tabview.get() == "Results":
            self._ptm_search_entry.focus_set()
            return "break"

    def _load_results(self) -> None:
        import pandas as pd

        wide_path = self._output_dir / "ptm_mutation_proximity_db.tsv"
        long_path = self._output_dir / "ptm_mutation_proximity_long.tsv"

        if not wide_path.exists():
            self._results_status.configure(
                text=f"No output found in {self._output_dir.name}/",
                text_color=_RED,
            )
            self._clear_treeview_fully(self._ptm_tv, self._ptm_tv_all_rows)
            self._clear_treeview_fully(self._mut_tv, self._mut_tv_all_rows)
            self._ptm_tv_all_rows = []
            self._mut_tv_all_rows = []
            self._results_df_wide = None
            self._results_df_long = None
            self._refresh_viz_selector(pd.DataFrame(columns=["gene", "ptm_site", "UniProt"]))
            return

        try:
            df_wide = pd.read_csv(wide_path, sep="\t", encoding="utf-16",
                                   dtype=str, keep_default_na=False)
        except Exception as exc:
            self._results_status.configure(text=f"Error loading file: {exc}",
                                            text_color=_RED)
            return

        df_long = None
        if long_path.exists():
            try:
                df_long = pd.read_csv(long_path, sep="\t", encoding="utf-16",
                                       dtype=str, keep_default_na=False)
            except Exception:
                pass

        self._results_df_wide = df_wide
        self._results_df_long = df_long

        n_sites = len(df_wide)
        n_proteins = df_wide["UniProt"].nunique() if "UniProt" in df_wide.columns else "?"
        long_note = (" · long format available" if df_long is not None
                     else " · enable long format for per-mutation detail")
        self._results_status.configure(
            text=f"{n_sites} PTM sites · {n_proteins} proteins{long_note}",
            text_color="gray60",
        )
        self._populate_ptm_tv(df_wide)
        self._clear_treeview_fully(self._mut_tv, self._mut_tv_all_rows)
        self._mut_tv_all_rows = []
        self._refresh_viz_selector(df_wide)

    def _populate_ptm_tv(self, df) -> None:
        tv = self._ptm_tv
        self._clear_treeview_fully(tv, self._ptm_tv_all_rows)
        for i, (_, row) in enumerate(df.iterrows(), 1):
            try:
                near_pts = int(float(row.get("nearby_muts_total_patient_count", "") or "0"))
                far_pts  = int(float(row.get("distant_muts_total_patient_count", "") or "0"))
                total_pts: int | str = near_pts + far_pts
            except ValueError:
                total_pts = ""
            try:
                near_unique = int(float(row.get("unique_mutation_position_count_within_5_positions", "") or "0"))
                far_unique  = int(float(row.get("unique_mutation_position_count_more_than_5_positions", "") or "0"))
                total_muts: int | str = near_unique + far_unique
            except ValueError:
                total_muts = ""
            disrupt  = "Yes" if row.get("confirmed_disrupting_mutations", "").strip() else ""
            diseases = row.get("ptm_diseases", "")
            if len(diseases) > 50:
                diseases = diseases[:47] + "…"
            tv.insert("", "end", iid=str(i), values=(
                i,
                row.get("gene", ""),
                row.get("ptm_site", ""),
                row.get("ptm_type", ""),
                row.get("mutation_count_within_5_positions", ""),
                row.get("mutation_count_more_than_5_positions", ""),
                total_muts,
                total_pts,
                row.get("total_cosmic_missense_patients", ""),
                row.get("mutation_at_ptm_site", ""),
                disrupt,
                row.get("1433pred_binding_site", ""),
                row.get("1433_confirmed_site", ""),
                diseases,
            ), tags=("odd" if i % 2 else "even",))
        self._ptm_tv_all_rows = self._capture_tv_rows(tv)
        self._filter_ptm_tv()

    def _on_ptm_select(self, *_) -> None:
        sel = self._ptm_tv.selection()
        if not sel or self._results_df_wide is None:
            return
        row = self._results_df_wide.iloc[int(sel[0]) - 1]
        if self._results_df_long is not None:
            uid  = row.get("UniProt", "")
            site = row.get("ptm_site", "")
            mask = (
                (self._results_df_long.get("uniprot_id", "") == uid) &
                (self._results_df_long.get("ptm_position", "") == site)
            )
            self._populate_mut_tv_long(self._results_df_long[mask])
        else:
            self._populate_mut_tv_wide(row)

    def _visualize_selected_ptm(self) -> None:
        """Jump to the Visualization tab and render the lollipop plot for the selected PTM row."""
        sel = self._ptm_tv.selection()
        if not sel or self._results_df_wide is None:
            self._results_status.configure(
                text="Select a PTM site in the table first.", text_color=_YELLOW,
            )
            return
        row = self._results_df_wide.iloc[int(sel[0]) - 1]
        label = f"{row.get('gene', '?')}  {row.get('ptm_site', '?')}  ({row.get('UniProt', '?')})"

        self._tabview.set("Visualization")
        self._viz_search_var.set("")
        if label in self._viz_ptm_rows:
            self._viz_combo.set(label)
            self._generate_lollipop_plot()
        else:
            self._viz_status.configure(
                text=f"Could not find '{label}' in the Visualization selector.", text_color=_RED,
            )

    def _populate_mut_tv_long(self, df) -> None:
        tv = self._mut_tv
        self._clear_treeview_fully(tv, self._mut_tv_all_rows)
        for i, (_, r) in enumerate(df.iterrows(), 1):
            tv.insert("", "end", iid=str(i), values=(
                i,
                r.get("mutation", ""),
                r.get("sequence_distance", ""),
                r.get("distance_angstrom", ""),
                r.get("mut_aiupred_binding", ""),
                r.get("mut_is_binding", ""),
                r.get("mut_aiupred_general", ""),
                r.get("mut_is_disordered", ""),
                r.get("polyphen_class", ""),
                r.get("polyphen_score", ""),
                r.get("mutation_plddt", ""),
                r.get("pair_pae", ""),
                r.get("patient_count", ""),
                r.get("confirmed_disrupting_mutation", ""),
            ), tags=("odd" if i % 2 else "even",))
        self._mut_tv_all_rows = self._capture_tv_rows(tv)
        self._filter_mut_tv()

    def _populate_mut_tv_wide(self, row) -> None:
        import re as _re
        tv = self._mut_tv
        self._clear_treeview_fully(tv, self._mut_tv_all_rows)
        ptm_m = _re.search(r"(\d+)", str(row.get("ptm_site", "")))
        ptm_pos = int(ptm_m.group(1)) if ptm_m else None
        i = 0
        for col_key, _ in [
            ("mutations_within_5_positions",  "≤5 pos"),
            ("mutations_more_than_5_positions", ">5 pos"),
        ]:
            for entry in (row.get(col_key, "") or "").split(", "):
                entry = entry.strip()
                if not entry:
                    continue
                m = _MUT_ENTRY_RE.match(entry)
                if not m:
                    continue
                i += 1
                pp_code = m.group(2) or ""
                mut_m   = _re.search(r"\d+", m.group(1))
                mut_pos = int(mut_m.group()) if mut_m else None
                seq_d   = abs(mut_pos - ptm_pos) if (mut_pos is not None and ptm_pos is not None) else ""
                tv.insert("", "end", iid=str(i), values=(
                    i,
                    m.group(1),
                    seq_d,
                    m.group(4),
                    "", "",  # binding score / is_binding (unavailable in wide format)
                    "", "",  # disorder score / is_disordered (unavailable in wide format)
                    _PP_LABEL.get(pp_code, ""),
                    m.group(3) or "",
                    "",
                    m.group(5) or "",
                    "",
                    "",
                ), tags=("odd" if i % 2 else "even",))
        self._mut_tv_all_rows = self._capture_tv_rows(tv)
        self._filter_mut_tv()

    # ── Cache management ─────────────────────────────────────────────────────

    def _manage_cache(self):
        import shutil

        dialog = ctk.CTkToplevel(self)
        dialog.title("Manage Cache")
        dialog.resizable(False, False)
        dialog.grab_set()

        ctk.CTkLabel(
            dialog, text="Manage Cache",
            font=ctk.CTkFont(size=15, weight="bold"),
        ).grid(row=0, column=0, columnspan=4, padx=20, pady=(16, 2), sticky="w")
        ctk.CTkLabel(
            dialog,
            text="Clear cached API results to force re-fetching on the next run.",
            text_color=_GRAY,
        ).grid(row=1, column=0, columnspan=4, padx=20, pady=(0, 10), sticky="w")

        for col, txt in enumerate(["Step", "Cache", "Entries", ""]):
            ctk.CTkLabel(
                dialog, text=txt, font=ctk.CTkFont(weight="bold"), text_color=_GRAY,
            ).grid(row=2, column=col, padx=(20 if col == 0 else 8, 8), pady=(0, 4), sticky="w")

        count_labels: list = []
        last_step = None

        def _clear(path: Path, is_dir: bool, lbl):
            if is_dir:
                shutil.rmtree(path, ignore_errors=True)
            else:
                path.unlink(missing_ok=True)
            lbl.configure(text="empty")

        def _clear_all():
            for _, _, path, is_dir in _CACHE_ITEMS:
                if is_dir:
                    shutil.rmtree(path, ignore_errors=True)
                else:
                    path.unlink(missing_ok=True)
            for lbl in count_labels:
                lbl.configure(text="empty")

        for i, (step, name, path, is_dir) in enumerate(_CACHE_ITEMS):
            row = 3 + i
            step_text = step if step != last_step else ""
            last_step = step

            ctk.CTkLabel(dialog, text=step_text, text_color=_GRAY).grid(
                row=row, column=0, padx=(20, 8), pady=3, sticky="w")
            ctk.CTkLabel(dialog, text=name).grid(
                row=row, column=1, padx=8, pady=3, sticky="w")

            count_lbl = ctk.CTkLabel(
                dialog, text=_cache_entry_count(path, is_dir),
                text_color=_GRAY, width=90, anchor="e",
            )
            count_lbl.grid(row=row, column=2, padx=8, pady=3, sticky="e")
            count_labels.append(count_lbl)

            ctk.CTkButton(
                dialog, text="Clear", width=64,
                fg_color="gray30", hover_color=_RED,
                command=lambda p=path, d=is_dir, lbl=count_lbl: _clear(p, d, lbl),
            ).grid(row=row, column=3, padx=(0, 20), pady=3)

        bottom_row = 3 + len(_CACHE_ITEMS)
        ctk.CTkButton(
            dialog, text="Clear All Caches",
            fg_color="gray30", hover_color=_RED,
            command=_clear_all,
        ).grid(row=bottom_row, column=0, columnspan=4, padx=20, pady=(10, 16), sticky="ew")

        dialog.after(50, dialog.lift)

    # ── Output folder ────────────────────────────────────────────────────────

    def _open_output_folder(self):
        out = self._output_dir
        out.mkdir(parents=True, exist_ok=True)
        if platform.system() == "Windows":
            subprocess.run(["explorer", str(out)], check=False)
        elif platform.system() == "Darwin":
            subprocess.run(["open", str(out)], check=False)
        else:
            subprocess.run(["xdg-open", str(out)], check=False)


if __name__ == "__main__":
    app = App()
    app.mainloop()
