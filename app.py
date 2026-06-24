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
        self.title("Mutation Cluster Proximity Pipeline")
        self.geometry("960x820")
        self.minsize(760, 560)

        self._queue: queue.Queue[tuple] = queue.Queue()
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
        self._tabview = ctk.CTkTabview(self)
        self._tabview.grid(row=0, column=0, padx=12, pady=12, sticky="nsew")
        pipeline_tab = self._tabview.add("Pipeline")
        help_tab = self._tabview.add("Help / Documentation")

        # ── Pipeline tab ──
        pipeline_tab.grid_columnconfigure(0, weight=1)
        pipeline_tab.grid_rowconfigure(0, weight=1)

        scroll = ctk.CTkScrollableFrame(pipeline_tab, fg_color="transparent")
        scroll.grid(row=0, column=0, sticky="nsew")
        scroll.grid_columnconfigure(0, weight=1)

        p = scroll  # all pipeline widgets go in the scrollable frame

        # Title
        ctk.CTkLabel(
            p,
            text="Mutation Cluster Proximity Pipeline",
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

        ctk.CTkLabel(settings_frame, text="Distance cutoff (A):").pack(side="left", padx=(8, 4), pady=8)
        self._cutoff_var = ctk.StringVar(value="10.0")
        ctk.CTkEntry(
            settings_frame, textvariable=self._cutoff_var, width=60,
        ).pack(side="left", padx=(0, 16), pady=8)

        ctk.CTkLabel(settings_frame, text="Min COSMIC samples:").pack(side="left", padx=(8, 4), pady=8)
        self._min_samples_var = ctk.StringVar(value="3")
        ctk.CTkEntry(
            settings_frame, textvariable=self._min_samples_var, width=60,
        ).pack(side="left", padx=(0, 16), pady=8)

        ctk.CTkLabel(settings_frame, text="Min pLDDT:").pack(side="left", padx=(8, 4), pady=8)
        self._min_plddt_var = ctk.StringVar(value="70")
        ctk.CTkEntry(
            settings_frame, textvariable=self._min_plddt_var, width=60,
        ).pack(side="left", padx=(0, 16), pady=8)

        ctk.CTkLabel(settings_frame, text="Max PAE:").pack(side="left", padx=(8, 4), pady=8)
        self._max_pae_var = ctk.StringVar(value="5.0")
        ctk.CTkEntry(
            settings_frame, textvariable=self._max_pae_var, width=60,
        ).pack(side="left", padx=(0, 12), pady=8)

        # Steps panel
        self._steps_outer = ctk.CTkFrame(p)
        self._steps_outer.grid(row=5, column=0, padx=24, pady=4, sticky="ew")
        self._steps_outer.grid_columnconfigure(1, weight=1)
        self._rebuild_step_rows()

        # Buttons
        btn_frame = ctk.CTkFrame(p, fg_color="transparent")
        btn_frame.grid(row=6, column=0, padx=24, pady=8, sticky="ew")

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

        self._timer_label = ctk.CTkLabel(
            btn_frame,
            text="",
            text_color=_GRAY,
            font=ctk.CTkFont(size=13),
        )
        self._timer_label.pack(side="right", padx=12)

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
        self._log_toggle.grid(row=7, column=0, padx=24, pady=(8, 0), sticky="w")

        self._log = ctk.CTkTextbox(
            p,
            font=ctk.CTkFont(family="Courier New", size=12),
            wrap="word",
            state="disabled",
            height=300,
        )

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

    # ── Timer ────────────────────────────────────────────────────────────────

    def _tick(self) -> None:
        if not self._running or self._pipeline_start is None:
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
        """Show or hide the raw log output panel."""
        if self._log_visible:
            self._log.grid_remove()
            self._log_toggle.configure(text="Show Details")
            self._log_visible = False
        else:
            self._log.grid(row=8, column=0, padx=24, pady=(4, 20), sticky="nsew")
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
        "mutation_cluster_db.tsv",
    ]

    @property
    def _output_files(self) -> list[Path]:
        return [self._output_dir / name for name in self._DEFAULT_OUTPUT_FILES]

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

        mode = self._mode.get()
        if mode == "single-protein":
            cif = getattr(self, "_single_cif_var", None)
            if not cif or not cif.get().strip():
                from tkinter import messagebox
                messagebox.showwarning("Missing input", "Please select a CIF file first.")
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

            self._q("log", f"Step 4: {cached_1433}/{n_proteins} 14-3-3 predictions cached, "
                    f"{cached_pp} PolyPhen pairs cached, {cached_kin} kinase windows cached")

            step4_est = (uncached_1433 * self._TIME_PER_1433_FETCH
                         + max(0, n_proteins * 2 - cached_pp) * self._TIME_PER_PP_FETCH
                         + max(0, n_proteins - cached_kin) * self._TIME_PER_KINASE_PREDICT
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

        python = [sys.executable, "-u"]
        input_tsv = PROJECT_ROOT / "data" / "steps" / "PTMD_TCGA_hotspots_by_protein.tsv"
        models_dir = PROJECT_ROOT / "cif_models"

        cutoff = self._cutoff_var.get().strip() or "10.0"
        min_samples = self._min_samples_var.get().strip() or "3"
        min_plddt = self._min_plddt_var.get().strip() or "0"
        max_pae = self._max_pae_var.get().strip() or "0"

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
             "--min-plddt", min_plddt, "--max-pae", max_pae],
        ]
        if mode == "ptm-proximity":
            cmds.append([*python, str(SCRIPTS_DIR / "4_annotate.py"),
                         "--output-dir", str(self._output_dir)])

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
                elif kind == "enable_open":
                    self._open_btn.configure(state="normal", fg_color=_BLUE, hover_color="#2563eb")
                elif kind == "finished":
                    was_cancelled = self._stop_requested
                    self._running = False
                    self._stop_requested = False
                    self._suspended = False
                    self._current_proc = None
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
