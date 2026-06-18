"""Desktop GUI for the Mutation Cluster Proximity pipeline.

Run with:  uv run app.py
"""
from __future__ import annotations

import json
import platform
import queue
import subprocess
import sys
import threading
import time
from pathlib import Path

import customtkinter as ctk

PROJECT_ROOT = Path(__file__).resolve().parent
SCRIPTS_DIR = PROJECT_ROOT / "scripts"
OUTPUT_DIR = PROJECT_ROOT / "Output"

_PTM_STEPS = [
    "Filter and merge PTMD + COSMIC data",
    "Download AlphaFold CIF models and PAE files",
    "Find nearby mutations and compute distances",
    "Annotate 14-3-3-Pred binding-site predictions",
    "Annotate mutations with PolyPhen-2 scores",
]
_CLUSTER_STEPS = [
    "Filter COSMIC hotspot mutations",
    "Download AlphaFold CIF models and PAE files",
    "Find mutation clusters in 3D space",
]

_REQUIRED_FILES: dict[str, Path] = {
    "PTMD": PROJECT_ROOT / "data" / "PTMD_disease_associated_ptms.tsv",
    "COSMIC": PROJECT_ROOT / "data" / "Cosmic_MutantCensus_v104_GRCh38.tsv",
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
        self.geometry("960x720")
        self.minsize(760, 560)

        self._queue: queue.Queue[tuple] = queue.Queue()
        self._running = False
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
        self.grid_rowconfigure(5, weight=1)

        # Title
        ctk.CTkLabel(
            self,
            text="Mutation Cluster Proximity Pipeline",
            font=ctk.CTkFont(size=22, weight="bold"),
        ).grid(row=0, column=0, padx=24, pady=(20, 4), sticky="w")

        # Data-file status bar
        self._file_frame = ctk.CTkFrame(self)
        self._file_frame.grid(row=1, column=0, padx=24, pady=4, sticky="ew")
        self._file_indicators: dict[str, ctk.CTkLabel] = {}

        ctk.CTkLabel(
            self._file_frame,
            text="Data files:",
            font=ctk.CTkFont(weight="bold"),
        ).pack(side="left", padx=(12, 8), pady=8)

        for name in _REQUIRED_FILES:
            lbl = ctk.CTkLabel(self._file_frame, text=f"{name} …")
            lbl.pack(side="left", padx=8, pady=8)
            self._file_indicators[name] = lbl

        # Mode selection
        mode_frame = ctk.CTkFrame(self)
        mode_frame.grid(row=2, column=0, padx=24, pady=4, sticky="ew")

        ctk.CTkLabel(
            mode_frame, text="Mode:", font=ctk.CTkFont(weight="bold")
        ).pack(side="left", padx=(12, 8), pady=10)

        self._mode = ctk.StringVar(value="ptm-proximity")
        ctk.CTkRadioButton(
            mode_frame,
            text="PTM Proximity",
            variable=self._mode,
            value="ptm-proximity",
            command=self._rebuild_step_rows,
        ).pack(side="left", padx=8, pady=10)
        ctk.CTkRadioButton(
            mode_frame,
            text="Mutation Clustering",
            variable=self._mode,
            value="mutation-clustering",
            command=self._rebuild_step_rows,
        ).pack(side="left", padx=8, pady=10)

        # Steps panel
        self._steps_outer = ctk.CTkFrame(self)
        self._steps_outer.grid(row=3, column=0, padx=24, pady=4, sticky="ew")
        self._steps_outer.grid_columnconfigure(1, weight=1)
        self._rebuild_step_rows()

        # Buttons
        btn_frame = ctk.CTkFrame(self, fg_color="transparent")
        btn_frame.grid(row=4, column=0, padx=24, pady=8, sticky="ew")

        self._run_btn = ctk.CTkButton(
            btn_frame,
            text="▶  Run Pipeline",
            command=self._start_pipeline,
            width=180,
            height=44,
            font=ctk.CTkFont(size=15, weight="bold"),
        )
        self._run_btn.pack(side="left", padx=(0, 12))

        self._open_btn = ctk.CTkButton(
            btn_frame,
            text="Open Output Folder",
            command=self._open_output_folder,
            width=180,
            height=44,
            state="disabled",
            fg_color="gray30",
            hover_color="gray40",
        )
        self._open_btn.pack(side="left")

        self._timer_label = ctk.CTkLabel(
            btn_frame,
            text="",
            text_color=_GRAY,
            font=ctk.CTkFont(size=13),
        )
        self._timer_label.pack(side="right", padx=12)

        # Log
        ctk.CTkLabel(
            self, text="Log", font=ctk.CTkFont(weight="bold")
        ).grid(row=5, column=0, padx=24, pady=(8, 0), sticky="w")

        self._log = ctk.CTkTextbox(
            self,
            font=ctk.CTkFont(family="Courier New", size=12),
            wrap="word",
            state="disabled",
        )
        self._log.grid(row=6, column=0, padx=24, pady=(4, 20), sticky="nsew")
        self.grid_rowconfigure(6, weight=1)

    def _rebuild_step_rows(self):
        for w in self._steps_outer.winfo_children():
            w.destroy()
        self._step_status_labels = []

        steps = _PTM_STEPS if self._mode.get() == "ptm-proximity" else _CLUSTER_STEPS

        ctk.CTkLabel(
            self._steps_outer,
            text="Steps",
            font=ctk.CTkFont(weight="bold"),
        ).grid(row=0, column=0, columnspan=3, padx=12, pady=(8, 2), sticky="w")

        for i, label in enumerate(steps, 1):
            ctk.CTkLabel(self._steps_outer, text=f"  {i}.", width=28).grid(
                row=i, column=0, padx=(12, 0), pady=5, sticky="w"
            )
            ctk.CTkLabel(self._steps_outer, text=label, anchor="w").grid(
                row=i, column=1, padx=6, pady=5, sticky="ew"
            )
            status = ctk.CTkLabel(
                self._steps_outer,
                text="●  Waiting",
                width=150,
                anchor="e",
                text_color=_GRAY,
            )
            status.grid(row=i, column=2, padx=12, pady=5, sticky="e")
            self._step_status_labels.append(status)

    # ── Timer ────────────────────────────────────────────────────────────────

    def _tick(self) -> None:
        if not self._running or self._pipeline_start is None:
            return
        elapsed = time.time() - self._pipeline_start

        text = f"Elapsed: {_fmt_time(elapsed)}"

        hist = self._historical_times
        if hist and len(hist) == self._total_steps:
            est_total = sum(hist)
            text += f"  |  Est. total: ~{_fmt_time(est_total)}"
        elif self._step_times:
            avg = sum(self._step_times) / len(self._step_times)
            est_total = avg * self._total_steps
            text += f"  |  Est. total: ~{_fmt_time(est_total)}"

        self._timer_label.configure(text=text)
        self.after(1000, self._tick)

    # ── File-status bar ──────────────────────────────────────────────────────

    def _refresh_file_status(self):
        for name, path in _REQUIRED_FILES.items():
            lbl = self._file_indicators[name]
            if path.exists():
                lbl.configure(text=f"✓  {name}", text_color=_GREEN)
            else:
                lbl.configure(text=f"✗  {name} missing", text_color=_RED)

    # ── Pipeline execution ───────────────────────────────────────────────────

    def _start_pipeline(self):
        if self._running:
            return

        self._running = True
        self._run_btn.configure(state="disabled")
        self._open_btn.configure(state="disabled")
        self._timer_label.configure(text="")

        self._log.configure(state="normal")
        self._log.delete("1.0", "end")
        self._log.configure(state="disabled")

        for lbl in self._step_status_labels:
            lbl.configure(text="●  Waiting", text_color=_GRAY)

        mode = self._mode.get()
        threading.Thread(target=self._run_pipeline, args=(mode,), daemon=True).start()

    def _run_pipeline(self, mode: str):
        python = sys.executable
        input_tsv = PROJECT_ROOT / "data" / "steps" / "PTMD_TCGA_hotspots_by_protein.tsv"
        models_dir = PROJECT_ROOT / "cif_models"

        cmds = [
            [python, str(SCRIPTS_DIR / "1_filter.py"), "--mode", mode],
            [
                python, str(SCRIPTS_DIR / "2_download_structures.py"),
                str(input_tsv),
                "--id_column", "uniprot_id",
                "--out_dir", str(models_dir),
                "--prefer", "cif",
                "--delay", "0.1",
                "--also_pae",
                "--logs_dir", str(OUTPUT_DIR / "logs"),
            ],
            [python, str(SCRIPTS_DIR / "3_find_nearby_mutations.py"), "--mode", mode],
        ]
        if mode == "ptm-proximity":
            cmds.append([python, str(SCRIPTS_DIR / "5_annotate_1433pred.py")])
            cmds.append([python, str(SCRIPTS_DIR / "6_annotate_polyphen.py")])

        steps = _PTM_STEPS if mode == "ptm-proximity" else _CLUSTER_STEPS
        run_type = _detect_run_type()
        self._q("pipeline_start", len(steps), mode, run_type)

        if run_type == "cold" and _load_runtimes(mode, "cold") is None:
            self._q("log", "Note: CIF files and API caches are not yet built.")
            self._q("log", "Step 2 (structure download) and Step 5 (14-3-3 predictions)")
            self._q("log", "may each take 20+ minutes on a first run. Subsequent runs")
            self._q("log", "will be much faster once these are cached.\n")

        all_ok = True
        step_times: list[float] = []
        for i, (label, cmd) in enumerate(zip(steps, cmds)):
            self._q("status", i, "▶  Running…", _BLUE)
            self._q("log", f"\n{'─' * 56}")
            self._q("log", f"Step {i + 1}/{len(steps)}: {label}")
            self._q("log", f"{'─' * 56}")
            self._q("step_start")

            t0 = time.time()
            ok = self._stream_cmd(cmd)
            elapsed = time.time() - t0

            if ok:
                step_times.append(elapsed)
                self._q("status", i, f"✓  {_fmt_time(elapsed)}", _GREEN)
                self._q("log", f"\n✓  Step {i + 1} completed in {_fmt_time(elapsed)}")
                self._q("step_complete", elapsed)
            else:
                self._q("status", i, "✗  Failed", _RED)
                self._q("log", f"\n✗  Step {i + 1} FAILED after {_fmt_time(elapsed)} — pipeline stopped")
                all_ok = False
                break

        if all_ok:
            self._q("save_runtimes", mode, run_type, step_times)
            self._q("log", f"\n{'═' * 56}")
            self._q("log", "Pipeline complete! Output saved to the Output/ folder.")
            self._q("log", f"{'═' * 56}")
            self._q("enable_open")

        self._q("finished")

    def _stream_cmd(self, cmd: list[str]) -> bool:
        proc = subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            encoding="utf-8",
            errors="replace",
        )
        for line in proc.stdout:
            self._q("log", line.rstrip())
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
                elif kind == "status":
                    _, idx, text, color = msg
                    if 0 <= idx < len(self._step_status_labels):
                        self._step_status_labels[idx].configure(text=text, text_color=color)
                elif kind == "pipeline_start":
                    _, total, mode, run_type = msg
                    self._pipeline_start = time.time()
                    self._step_start = None
                    self._total_steps = total
                    self._steps_done = 0
                    self._step_times = []
                    self._historical_times = _load_runtimes(mode, run_type)
                    self._tick()
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
                elif kind == "enable_open":
                    self._open_btn.configure(state="normal", fg_color=_BLUE, hover_color="#2563eb")
                elif kind == "finished":
                    self._running = False
                    self._run_btn.configure(state="normal")
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

    # ── Output folder ────────────────────────────────────────────────────────

    def _open_output_folder(self):
        OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
        if platform.system() == "Windows":
            subprocess.run(["explorer", str(OUTPUT_DIR)], check=False)
        elif platform.system() == "Darwin":
            subprocess.run(["open", str(OUTPUT_DIR)], check=False)
        else:
            subprocess.run(["xdg-open", str(OUTPUT_DIR)], check=False)


if __name__ == "__main__":
    app = App()
    app.mainloop()
