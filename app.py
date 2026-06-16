"""Desktop GUI for the Mutation Cluster Proximity pipeline.

Run with:  uv run app.py
"""
from __future__ import annotations

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
    "Merge HTP/LTP scores into proximity database",
    "Annotate 14-3-3-Pred binding-site predictions",
]
_CLUSTER_STEPS = [
    "Filter COSMIC hotspot mutations",
    "Download AlphaFold CIF models and PAE files",
    "Find mutation clusters in 3D space",
]

_REQUIRED_FILES: dict[str, Path] = {
    "PTMD": PROJECT_ROOT / "data" / "PTMD_disease_associated_ptms.tsv",
    "COSMIC": PROJECT_ROOT / "data" / "Cosmic_MutantCensus_v104_GRCh38.tsv",
    "HTP/LTP": PROJECT_ROOT / "data" / "htp_ltp_scores.tsv",
}

ctk.set_appearance_mode("dark")
ctk.set_default_color_theme("blue")

_GRAY = "gray"
_BLUE = "#3a86ff"
_GREEN = "#2ecc71"
_RED = "#e74c3c"
_YELLOW = "#f1c40f"


class App(ctk.CTk):
    def __init__(self):
        super().__init__()
        self.title("Mutation Cluster Proximity Pipeline")
        self.geometry("960x720")
        self.minsize(760, 560)

        self._queue: queue.Queue[tuple] = queue.Queue()
        self._running = False
        self._step_status_labels: list[ctk.CTkLabel] = []

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
            cmds.append([python, str(SCRIPTS_DIR / "4_merge_htp_ltp.py")])
            cmds.append([python, str(SCRIPTS_DIR / "5_annotate_1433pred.py")])

        steps = _PTM_STEPS if mode == "ptm-proximity" else _CLUSTER_STEPS

        all_ok = True
        for i, (label, cmd) in enumerate(zip(steps, cmds)):
            self._q("status", i, "▶  Running…", _BLUE)
            self._q("log", f"\n{'─' * 56}")
            self._q("log", f"Step {i + 1}/{len(steps)}: {label}")
            self._q("log", f"{'─' * 56}")

            t0 = time.time()
            ok = self._stream_cmd(cmd)
            elapsed = time.time() - t0

            if ok:
                self._q("status", i, f"✓  {elapsed:.0f}s", _GREEN)
                self._q("log", f"\n✓  Step {i + 1} completed in {elapsed:.1f}s")
            else:
                self._q("status", i, "✗  Failed", _RED)
                self._q("log", f"\n✗  Step {i + 1} FAILED after {elapsed:.1f}s — pipeline stopped")
                all_ok = False
                break

        if all_ok:
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
                elif kind == "enable_open":
                    self._open_btn.configure(state="normal", fg_color=_BLUE, hover_color="#2563eb")
                elif kind == "finished":
                    self._running = False
                    self._run_btn.configure(state="normal")
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
