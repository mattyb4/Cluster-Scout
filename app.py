"""Desktop GUI for the Mutation Cluster Proximity pipeline.

Run with:  uv run app.py

The App class itself just owns window setup and the top-level tab assembly —
each tab's widgets/logic live in their own ui/*.py mixin module (see ui/common.py
for the shared constants/helpers those mixins import).
"""
from __future__ import annotations

import queue
import subprocess

import customtkinter as ctk

from ui.common import PROJECT_ROOT
from ui.pipeline_panels import PipelineTabMixin
from ui.pipeline_runner import PipelineRunnerMixin
from ui.results_tab import ResultsTabMixin
from ui.visualization_tab import VisualizationTabMixin
from ui.analysis_tools_tab import AnalysisToolsTabMixin
from ui.help_tab import HelpTabMixin

ctk.set_appearance_mode("dark")
ctk.set_default_color_theme("blue")


class App(
    PipelineTabMixin,
    PipelineRunnerMixin,
    ResultsTabMixin,
    VisualizationTabMixin,
    AnalysisToolsTabMixin,
    HelpTabMixin,
    ctk.CTk,
):
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
        analysis_tools_tab = self._tabview.add("Analysis Tools")
        help_tab = self._tabview.add("Help / Documentation")

        self._build_pipeline_tab(pipeline_tab)
        self._build_results_tab(results_tab)
        self._build_viz_tab(viz_tab)
        self._build_analysis_tools_tab(analysis_tools_tab)
        self._build_help_tab(help_tab)

    def _on_tab_change(self) -> None:
        tab = self._tabview.get()
        if tab == "Results":
            self._load_results()
        elif tab == "Visualization" and self._results_df_wide is None:
            self._load_results()


if __name__ == "__main__":
    app = App()
    app.mainloop()
