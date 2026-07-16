"""Desktop GUI for the Mutation Cluster Proximity pipeline.

Run with:  uv run app.py

The App class itself just owns window setup and the top-level tab assembly —
each tab's widgets/logic live in their own ui/*.py mixin module (see ui/common.py
for the shared constants/helpers those mixins import).
"""
from __future__ import annotations

import queue
import subprocess
import sys

import customtkinter as ctk
from PIL import Image, ImageTk

from ui.common import PROJECT_ROOT, MIN_UI_SCALE, MAX_UI_SCALE, UI_SCALE_STEP, _BLUE
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
        self.title("Cluster-Scout (Beta)")
        self.geometry("1100x820")
        self.minsize(900, 560)
        ico = PROJECT_ROOT / "cluster_scout.ico"
        if ico.exists():
            if sys.platform == "win32":
                self.iconbitmap(str(ico))
            else:
                # .ico isn't a valid bitmap format for iconbitmap() outside
                # Windows (raises TclError on macOS/Linux); iconphoto() works
                # everywhere, and Pillow can read .ico regardless of platform.
                self._icon_image = ImageTk.PhotoImage(Image.open(ico))
                self.iconphoto(True, self._icon_image)

        self._queue: queue.Queue[tuple] = queue.Queue()
        self._radius_sweep_result = None
        self._cif_variance_result = None
        self._running = False
        self._at_run_active = False
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
        self._ui_scale = 1.0
        self._zoom_indicator_hide_job: str | None = None

        self._build_ui()
        self._refresh_file_status()
        self._poll_queue()
        self._bind_zoom_shortcuts()

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

        # Transient zoom-percentage toast (see _show_zoom_indicator) — a sibling
        # of the tabview, floated on top via place() so it stays visible
        # regardless of which tab is active.
        self._zoom_indicator = ctk.CTkLabel(
            self, font=ctk.CTkFont(size=14, weight="bold"),
            fg_color="#242424", text_color=_BLUE, corner_radius=8,
            padx=16, pady=8,
        )

    def _on_tab_change(self) -> None:
        tab = self._tabview.get()
        if tab == "Results":
            self._load_results()
        elif tab == "Visualization" and self._results_df_wide is None:
            self._load_results()

    # ── Ctrl+scroll zoom ─────────────────────────────────────────────────────

    def _bind_zoom_shortcuts(self) -> None:
        """App-wide Ctrl+scroll zoom (Windows/macOS wheel + X11 button-4/5)."""
        self.bind_all("<Control-MouseWheel>", self._on_ctrl_scroll_zoom)
        self.bind_all("<Control-Button-4>", self._on_ctrl_scroll_zoom)
        self.bind_all("<Control-Button-5>", self._on_ctrl_scroll_zoom)

    def _on_ctrl_scroll_zoom(self, event) -> str:
        """Nudge the app-wide CTk widget scaling up/down a step and re-clamp.

        Bound both at bind_all (app-wide) and directly on the Results-tab
        treeviews (see ResultsTabMixin) — Tk's unmodified `<MouseWheel>`
        class-binding on Treeview still matches Ctrl+wheel events (extra
        modifiers don't block a pattern that doesn't specify them), so
        without the widget-level override, scrolling those tables would
        also fire alongside the zoom.
        """
        if getattr(event, "num", None) == 4:
            direction = 1
        elif getattr(event, "num", None) == 5:
            direction = -1
        else:
            direction = 1 if event.delta > 0 else -1

        new_scale = round(max(MIN_UI_SCALE, min(MAX_UI_SCALE, self._ui_scale + direction * UI_SCALE_STEP)), 2)
        if new_scale != self._ui_scale:
            self._ui_scale = new_scale
            ctk.set_widget_scaling(new_scale)
        self._show_zoom_indicator()
        return "break"

    _ZOOM_INDICATOR_HIDE_MS = 1200

    def _show_zoom_indicator(self) -> None:
        """Flash the current zoom % in a corner toast, browser-style, then
        auto-hide after a short delay — rescheduled on every scroll so it
        stays up while the user keeps zooming.
        """
        self._zoom_indicator.configure(text=f"{round(self._ui_scale * 100)}%")
        self._zoom_indicator.place(relx=0.98, rely=0.02, anchor="ne")
        self._zoom_indicator.lift()
        if self._zoom_indicator_hide_job is not None:
            self.after_cancel(self._zoom_indicator_hide_job)
        self._zoom_indicator_hide_job = self.after(self._ZOOM_INDICATOR_HIDE_MS, self._hide_zoom_indicator)

    def _hide_zoom_indicator(self) -> None:
        self._zoom_indicator.place_forget()
        self._zoom_indicator_hide_job = None


if __name__ == "__main__":
    app = App()
    app.mainloop()
