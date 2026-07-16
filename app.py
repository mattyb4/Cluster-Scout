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
import tkinter as tk

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


def _patch_customtkinter_textbox_scroll_callback() -> None:
    """CTkTextbox schedules a recurring self.after() poll to auto-show/hide its
    scrollbars, guarded by a winfo_exists() check -- but only *before*
    rescheduling the *next* call, not before touching the widget on the
    *current* one. If the widget was destroyed since the poll was queued (e.g.
    a Pipeline-tab mode switch tearing down and rebuilding its steps frame),
    the current call's own xview()/yview() calls raise TclError first, which
    aborts the function before it reaches that guard -- so the loop still
    stops correctly, but the exception escapes uncaught in the meantime.
    """
    try:
        from customtkinter.windows.widgets import ctk_textbox
    except ImportError:
        return

    original = ctk_textbox.CTkTextbox._check_if_scrollbars_needed

    def _safe_check_if_scrollbars_needed(self, event=None, continue_loop: bool = False):
        try:
            original(self, event, continue_loop=continue_loop)
        except tk.TclError:
            return

    ctk_textbox.CTkTextbox._check_if_scrollbars_needed = _safe_check_if_scrollbars_needed


def _patch_customtkinter_scaling_tracker() -> None:
    """ScalingTracker polls every registered window every 100ms, forever, to
    detect OS-level per-monitor DPI changes -- a real feature on Windows, but
    ScalingTracker.get_window_dpi_scaling() hard-codes a return of 1 on macOS
    and Linux ("scaling works automatically on macOS" / "not implemented" on
    Linux), so on those platforms the loop can never detect a change: it's
    pure overhead for the app's entire lifetime, and its winfo_exists() check
    on each window isn't enough to stop it throwing TclError once the app
    starts closing (the interpreter can be torn down between the check and
    the following calls in the same pass). Since it provably does nothing
    outside Windows, skip it there entirely instead of just guarding it --
    that fixes the wasted CPU and removes the crash source in one place.
    """
    if sys.platform == "win32":
        return

    try:
        from customtkinter.windows.widgets.scaling import scaling_tracker
    except ImportError:
        return

    @classmethod
    def _noop_check_dpi_scaling(cls):
        pass

    scaling_tracker.ScalingTracker.check_dpi_scaling = _noop_check_dpi_scaling


def _patch_customtkinter_appearance_mode_tracker() -> None:
    """AppearanceModeTracker polls every 30ms, forever -- three times more
    often than the ScalingTracker loop above -- to detect a live OS-level
    light/dark mode change via the darkdetect package. This app hardcodes
    ctk.set_appearance_mode("dark") once at import and never calls "system"
    or offers any way to change it at runtime, so the live-detection this
    loop exists for never applies here: every widget already gets its
    correct color at creation time from the mode set once up front. Skip
    the loop entirely on every platform (not just macOS) since it's dead
    weight for this app regardless of OS -- and at 30ms it's the more
    likely of the two loops to be felt as UI lag, since it fires while
    genuinely idle far more often than the 100ms scaling loop did.
    """
    try:
        from customtkinter.windows.widgets.appearance_mode import appearance_mode_tracker
    except ImportError:
        return

    @classmethod
    def _noop_update(cls):
        pass

    appearance_mode_tracker.AppearanceModeTracker.update = _noop_update


_patch_customtkinter_textbox_scroll_callback()
_patch_customtkinter_scaling_tracker()
_patch_customtkinter_appearance_mode_tracker()


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
