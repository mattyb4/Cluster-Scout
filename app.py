"""Desktop GUI for the Mutation Cluster Proximity pipeline.

Run with:  uv run app.py

The App class itself just owns window setup and the top-level tab assembly —
each tab's widgets/logic live in their own ui/*.py mixin module (see ui/common.py
for the shared constants/helpers those mixins import).
"""
from __future__ import annotations

import os
import queue
import subprocess
import sys
from pathlib import Path


def _fix_tcl_tk_library_paths() -> None:
    """uv-managed (python-build-standalone) macOS Pythons can fail to locate
    their bundled Tcl/Tk init.tcl when run from inside a venv (TclError:
    "Can't find a usable init.tcl"), since the interpreter's relocation logic
    doesn't account for a venv's symlinked layout. sys.base_prefix still
    points at the real base install, so pointing TCL_LIBRARY/TK_LIBRARY there
    sidesteps the bug. Must run before `import tkinter`, which loads the
    _tkinter C extension and triggers the lookup. No-op if TCL_LIBRARY is
    already set, or the base install doesn't have this lib/tclX.Y layout.
    """
    if sys.platform != "darwin" or os.environ.get("TCL_LIBRARY"):
        return
    lib_dir = Path(sys.base_prefix) / "lib"
    tcl_dirs = sorted(lib_dir.glob("tcl8.*"))
    tk_dirs = sorted(lib_dir.glob("tk8.*"))
    if tcl_dirs and (tcl_dirs[-1] / "init.tcl").exists():
        os.environ["TCL_LIBRARY"] = str(tcl_dirs[-1])
    if tk_dirs:
        os.environ["TK_LIBRARY"] = str(tk_dirs[-1])


_fix_tcl_tk_library_paths()

import tkinter as tk  # noqa: E402

import customtkinter as ctk  # noqa: E402
from PIL import Image, ImageTk  # noqa: E402

from ui.common import PROJECT_ROOT, MIN_UI_SCALE, MAX_UI_SCALE, UI_SCALE_STEP, _BLUE  # noqa: E402
from ui.pipeline_panels import PipelineTabMixin  # noqa: E402
from ui.pipeline_runner import PipelineRunnerMixin  # noqa: E402
from ui.results_tab import ResultsTabMixin  # noqa: E402
from ui.visualization_tab import VisualizationTabMixin  # noqa: E402
from ui.analysis_tools_tab import AnalysisToolsTabMixin  # noqa: E402
from ui.help_tab import HelpTabMixin  # noqa: E402

ctk.set_appearance_mode("dark")
ctk.set_default_color_theme("blue")


def _patch_customtkinter_textbox_scroll_callback() -> None:
    """CTkTextbox's recurring self.after() poll for auto-show/hide scrollbars
    only checks winfo_exists() before scheduling the *next* call, not before
    touching the widget on the *current* one -- so if the widget was destroyed
    since the poll was queued, its own xview()/yview() calls raise TclError
    uncaught before reaching that guard.
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
    """ScalingTracker polls every window every 100ms forever to detect OS-level
    per-monitor DPI changes -- a real feature on Windows, but
    get_window_dpi_scaling() hard-codes a return of 1 on macOS/Linux, so the
    loop can never detect a change there: pure overhead, and its
    winfo_exists() check isn't enough to stop TclError once the app starts
    closing. Since it does nothing outside Windows, skip it there entirely.
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
    """AppearanceModeTracker polls every 30ms forever to detect a live OS-level
    light/dark mode change via darkdetect. This app hardcodes
    ctk.set_appearance_mode("dark") once at import and never offers a way to
    change it at runtime, so the loop this exists for never applies -- skip
    it on every platform.
    """
    try:
        from customtkinter.windows.widgets.appearance_mode import appearance_mode_tracker
    except ImportError:
        return

    @classmethod
    def _noop_update(cls):
        pass

    appearance_mode_tracker.AppearanceModeTracker.update = _noop_update


def _patch_customtkinter_tabview_switch() -> None:
    """CTkTabview switches tabs by grid_forget()-ing the outgoing frame and
    grid()-ing the incoming one -- on Tk's macOS Aqua backend each widget is a
    real Cocoa NSView, so unmapping/mapping a tab's whole subtree is a genuine
    per-widget Cocoa operation, causing a 50-211ms main-thread stall per switch.

    Once a tab has been shown at least once, keeping it permanently gridded
    and switching only which one is on top via tkraise() avoids that cost:
    grid() on an already-managed widget is a cheap no-op, and tkraise() only
    reorders stacking. Safe here since no code relies on hidden tabs being
    actually unmapped, and tabs are only ever added at startup.
    """
    try:
        from customtkinter.windows.widgets import ctk_tabview
    except ImportError:
        return

    def _fast_segmented_button_callback(self, selected_name):
        self._current_name = selected_name
        self._set_grid_current_tab()
        self._tab_dict[self._current_name].tkraise()
        if self._command is not None:
            self._command()

    def _fast_set(self, name: str):
        if name in self._tab_dict:
            self._current_name = name
            self._segmented_button.set(name)
            self._set_grid_current_tab()
            self._tab_dict[name].tkraise()
        else:
            raise ValueError(f"CTkTabview has no tab named '{name}'")

    ctk_tabview.CTkTabview._segmented_button_callback = _fast_segmented_button_callback
    ctk_tabview.CTkTabview.set = _fast_set


_patch_customtkinter_textbox_scroll_callback()
_patch_customtkinter_scaling_tracker()
_patch_customtkinter_appearance_mode_tracker()
_patch_customtkinter_tabview_switch()


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

        # Transient zoom-percentage toast, floated over the tabview via place()
        # so it stays visible regardless of which tab is active
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

        Also bound directly on the Results-tab treeviews (ResultsTabMixin):
        Tk's unmodified `<MouseWheel>` class-binding on Treeview still matches
        Ctrl+wheel, so without that override those tables would scroll too.
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
