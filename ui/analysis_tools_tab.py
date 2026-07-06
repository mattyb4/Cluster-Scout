"""Analysis Tools tab: embedded Radius Sweep / CIF Variance plots, shared toggle.

`self._analysis_subtool_var` (the segmented-button toggle) is shared with
PipelineTabMixin's Pipeline-tab panel (see _build_analysis_tools_panel there)
— both use a `hasattr` guard so whichever one builds first creates the
StringVar and the other just reuses it. `_style_dark_figure` calls
`self._style_lollipop_axis`, a VisualizationTabMixin method.
"""
from __future__ import annotations

import customtkinter as ctk

from ui.common import _GREEN, _RED


class AnalysisToolsTabMixin:
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

    # ── Analysis Tools tab (Radius Sweep + CIF Variance, shared toggle) ─────

    def _build_analysis_tools_tab(self, tab) -> None:
        import matplotlib
        matplotlib.use("TkAgg")
        from matplotlib.figure import Figure
        from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

        if not hasattr(self, "_analysis_subtool_var"):
            self._analysis_subtool_var = ctk.StringVar(value="Radius Sweep")

        tab.grid_columnconfigure(0, weight=1)
        tab.grid_rowconfigure(1, weight=1)

        controls = ctk.CTkFrame(tab)
        controls.grid(row=0, column=0, padx=12, pady=(12, 6), sticky="ew")

        ctk.CTkSegmentedButton(
            controls, values=["Radius Sweep", "CIF Variance"],
            variable=self._analysis_subtool_var,
            command=lambda _v: self._show_active_analysis_plot(),
        ).pack(side="left", padx=(12, 12), pady=10)

        ctk.CTkButton(
            controls, text="Save PNG", width=100,
            fg_color="gray30", hover_color="gray40",
            command=self._save_active_analysis_plot,
        ).pack(side="left", padx=(0, 12), pady=10)

        self._analysis_tools_status = ctk.CTkLabel(
            controls, text="Run Radius Sweep or CIF Variance mode from the Pipeline tab to see results here.",
            text_color="gray60", font=ctk.CTkFont(size=11),
        )
        self._analysis_tools_status.pack(side="left", padx=(0, 12), pady=10)

        canvas_frame = ctk.CTkFrame(tab, fg_color="#2b2b2b")
        canvas_frame.grid(row=1, column=0, padx=12, pady=(0, 12), sticky="nsew")
        canvas_frame.grid_columnconfigure(0, weight=1)
        canvas_frame.grid_rowconfigure(0, weight=1)

        # Both figures live in the same canvas_frame; only one is gridded at a time
        # (toggled by the segmented button above) so switching is instant and each
        # tool keeps its own plot in memory without re-running anything.
        self._radius_sweep_fig = Figure(figsize=(14, 9), dpi=100, facecolor="#2b2b2b")
        self._radius_sweep_canvas = FigureCanvasTkAgg(self._radius_sweep_fig, master=canvas_frame)
        self._radius_sweep_canvas.get_tk_widget().grid(row=0, column=0, sticky="nsew")

        self._cif_variance_fig = Figure(figsize=(14, 8), dpi=100, facecolor="#2b2b2b")
        self._cif_variance_canvas = FigureCanvasTkAgg(self._cif_variance_fig, master=canvas_frame)
        self._cif_variance_canvas.get_tk_widget().grid(row=0, column=0, sticky="nsew")

        self._show_active_analysis_plot()

    def _show_active_analysis_plot(self) -> None:
        """Show whichever canvas matches the current segmented-button selection."""
        if self._analysis_subtool_var.get() == "Radius Sweep":
            self._cif_variance_canvas.get_tk_widget().grid_remove()
            self._radius_sweep_canvas.get_tk_widget().grid()
        else:
            self._radius_sweep_canvas.get_tk_widget().grid_remove()
            self._cif_variance_canvas.get_tk_widget().grid()

    def _save_active_analysis_plot(self) -> None:
        if self._analysis_subtool_var.get() == "Radius Sweep":
            self._save_radius_sweep_plot()
        else:
            self._save_cif_variance_plot()

    def _draw_radius_sweep_plot(self) -> None:
        from radius_sweep import build_sweep_figure

        self._radius_sweep_fig.clf()
        build_sweep_figure(self._radius_sweep_result, fig=self._radius_sweep_fig)
        self._style_dark_figure(self._radius_sweep_fig)
        self._radius_sweep_canvas.draw()
        self._show_active_analysis_plot()
        self._analysis_tools_status.configure(text="Sweep complete.", text_color=_GREEN)

    def _save_radius_sweep_plot(self) -> None:
        if not self._radius_sweep_fig.axes:
            self._analysis_tools_status.configure(text="Run a sweep before saving.", text_color=_RED)
            return
        out_dir = self._output_dir
        out_dir.mkdir(parents=True, exist_ok=True)
        out_path = out_dir / "radius_sweep_plot.png"
        self._radius_sweep_fig.savefig(
            out_path, dpi=200, facecolor=self._radius_sweep_fig.get_facecolor(), bbox_inches="tight",
        )
        self._analysis_tools_status.configure(text=f"Saved to {out_path}", text_color=_GREEN)

    def _draw_cif_variance_plot(self) -> None:
        from cif_variance import build_variance_figure

        self._cif_variance_fig.clf()
        build_variance_figure(self._cif_variance_result, fig=self._cif_variance_fig)
        self._style_dark_figure(self._cif_variance_fig)
        self._cif_variance_canvas.draw()
        self._show_active_analysis_plot()
        self._analysis_tools_status.configure(text="Analysis complete.", text_color=_GREEN)

    def _save_cif_variance_plot(self) -> None:
        if not self._cif_variance_fig.axes:
            self._analysis_tools_status.configure(text="Run an analysis before saving.", text_color=_RED)
            return
        out_dir = self._output_dir
        out_dir.mkdir(parents=True, exist_ok=True)
        out_path = out_dir / "cif_variance_plot.png"
        self._cif_variance_fig.savefig(
            out_path, dpi=200, facecolor=self._cif_variance_fig.get_facecolor(), bbox_inches="tight",
        )
        self._analysis_tools_status.configure(text=f"Saved to {out_path}", text_color=_GREEN)
