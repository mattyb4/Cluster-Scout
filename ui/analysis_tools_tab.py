"""Analysis Tools tab: self-contained Radius Sweep / CIF Variance UI —
parameters, a Run trigger, status/progress/log feedback, and the result plots.

Shares the app's single background-run slot with the Pipeline tab
(`self._running`); `_poll_queue` (pipeline_runner.py) routes queue messages
here instead of the Pipeline tab's widgets whenever `self._at_run_active` is
set. `_style_themed_figure` calls `self._style_lollipop_axis` (VisualizationTabMixin).
"""
from __future__ import annotations

import queue
import threading
from pathlib import Path
from tkinter import filedialog

import customtkinter as ctk

from ui import theme
from ui.common import (
    _BTN2,
    _BTN2_HOVER,
    _CHIP,
    _CHIP_HOVER,
    _CIF_VARIANCE_HELP,
    _GRAY,
    _GREEN,
    _PLOT_BG,
    _RADIUS_SWEEP_HELP,
    _RED,
    _YELLOW,
    OUTPUT_DIR,
    add_resize_grip,
    help_icon,
    isolate_textbox_scroll,
    section_header,
)

_SIDE_W = 340  # px; fixed width of the parameters panel


class AnalysisToolsTabMixin:
    def _style_themed_figure(self, fig) -> None:
        """Apply the current light/dark theme to every axes of a Figure built by an external script."""
        fig.patch.set_facecolor(theme.c("plot_bg"))
        for ax in fig.axes:
            self._style_lollipop_axis(ax)
            if ax.get_title():
                ax.title.set_color(theme.c("plot_fg"))
            legend = ax.get_legend()
            if legend is not None:
                legend.get_frame().set_facecolor(theme.c("plot_legend_bg"))
                legend.get_frame().set_edgecolor(theme.c("plot_legend_edge"))
                for text in legend.get_texts():
                    text.set_color(theme.c("plot_fg"))

    def _retheme_analysis(self) -> None:
        """Redraw any finished Radius Sweep / CIF Variance plot in the new theme."""
        for fig, canvas, result, draw in (
            (self._radius_sweep_fig, self._radius_sweep_canvas,
             self._radius_sweep_result, self._draw_radius_sweep_plot),
            (self._cif_variance_fig, self._cif_variance_canvas,
             self._cif_variance_result, self._draw_cif_variance_plot),
        ):
            fig.patch.set_facecolor(theme.c("plot_bg"))
            if fig.axes and result is not None:
                status = self._at_status_label.cget("text"), self._at_status_label.cget("text_color")
                draw()
                self._at_status_label.configure(text=status[0], text_color=status[1])
            else:
                if fig.axes:
                    self._style_themed_figure(fig)
                canvas.draw()

    # ── Analysis Tools tab (Radius Sweep + CIF Variance, shared toggle) ─────

    def _build_analysis_tools_tab(self, tab) -> None:
        import matplotlib
        matplotlib.use("TkAgg")
        from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
        from matplotlib.figure import Figure

        self._analysis_subtool_var = ctk.StringVar(value="Radius Sweep")

        # Parameters panel on the left (fixed width, Run pinned at its bottom), a divider,
        # and the plot workspace filling everything else.
        tab.grid_columnconfigure(2, weight=1)
        tab.grid_rowconfigure(0, weight=1)

        side = ctk.CTkFrame(tab, fg_color="transparent", corner_radius=0, border_width=0, width=_SIDE_W)
        side.grid(row=0, column=0, sticky="ns", padx=(14, 10), pady=12)
        side.grid_propagate(False)
        side.grid_columnconfigure(0, weight=1)
        side.grid_rowconfigure(2, weight=1)

        ctk.CTkSegmentedButton(
            side, values=["Radius Sweep", "CIF Variance"],
            variable=self._analysis_subtool_var,
            command=lambda _v: (self._rebuild_analysis_tool_fields(), self._show_active_analysis_plot()),
        ).grid(row=0, column=0, sticky="ew")
        section_header(side, "Parameters", row=1, pady=(14, 2))

        # Explicit fg_color (not "transparent"): children of a transparent CTkScrollableFrame
        # capture a static bg color at construction and show stale corners after a light/dark switch.
        params_scroll = ctk.CTkScrollableFrame(side, fg_color=theme.pair("bg"))
        params_scroll.grid(row=2, column=0, sticky="nsew")
        params_scroll.grid_columnconfigure(0, weight=1)
        # Cleared and rebuilt by _rebuild_analysis_tool_fields on Radius Sweep <-> CIF Variance toggle
        self._at_params_frame = ctk.CTkFrame(params_scroll, fg_color="transparent")
        self._at_params_frame.grid(row=0, column=0, sticky="ew")
        self._at_params_frame.grid_columnconfigure(0, weight=1)

        # Persistent, built once -- unlike the params frame, this is never torn down on toggle
        run_area = ctk.CTkFrame(side, fg_color="transparent")
        run_area.grid(row=3, column=0, sticky="ew", pady=(8, 0))
        run_area.grid_columnconfigure(0, weight=1)
        self._at_progress_bar = ctk.CTkProgressBar(run_area, height=10)
        self._at_progress_bar.set(0)
        self._at_progress_bar.grid(row=0, column=0, sticky="ew", pady=(0, 8))
        self._at_progress_bar.grid_remove()
        self._at_run_btn = ctk.CTkButton(
            run_area, text="▶  Run", height=36,
            font=ctk.CTkFont(size=14, weight="bold"),
            command=lambda: self._start_analysis_tool_run(
                "radius-sweep" if self._analysis_subtool_var.get() == "Radius Sweep" else "cif-variance"
            ),
        )
        self._at_run_btn.grid(row=1, column=0, sticky="ew")

        ctk.CTkFrame(
            tab, width=1, corner_radius=0, border_width=0, fg_color=theme.pair("line"),
        ).grid(row=0, column=1, sticky="ns", pady=12)

        main = ctk.CTkFrame(tab, fg_color="transparent", corner_radius=0, border_width=0)
        main.grid(row=0, column=2, sticky="nsew", padx=14, pady=12)
        main.grid_columnconfigure(0, weight=1)
        main.grid_rowconfigure(1, weight=1)

        # Results heading: title, live status, hairline, Save PNG
        head = ctk.CTkFrame(main, fg_color="transparent")
        head.grid(row=0, column=0, sticky="ew", pady=(0, 6))
        head.grid_columnconfigure(2, weight=1)
        title_frame = ctk.CTkFrame(head, fg_color="transparent", width=0, height=0)
        title_frame.grid(row=0, column=0, padx=(0, 10))
        ctk.CTkLabel(
            title_frame, text="RESULTS", text_color=_GRAY, font=ctk.CTkFont(size=11, weight="bold"),
        ).pack(side="left")
        # Explains the Radius Sweep panels; only shown while that tool is selected
        self._at_plots_help = help_icon(title_frame, _RADIUS_SWEEP_HELP["plots"])
        self._at_status_label = ctk.CTkLabel(
            head, text="Set parameters on the left, then click Run.",
            text_color=_GRAY, font=ctk.CTkFont(size=11), anchor="w",
        )
        self._at_status_label.grid(row=0, column=1, padx=(0, 10))
        ctk.CTkFrame(
            head, width=0, height=1, corner_radius=0, border_width=0, fg_color=theme.pair("line"),
        ).grid(row=0, column=2, sticky="ew")
        ctk.CTkButton(
            head, text="Save PNG", width=90, height=26, fg_color=_BTN2, hover_color=_BTN2_HOVER,
            font=ctk.CTkFont(size=12), command=self._save_active_analysis_plot,
        ).grid(row=0, column=3, padx=(10, 0))

        canvas_frame = ctk.CTkFrame(main, fg_color=_PLOT_BG)
        canvas_frame.grid(row=1, column=0, sticky="nsew")
        canvas_frame.grid_columnconfigure(0, weight=1)
        canvas_frame.grid_rowconfigure(0, weight=1)
        # Prevents the canvas from resizing this frame back (avoids a resize loop).
        canvas_frame.grid_propagate(False)

        # Both figures live in canvas_frame; only one is gridded at a time so
        # switching is instant and each tool keeps its plot without re-running.
        # Each follows the frame's size (see _bind_analysis_plot_resize).
        self._radius_sweep_fig = Figure(figsize=(10, 6), dpi=100, facecolor=theme.c("plot_bg"))
        self._radius_sweep_canvas = FigureCanvasTkAgg(self._radius_sweep_fig, master=canvas_frame)
        self._radius_sweep_canvas.get_tk_widget().grid(row=0, column=0, sticky="nsew")
        self._bind_analysis_plot_resize(canvas_frame, self._radius_sweep_canvas, self._radius_sweep_fig)

        self._cif_variance_fig = Figure(figsize=(10, 6), dpi=100, facecolor=theme.c("plot_bg"))
        self._cif_variance_canvas = FigureCanvasTkAgg(self._cif_variance_fig, master=canvas_frame)
        self._cif_variance_canvas.get_tk_widget().grid(row=0, column=0, sticky="nsew")
        self._bind_analysis_plot_resize(canvas_frame, self._cif_variance_canvas, self._cif_variance_fig)

        # Shown over the plot area until the active tool has a result to draw
        self._at_empty_label = ctk.CTkLabel(
            canvas_frame, text="", text_color=_GRAY, font=ctk.CTkFont(size=13),
        )

        # Collapsed log drawer under the plot (the plot is this tab's primary content)
        self._at_log_visible = False
        self._at_log_frame = ctk.CTkFrame(main, fg_color="transparent")
        self._at_log = ctk.CTkTextbox(
            self._at_log_frame, height=140, font=ctk.CTkFont(family="Courier New", size=12),
            wrap="word", state="disabled", border_width=1, border_color=theme.pair("line"),
        )
        self._at_log.pack(fill="both", expand=True)
        isolate_textbox_scroll(self._at_log)
        add_resize_grip(self._at_log).pack(fill="x")
        self._at_log_toggle = ctk.CTkButton(
            main, text="▸  Show details", height=26, anchor="w", width=0,
            fg_color="transparent", hover_color=_CHIP_HOVER,
            font=ctk.CTkFont(size=12), command=self._at_toggle_log,
        )
        self._at_log_toggle.grid(row=3, column=0, sticky="w", pady=(6, 0))

        self._rebuild_analysis_tool_fields()
        self._show_active_analysis_plot()

    def _bind_analysis_plot_resize(self, frame, canvas, fig) -> None:
        """Keep *fig* sized to fill *frame* as the window resizes. Debounced; the
        layout is re-fitted after each resize because the figure's own
        tight_layout ran at the old size.
        """
        state = {"after_id": None}

        def _apply(width_px: int, height_px: int) -> None:
            state["after_id"] = None
            if width_px < 200 or height_px < 150:
                return
            dpi = fig.get_dpi()
            new_w, new_h = width_px / dpi, height_px / dpi
            if abs(new_w - fig.get_figwidth()) < 0.05 and abs(new_h - fig.get_figheight()) < 0.05:
                return
            fig.set_size_inches(new_w, new_h)
            if fig.axes:
                try:
                    fig.tight_layout()
                except Exception:
                    pass
                canvas.draw()

        def _on_configure(event) -> None:
            if state["after_id"] is not None:
                self.after_cancel(state["after_id"])
            state["after_id"] = self.after(150, lambda: _apply(event.width, event.height))

        frame.bind("<Configure>", _on_configure, add="+")

    def _at_toggle_log(self) -> None:
        """Show or hide this tab's own collapsible log panel."""
        if self._at_log_visible:
            self._at_log_frame.grid_remove()
            self._at_log_toggle.configure(text="▸  Show details")
            self._at_log_visible = False
        else:
            self._at_log_frame.grid(row=2, column=0, sticky="ew", pady=(8, 0))
            self._at_log_toggle.configure(text="▾  Hide details")
            self._at_log_visible = True

    def _rebuild_analysis_tool_fields(self) -> None:
        """Clear and rebuild the params frame for whichever sub-tool is selected."""
        for w in self._at_params_frame.winfo_children():
            w.destroy()
        self._at_params_frame.grid_columnconfigure(0, weight=1)
        if self._analysis_subtool_var.get() == "Radius Sweep":
            self._build_radius_sweep_fields()
        else:
            self._build_cif_variance_fields()

    def _at_field_label(self, row: int, text: str, help_text: str | None = None) -> None:
        """Small muted label (plus optional help icon) above a field in the params panel."""
        frame = ctk.CTkFrame(self._at_params_frame, fg_color="transparent", width=0, height=0)
        frame.grid(row=row, column=0, sticky="w", pady=(12, 3))
        ctk.CTkLabel(
            frame, text=text, anchor="w", text_color=_GRAY, font=ctk.CTkFont(size=12),
        ).pack(side="left")
        if help_text:
            help_icon(frame, help_text).pack(side="left", padx=(5, 0))

    def _build_radius_sweep_fields(self):
        """Radius Sweep parameter fields, built into self._at_params_frame."""
        f = self._at_params_frame
        wrap = _SIDE_W - 50

        # Genes. self._radius_genes is the source of truth, persisted across sub-tool toggles
        self._at_field_label(0, "Genes / UniProt IDs", _RADIUS_SWEEP_HELP["genes"])
        if not hasattr(self, "_radius_genes"):
            self._radius_genes: list[str] = []
        gene_input_frame = ctk.CTkFrame(f, fg_color="transparent")
        gene_input_frame.grid(row=1, column=0, sticky="ew")
        gene_input_frame.grid_columnconfigure(0, weight=1)
        self._radius_gene_input_var = ctk.StringVar(value="")
        gene_entry = ctk.CTkEntry(
            gene_input_frame, textvariable=self._radius_gene_input_var,
            placeholder_text="Gene symbol or UniProt ID",
        )
        gene_entry.grid(row=0, column=0, sticky="ew", padx=(0, 6))
        gene_entry.bind("<Return>", lambda _e: self._add_radius_gene())
        ctk.CTkButton(
            gene_input_frame, text="+ Add", width=64, height=28,
            command=self._add_radius_gene,
        ).grid(row=0, column=1)

        # Error feedback for _add_radius_gene — hidden until a lookup fails.
        self._radius_gene_error_label = ctk.CTkLabel(
            f, text="", text_color=_RED, font=ctk.CTkFont(size=11), anchor="w",
            justify="left", wraplength=wrap,
        )
        self._radius_gene_error_label.grid(row=2, column=0, sticky="w", pady=(4, 0))
        self._radius_gene_error_label.grid_remove()

        self._radius_genes_list_frame = ctk.CTkFrame(f, fg_color="transparent")
        self._radius_genes_list_frame.grid(row=3, column=0, sticky="ew", pady=(6, 0))
        self._refresh_radius_gene_chips()

        # Radius range
        self._at_field_label(4, "Radius range (Å)", _RADIUS_SWEEP_HELP["radius_range"])
        range_frame = ctk.CTkFrame(f, fg_color="transparent")
        range_frame.grid(row=5, column=0, sticky="w")
        if not hasattr(self, "_radius_start_var"):
            self._radius_start_var = ctk.StringVar(value="4")
            self._radius_stop_var = ctk.StringVar(value="20")
            self._radius_step_var = ctk.StringVar(value="1")
        for label, var in [
            ("start", self._radius_start_var),
            ("stop", self._radius_stop_var),
            ("step", self._radius_step_var),
        ]:
            ctk.CTkLabel(range_frame, text=label, text_color=_GRAY, font=ctk.CTkFont(size=11)).pack(
                side="left", padx=(0, 4),
            )
            ctk.CTkEntry(range_frame, textvariable=var, width=52).pack(side="left", padx=(0, 12))

        # Hotspot recurrence threshold, computed live from raw COSMIC -- independent
        # of the "Min samples" value the Pipeline tab used to build the intermediate TSV
        self._at_field_label(6, "Min samples")
        if not hasattr(self, "_radius_min_cases_var"):
            from radius_sweep import DEFAULT_MIN_CASES
            self._radius_min_cases_var = ctk.StringVar(value=str(DEFAULT_MIN_CASES))
        ctk.CTkEntry(f, textvariable=self._radius_min_cases_var, width=64).grid(row=7, column=0, sticky="w")
        ctk.CTkLabel(
            f, text="Minimum distinct COSMIC samples for a mutation to count as a hotspot",
            text_color=_GRAY, font=ctk.CTkFont(size=11), anchor="w", justify="left", wraplength=wrap,
        ).grid(row=8, column=0, sticky="w", pady=(3, 0))

        # Unfiltered comparison
        if not hasattr(self, "_radius_unfiltered_var"):
            self._radius_unfiltered_var = ctk.BooleanVar(value=False)
        unfiltered_frame = ctk.CTkFrame(f, fg_color="transparent")
        unfiltered_frame.grid(row=9, column=0, sticky="w", pady=(14, 0))
        ctk.CTkCheckBox(
            unfiltered_frame, text="Include unfiltered COSMIC comparison",
            variable=self._radius_unfiltered_var,
            checkbox_width=18, checkbox_height=18,
        ).pack(side="left")
        help_icon(unfiltered_frame, _RADIUS_SWEEP_HELP["unfiltered"]).pack(side="left", padx=(4, 0))

        # Output folder
        self._at_field_label(
            10, "Output folder",
            "Where radius_sweep.tsv (and radius_sweep_plot.png, if you use "
            "Save PNG) are written. Defaults to the app's main Output folder.",
        )
        if not hasattr(self, "_radius_output_dir_var"):
            self._radius_output_dir_var = ctk.StringVar(value=str(OUTPUT_DIR))
        ctk.CTkEntry(
            f, textvariable=self._radius_output_dir_var, state="readonly",
        ).grid(row=11, column=0, sticky="ew")
        buttons = ctk.CTkFrame(f, fg_color="transparent")
        buttons.grid(row=12, column=0, sticky="w", pady=(6, 0))
        ctk.CTkButton(
            buttons, text="Browse", width=70, height=26, font=ctk.CTkFont(size=12),
            command=self._browse_radius_output_dir,
        ).pack(side="left", padx=(0, 6))
        ctk.CTkButton(
            buttons, text="Open Folder", width=90, height=26, font=ctk.CTkFont(size=12),
            fg_color=_BTN2, hover_color=_BTN2_HOVER,
            command=lambda: self._open_folder(Path(self._radius_output_dir_var.get())),
        ).pack(side="left")

    def _browse_radius_output_dir(self) -> None:
        path = filedialog.askdirectory(
            title="Select Radius Sweep output folder",
            initialdir=self._radius_output_dir_var.get(),
        )
        if path:
            self._radius_output_dir_var.set(path)

    def _add_radius_gene(self) -> None:
        """Add a gene (symbol or UniProt accession, auto-detected). Validates
        upfront that the dataset has an entry and its AlphaFold CIF is downloaded,
        rather than silently skipping it later at run time.
        """
        token = self._radius_gene_input_var.get().strip()
        if not token:
            return

        from radius_sweep import PTM_TSV, has_cif, has_multiple_fragments, resolve_gene_token

        if not PTM_TSV.exists():
            # Can't validate yet -- the Run button's own check reports this clearly
            gene = token.upper()
            if gene not in self._radius_genes:
                self._radius_gene_error_label.grid_remove()
                self._radius_genes.append(gene)
                self._refresh_radius_gene_chips()
            self._radius_gene_input_var.set("")
            return

        resolved = resolve_gene_token(token)
        if resolved is None:
            self._radius_gene_error_label.configure(
                text=f"⚠  No data found for '{token}'. Check the spelling, or run "
                     f"the PTM Proximity/Mutation Clustering pipeline (step 1) first.",
                text_color=_RED,
            )
            self._radius_gene_error_label.grid()
            return

        gene, uid = resolved

        if gene in self._radius_genes:
            self._radius_gene_error_label.grid_remove()
            self._radius_gene_input_var.set("")
            return

        if not has_cif(uid):
            self._radius_gene_error_label.configure(
                text=f"⚠  No AlphaFold structure (CIF file) downloaded for {gene} ({uid}). "
                     f"Download it via the Pipeline tab (step 2) first.",
                text_color=_RED,
            )
            self._radius_gene_error_label.grid()
            return

        self._radius_genes.append(gene)
        self._refresh_radius_gene_chips()
        self._radius_gene_input_var.set("")

        # Non-blocking warning: only fragment 1 is analyzed for multi-fragment proteins
        if has_multiple_fragments(uid):
            self._radius_gene_error_label.configure(
                text=f"⚠  {gene} ({uid}) spans multiple AlphaFold fragments — only "
                     f"fragment 1 is analyzed, so results may be incomplete for this protein.",
                text_color=_YELLOW,
            )
            self._radius_gene_error_label.grid()
        else:
            self._radius_gene_error_label.grid_remove()

    def _remove_radius_gene(self, gene: str) -> None:
        if gene in self._radius_genes:
            self._radius_genes.remove(gene)
            self._refresh_radius_gene_chips()

    def _refresh_radius_gene_chips(self) -> None:
        """Redraw the added-genes chip list, each removable via its own ✕ button."""
        for w in self._radius_genes_list_frame.winfo_children():
            w.destroy()

        if not self._radius_genes:
            ctk.CTkLabel(
                self._radius_genes_list_frame, text="None added yet.",
                text_color=_GRAY, font=ctk.CTkFont(size=11),
            ).pack(anchor="w")
            return

        row = None
        for i, gene in enumerate(self._radius_genes):
            if i % 3 == 0:
                row = ctk.CTkFrame(self._radius_genes_list_frame, fg_color="transparent")
                row.pack(anchor="w", pady=(0, 4))
            chip = ctk.CTkFrame(row, fg_color=_CHIP, corner_radius=6)
            chip.pack(side="left", padx=(0, 6))
            ctk.CTkLabel(chip, text=gene, font=ctk.CTkFont(size=12)).pack(
                side="left", padx=(8, 4), pady=4,
            )
            ctk.CTkButton(
                chip, text="✕", width=20, height=20, fg_color="transparent",
                hover_color=_CHIP_HOVER, font=ctk.CTkFont(size=11),
                command=lambda g=gene: self._remove_radius_gene(g),
            ).pack(side="left", padx=(0, 6), pady=4)

    def _build_cif_variance_fields(self):
        """CIF Variance parameter fields, built into self._at_params_frame."""
        f = self._at_params_frame
        wrap = _SIDE_W - 50

        # Input folder
        self._at_field_label(0, "Input folder", _CIF_VARIANCE_HELP["input_dir"])
        if not hasattr(self, "_variance_input_dir_var"):
            from cif_variance import DEFAULT_INPUT_DIR
            self._variance_input_dir_var = ctk.StringVar(value=str(DEFAULT_INPUT_DIR))
            self._variance_input_dir_var.trace_add("write", self._update_variance_cif_count)
        input_row = ctk.CTkFrame(f, fg_color="transparent")
        input_row.grid(row=1, column=0, sticky="ew")
        input_row.grid_columnconfigure(0, weight=1)
        ctk.CTkEntry(input_row, textvariable=self._variance_input_dir_var).grid(
            row=0, column=0, sticky="ew", padx=(0, 6),
        )
        ctk.CTkButton(
            input_row, text="Browse", width=70, height=26, font=ctk.CTkFont(size=12),
            command=self._browse_variance_input_dir,
        ).grid(row=0, column=1)
        self._variance_cif_count_label = ctk.CTkLabel(
            f, text="", anchor="w", font=ctk.CTkFont(size=11),
        )
        self._variance_cif_count_label.grid(row=2, column=0, sticky="w", pady=(3, 0))
        self._update_variance_cif_count()

        # Top N
        self._at_field_label(3, "Top N residues", _CIF_VARIANCE_HELP["top_n"])
        if not hasattr(self, "_variance_top_var"):
            self._variance_top_var = ctk.StringVar(value="10")
        ctk.CTkEntry(f, textvariable=self._variance_top_var, width=64).grid(row=4, column=0, sticky="w")

        # Report / align ranges
        # (Placeholder text doesn't show on entries bound to a StringVar, so each range
        # spells out "from ... to ..." and what a blank end means.)
        def _range_row(row: int, start_var, end_var, blank_hint: str) -> None:
            frame = ctk.CTkFrame(f, fg_color="transparent")
            frame.grid(row=row, column=0, sticky="w")
            small = ctk.CTkFont(size=11)
            ctk.CTkLabel(frame, text="from", text_color=_GRAY, font=small).pack(side="left", padx=(0, 4))
            ctk.CTkEntry(frame, textvariable=start_var, width=70).pack(side="left", padx=(0, 10))
            ctk.CTkLabel(frame, text="to", text_color=_GRAY, font=small).pack(side="left", padx=(0, 4))
            ctk.CTkEntry(frame, textvariable=end_var, width=70).pack(side="left")
            ctk.CTkLabel(
                f, text=blank_hint, text_color=_GRAY, font=small, anchor="w",
            ).grid(row=row + 1, column=0, sticky="w", pady=(3, 0))

        self._at_field_label(5, "Report range", _CIF_VARIANCE_HELP["report_range"])
        if not hasattr(self, "_variance_range_start_var"):
            self._variance_range_start_var = ctk.StringVar(value="")
            self._variance_range_end_var = ctk.StringVar(value="")
        _range_row(6, self._variance_range_start_var, self._variance_range_end_var, "Blank end = all residues")

        self._at_field_label(8, "Align range", _CIF_VARIANCE_HELP["align_range"])
        if not hasattr(self, "_variance_align_start_var"):
            self._variance_align_start_var = ctk.StringVar(value="")
            self._variance_align_end_var = ctk.StringVar(value="")
        _range_row(
            9, self._variance_align_start_var, self._variance_align_end_var, "Blank = same as the report range",
        )

        # UniProt / gene overrides
        self._at_field_label(11, "UniProt override", _CIF_VARIANCE_HELP["uniprot_override"])
        if not hasattr(self, "_variance_uniprot_var"):
            self._variance_uniprot_var = ctk.StringVar(value="")
        ctk.CTkEntry(
            f, textvariable=self._variance_uniprot_var, placeholder_text="auto-detected from CIF",
        ).grid(row=12, column=0, sticky="ew")

        self._at_field_label(13, "Gene override", _CIF_VARIANCE_HELP["gene_override"])
        if not hasattr(self, "_variance_gene_var"):
            self._variance_gene_var = ctk.StringVar(value="")
        ctk.CTkEntry(
            f, textvariable=self._variance_gene_var, placeholder_text="optional, for UniProt lookup",
        ).grid(row=14, column=0, sticky="ew")

        # AlphaFold Server seed-JSON export
        from cif_variance import DEFAULT_SEEDS
        self._at_field_label(15, "AlphaFold seeds", _CIF_VARIANCE_HELP["seed_json"])
        seed_json_frame = ctk.CTkFrame(f, fg_color="transparent")
        seed_json_frame.grid(row=16, column=0, sticky="w")
        if not hasattr(self, "_variance_num_seeds_var"):
            self._variance_num_seeds_var = ctk.StringVar(value=str(len(DEFAULT_SEEDS)))
        ctk.CTkEntry(seed_json_frame, textvariable=self._variance_num_seeds_var, width=44).pack(
            side="left", padx=(0, 8),
        )
        self._variance_seed_json_btn = ctk.CTkButton(
            seed_json_frame, text="⬇  Generate seeds JSON", width=180, height=28,
            font=ctk.CTkFont(size=12), fg_color=_BTN2, hover_color=_BTN2_HOVER,
            command=self._generate_alphafold_seed_json,
        )
        self._variance_seed_json_btn.pack(side="left")
        self._variance_seed_json_status = ctk.CTkLabel(
            f, text="", anchor="w", justify="left", wraplength=wrap, font=ctk.CTkFont(size=11),
        )
        self._variance_seed_json_status.grid(row=17, column=0, sticky="w", pady=(3, 0))

        # Output folder
        self._at_field_label(
            18, "Output folder",
            "Where the variance plot/data, pairwise RMSD, and any generated "
            "AlphaFold Seeds JSON are written. Defaults to Output/cif_variance.",
        )
        if not hasattr(self, "_variance_output_dir_var"):
            self._variance_output_dir_var = ctk.StringVar(value=str(OUTPUT_DIR / "cif_variance"))
        ctk.CTkEntry(
            f, textvariable=self._variance_output_dir_var, state="readonly",
        ).grid(row=19, column=0, sticky="ew")
        buttons = ctk.CTkFrame(f, fg_color="transparent")
        buttons.grid(row=20, column=0, sticky="w", pady=(6, 0))
        ctk.CTkButton(
            buttons, text="Browse", width=70, height=26, font=ctk.CTkFont(size=12),
            command=self._browse_variance_output_dir,
        ).pack(side="left", padx=(0, 6))
        ctk.CTkButton(
            buttons, text="Open Folder", width=90, height=26, font=ctk.CTkFont(size=12),
            fg_color=_BTN2, hover_color=_BTN2_HOVER,
            command=lambda: self._open_folder(Path(self._variance_output_dir_var.get())),
        ).pack(side="left")

    def _browse_variance_output_dir(self) -> None:
        path = filedialog.askdirectory(
            title="Select CIF Variance output folder",
            initialdir=self._variance_output_dir_var.get(),
        )
        if path:
            self._variance_output_dir_var.set(path)

    def _generate_alphafold_seed_json(self) -> None:
        """Resolve the CIF Variance tool's target protein and write an
        AlphaFold Server batch JSON of one separate job per seed (1-N, see
        the Seeds field), so the researcher can upload it at
        alphafoldserver.com and drop the resulting CIFs straight into the
        input folder above for comparison.

        The UniProt-sequence fetch is network I/O, so it runs off the main
        thread like the app's other network lookups (e.g. CA-coordinate
        export). Tkinter itself isn't thread-safe -- the worker thread must
        never touch a Tk widget or call self.after() directly, so it only
        puts its result on a plain queue.Queue; a self.after()-scheduled poll
        (always running on the main thread) drains it and does the actual
        widget updates.
        """
        num_seeds_text = self._variance_num_seeds_var.get().strip()
        if not num_seeds_text.isdigit() or int(num_seeds_text) < 1:
            self._variance_seed_json_status.configure(
                text="⚠  Seeds must be a whole number of 1 or more.", text_color=_RED,
            )
            return
        num_seeds = int(num_seeds_text)

        input_dir = Path(self._variance_input_dir_var.get().strip())
        uniprot = self._variance_uniprot_var.get().strip() or None
        gene = self._variance_gene_var.get().strip() or None
        output_dir = Path(self._variance_output_dir_var.get())

        self._variance_seed_json_btn.configure(state="disabled")
        self._variance_seed_json_status.configure(
            text="Resolving protein and fetching its sequence from UniProt…", text_color=_GRAY,
        )

        result_q: queue.Queue = queue.Queue()

        def _worker():
            from cif_variance import generate_alphafold_seed_json
            try:
                path = generate_alphafold_seed_json(
                    input_dir=input_dir, output_dir=output_dir, uniprot=uniprot, gene=gene,
                    seeds=list(range(1, num_seeds + 1)),
                )
            except Exception as exc:
                result_q.put(("error", str(exc)))
            else:
                result_q.put(("ok", path))

        threading.Thread(target=_worker, daemon=True).start()
        self.after(100, lambda: self._poll_seed_json_result(result_q))

    def _poll_seed_json_result(self, result_q: "queue.Queue") -> None:
        try:
            kind, value = result_q.get_nowait()
        except queue.Empty:
            self.after(100, lambda: self._poll_seed_json_result(result_q))
            return
        if kind == "error":
            self._finish_seed_json(error=value)
        else:
            self._finish_seed_json(path=value)

    def _finish_seed_json(self, path: Path | None = None, error: str | None = None) -> None:
        self._variance_seed_json_btn.configure(state="normal")
        if error:
            self._variance_seed_json_status.configure(text=f"⚠  {error}", text_color=_RED)
        else:
            self._variance_seed_json_status.configure(
                text=f"✓  Wrote {path.name} to {path.parent}", text_color=_GREEN,
            )

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

    def _show_active_analysis_plot(self) -> None:
        """Show whichever canvas matches the current segmented-button selection."""
        if self._analysis_subtool_var.get() == "Radius Sweep":
            self._cif_variance_canvas.get_tk_widget().grid_remove()
            self._radius_sweep_canvas.get_tk_widget().grid()
            fig, hint = self._radius_sweep_fig, "Add genes and click Run to plot hotspot clusters by radius."
            self._at_plots_help.pack(side="left", padx=(6, 0))
        else:
            self._radius_sweep_canvas.get_tk_widget().grid_remove()
            self._cif_variance_canvas.get_tk_widget().grid()
            fig, hint = self._cif_variance_fig, "Choose a folder of .cif models and click Run to compare them."
            self._at_plots_help.pack_forget()
        if fig.axes:
            self._at_empty_label.place_forget()
        else:
            self._at_empty_label.configure(text=hint)
            self._at_empty_label.place(relx=0.5, rely=0.5, anchor="center")
            self._at_empty_label.lift()

    def _save_active_analysis_plot(self) -> None:
        if self._analysis_subtool_var.get() == "Radius Sweep":
            self._save_radius_sweep_plot()
        else:
            self._save_cif_variance_plot()

    def _draw_radius_sweep_plot(self) -> None:
        from radius_sweep import build_sweep_figure

        self._radius_sweep_fig.clf()
        build_sweep_figure(self._radius_sweep_result, fig=self._radius_sweep_fig)
        self._style_themed_figure(self._radius_sweep_fig)
        self._radius_sweep_canvas.draw()
        self._show_active_analysis_plot()
        self._at_status_label.configure(text="Sweep complete.", text_color=_GREEN)

    def _save_radius_sweep_plot(self) -> None:
        if not self._radius_sweep_fig.axes:
            self._at_status_label.configure(text="Run a sweep before saving.", text_color=_RED)
            return
        out_dir = Path(self._radius_output_dir_var.get())
        out_dir.mkdir(parents=True, exist_ok=True)
        out_path = out_dir / "radius_sweep_plot.png"
        self._radius_sweep_fig.savefig(
            out_path, dpi=200, facecolor=self._radius_sweep_fig.get_facecolor(), bbox_inches="tight",
        )
        self._at_status_label.configure(text=f"Saved to {out_path}", text_color=_GREEN)

    def _draw_cif_variance_plot(self) -> None:
        from cif_variance import build_variance_figure

        self._cif_variance_fig.clf()
        build_variance_figure(self._cif_variance_result, fig=self._cif_variance_fig)
        self._style_themed_figure(self._cif_variance_fig)
        self._cif_variance_canvas.draw()
        self._show_active_analysis_plot()
        self._at_status_label.configure(text="Analysis complete.", text_color=_GREEN)

    def _save_cif_variance_plot(self) -> None:
        if not self._cif_variance_fig.axes:
            self._at_status_label.configure(text="Run an analysis before saving.", text_color=_RED)
            return
        out_dir = Path(self._variance_output_dir_var.get())
        out_dir.mkdir(parents=True, exist_ok=True)
        out_path = out_dir / "cif_variance_plot.png"
        self._cif_variance_fig.savefig(
            out_path, dpi=200, facecolor=self._cif_variance_fig.get_facecolor(), bbox_inches="tight",
        )
        self._at_status_label.configure(text=f"Saved to {out_path}", text_color=_GREEN)
