"""Analysis Tools tab: self-contained Radius Sweep / CIF Variance UI —
parameters, a Run trigger, status/progress/log feedback, and the result plots.

Running a tool here shares the app's single background-run slot with the
Pipeline tab (`self._running`, guarded in PipelineRunnerMixin._start_analysis_tool_run)
since both use the same thread+queue execution engine — only one run (of any
kind) can be active at a time. `_poll_queue` (pipeline_runner.py) routes the
generic "status"/"progress"/"log"/etc. queue messages to this tab's widgets
instead of the Pipeline tab's step-indexed ones whenever `self._at_run_active`
is set, via `_status_target`/`_progress_target`. `_style_dark_figure` calls
`self._style_lollipop_axis`, a VisualizationTabMixin method.
"""
from __future__ import annotations

from pathlib import Path
from tkinter import filedialog

import customtkinter as ctk

from ui.common import (
    _GRAY, _GREEN, _RED, _YELLOW, add_resize_grip, isolate_textbox_scroll,
    help_icon, _RADIUS_SWEEP_HELP, _CIF_VARIANCE_HELP,
)


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

        # Sole creation site now that the Pipeline tab no longer has its own
        # Analysis Tools mode/panel.
        self._analysis_subtool_var = ctk.StringVar(value="Radius Sweep")

        # One scrollable section for the whole tab (controls, params, run row,
        # log, and the plot all together) — mirrors the Pipeline tab's own
        # _build_pipeline_tab exactly, so there's a single scrollbar for the
        # tab instead of the plot having its own separate scrolling area.
        tab.grid_columnconfigure(0, weight=1)
        tab.grid_rowconfigure(0, weight=1)

        scroll = ctk.CTkScrollableFrame(tab, fg_color="transparent")
        scroll.grid(row=0, column=0, sticky="nsew")
        scroll.grid_columnconfigure(0, weight=1)

        # CTkScrollableFrame's vertical mode only ever shows a vertical
        # scrollbar, and forces content width to exactly match the viewport —
        # fine for narrow content, but it clips anything wider (like a big
        # plot) with no way to reach the rest. It already has working
        # horizontal-scroll support internally (Shift+wheel), so add a visible
        # horizontal scrollbar and replace the forced-width handler with one
        # that only stretches narrow content up to the viewport width, never
        # shrinks wide content down to it.
        scroll._parent_canvas.unbind("<Configure>")

        def _fit_width_at_least_viewport(event):
            natural_w = scroll.winfo_reqwidth()
            scroll._parent_canvas.itemconfigure(
                scroll._create_window_id, width=max(natural_w, event.width),
            )

        scroll._parent_canvas.bind("<Configure>", _fit_width_at_least_viewport)

        at_h_scrollbar = ctk.CTkScrollbar(
            scroll._parent_frame, orientation="horizontal", command=scroll._parent_canvas.xview,
        )
        at_h_scrollbar.grid(row=2, column=0, sticky="ew")
        scroll._parent_canvas.configure(xscrollcommand=at_h_scrollbar.set)

        p = scroll  # everything below goes in the scrollable frame

        controls = ctk.CTkFrame(p)
        controls.grid(row=0, column=0, padx=12, pady=(12, 6), sticky="ew")

        ctk.CTkSegmentedButton(
            controls, values=["Radius Sweep", "CIF Variance"],
            variable=self._analysis_subtool_var,
            command=lambda _v: (self._rebuild_analysis_tool_fields(), self._show_active_analysis_plot()),
        ).pack(side="left", padx=(12, 12), pady=10)

        ctk.CTkButton(
            controls, text="Save PNG", width=100,
            fg_color="gray30", hover_color="gray40",
            command=self._save_active_analysis_plot,
        ).pack(side="left", padx=(0, 12), pady=10)

        # Parameter fields — cleared and rebuilt by _rebuild_analysis_tool_fields
        # whenever the segmented button above toggles Radius Sweep <-> CIF Variance.
        self._at_params_frame = ctk.CTkFrame(p)
        # sticky="w", not "ew": this frame shares its grid column with the much
        # wider plot canvas below it (canvas_frame, row 5 -- the matplotlib
        # figure alone renders at 1400px). With "ew" this frame was forced to
        # stretch to match that column width, dragging its own Browse button
        # far past the edge of any normal-sized window. "w" keeps it sized to
        # its own natural content instead of the unrelated plot's width.
        self._at_params_frame.grid(row=1, column=0, padx=12, pady=(0, 6), sticky="w")
        self._at_params_frame.grid_columnconfigure(1, weight=1)

        # Run trigger + status/progress — persistent, built once (unlike the
        # params frame above, this tab is never torn down on toggle).
        run_row = ctk.CTkFrame(p, fg_color="transparent")
        run_row.grid(row=2, column=0, padx=12, pady=(0, 6), sticky="ew")
        run_row.grid_columnconfigure(2, weight=1)

        self._at_run_btn = ctk.CTkButton(
            run_row, text="▶  Run", width=100, height=32,
            font=ctk.CTkFont(size=13, weight="bold"),
            command=lambda: self._start_analysis_tool_run(
                "radius-sweep" if self._analysis_subtool_var.get() == "Radius Sweep" else "cif-variance"
            ),
        )
        self._at_run_btn.grid(row=0, column=0, padx=(12, 12), pady=8)

        self._at_progress_bar = ctk.CTkProgressBar(run_row, width=140, height=14)
        self._at_progress_bar.set(0)
        self._at_progress_bar.grid(row=0, column=1, padx=(0, 12), pady=8, sticky="w")
        self._at_progress_bar.grid_remove()

        self._at_status_label = ctk.CTkLabel(
            run_row, text="Set parameters above, then click Run.",
            text_color=_GRAY, font=ctk.CTkFont(size=11), anchor="w",
        )
        self._at_status_label.grid(row=0, column=2, padx=(0, 12), pady=8, sticky="w")

        # Log (collapsible, mirrors the Pipeline tab's _toggle_log — collapsed by
        # default here since the plot is this tab's primary content).
        self._at_log_visible = False
        self._at_log_toggle = ctk.CTkButton(
            p, text="Show Details", width=120, height=28,
            font=ctk.CTkFont(size=12), fg_color="gray30", hover_color="gray40",
            command=self._at_toggle_log,
        )
        self._at_log_toggle.grid(row=3, column=0, padx=12, pady=(0, 4), sticky="w")

        self._at_log_frame = ctk.CTkFrame(p, fg_color="transparent")
        self._at_log = ctk.CTkTextbox(
            self._at_log_frame, height=120, font=ctk.CTkFont(family="Courier New", size=12),
            wrap="word", state="disabled",
        )
        self._at_log.pack(fill="both", expand=True)
        isolate_textbox_scroll(self._at_log)
        add_resize_grip(self._at_log).pack(fill="x")

        # Plot area: a plain frame, part of the same single scrollable section
        # as everything above it — no independent scrollbar of its own.
        canvas_frame = ctk.CTkFrame(p, fg_color="#2b2b2b")
        canvas_frame.grid(row=5, column=0, padx=12, pady=(0, 12), sticky="nsew")

        # Both figures live in the same canvas_frame; only one is gridded at a time
        # (toggled by the segmented button above) so switching is instant and each
        # tool keeps its own plot in memory without re-running anything.
        self._radius_sweep_fig = Figure(figsize=(14, 9), dpi=100, facecolor="#2b2b2b")
        self._radius_sweep_canvas = FigureCanvasTkAgg(self._radius_sweep_fig, master=canvas_frame)
        self._radius_sweep_canvas.get_tk_widget().grid(row=0, column=0)

        self._cif_variance_fig = Figure(figsize=(14, 8), dpi=100, facecolor="#2b2b2b")
        self._cif_variance_canvas = FigureCanvasTkAgg(self._cif_variance_fig, master=canvas_frame)
        self._cif_variance_canvas.get_tk_widget().grid(row=0, column=0)

        self._rebuild_analysis_tool_fields()
        self._show_active_analysis_plot()

        # Hide each scrollbar when content already fits that axis; show it
        # only when scrolling is needed — same convention as the Pipeline
        # tab's own _update_scrollbar, extended to the horizontal axis too.
        def _update_scrollbar(*_):
            try:
                canvas = scroll._parent_canvas
                sr = canvas.cget("scrollregion")
                if not sr:
                    scroll._scrollbar.grid_remove()
                    at_h_scrollbar.grid_remove()
                    return
                _x0, _y0, x1, y1 = (int(float(v)) for v in sr.split())
                if y1 > canvas.winfo_height():
                    scroll._scrollbar.grid()
                else:
                    scroll._scrollbar.grid_remove()
                if x1 > canvas.winfo_width():
                    at_h_scrollbar.grid()
                else:
                    at_h_scrollbar.grid_remove()
            except Exception:
                pass

        scroll.bind("<Configure>", _update_scrollbar, add="+")
        tab.bind("<Configure>", _update_scrollbar, add="+")
        self.after(200, _update_scrollbar)

    def _at_toggle_log(self) -> None:
        """Show or hide this tab's own collapsible log panel."""
        if self._at_log_visible:
            self._at_log_frame.grid_remove()
            self._at_log_toggle.configure(text="Show Details")
            self._at_log_visible = False
        else:
            self._at_log_frame.grid(row=4, column=0, padx=12, pady=(0, 8), sticky="ew")
            self._at_log_toggle.configure(text="Hide Details")
            self._at_log_visible = True

    def _rebuild_analysis_tool_fields(self) -> None:
        """Clear and rebuild the params frame for whichever sub-tool is selected."""
        for w in self._at_params_frame.winfo_children():
            w.destroy()
        self._at_params_frame.grid_columnconfigure(1, weight=1)
        if self._analysis_subtool_var.get() == "Radius Sweep":
            self._build_radius_sweep_fields()
        else:
            self._build_cif_variance_fields()

    def _build_radius_sweep_fields(self):
        """Radius Sweep parameter fields, built into self._at_params_frame."""
        # Genes — added one at a time (no pre-filled defaults); self._radius_genes
        # is the source of truth, persisted across sub-tool toggles via the guard.
        genes_label_frame = ctk.CTkFrame(self._at_params_frame, fg_color="transparent")
        genes_label_frame.grid(row=0, column=0, padx=(12, 6), pady=6, sticky="nw")
        ctk.CTkLabel(genes_label_frame, text="Genes / UniProt IDs:", anchor="w").pack(side="left")
        help_icon(genes_label_frame, _RADIUS_SWEEP_HELP["genes"]).pack(side="left", padx=(4, 0))
        if not hasattr(self, "_radius_genes"):
            self._radius_genes: list[str] = []

        gene_input_frame = ctk.CTkFrame(self._at_params_frame, fg_color="transparent")
        gene_input_frame.grid(row=0, column=1, columnspan=3, padx=6, pady=6, sticky="w")

        self._radius_gene_input_var = ctk.StringVar(value="")
        gene_entry = ctk.CTkEntry(
            gene_input_frame, textvariable=self._radius_gene_input_var, width=180,
            placeholder_text="Gene symbol or UniProt ID",
        )
        gene_entry.pack(side="left", padx=(0, 6))
        gene_entry.bind("<Return>", lambda _e: self._add_radius_gene())

        ctk.CTkButton(
            gene_input_frame, text="+ Add", width=70, height=28,
            command=self._add_radius_gene,
        ).pack(side="left")

        # Error feedback for _add_radius_gene — hidden until a lookup fails.
        self._radius_gene_error_label = ctk.CTkLabel(
            self._at_params_frame, text="", text_color=_RED,
            font=ctk.CTkFont(size=11), anchor="w",
        )
        self._radius_gene_error_label.grid(row=1, column=1, columnspan=3, padx=6, pady=(0, 4), sticky="w")
        self._radius_gene_error_label.grid_remove()

        self._radius_genes_list_frame = ctk.CTkFrame(self._at_params_frame, fg_color="transparent")
        self._radius_genes_list_frame.grid(row=2, column=1, columnspan=3, padx=6, pady=(0, 6), sticky="ew")
        self._refresh_radius_gene_chips()

        # Radius range
        range_label_frame = ctk.CTkFrame(self._at_params_frame, fg_color="transparent")
        range_label_frame.grid(row=3, column=0, padx=(12, 6), pady=6, sticky="w")
        ctk.CTkLabel(range_label_frame, text="Radius range (Å):", anchor="w").pack(side="left")
        help_icon(range_label_frame, _RADIUS_SWEEP_HELP["radius_range"]).pack(side="left", padx=(4, 0))
        range_frame = ctk.CTkFrame(self._at_params_frame, fg_color="transparent")
        range_frame.grid(row=3, column=1, columnspan=3, padx=6, pady=6, sticky="w")
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

        # Min samples — the hotspot recurrence threshold, computed live from raw
        # COSMIC by Radius Sweep itself, independent of whatever "Min samples"
        # value the main Pipeline tab used to build the intermediate TSV.
        ctk.CTkLabel(self._at_params_frame, text="Min samples:", anchor="w").grid(
            row=4, column=0, padx=(12, 6), pady=6, sticky="w"
        )
        if not hasattr(self, "_radius_min_cases_var"):
            from radius_sweep import DEFAULT_MIN_CASES
            self._radius_min_cases_var = ctk.StringVar(value=str(DEFAULT_MIN_CASES))
        ctk.CTkEntry(
            self._at_params_frame, textvariable=self._radius_min_cases_var, width=60,
        ).grid(row=4, column=1, padx=6, pady=6, sticky="w")
        ctk.CTkLabel(
            self._at_params_frame,
            text="Minimum distinct COSMIC samples for a mutation to count as a hotspot",
            text_color=_GRAY, font=ctk.CTkFont(size=11), anchor="w",
        ).grid(row=4, column=2, columnspan=2, padx=6, pady=6, sticky="w")

        # Unfiltered comparison
        if not hasattr(self, "_radius_unfiltered_var"):
            self._radius_unfiltered_var = ctk.BooleanVar(value=False)
        unfiltered_frame = ctk.CTkFrame(self._at_params_frame, fg_color="transparent")
        unfiltered_frame.grid(row=5, column=0, columnspan=4, padx=12, pady=(6, 10), sticky="w")
        ctk.CTkCheckBox(
            unfiltered_frame, text="Include unfiltered COSMIC comparison",
            variable=self._radius_unfiltered_var,
            checkbox_width=18, checkbox_height=18,
        ).pack(side="left")
        help_icon(unfiltered_frame, _RADIUS_SWEEP_HELP["unfiltered"]).pack(side="left", padx=(4, 0))

    def _add_radius_gene(self) -> None:
        """Add a gene, accepting either a gene symbol or a UniProt accession
        (auto-detected by format). Validates upfront — both that the dataset
        has an entry for it, and that its AlphaFold CIF is actually
        downloaded — rather than silently skipping it later at run time.
        """
        token = self._radius_gene_input_var.get().strip()
        if not token:
            return

        from radius_sweep import PTM_TSV, resolve_gene_token, has_cif, has_multiple_fragments

        if not PTM_TSV.exists():
            # Can't validate yet — the Run button's own check already reports
            # this clearly ("Run the pipeline (step 1) first").
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

        # Non-blocking: the gene is still added, but large proteins AlphaFold
        # split into multiple fragments only get fragment 1 analyzed here.
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
            if i % 6 == 0:
                row = ctk.CTkFrame(self._radius_genes_list_frame, fg_color="transparent")
                row.pack(anchor="w", pady=(0, 4))
            chip = ctk.CTkFrame(row, fg_color="#3a3a3a", corner_radius=6)
            chip.pack(side="left", padx=(0, 6))
            ctk.CTkLabel(chip, text=gene, font=ctk.CTkFont(size=12)).pack(
                side="left", padx=(8, 4), pady=4,
            )
            ctk.CTkButton(
                chip, text="✕", width=20, height=20, fg_color="transparent",
                hover_color="#4a4a4a", font=ctk.CTkFont(size=11),
                command=lambda g=gene: self._remove_radius_gene(g),
            ).pack(side="left", padx=(0, 6), pady=4)

    def _build_cif_variance_fields(self):
        """CIF Variance parameter fields, built into self._at_params_frame."""
        # Input folder
        input_label_frame = ctk.CTkFrame(self._at_params_frame, fg_color="transparent")
        input_label_frame.grid(row=0, column=0, padx=(12, 6), pady=6, sticky="w")
        ctk.CTkLabel(input_label_frame, text="Input folder:", anchor="w").pack(side="left")
        help_icon(input_label_frame, _CIF_VARIANCE_HELP["input_dir"]).pack(side="left", padx=(4, 0))
        if not hasattr(self, "_variance_input_dir_var"):
            from cif_variance import DEFAULT_INPUT_DIR
            self._variance_input_dir_var = ctk.StringVar(value=str(DEFAULT_INPUT_DIR))
            self._variance_input_dir_var.trace_add("write", self._update_variance_cif_count)
        ctk.CTkEntry(
            self._at_params_frame, textvariable=self._variance_input_dir_var, width=280,
        ).grid(row=0, column=1, padx=6, pady=6, sticky="ew")
        ctk.CTkButton(
            self._at_params_frame, text="Browse", width=70, height=26,
            font=ctk.CTkFont(size=12),
            command=self._browse_variance_input_dir,
        ).grid(row=0, column=2, padx=12, pady=6, sticky="e")

        self._variance_cif_count_label = ctk.CTkLabel(
            self._at_params_frame, text="", anchor="w", font=ctk.CTkFont(size=11),
        )
        self._variance_cif_count_label.grid(row=1, column=1, columnspan=2, padx=6, pady=(0, 6), sticky="w")
        self._update_variance_cif_count()

        # Top N
        top_label_frame = ctk.CTkFrame(self._at_params_frame, fg_color="transparent")
        top_label_frame.grid(row=2, column=0, padx=(12, 6), pady=6, sticky="w")
        ctk.CTkLabel(top_label_frame, text="Top N residues:", anchor="w").pack(side="left")
        help_icon(top_label_frame, _CIF_VARIANCE_HELP["top_n"]).pack(side="left", padx=(4, 0))
        if not hasattr(self, "_variance_top_var"):
            self._variance_top_var = ctk.StringVar(value="10")
        ctk.CTkEntry(
            self._at_params_frame, textvariable=self._variance_top_var, width=60,
        ).grid(row=2, column=1, padx=6, pady=6, sticky="w")

        # Report range
        report_label_frame = ctk.CTkFrame(self._at_params_frame, fg_color="transparent")
        report_label_frame.grid(row=3, column=0, padx=(12, 6), pady=6, sticky="w")
        ctk.CTkLabel(report_label_frame, text="Report range:", anchor="w").pack(side="left")
        help_icon(report_label_frame, _CIF_VARIANCE_HELP["report_range"]).pack(side="left", padx=(4, 0))
        report_frame = ctk.CTkFrame(self._at_params_frame, fg_color="transparent")
        report_frame.grid(row=3, column=1, columnspan=2, padx=6, pady=6, sticky="w")
        if not hasattr(self, "_variance_range_start_var"):
            self._variance_range_start_var = ctk.StringVar(value="")
            self._variance_range_end_var = ctk.StringVar(value="")
        ctk.CTkEntry(report_frame, textvariable=self._variance_range_start_var, width=70,
                     placeholder_text="start").pack(side="left", padx=(0, 6))
        ctk.CTkEntry(report_frame, textvariable=self._variance_range_end_var, width=70,
                     placeholder_text="end (blank = all)").pack(side="left")

        # Align range
        align_label_frame = ctk.CTkFrame(self._at_params_frame, fg_color="transparent")
        align_label_frame.grid(row=4, column=0, padx=(12, 6), pady=6, sticky="w")
        ctk.CTkLabel(align_label_frame, text="Align range:", anchor="w").pack(side="left")
        help_icon(align_label_frame, _CIF_VARIANCE_HELP["align_range"]).pack(side="left", padx=(4, 0))
        align_frame = ctk.CTkFrame(self._at_params_frame, fg_color="transparent")
        align_frame.grid(row=4, column=1, columnspan=2, padx=6, pady=6, sticky="w")
        if not hasattr(self, "_variance_align_start_var"):
            self._variance_align_start_var = ctk.StringVar(value="")
            self._variance_align_end_var = ctk.StringVar(value="")
        ctk.CTkEntry(align_frame, textvariable=self._variance_align_start_var, width=70,
                     placeholder_text="start").pack(side="left", padx=(0, 6))
        ctk.CTkEntry(align_frame, textvariable=self._variance_align_end_var, width=70,
                     placeholder_text="end (blank = same as report range)").pack(side="left")

        # UniProt / gene overrides
        uniprot_label_frame = ctk.CTkFrame(self._at_params_frame, fg_color="transparent")
        uniprot_label_frame.grid(row=5, column=0, padx=(12, 6), pady=6, sticky="w")
        ctk.CTkLabel(uniprot_label_frame, text="UniProt override:", anchor="w").pack(side="left")
        help_icon(uniprot_label_frame, _CIF_VARIANCE_HELP["uniprot_override"]).pack(side="left", padx=(4, 0))
        if not hasattr(self, "_variance_uniprot_var"):
            self._variance_uniprot_var = ctk.StringVar(value="")
        ctk.CTkEntry(
            self._at_params_frame, textvariable=self._variance_uniprot_var, width=150,
            placeholder_text="auto-detected from CIF",
        ).grid(row=5, column=1, padx=6, pady=6, sticky="w")

        gene_label_frame = ctk.CTkFrame(self._at_params_frame, fg_color="transparent")
        gene_label_frame.grid(row=6, column=0, padx=(12, 6), pady=6, sticky="w")
        ctk.CTkLabel(gene_label_frame, text="Gene override:", anchor="w").pack(side="left")
        help_icon(gene_label_frame, _CIF_VARIANCE_HELP["gene_override"]).pack(side="left", padx=(4, 0))
        if not hasattr(self, "_variance_gene_var"):
            self._variance_gene_var = ctk.StringVar(value="")
        ctk.CTkEntry(
            self._at_params_frame, textvariable=self._variance_gene_var, width=150,
            placeholder_text="optional, for UniProt lookup",
        ).grid(row=6, column=1, padx=6, pady=(6, 10), sticky="w")

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
        self._at_status_label.configure(text="Sweep complete.", text_color=_GREEN)

    def _save_radius_sweep_plot(self) -> None:
        if not self._radius_sweep_fig.axes:
            self._at_status_label.configure(text="Run a sweep before saving.", text_color=_RED)
            return
        out_dir = self._output_dir
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
        self._style_dark_figure(self._cif_variance_fig)
        self._cif_variance_canvas.draw()
        self._show_active_analysis_plot()
        self._at_status_label.configure(text="Analysis complete.", text_color=_GREEN)

    def _save_cif_variance_plot(self) -> None:
        if not self._cif_variance_fig.axes:
            self._at_status_label.configure(text="Run an analysis before saving.", text_color=_RED)
            return
        out_dir = self._output_dir
        out_dir.mkdir(parents=True, exist_ok=True)
        out_path = out_dir / "cif_variance_plot.png"
        self._cif_variance_fig.savefig(
            out_path, dpi=200, facecolor=self._cif_variance_fig.get_facecolor(), bbox_inches="tight",
        )
        self._at_status_label.configure(text=f"Saved to {out_path}", text_color=_GREEN)
