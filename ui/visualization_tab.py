"""Visualization tab: lollipop (needle) plot of mutations near a PTM site.

Reads ResultsTabMixin's `_results_df_wide`/`_results_df_long` state directly.
`_style_lollipop_axis` is also reused by AnalysisToolsTabMixin.
"""
from __future__ import annotations

import re

import customtkinter as ctk

from ui.common import (
    _CIF_DIR,
    _DOMAIN_TYPE_COLORS,
    _DOMAIN_TYPE_FALLBACK_COLOR,
    _DOMAIN_TYPE_FALLBACK_LANE,
    _DOMAIN_TYPE_LANES,
    _GREEN,
    _MUT_ENTRY_RE,
    _NEEDLE_DEFAULT_COLOR,
    _PP_COLORS,
    _PP_LABEL,
    _PTM_MARKER_COLOR,
    _RED,
    _YELLOW,
    _load_interpro_entries,
    get_protein_length,
)


class VisualizationTabMixin:
    def _build_viz_tab(self, tab) -> None:
        import matplotlib
        matplotlib.use("TkAgg")
        from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
        from matplotlib.figure import Figure

        tab.grid_columnconfigure(0, weight=1)
        tab.grid_rowconfigure(3, weight=1)

        self._viz_ptm_rows: dict[str, int] = {}
        self._viz_all_labels: list[str] = []
        self._viz_protein_rows: dict[str, list[int]] = {}
        self._viz_all_protein_labels: list[str] = []
        self._viz_cluster_rows: dict[str, int] = {}
        self._viz_all_cluster_labels: list[str] = []
        self._viz_cluster_protein_rows: dict[str, list[int]] = {}
        self._viz_all_cluster_protein_labels: list[str] = []

        # ── Controls row 0: data source + view mode ──
        mode_row = ctk.CTkFrame(tab)
        mode_row.grid(row=0, column=0, padx=12, pady=(12, 6), sticky="ew")

        ctk.CTkLabel(mode_row, text="Data:", font=ctk.CTkFont(weight="bold")).pack(
            side="left", padx=(12, 6), pady=10,
        )
        self._viz_data_source_var = ctk.StringVar(value="PTM Proximity")
        ctk.CTkSegmentedButton(
            mode_row, values=["PTM Proximity", "Mutation Clusters"],
            variable=self._viz_data_source_var,
            command=self._on_viz_data_source_change,
        ).pack(side="left", padx=(0, 18), pady=10)

        ctk.CTkLabel(mode_row, text="View:", font=ctk.CTkFont(weight="bold")).pack(
            side="left", padx=(0, 6), pady=10,
        )
        self._viz_view_mode_var = ctk.StringVar(value="Single PTM")
        ctk.CTkSegmentedButton(
            mode_row, values=["Single PTM", "Whole protein"],
            variable=self._viz_view_mode_var,
            command=self._on_viz_view_mode_change,
        ).pack(side="left", padx=(0, 12), pady=10)

        # ── Controls row 1a: PTM/anchor selection (Single PTM mode) ──
        self._viz_ptm_controls = ctk.CTkFrame(tab)
        self._viz_ptm_controls.grid(row=1, column=0, padx=12, pady=(0, 6), sticky="ew")
        controls = self._viz_ptm_controls

        self._viz_anchor_label = ctk.CTkLabel(
            controls, text="PTM site:", font=ctk.CTkFont(weight="bold"),
        )
        self._viz_anchor_label.pack(side="left", padx=(12, 6), pady=10)

        self._viz_search_var = ctk.StringVar(value="")
        search_entry = ctk.CTkEntry(
            controls, textvariable=self._viz_search_var, width=160,
            placeholder_text="Search gene or site…",
        )
        search_entry.pack(side="left", padx=(0, 6), pady=10)
        self._viz_search_var.trace_add("write", self._on_viz_search)

        self._viz_combo = ctk.CTkComboBox(
            controls, width=220, values=[], command=lambda _v: self._generate_current_view(),
        )
        self._viz_combo.pack(side="left", padx=(0, 12), pady=10)

        ctk.CTkLabel(controls, text="±aa window:").pack(side="left", padx=(8, 4), pady=10)
        self._viz_window_var = ctk.StringVar(value="15")
        self._viz_window_var.trace_add("write", self._on_viz_window_changed)
        ctk.CTkEntry(
            controls, textvariable=self._viz_window_var, width=50,
        ).pack(side="left", padx=(0, 12), pady=10)

        # ── Controls row 1b: protein selection (Whole protein mode) ──
        self._viz_protein_controls = ctk.CTkFrame(tab)
        self._viz_protein_controls.grid(row=1, column=0, padx=12, pady=(0, 6), sticky="ew")
        self._viz_protein_controls.grid_remove()
        pcontrols = self._viz_protein_controls

        ctk.CTkLabel(
            pcontrols, text="Protein:", font=ctk.CTkFont(weight="bold"),
        ).pack(side="left", padx=(12, 6), pady=10)

        self._viz_protein_search_var = ctk.StringVar(value="")
        ctk.CTkEntry(
            pcontrols, textvariable=self._viz_protein_search_var, width=160,
            placeholder_text="Search gene or UniProt…",
        ).pack(side="left", padx=(0, 6), pady=10)
        self._viz_protein_search_var.trace_add("write", self._on_viz_protein_search)

        self._viz_protein_combo = ctk.CTkComboBox(
            pcontrols, width=220, values=[], command=lambda _v: self._generate_current_view(),
        )
        self._viz_protein_combo.pack(side="left", padx=(0, 12), pady=10)

        ctk.CTkLabel(pcontrols, text="±aa window:").pack(side="left", padx=(8, 4), pady=10)
        ctk.CTkEntry(
            pcontrols, textvariable=self._viz_window_var, width=50,
        ).pack(side="left", padx=(0, 12), pady=10)

        # ── Controls row 2: display mode + actions (shared) ──
        controls2 = ctk.CTkFrame(tab)
        controls2.grid(row=2, column=0, padx=12, pady=(0, 6), sticky="ew")

        ctk.CTkLabel(controls2, text="Show:", font=ctk.CTkFont(weight="bold")).pack(
            side="left", padx=(12, 6), pady=10,
        )
        self._viz_mode_var = ctk.StringVar(value="All mutations")
        ctk.CTkSegmentedButton(
            controls2, values=["All mutations", "Unique per position"],
            variable=self._viz_mode_var,
            command=lambda _v: self._generate_current_view(),
        ).pack(side="left", padx=(0, 12), pady=10)

        ctk.CTkButton(
            controls2, text="Save PNG", width=100,
            fg_color="gray30", hover_color="gray40",
            command=self._save_viz_plot,
        ).pack(side="left", padx=(0, 12), pady=10)

        self._viz_status = ctk.CTkLabel(
            controls2, text="Load results to select a PTM site.",
            text_color="gray60", font=ctk.CTkFont(size=11),
        )
        self._viz_status.pack(side="left", padx=(4, 12), pady=10)

        # ── Plot area: single-PTM canvas (row 3, Single PTM mode) ──
        self._viz_canvas_frame = ctk.CTkFrame(tab, fg_color="#2b2b2b")
        self._viz_canvas_frame.grid(row=3, column=0, padx=12, pady=(0, 12), sticky="nsew")
        self._viz_canvas_frame.grid_columnconfigure(0, weight=1)
        self._viz_canvas_frame.grid_rowconfigure(0, weight=1)
        # Prevents the canvas from resizing this frame back (avoids a resize loop).
        self._viz_canvas_frame.grid_propagate(False)

        self._viz_fig = Figure(figsize=(10, 6), dpi=100, facecolor="#2b2b2b")
        self._viz_canvas = FigureCanvasTkAgg(self._viz_fig, master=self._viz_canvas_frame)
        self._viz_canvas.get_tk_widget().grid(row=0, column=0, sticky="nsew")
        self._bind_viz_canvas_resize(self._viz_canvas_frame, self._viz_canvas, self._viz_fig)

        # ── Plot area: domain map + scrollable lollipop stack (row 3, Whole protein mode) ──
        self._build_viz_whole_protein_area(tab)

    def _on_viz_view_mode_change(self, _value: str = "") -> None:
        whole = self._viz_view_mode_var.get() == "Whole protein"
        if whole:
            self._viz_ptm_controls.grid_remove()
            self._viz_protein_controls.grid()
            self._viz_canvas_frame.grid_remove()
            self._viz_whole_protein_frame.grid()
            self._carry_ptm_selection_to_protein()
        else:
            self._viz_protein_controls.grid_remove()
            self._viz_ptm_controls.grid()
            self._viz_whole_protein_frame.grid_remove()
            self._viz_canvas_frame.grid()

    # ── Visualization: PTM Proximity / Mutation Clusters data source ────────

    def _viz_is_cluster(self) -> bool:
        return self._viz_data_source_var.get() == "Mutation Clusters"

    def _viz_active_wide_df(self):
        return self._cluster_df_wide if self._viz_is_cluster() else self._results_df_wide

    def _viz_active_rows(self) -> dict:
        return self._viz_cluster_rows if self._viz_is_cluster() else self._viz_ptm_rows

    def _viz_active_labels(self) -> list:
        return self._viz_all_cluster_labels if self._viz_is_cluster() else self._viz_all_labels

    def _viz_active_protein_rows(self) -> dict:
        return self._viz_cluster_protein_rows if self._viz_is_cluster() else self._viz_protein_rows

    def _viz_active_protein_labels(self) -> list:
        return (self._viz_all_cluster_protein_labels if self._viz_is_cluster()
                else self._viz_all_protein_labels)

    def _viz_active_anchor_col(self) -> str:
        return "anchor_mutation" if self._viz_is_cluster() else "ptm_site"

    def _on_viz_data_source_change(self, _value: str = "") -> None:
        cluster = self._viz_is_cluster()
        self._viz_anchor_label.configure(text="Anchor mutation:" if cluster else "PTM site:")
        self._viz_stack_legend_labels["domain_marker"].configure(
            text="▼ Anchor mutation (domain map)" if cluster else "▼ PTM position (domain map)",
        )
        self._viz_stack_legend_labels["site_marker"].configure(
            text="★ Anchor mutation" if cluster else "★ PTM site",
        )
        # Loads the target data source if not already loaded.
        self._load_results()
        self._viz_search_var.set("")
        labels = self._viz_active_labels()
        self._viz_combo.configure(values=labels)
        self._viz_combo.set(labels[0] if labels else "")
        self._viz_protein_search_var.set("")
        plabels = self._viz_active_protein_labels()
        self._viz_protein_combo.configure(values=plabels)
        self._viz_protein_combo.set(plabels[0] if plabels else "")
        self._generate_current_view()

    def _carry_ptm_selection_to_protein(self) -> None:
        """Default Whole Protein mode to the PTM/anchor's own protein and render it."""
        df = self._viz_active_wide_df()
        if df is None:
            return
        idx = self._viz_active_rows().get(self._viz_combo.get())
        if idx is None:
            return
        row = df.loc[idx]
        plabel = f"{row.get('gene', '?')} ({row.get('UniProt', '?')})"
        if plabel not in self._viz_active_protein_rows():
            return
        self._viz_protein_search_var.set("")
        self._viz_protein_combo.set(plabel)
        self._generate_whole_protein_view()

    def _generate_current_view(self) -> None:
        if self._viz_view_mode_var.get() == "Whole protein":
            self._generate_whole_protein_view()
        else:
            self._generate_lollipop_plot()

    def _on_viz_window_changed(self, *_args) -> None:
        """Regenerate the plot once the ±aa window value settles (debounced)."""
        if getattr(self, "_viz_window_after_id", None) is not None:
            self.after_cancel(self._viz_window_after_id)
        self._viz_window_after_id = self.after(600, self._generate_current_view)

    # ── Visualization: PTM / protein selection ──────────────────────────────

    def _refresh_viz_selector(self, df) -> None:
        """(Re)populate the PTM-site and protein selectors from loaded results."""
        self._viz_ptm_rows = {}
        labels: list[str] = []
        self._viz_protein_rows = {}
        protein_labels: list[str] = []
        for idx, row in df.iterrows():
            gene = row.get("gene", "?")
            site = row.get("ptm_site", "?")
            uid = row.get("UniProt", "?")
            label = f"{gene}  {site}  ({uid})"
            labels.append(label)
            self._viz_ptm_rows[label] = idx

            plabel = f"{gene} ({uid})"
            if plabel not in self._viz_protein_rows:
                protein_labels.append(plabel)
                self._viz_protein_rows[plabel] = []
            self._viz_protein_rows[plabel].append(idx)

        self._viz_all_labels = labels
        current = self._viz_combo.get()
        self._viz_combo.configure(values=labels)
        if labels and current not in labels:
            self._viz_combo.set(labels[0])
        elif not labels:
            self._viz_combo.set("")

        self._viz_all_protein_labels = protein_labels
        pcurrent = self._viz_protein_combo.get()
        self._viz_protein_combo.configure(values=protein_labels)
        if protein_labels and pcurrent not in protein_labels:
            self._viz_protein_combo.set(protein_labels[0])
        elif not protein_labels:
            self._viz_protein_combo.set("")

    def _refresh_cluster_viz_selector(self, df) -> None:
        """Cluster-mode counterpart of _refresh_viz_selector."""
        self._viz_cluster_rows = {}
        labels: list[str] = []
        self._viz_cluster_protein_rows = {}
        protein_labels: list[str] = []
        for idx, row in df.iterrows():
            gene = row.get("gene", "?")
            anchor = row.get("anchor_mutation", "?")
            uid = row.get("UniProt", "?")
            label = f"{gene}  {anchor}  ({uid})"
            labels.append(label)
            self._viz_cluster_rows[label] = idx

            plabel = f"{gene} ({uid})"
            if plabel not in self._viz_cluster_protein_rows:
                protein_labels.append(plabel)
                self._viz_cluster_protein_rows[plabel] = []
            self._viz_cluster_protein_rows[plabel].append(idx)

        self._viz_all_cluster_labels = labels
        self._viz_all_cluster_protein_labels = protein_labels
        if self._viz_is_cluster():
            current = self._viz_combo.get()
            self._viz_combo.configure(values=labels)
            if labels and current not in labels:
                self._viz_combo.set(labels[0])
            elif not labels:
                self._viz_combo.set("")

            pcurrent = self._viz_protein_combo.get()
            self._viz_protein_combo.configure(values=protein_labels)
            if protein_labels and pcurrent not in protein_labels:
                self._viz_protein_combo.set(protein_labels[0])
            elif not protein_labels:
                self._viz_protein_combo.set("")

    def _on_viz_search(self, *_args) -> None:
        query = self._viz_search_var.get().strip().lower()
        all_labels = self._viz_active_labels()
        filtered = (
            [label for label in all_labels if query in label.lower()]
            if query else all_labels
        )
        self._viz_combo.configure(values=filtered)
        if filtered and self._viz_combo.get() not in filtered:
            self._viz_combo.set(filtered[0])

    def _on_viz_protein_search(self, *_args) -> None:
        query = self._viz_protein_search_var.get().strip().lower()
        all_labels = self._viz_active_protein_labels()
        filtered = (
            [label for label in all_labels if query in label.lower()]
            if query else all_labels
        )
        self._viz_protein_combo.configure(values=filtered)
        if filtered and self._viz_protein_combo.get() not in filtered:
            self._viz_protein_combo.set(filtered[0])

    # ── Visualization: whole-protein plot area ──────────────────────────────

    def _build_viz_whole_protein_area(self, tab) -> None:
        from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
        from matplotlib.figure import Figure

        self._viz_whole_protein_frame = ctk.CTkFrame(tab, fg_color="transparent")
        self._viz_whole_protein_frame.grid(row=3, column=0, padx=12, pady=(0, 12), sticky="nsew")
        self._viz_whole_protein_frame.grid_remove()
        self._viz_whole_protein_frame.grid_columnconfigure(0, weight=1)
        self._viz_whole_protein_frame.grid_rowconfigure(1, weight=1)

        # Plain Tk title/legend, not matplotlib artists -- fig.suptitle()/
        # fig.legend() break down on a figure that can be 100+ inches tall.
        stack_header = ctk.CTkFrame(self._viz_whole_protein_frame, fg_color="transparent")
        stack_header.grid(row=0, column=0, sticky="ew", pady=(0, 4))

        self._viz_stack_title_label = ctk.CTkLabel(
            stack_header, text="", font=ctk.CTkFont(size=14, weight="bold"),
        )
        self._viz_stack_title_label.pack(side="left", padx=(4, 16))

        legend_frame = ctk.CTkFrame(stack_header, fg_color="transparent")
        legend_frame.pack(side="left")
        self._viz_stack_legend_labels: dict[str, ctk.CTkLabel] = {}
        for key, text, color in [
            ("domain_marker", "▼ PTM position (domain map)", _PTM_MARKER_COLOR),
            ("site_marker", "★ PTM site", _PTM_MARKER_COLOR),
            (None, "Benign", _PP_COLORS["benign"]),
            (None, "Possibly damaging", _PP_COLORS["possibly_damaging"]),
            (None, "Probably damaging", _PP_COLORS["probably_damaging"]),
            (None, "Unknown", _NEEDLE_DEFAULT_COLOR),
        ]:
            swatch = ctk.CTkFrame(legend_frame, fg_color=color, width=14, height=14, corner_radius=3)
            swatch.pack(side="left", padx=(10, 3))
            swatch.pack_propagate(False)
            label = ctk.CTkLabel(legend_frame, text=text, font=ctk.CTkFont(size=11))
            label.pack(side="left")
            if key:
                self._viz_stack_legend_labels[key] = label

        # ── Scrollable lollipop stack ──
        stack_outer = ctk.CTkFrame(self._viz_whole_protein_frame, fg_color="#2b2b2b")
        stack_outer.grid(row=1, column=0, sticky="nsew")
        stack_outer.grid_columnconfigure(0, weight=1)
        stack_outer.grid_rowconfigure(0, weight=1)
        stack_outer.grid_propagate(False)  # same reasoning as _viz_canvas_frame above
        self._bind_viz_stack_resize(stack_outer)

        scroll = ctk.CTkScrollableFrame(stack_outer, fg_color="transparent")
        scroll.grid(row=0, column=0, sticky="nsew")
        scroll.grid_columnconfigure(0, weight=1)

        # Width-clipping workaround (see analysis_tools_tab.py) so wide
        # content scrolls instead of being clipped.
        scroll._parent_canvas.unbind("<Configure>")

        def _fit_width_at_least_viewport(event):
            natural_w = scroll.winfo_reqwidth()
            scroll._parent_canvas.itemconfigure(
                scroll._create_window_id, width=max(natural_w, event.width),
            )

        scroll._parent_canvas.bind("<Configure>", _fit_width_at_least_viewport)

        stack_h_scrollbar = ctk.CTkScrollbar(
            scroll._parent_frame, orientation="horizontal", command=scroll._parent_canvas.xview,
        )
        stack_h_scrollbar.grid(row=2, column=0, sticky="ew")
        scroll._parent_canvas.configure(xscrollcommand=stack_h_scrollbar.set)

        self._viz_stack_scroll = scroll
        self._viz_stack_fig = Figure(figsize=(10, 2), dpi=100, facecolor="#2b2b2b")
        self._viz_stack_canvas = FigureCanvasTkAgg(self._viz_stack_fig, master=scroll)
        self._viz_stack_canvas.get_tk_widget().grid(row=0, column=0, sticky="nw")

        self._viz_stack_build_token = None

    def _current_viz_canvas_size_in(self) -> tuple[float, float]:
        """Current single-PTM canvas size in inches, so the figure fills the
        window. update_idletasks() first avoids reading a stale size right
        after a mode switch reveals this frame.
        """
        self.update_idletasks()
        dpi = self._viz_fig.get_dpi()
        w = self._viz_canvas_frame.winfo_width()
        h = self._viz_canvas_frame.winfo_height()
        if w <= 1 or h <= 1:
            return (10.0, 6.0)
        return (max(w / dpi, 6.0), max(h / dpi, 4.0))

    def _bind_viz_canvas_resize(self, frame, canvas, fig) -> None:
        """Keep *fig* sized to fill *frame* live as the window is resized.

        Binds on *frame*, not the canvas -- binding the canvas itself causes
        a resize feedback loop. Debounced to avoid redrawing on every pixel
        of a drag.
        """
        state = {"after_id": None}

        def _apply_resize(width_px: int, height_px: int) -> None:
            state["after_id"] = None
            if width_px < 100 or height_px < 100:
                return
            dpi = fig.get_dpi()
            new_w_in, new_h_in = width_px / dpi, height_px / dpi
            if (abs(new_w_in - fig.get_figwidth()) < 0.05
                    and abs(new_h_in - fig.get_figheight()) < 0.05):
                return
            fig.set_size_inches(new_w_in, new_h_in)
            if fig.axes:
                canvas.draw()

        def _on_configure(event) -> None:
            if state["after_id"] is not None:
                self.after_cancel(state["after_id"])
            state["after_id"] = self.after(150, lambda: _apply_resize(event.width, event.height))

        frame.bind("<Configure>", _on_configure)

    def _bind_viz_stack_resize(self, frame) -> None:
        """Regenerate the whole-protein stack on resize (full rebuild, not a
        resize-in-place, since layout depends on width at build time).
        """
        state = {"after_id": None, "last_width": None}

        def _apply(width_px: int) -> None:
            state["after_id"] = None
            if width_px < 200:
                return
            if state["last_width"] is not None and abs(width_px - state["last_width"]) < 20:
                return
            state["last_width"] = width_px
            if self._viz_view_mode_var.get() == "Whole protein" and self._viz_stack_title_label.cget("text"):
                self._generate_whole_protein_view()

        def _on_configure(event) -> None:
            if state["after_id"] is not None:
                self.after_cancel(state["after_id"])
            state["after_id"] = self.after(500, lambda: _apply(event.width))

        frame.bind("<Configure>", _on_configure)

    # ── Visualization: data prep ─────────────────────────────────────────────

    def _get_viz_mutation_df(self, uid: str, ptm_site: str, wide_row):
        """Return (DataFrame, note) of nearby mutations for a PTM site.

        Prefers the long-format table (has per-mutation patient counts);
        falls back to parsing the wide-format summary columns otherwise.
        """
        import pandas as pd

        if self._results_df_long is not None:
            df = self._results_df_long
            mask = (
                (df.get("uniprot_id", "") == uid) &
                (df.get("ptm_position", "") == ptm_site)
            )
            sub = df[mask].copy()
            if not sub.empty:
                def _pos(m):
                    match = re.search(r"\d+", str(m))
                    return int(match.group()) if match else None

                sub["mutation_position"] = sub["mutation"].apply(_pos)
                sub = sub.dropna(subset=["mutation_position"])
                sub["patient_count"] = pd.to_numeric(
                    sub.get("patient_count", 0), errors="coerce"
                ).fillna(0)
                sub["distance_angstrom"] = pd.to_numeric(
                    sub.get("distance_angstrom", 0), errors="coerce"
                ).fillna(0)
                if "polyphen_class" not in sub.columns:
                    sub["polyphen_class"] = ""
                return sub[[
                    "mutation", "mutation_position", "patient_count",
                    "distance_angstrom", "polyphen_class",
                ]], None

        # Fallback: wide-format columns only (no per-mutation patient counts)
        rows = []
        for col_key in ("mutations_within_5_positions", "mutations_more_than_5_positions"):
            for entry in (wide_row.get(col_key, "") or "").split(", "):
                entry = entry.strip()
                if not entry:
                    continue
                m = _MUT_ENTRY_RE.match(entry)
                if not m:
                    continue
                pos_m = re.search(r"\d+", m.group(1))
                if not pos_m:
                    continue
                rows.append({
                    "mutation": m.group(1),
                    "mutation_position": int(pos_m.group()),
                    "patient_count": 1,
                    "distance_angstrom": float(m.group(4)),
                    "polyphen_class": _PP_LABEL.get(m.group(2) or "", ""),
                })
        df_fallback = pd.DataFrame(rows)
        note = None
        if not df_fallback.empty:
            note = "patient counts unavailable — re-run the pipeline to generate the long-format table"
        return df_fallback, note

    def _get_viz_cluster_mutation_df(self, uid: str, anchor_mutation: str, wide_row):
        """Cluster-mode counterpart of _get_viz_mutation_df."""
        import pandas as pd

        if self._cluster_df_long is not None:
            df = self._cluster_df_long
            mask = (
                (df.get("UniProt", "") == uid) &
                (df.get("anchor_mutation", "") == anchor_mutation)
            )
            sub = df[mask].copy()
            if not sub.empty:
                sub["mutation_position"] = pd.to_numeric(
                    sub.get("mutation_position", 0), errors="coerce"
                )
                sub = sub.dropna(subset=["mutation_position"])
                sub["patient_count"] = pd.to_numeric(
                    sub.get("patient_count", 0), errors="coerce"
                ).fillna(0)
                sub["distance_angstrom"] = pd.to_numeric(
                    sub.get("distance_angstrom", 0), errors="coerce"
                ).fillna(0)
                if "polyphen_class" not in sub.columns:
                    sub["polyphen_class"] = ""
                return sub[[
                    "mutation", "mutation_position", "patient_count",
                    "distance_angstrom", "polyphen_class",
                ]], None

        # Fallback: wide-format "nearby_mutations" column (no patient counts).
        rows = []
        for entry in (wide_row.get("nearby_mutations", "") or "").split(", "):
            entry = entry.strip()
            if not entry:
                continue
            m = _MUT_ENTRY_RE.match(entry)
            if not m:
                continue
            pos_m = re.search(r"\d+", m.group(1))
            if not pos_m:
                continue
            rows.append({
                "mutation": m.group(1),
                "mutation_position": int(pos_m.group()),
                "patient_count": 1,
                "distance_angstrom": float(m.group(4)),
                "polyphen_class": _PP_LABEL.get(m.group(2) or "", ""),
            })
        df_fallback = pd.DataFrame(rows)
        note = None
        if not df_fallback.empty:
            note = "patient counts unavailable — re-run Mutation Clustering mode to generate the long-format table"
        return df_fallback, note

    # ── Visualization: plotting ──────────────────────────────────────────────

    def _style_lollipop_axis(self, ax) -> None:
        ax.set_facecolor("#2b2b2b")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["left"].set_color("#555555")
        ax.spines["bottom"].set_color("#555555")
        ax.tick_params(colors="#dcdcdc")
        ax.yaxis.label.set_color("#dcdcdc")
        ax.xaxis.label.set_color("#dcdcdc")

    _MUT_LABEL_FONTSIZE = 8
    _MUT_LABEL_BASE_OFFSET_PTS = (0, 10)  # directly above the marker (label is rotated fully vertical)
    _MUT_LABEL_PUSH_STEP_PTS = 8          # increment used to push a colliding label clear
    _MUT_LABEL_PAD_PX = 3                 # breathing room added around each label's box

    def _layout_mutation_label_offsets(self, ax, mutations: list[str],
                                        x_plot: list[float], y_plot: list[float]) -> dict[int, float]:
        """Push mutation labels straight up, one at a time, until each one's
        actual rendered footprint (in real screen pixels) no longer overlaps
        any label already placed.

        Needed because a label's on-screen position depends on *both* its
        residue position (x) and its marker height (y = patient count) --
        two mutations a residue apart can still land far apart vertically if
        their patient counts differ, while two at very different positions
        can still collide if one has a tall needle. A purely x-order row
        assignment (as used for the domain map, where every label sits on
        the same baseline) can't account for that, so this checks true
        pixel-space bounding boxes instead. Returns {df_row_index: extra
        y-offset in points} to add on top of _MUT_LABEL_BASE_OFFSET_PTS;
        entries with no collision map to 0.0, reproducing the original fixed
        offset.
        """
        if not mutations:
            return {}

        fig = ax.figure
        px_per_pt = fig.dpi / 72.0
        base_x_px, base_y_px = (v * px_per_pt for v in self._MUT_LABEL_BASE_OFFSET_PTS)
        push_step_px = self._MUT_LABEL_PUSH_STEP_PTS * px_per_pt
        pad = self._MUT_LABEL_PAD_PX

        px_lengths = self._measure_text_widths_px(mutations, self._MUT_LABEL_FONTSIZE)
        text_h_px = self._MUT_LABEL_FONTSIZE * px_per_pt * 1.2  # approx line height incl. leading
        # Labels are rotated fully vertical (90 degrees): the string's own
        # rendered length becomes its vertical footprint, and the font's
        # line height becomes its (much narrower) horizontal footprint.
        box_sizes = [(text_h_px + pad, length + pad) for length in px_lengths]

        anchors_px = ax.transData.transform(list(zip(x_plot, y_plot)))

        order = sorted(range(len(mutations)), key=lambda i: x_plot[i])
        placed: list[tuple[float, float, float, float]] = []
        offsets: dict[int, float] = {}
        for i in order:
            ax_px, ay_px = anchors_px[i]
            w, h = box_sizes[i]
            x0 = ax_px + base_x_px - w / 2  # centered horizontally on the marker
            extra_px = 0.0
            while True:
                y0 = ay_px + base_y_px + extra_px
                box = (x0, y0, x0 + w, y0 + h)
                conflict = any(
                    box[0] < o[2] and box[2] > o[0] and box[1] < o[3] and box[3] > o[1]
                    for o in placed
                )
                if not conflict:
                    placed.append(box)
                    break
                extra_px += push_step_px
            offsets[i] = extra_px / px_per_pt
        return offsets

    def _draw_mutation_needles(self, ax, df, jitter: bool, x_col: str = "mutation_position") -> None:
        from collections import defaultdict

        positions = df[x_col].tolist()
        groups: dict = defaultdict(list)
        for i, pos in enumerate(positions):
            groups[round(pos)].append(i)

        x_plot = list(positions)
        if jitter:
            for pos, idxs in groups.items():
                n = len(idxs)
                if n > 1:
                    for k, i in enumerate(idxs):
                        x_plot[i] = pos + (k - (n - 1) / 2) * 0.6

        label_offsets = (
            self._layout_mutation_label_offsets(
                ax, df["mutation"].astype(str).tolist(), x_plot, df["patient_count"].tolist(),
            ) if jitter else {}
        )

        for i, (_, r) in enumerate(df.iterrows()):
            x = x_plot[i]
            count = r["patient_count"]
            color = _PP_COLORS.get(str(r.get("polyphen_class", "")).strip(), _NEEDLE_DEFAULT_COLOR)
            ax.vlines(x, 0, count, color=color, linewidth=2, zorder=2)
            ax.scatter([x], [count], s=90, color=color, edgecolor="black", zorder=3)
            # The far (categorical) panel already names each mutation via its x-tick
            # label, so only annotate above the marker for the true-scale local panel.
            if jitter:
                base_x, base_y = self._MUT_LABEL_BASE_OFFSET_PTS
                ax.annotate(
                    str(r["mutation"]), (x, count),
                    xytext=(base_x, base_y + label_offsets.get(i, 0.0)),
                    textcoords="offset points", ha="center", va="bottom",
                    fontsize=self._MUT_LABEL_FONTSIZE, color="#dcdcdc", rotation=90,
                )

    _DOMAIN_LABEL_FONTSIZE = 7
    _DOMAIN_LANE_HEIGHT_IN = 0.16  # physical height of one lane row (boxes)
    _DOMAIN_BACKBONE_MARGIN_IN = 0.12  # space below the lowest lane, for the backbone/PTM marker
    _DOMAIN_LABEL_TOP_PADDING_IN = 0.14  # gap between the top lane and the first label row
    _DOMAIN_LABEL_ROW_INCHES = 0.26  # physical row pitch -- see _layout_domain_labels
    _DOMAIN_STRIP_MARGIN_FRACTION = 0.89  # matches _draw_lollipop_group's left=0.08/right=0.97

    def _measure_text_widths_px(self, texts: list[str], fontsize: float) -> list[float]:
        """Real rendered pixel widths for *texts*, via a scratch figure."""
        from matplotlib.backends.backend_agg import FigureCanvasAgg
        from matplotlib.figure import Figure as _Figure

        scratch_fig = _Figure(dpi=100)
        canvas = FigureCanvasAgg(scratch_fig)
        ax = scratch_fig.add_subplot(111)
        renderer = canvas.get_renderer()
        widths = []
        for text in texts:
            t = ax.text(0, 0, text, fontsize=fontsize, alpha=0)
            widths.append(t.get_window_extent(renderer=renderer).width)
            t.remove()
        return widths

    _DOMAIN_ARROW_FAN_PX = 6  # min horizontal separation between fanned-out arrow tips

    def _layout_domain_labels(self, entries: list[dict], length: float,
                               width_in: float) -> tuple[list[tuple[dict, int, float]], int]:
        """Assign each label to the lowest row where it doesn't overlap a
        neighbor's text or arrow anchor (never shifts a label sideways).
        Returns (assignments, n_rows_used); each assignment is
        (entry, row, arrow_x_offset).
        """
        if not entries:
            return [], 0

        axis_width_px = max(width_in * self._DOMAIN_STRIP_MARGIN_FRACTION * 100, 1.0)  # dpi=100
        data_units_per_px = length / axis_width_px
        gap = 10 * data_units_per_px
        min_arrow_gap = 24 * data_units_per_px

        sorted_entries = sorted(entries, key=lambda e: e["start"])
        px_widths = self._measure_text_widths_px(
            [e["name"] for e in sorted_entries], self._DOMAIN_LABEL_FONTSIZE,
        )

        row_occupied: list[list[tuple[float, float, float]]] = []
        assignments = []
        for e, px_w in zip(sorted_entries, px_widths):
            label_w = px_w * data_units_per_px
            center_x = e["start"] + (e["end"] - e["start"]) / 2
            label_start, label_end = center_x - label_w / 2, center_x + label_w / 2

            row = 0
            while True:
                if row == len(row_occupied):
                    row_occupied.append([])
                conflict = any(
                    (label_start - gap < o_end and label_end + gap > o_start)
                    or abs(center_x - o_center) < min_arrow_gap
                    for o_start, o_end, o_center in row_occupied[row]
                )
                if not conflict:
                    row_occupied[row].append((label_start, label_end, center_x))
                    assignments.append((e, row))
                    break
                row += 1

        # Entries whose box positions are so close together that their
        # (otherwise perfectly vertical) arrows would land on top of each
        # other -- stacking them across rows already keeps their *text* from
        # overlapping, but the arrow lines themselves still converge on
        # nearly the same point, making it impossible to tell which arrow
        # belongs to which label. Fan each cluster's arrow tips out by a
        # small, fixed pixel step so every arrow traces a visually distinct
        # path; chained via nearest-neighbor distance, same threshold used
        # above to decide these arrows were too close in the first place.
        fan_step = self._DOMAIN_ARROW_FAN_PX * data_units_per_px
        by_center = sorted(assignments, key=lambda a: a[0]["start"] + (a[0]["end"] - a[0]["start"]) / 2)
        offsets: dict[int, float] = {}
        cluster: list[int] = []
        prev_center = None
        for i, (e, _row) in enumerate(by_center):
            center_x = e["start"] + (e["end"] - e["start"]) / 2
            if prev_center is not None and center_x - prev_center >= min_arrow_gap and cluster:
                for k, j in enumerate(cluster):
                    offsets[j] = (k - (len(cluster) - 1) / 2) * fan_step
                cluster = []
            cluster.append(i)
            prev_center = center_x
        if cluster:
            for k, j in enumerate(cluster):
                offsets[j] = (k - (len(cluster) - 1) / 2) * fan_step

        assignments = [(e, row, offsets.get(i, 0.0)) for i, (e, row) in enumerate(by_center)]
        return assignments, len(row_occupied)

    def _domain_strip_height_in(self, n_label_rows: int) -> float:
        """Physical height (inches) needed for the domain strip's backbone,
        boxes, and label rows. 1 data-unit-of-y == 1 real inch.
        """
        n_lanes = max(_DOMAIN_TYPE_LANES.values(), default=0) + 1
        return (self._DOMAIN_BACKBONE_MARGIN_IN + n_lanes * self._DOMAIN_LANE_HEIGHT_IN
                + self._DOMAIN_LABEL_TOP_PADDING_IN + max(n_label_rows, 1) * self._DOMAIN_LABEL_ROW_INCHES)

    def _draw_domain_strip(self, ax, entries: list[dict], length: float, ptm_pos: int,
                            domain_layout: tuple[list[tuple[dict, int, float]], int]) -> None:
        """Draw a linear domain-map strip on *ax*: backbone, InterPro boxes
        by lane, labels with vertical arrows, and a marker at *ptm_pos*.
        Y-axis is inches (see _domain_strip_height_in). *domain_layout* is
        pre-computed once per protein, not per row.
        """
        from matplotlib.patches import Rectangle

        self._style_lollipop_axis(ax)

        assignments, n_rows = domain_layout
        strip_height_in = self._domain_strip_height_in(n_rows)
        n_lanes = max(_DOMAIN_TYPE_LANES.values(), default=0) + 1
        backbone_y = 0.0
        lanes_top_y = backbone_y + self._DOMAIN_BACKBONE_MARGIN_IN + n_lanes * self._DOMAIN_LANE_HEIGHT_IN

        ax.hlines(backbone_y, 0, length, color=_NEEDLE_DEFAULT_COLOR, linewidth=2.5, zorder=1)

        def _box_position(e: dict) -> tuple[float, float, float]:
            lane = _DOMAIN_TYPE_LANES.get(e["type"], _DOMAIN_TYPE_FALLBACK_LANE)
            y0 = backbone_y + self._DOMAIN_BACKBONE_MARGIN_IN + lane * self._DOMAIN_LANE_HEIGHT_IN
            return y0, y0 + self._DOMAIN_LANE_HEIGHT_IN * 0.8, e["start"] + (e["end"] - e["start"]) / 2

        for e in entries:
            y0, _box_top, _center = _box_position(e)
            color = _DOMAIN_TYPE_COLORS.get(e["type"], _DOMAIN_TYPE_FALLBACK_COLOR)
            ax.add_patch(Rectangle(
                (e["start"], y0), e["end"] - e["start"], self._DOMAIN_LANE_HEIGHT_IN * 0.8,
                facecolor=color, edgecolor="black", linewidth=0.4, alpha=0.85, zorder=2,
            ))

        for e, row, arrow_x_offset in assignments:
            _y0, box_top_y, box_center_x = _box_position(e)
            label_y = lanes_top_y + self._DOMAIN_LABEL_TOP_PADDING_IN + row * self._DOMAIN_LABEL_ROW_INCHES
            # The label itself stays centered on the true box position (text
            # collisions are already resolved by row alone); only the arrow's
            # tip is nudged by arrow_x_offset, so tightly-clustered arrows fan
            # out into distinguishable paths instead of overlapping straight
            # vertical lines.
            ax.annotate(
                e["name"], xy=(box_center_x + arrow_x_offset, box_top_y), xytext=(box_center_x, label_y),
                fontsize=self._DOMAIN_LABEL_FONTSIZE, color="#dcdcdc", ha="center", va="bottom",
                arrowprops=dict(arrowstyle="->", color="#888888", lw=0.7, shrinkA=0, shrinkB=2),
                clip_on=False, zorder=3,
            )

        ax.plot(
            ptm_pos, backbone_y, marker="v", markersize=8,
            color=_PTM_MARKER_COLOR, markeredgecolor="white", markeredgewidth=0.6, zorder=4,
        )

        ax.set_xlim(-length * 0.02, length * 1.02)
        ax.set_ylim(-self._DOMAIN_BACKBONE_MARGIN_IN, strip_height_in - self._DOMAIN_BACKBONE_MARGIN_IN)
        ax.set_yticks([])
        ax.tick_params(labelsize=7)
        # No x-label here -- the lollipop panel below already has one on
        # the same residue-number scale.
        for spine in ("top", "left", "right"):
            ax.spines[spine].set_visible(False)

        if not entries:
            ax.annotate(
                "No InterPro domain data available for this protein", (0.5, 0.5),
                xycoords="axes fraction", ha="center", va="center",
                color="#888888", fontsize=7,
            )

    def _draw_lollipop_group(self, fig, subplotspec, gene: str, ptm_site: str, ptm_pos: int,
                              mut_df, local_window: float, domain_entries: list[dict],
                              length: float, domain_layout: tuple[list[tuple[dict, int, float]], int],
                              available_height_in: float, mutation_count: int, unique_count: int):
        """Draw one PTM site's domain strip + local/far lollipop panels into *fig*.

        *subplotspec*=None lays out on the whole figure (single-PTM view);
        otherwise nests via subgridspec (whole-protein view) -- nested
        gridspecs can't take margin kwargs. *domain_entries*/*length*/
        *domain_layout* are pre-computed once per protein. *available_height_in*
        splits domain-strip vs. lollipop height. *mutation_count*/*unique_count*
        are the pre-collapse totals for this cluster (mut_df itself may already
        be collapsed to one row per position by the "Unique per position" view
        mode, so they can't be re-derived from mut_df here). Returns
        (ax_local, domain_ax).
        """
        local_df = mut_df[(mut_df["mutation_position"] - ptm_pos).abs() <= local_window] \
            .sort_values("mutation_position")
        far_df = mut_df[(mut_df["mutation_position"] - ptm_pos).abs() > local_window] \
            .sort_values("mutation_position").reset_index(drop=True)
        has_far = not far_df.empty

        max_count = max(float(mut_df["patient_count"].max()), 1.0)
        ptm_height = max_count * 1.15
        y_top = ptm_height * 1.45

        n_label_rows = domain_layout[1]
        domain_h = min(self._domain_strip_height_in(n_label_rows), available_height_in * 0.7)
        lollipop_h = max(available_height_in - domain_h, 0.5)

        if subplotspec is None:
            outer_gs = fig.add_gridspec(
                2, 1, height_ratios=[domain_h, lollipop_h], hspace=0.6,
                left=0.08, right=0.97, top=0.86, bottom=0.14,
            )
        else:
            outer_gs = subplotspec.subgridspec(2, 1, height_ratios=[domain_h, lollipop_h], hspace=0.6)

        domain_ax = fig.add_subplot(outer_gs[0])
        self._draw_domain_strip(domain_ax, domain_entries, length, ptm_pos, domain_layout)

        lollipop_cell = outer_gs[1]
        gs = (lollipop_cell.subgridspec(1, 3, width_ratios=[0.68, 0.05, 0.27], wspace=0.05)
              if has_far else lollipop_cell.subgridspec(1, 1))

        if has_far:
            ax_local = fig.add_subplot(gs[0])
            ax_gap = fig.add_subplot(gs[1])
            ax_far = fig.add_subplot(gs[2])
        else:
            ax_local = fig.add_subplot(gs[0])
            ax_gap = None
            ax_far = None

        # ── Local panel (true sequence scale) ──
        self._style_lollipop_axis(ax_local)
        x_min, x_max = ptm_pos - local_window, ptm_pos + local_window
        ax_local.hlines(0, x_min, x_max, color=_NEEDLE_DEFAULT_COLOR, linewidth=2, zorder=1)
        ax_local.set_xlim(x_min - 0.5, x_max + 0.5)
        ax_local.set_ylim(0, y_top)
        ax_local.set_ylabel("Patient count")
        ax_local.set_xlabel("Residue position")
        if has_far:
            ax_local.spines["right"].set_visible(False)

        ax_local.vlines(ptm_pos, 0, ptm_height, color=_PTM_MARKER_COLOR, linewidth=2, zorder=2)
        ax_local.scatter(
            [ptm_pos], [ptm_height], marker="*", s=420,
            color=_PTM_MARKER_COLOR, edgecolor="white", linewidth=0.8, zorder=4,
        )
        ax_local.annotate(
            str(ptm_site), (ptm_pos, ptm_height), xytext=(0, 10),
            textcoords="offset points", ha="center", color="white",
            fontsize=11, fontweight="bold",
        )

        self._draw_mutation_needles(ax_local, local_df, jitter=True)

        # ── Broken-axis gap + far panel (categorical scale) ──
        if has_far:
            ax_gap.axis("off")
            ax_gap.set_xlim(0, 1)
            ax_gap.set_ylim(0, 1)
            for xc in (0.3, 0.7):
                ax_gap.plot([xc - 0.12, xc + 0.12], [0.42, 0.52],
                            color="#888888", lw=1.2, transform=ax_gap.transAxes, clip_on=False)
                ax_gap.plot([xc - 0.12, xc + 0.12], [0.48, 0.58],
                            color="#888888", lw=1.2, transform=ax_gap.transAxes, clip_on=False)

            self._style_lollipop_axis(ax_far)
            far_df["_x"] = range(len(far_df))
            n_far = len(far_df)
            ax_far.hlines(0, -0.5, n_far - 0.5, color=_NEEDLE_DEFAULT_COLOR, linewidth=2, zorder=1)
            ax_far.set_xlim(-0.5, n_far - 0.5)
            ax_far.set_ylim(0, y_top)
            ax_far.set_yticklabels([])
            ax_far.spines["left"].set_visible(False)

            self._draw_mutation_needles(ax_far, far_df, jitter=False, x_col="_x")
            ax_far.set_xticks(list(range(n_far)))
            ax_far.set_xticklabels(far_df["mutation"], rotation=45, ha="right", fontsize=8)
            ax_far.set_title(f"> {local_window:g} aa away", fontsize=9, color="gray")

        # Quick-reference counts, anchored to the figure (not an axes) in the
        # left margin, vertically centered on the domain strip specifically
        # (not the whole domain+lollipop row). domain_ax has no y-ticks/ylabel
        # of its own and a guaranteed minimum height regardless of window
        # size, so this placement can't collide with anything -- unlike
        # centering over the full row, which put it in ax_local's territory
        # on short/sparse-domain renders, where ax_local's "Patient count"
        # ylabel (a fixed-length rotated string) overflows past its own tiny
        # axes bounds and collides with it.
        domain_pos = domain_ax.get_position()
        left_x = domain_pos.x0
        center_y = (domain_pos.y0 + domain_pos.y1) / 2
        fig.text(
            left_x - 0.01, center_y,
            f"{mutation_count} mutation{'s' if mutation_count != 1 else ''}\n"
            f"{unique_count} unique",
            transform=fig.transFigure, ha="right", va="center",
            fontsize=8, color="#dcdcdc", linespacing=1.7,
        )

        return ax_local, domain_ax

    def _domain_map_legend_handles(self):
        from matplotlib.lines import Line2D

        label = "Anchor mutation" if self._viz_is_cluster() else "PTM position"
        return [
            Line2D([0], [0], marker="v", color="none", markerfacecolor=_PTM_MARKER_COLOR,
                   markeredgecolor="white", markersize=10, label=label),
        ]

    def _lollipop_legend_handles(self):
        from matplotlib.lines import Line2D
        from matplotlib.patches import Patch

        marker_label = "Anchor mutation" if self._viz_is_cluster() else "PTM site"
        return [
            Line2D([0], [0], marker="*", color="none", markerfacecolor=_PTM_MARKER_COLOR,
                   markeredgecolor="white", markersize=14, label=marker_label),
            Patch(facecolor=_PP_COLORS["benign"], label="Benign"),
            Patch(facecolor=_PP_COLORS["possibly_damaging"], label="Possibly damaging"),
            Patch(facecolor=_PP_COLORS["probably_damaging"], label="Probably damaging"),
            Patch(facecolor=_NEEDLE_DEFAULT_COLOR, label="Unknown / not scored"),
        ]

    def _draw_lollipop(self, gene: str, ptm_site: str, ptm_pos: int, mut_df, local_window: float,
                        uid: str, mutation_count: int, unique_count: int) -> None:
        fig = self._viz_fig
        fig.clf()
        fig.patch.set_facecolor("#2b2b2b")
        fig.set_size_inches(*self._current_viz_canvas_size_in())

        domain_entries = _load_interpro_entries(uid)
        protein_length = get_protein_length(_CIF_DIR / uid)
        max_end = max((e["end"] for e in domain_entries), default=0)
        length = max(protein_length or 0, max_end, ptm_pos, 1)
        domain_layout = self._layout_domain_labels(domain_entries, length, fig.get_figwidth())
        available_height_in = fig.get_figheight() * (0.86 - 0.14)

        ax_local, domain_ax = self._draw_lollipop_group(
            fig, None, gene, ptm_site, ptm_pos, mut_df, local_window,
            domain_entries, length, domain_layout, available_height_in,
            mutation_count, unique_count,
        )

        ax_local.legend(
            handles=self._lollipop_legend_handles(), loc="upper left", fontsize=8,
            facecolor="#3a3a3a", edgecolor="#555555", labelcolor="#dcdcdc",
        )
        # Anchored to the figure (not domain_ax) so it can't overlap domain
        # boxes/labels near position 0.
        fig.legend(
            handles=self._domain_map_legend_handles(), loc="lower left",
            bbox_to_anchor=(0.08, 0.87), fontsize=7,
            facecolor="#3a3a3a", edgecolor="#555555", labelcolor="#dcdcdc",
        )

        fig.suptitle(f"{gene} — {ptm_site} and nearby mutations", color="white",
                     fontsize=14, fontweight="bold")

        self._viz_canvas.draw()

    def _generate_lollipop_plot(self) -> None:
        cluster = self._viz_is_cluster()
        df_wide = self._viz_active_wide_df()
        if df_wide is None:
            missing = "Mutation Clustering" if cluster else "PTM Proximity"
            self._viz_status.configure(
                text=f"No {missing} results loaded — run that pipeline mode or open the Results tab first.",
                text_color=_RED,
            )
            return

        label = self._viz_combo.get()
        idx = self._viz_active_rows().get(label)
        if idx is None:
            what = "an anchor mutation" if cluster else "a PTM site"
            self._viz_status.configure(text=f"Select {what} to plot.", text_color=_RED)
            return

        row = df_wide.loc[idx]
        gene = row.get("gene", "?")
        ptm_site = row.get(self._viz_active_anchor_col(), "?")
        uid = row.get("UniProt", "?")

        ptm_m = re.search(r"\d+", str(ptm_site))
        if not ptm_m:
            self._viz_status.configure(
                text=f"Could not parse a residue position from '{ptm_site}'.", text_color=_RED,
            )
            return
        ptm_pos = int(ptm_m.group())

        mut_df, note = (
            self._get_viz_cluster_mutation_df(uid, ptm_site, row) if cluster
            else self._get_viz_mutation_df(uid, ptm_site, row)
        )
        if mut_df.empty:
            self._viz_fig.clf()
            self._viz_canvas.draw()
            self._viz_status.configure(
                text=f"No nearby mutations found for {gene} {ptm_site}.", text_color=_YELLOW,
            )
            return

        mutation_count = len(mut_df)
        unique_count = mut_df["mutation_position"].nunique()

        unique_only = self._viz_mode_var.get() == "Unique per position"
        if unique_only:
            mut_df = self._collapse_to_unique_positions(mut_df)

        try:
            local_window = float(self._viz_window_var.get().strip() or "15")
        except ValueError:
            local_window = 15.0

        self._draw_lollipop(gene, ptm_site, ptm_pos, mut_df, local_window, uid,
                             mutation_count, unique_count)

        kind = "unique position(s)" if unique_only else "nearby mutation(s)"
        status = f"{len(mut_df)} {kind} for {gene} {ptm_site}"
        if note:
            status += f" — {note}"
        self._viz_status.configure(text=status, text_color="gray60")

    def _collapse_to_unique_positions(self, df):
        """Merge same-position mutations into one entry (summed patient count,
        worst PolyPhen class for color).
        """
        import pandas as pd

        _SEVERITY = {"probably_damaging": 3, "possibly_damaging": 2, "benign": 1}

        def _worst_class(classes) -> str:
            best, best_rank = "", -1
            for c in classes:
                c = str(c).strip()
                rank = _SEVERITY.get(c, 0)
                if rank > best_rank:
                    best, best_rank = c, rank
            return best

        def _position_label(mutation: str) -> str:
            m = re.match(r"^([A-Za-z]+\d+)", str(mutation))
            return m.group(1) if m else str(mutation)

        rows = []
        for pos, grp in df.groupby("mutation_position"):
            label = _position_label(grp["mutation"].iloc[0])
            if len(grp) > 1:
                label = f"{label} (×{len(grp)})"
            rows.append({
                "mutation": label,
                "mutation_position": pos,
                "patient_count": grp["patient_count"].sum(),
                "distance_angstrom": grp["distance_angstrom"].min(),
                "polyphen_class": _worst_class(grp["polyphen_class"]),
            })
        return pd.DataFrame(rows)

    # ── Visualization: whole-protein view ────────────────────────────────────

    def _generate_whole_protein_view(self) -> None:
        cluster = self._viz_is_cluster()
        df_wide = self._viz_active_wide_df()
        if df_wide is None:
            missing = "Mutation Clustering" if cluster else "PTM Proximity"
            self._viz_status.configure(
                text=f"No {missing} results loaded — run that pipeline mode or open the Results tab first.",
                text_color=_RED,
            )
            return

        label = self._viz_protein_combo.get()
        indices = self._viz_active_protein_rows().get(label)
        if not indices:
            self._viz_status.configure(text="Select a protein to plot.", text_color=_RED)
            return

        rows = df_wide.loc[indices]
        first = rows.iloc[0]
        gene = first.get("gene", "?")
        uid = first.get("UniProt", "?")
        anchor_col = self._viz_active_anchor_col()

        try:
            local_window = float(self._viz_window_var.get().strip() or "15")
        except ValueError:
            local_window = 15.0
        unique_only = self._viz_mode_var.get() == "Unique per position"

        ptm_entries = []
        for _, row in rows.iterrows():
            ptm_site = row.get(anchor_col, "?")
            ptm_m = re.search(r"\d+", str(ptm_site))
            if not ptm_m:
                continue
            ptm_pos = int(ptm_m.group())
            mut_df, _note = (
                self._get_viz_cluster_mutation_df(uid, ptm_site, row) if cluster
                else self._get_viz_mutation_df(uid, ptm_site, row)
            )
            if mut_df.empty:
                continue
            mutation_count = len(mut_df)
            unique_count = mut_df["mutation_position"].nunique()
            if unique_only:
                mut_df = self._collapse_to_unique_positions(mut_df)
            ptm_entries.append((ptm_site, ptm_pos, mut_df, mutation_count, unique_count))

        kind_label = "anchor mutation(s)" if cluster else "PTM site(s)"
        if not ptm_entries:
            self._viz_stack_fig.clf()
            self._viz_stack_canvas.draw()
            self._viz_status.configure(
                text=f"No nearby mutations found for any {kind_label} on {gene} ({uid}).",
                text_color=_YELLOW,
            )
            return

        self._viz_status.configure(
            text=f"Rendering {len(ptm_entries)} {kind_label} for {gene} ({uid})…", text_color="gray60",
        )
        domain_entries = _load_interpro_entries(uid)
        protein_length = get_protein_length(_CIF_DIR / uid)
        self._draw_whole_protein_view(uid, gene, ptm_entries, local_window, domain_entries, protein_length)

    def _current_viz_stack_width_in(self) -> float:
        """Current viewport width in inches for the whole-protein stack.
        update_idletasks() first avoids reading a stale pre-reveal size
        after a mode switch.
        """
        self.update_idletasks()
        viewport_w = self._viz_stack_scroll._parent_canvas.winfo_width()
        if viewport_w <= 1:
            return 10.0
        return max(viewport_w / self._viz_stack_fig.get_dpi(), 6.0)

    _WHOLE_PROTEIN_LOLLIPOP_HEIGHT_IN = 2.2  # target lollipop-panel height per row

    def _draw_whole_protein_view(self, uid: str, gene: str, ptm_entries: list, local_window: float,
                                  domain_entries: list[dict], protein_length: int | None) -> None:
        """Draw one lollipop row per PTM site into self._viz_stack_fig.

        Built in batches via self.after() to avoid freezing the UI; a token
        guards against a stale build after the user switches proteins.
        Domain layout is computed once for the whole protein, not per row.
        """
        fig = self._viz_stack_fig
        fig.clf()
        fig.patch.set_facecolor("#2b2b2b")

        n = len(ptm_entries)
        width_in = self._current_viz_stack_width_in()

        max_ptm_pos = max((pos for _, pos, _, _, _ in ptm_entries), default=0)
        max_end = max((e["end"] for e in domain_entries), default=0)
        length = max(protein_length or 0, max_end, max_ptm_pos, 1)
        domain_layout = self._layout_domain_labels(domain_entries, length, width_in)
        row_height_in = self._domain_strip_height_in(domain_layout[1]) + self._WHOLE_PROTEIN_LOLLIPOP_HEIGHT_IN

        total_height_in = 0.6 + n * row_height_in
        dpi = fig.get_dpi()
        fig.set_size_inches(width_in, total_height_in)

        margin_in = 0.3
        outer_gs = fig.add_gridspec(
            n, 1, hspace=0.9,
            left=0.08, right=0.97,
            top=1 - margin_in / total_height_in,
            bottom=margin_in / total_height_in,
        )

        token = object()
        self._viz_stack_build_token = token
        batch_size = 3

        def _build_batch(start: int) -> None:
            if self._viz_stack_build_token is not token:
                return  # superseded by a newer render (protein/mode changed mid-build)
            end = min(start + batch_size, n)
            for i in range(start, end):
                ptm_site, ptm_pos, mut_df, mutation_count, unique_count = ptm_entries[i]
                ax_local, _domain_ax = self._draw_lollipop_group(
                    fig, outer_gs[i, 0], gene, ptm_site, ptm_pos, mut_df, local_window,
                    domain_entries, length, domain_layout, row_height_in,
                    mutation_count, unique_count,
                )
                ax_local.set_title(str(ptm_site), fontsize=10, color="#dcdcdc")
            if end < n:
                self.after(1, _build_batch, end)
            else:
                self._finish_whole_protein_stack(gene, uid, n, dpi, width_in, total_height_in)

        _build_batch(0)

    def _finish_whole_protein_stack(self, gene: str, uid: str, n: int, dpi: float,
                                     width_in: float, total_height_in: float) -> None:
        kind_label = "anchor mutation(s)" if self._viz_is_cluster() else "PTM site(s)"
        self._viz_stack_title_label.configure(text=f"{gene} — {uid} — {n} {kind_label}")

        # matplotlib doesn't resize the Tk widget when the Figure's size
        # changes -- set it explicitly so the scroll region grows too.
        self._viz_stack_canvas.get_tk_widget().configure(
            width=int(width_in * dpi), height=int(total_height_in * dpi),
        )
        self._viz_stack_canvas.draw()

        self._viz_status.configure(
            text=f"{n} {kind_label} for {gene} ({uid})", text_color="gray60",
        )

    def _save_viz_plot(self) -> None:
        whole = self._viz_view_mode_var.get() == "Whole protein"
        fig = self._viz_stack_fig if whole else self._viz_fig
        label = self._viz_protein_combo.get() if whole else self._viz_combo.get()

        if not fig.axes:
            self._viz_status.configure(text="Select something to plot before saving.", text_color=_RED)
            return
        safe = re.sub(r"[^\w-]+", "_", label).strip("_") or "lollipop"
        out_dir = self._output_dir
        out_dir.mkdir(parents=True, exist_ok=True)
        out_path = out_dir / f"lollipop_{safe}.png"
        fig.savefig(out_path, dpi=200, facecolor=fig.get_facecolor(), bbox_inches="tight")
        self._viz_status.configure(text=f"Saved to {out_path}", text_color=_GREEN)
