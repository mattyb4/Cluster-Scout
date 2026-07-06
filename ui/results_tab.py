"""Results tab: PTM Sites / Mutation Details treeviews, search/filter, loading.

`_visualize_selected_ptm` reaches directly into VisualizationTabMixin state
(`self._viz_search_var`, `self._viz_combo`, `self._viz_ptm_rows`, `self._viz_status`)
and calls `self._generate_lollipop_plot()` — the heaviest cross-mixin coupling
in the app. `_load_results` also calls `self._refresh_viz_selector(...)`
(VisualizationTabMixin). Both work correctly via the shared `self` all mixins
compose into.
"""
from __future__ import annotations

import customtkinter as ctk

from ui.common import _MUT_ENTRY_RE, _PP_LABEL, _PTM_TV_COLS, _MUT_TV_COLS, _RED, _GREEN, _YELLOW, _BLUE


class ResultsTabMixin:
    def _setup_treeview_style(self) -> None:
        import tkinter.ttk as ttk
        style = ttk.Style()
        style.theme_use("default")
        bg, fg = "#2b2b2b", "#dcdcdc"
        heading_bg = "#3a3a3a"
        for name in ("Results.Treeview",):
            style.configure(name,
                background=bg, foreground=fg, rowheight=24,
                fieldbackground=bg, borderwidth=0, relief="flat",
                font=("Segoe UI", 10),
            )
            style.configure(f"{name}.Heading",
                background=heading_bg, foreground=fg,
                relief="flat", borderwidth=1,
                font=("Segoe UI", 10, "bold"),
            )
            style.map(name,
                background=[("selected", _BLUE)],
                foreground=[("selected", "white")],
            )
            style.map(f"{name}.Heading",
                background=[("active", "#4a4a4a")],
            )

    def _make_treeview(self, parent, col_defs: list):
        import tkinter as tk
        import tkinter.ttk as ttk

        frame = tk.Frame(parent, bg="#2b2b2b")
        frame.pack(fill=tk.BOTH, expand=True, padx=8, pady=(0, 8))
        frame.grid_rowconfigure(0, weight=1)
        frame.grid_columnconfigure(0, weight=1)

        tv = ttk.Treeview(frame, style="Results.Treeview", show="headings",
                           selectmode="browse")
        col_ids = [c[1] for c in col_defs]
        tv["columns"] = col_ids

        numeric_cols = {"#col", "near", "far", "total", "pts", "cosmic", "seqd", "dist", "pae", "mpld", "maxlin"}
        for display, col_id, width in col_defs:
            anchor = "e" if col_id in numeric_cols else "w"
            tv.heading(col_id, text=display,
                       command=lambda c=col_id: self._sort_tv(tv, c, False))
            stretch = col_id == col_ids[-1]
            tv.column(col_id, width=width, minwidth=30, stretch=stretch, anchor=anchor)

        vsb = ttk.Scrollbar(frame, orient="vertical", command=tv.yview)
        hsb = ttk.Scrollbar(frame, orient="horizontal", command=tv.xview)
        tv.configure(yscrollcommand=vsb.set, xscrollcommand=hsb.set)

        tv.tag_configure("odd",  background="#2b2b2b")
        tv.tag_configure("even", background="#313131")

        tv.grid(row=0, column=0, sticky="nsew")
        vsb.grid(row=0, column=1, sticky="ns")
        hsb.grid(row=1, column=0, sticky="ew")

        return tv

    def _sort_tv(self, tv, col: str, reverse: bool) -> None:
        def _key(v: str):
            try:
                return (0, float(v))
            except (ValueError, TypeError):
                return (1, str(v).lower())

        rows = [(tv.set(k, col), k) for k in tv.get_children("")]
        rows.sort(key=lambda x: _key(x[0]), reverse=reverse)
        for idx, (_, k) in enumerate(rows):
            tv.move(k, "", idx)
        if "#col" in tv["columns"]:
            for idx, k in enumerate(tv.get_children(""), 1):
                tv.set(k, "#col", idx)
        tv.heading(col, command=lambda: self._sort_tv(tv, col, not reverse))

    def _clear_treeview_fully(self, tv, all_rows: list) -> None:
        """Delete every row a treeview has ever held, including ones currently
        hidden by a search filter (plain tv.delete(*tv.get_children()) would
        miss detached rows, leaking them and causing iid collisions on reinsert).
        """
        all_iids = set(tv.get_children("")) | {iid for iid, _, _ in all_rows}
        if all_iids:
            tv.delete(*all_iids)

    def _capture_tv_rows(self, tv) -> list:
        """Snapshot (iid, values, tags) for every row, to filter against later.

        Must be called immediately after a full (unfiltered) populate, while
        every row is still attached.
        """
        return [(iid, tv.item(iid, "values"), tv.item(iid, "tags")) for iid in tv.get_children("")]

    def _filter_treeview(self, tv, all_rows: list, query: str) -> None:
        """Show only rows whose displayed values contain *query* (case-insensitive).

        Uses detach()/move() rather than delete(), so hidden rows keep their
        iid and can reappear — this preserves the iid-as-dataframe-position
        scheme that selection handlers rely on.
        """
        query = query.strip().lower()
        attached = set(tv.get_children(""))
        shown = 0
        for iid, values, _tags in all_rows:
            haystack = " ".join(str(v) for v in values).lower()
            if not query or query in haystack:
                tv.move(iid, "", shown)
                shown += 1
            elif iid in attached:
                tv.detach(iid)

    def _build_results_tab(self, tab) -> None:
        import tkinter as tk

        self._results_df_wide = None
        self._results_df_long = None
        self._ptm_tv_all_rows: list = []
        self._mut_tv_all_rows: list = []
        self._setup_treeview_style()

        tab.grid_columnconfigure(0, weight=1)
        tab.grid_rowconfigure(0, weight=1)

        outer = ctk.CTkFrame(tab, fg_color="transparent")
        outer.grid(row=0, column=0, sticky="nsew", padx=4, pady=4)
        outer.grid_columnconfigure(0, weight=1)
        outer.grid_rowconfigure(1, weight=3)
        outer.grid_rowconfigure(3, weight=2)

        # ── Top panel header ──
        top_header = ctk.CTkFrame(outer, fg_color="transparent")
        top_header.grid(row=0, column=0, sticky="ew", padx=8, pady=(8, 2))
        ctk.CTkLabel(top_header, text="PTM Sites",
                     font=ctk.CTkFont(size=14, weight="bold")).pack(side=tk.LEFT)
        self._results_status = ctk.CTkLabel(top_header, text="",
                                             text_color="gray60",
                                             font=ctk.CTkFont(size=11))
        self._results_status.pack(side=tk.LEFT, padx=(12, 0))
        ctk.CTkButton(top_header, text="↺  Refresh", width=90, height=28,
                       font=ctk.CTkFont(size=12),
                       command=self._load_results).pack(side=tk.RIGHT)
        ctk.CTkButton(top_header, text="📈  Visualize", width=100, height=28,
                       font=ctk.CTkFont(size=12),
                       command=self._visualize_selected_ptm).pack(side=tk.RIGHT, padx=(0, 6))
        self._results_ptm_search_var = ctk.StringVar(value="")
        self._ptm_search_entry = ctk.CTkEntry(
            top_header, textvariable=self._results_ptm_search_var,
            width=220, placeholder_text="🔍 Search PTM sites… (Ctrl+F)",
        )
        self._ptm_search_entry.pack(side=tk.RIGHT, padx=(0, 12))
        self._results_ptm_search_var.trace_add("write", self._filter_ptm_tv)

        # ── PTM site treeview ──
        top_frame = ctk.CTkFrame(outer)
        top_frame.grid(row=1, column=0, sticky="nsew", padx=8, pady=(0, 4))
        top_frame.grid_rowconfigure(0, weight=1)
        top_frame.grid_columnconfigure(0, weight=1)
        self._ptm_tv = self._make_treeview(top_frame, _PTM_TV_COLS)
        self._ptm_tv.bind("<<TreeviewSelect>>", self._on_ptm_select)
        self._ptm_tv.bind("<Double-Button-1>", lambda _e: self._visualize_selected_ptm())

        # ── Detail panel header ──
        bot_header = ctk.CTkFrame(outer, fg_color="transparent")
        bot_header.grid(row=2, column=0, sticky="ew", padx=8, pady=(4, 2))
        ctk.CTkLabel(bot_header, text="Mutation Details",
                     font=ctk.CTkFont(size=14, weight="bold")).pack(side=tk.LEFT, padx=(8, 0))
        self._results_mut_search_var = ctk.StringVar(value="")
        self._mut_search_entry = ctk.CTkEntry(
            bot_header, textvariable=self._results_mut_search_var,
            width=220, placeholder_text="🔍 Search mutations…",
        )
        self._mut_search_entry.pack(side=tk.RIGHT, padx=(0, 6))
        self._results_mut_search_var.trace_add("write", self._filter_mut_tv)

        # ── Mutation detail treeview ──
        bot_frame = ctk.CTkFrame(outer)
        bot_frame.grid(row=3, column=0, sticky="nsew", padx=8, pady=(0, 8))
        bot_frame.grid_rowconfigure(0, weight=1)
        bot_frame.grid_columnconfigure(0, weight=1)
        self._mut_tv = self._make_treeview(bot_frame, _MUT_TV_COLS)

        # Ctrl+F focuses the PTM search box whenever the Results tab is active
        self.bind_all("<Control-f>", self._focus_results_search)

    def _filter_ptm_tv(self, *_args) -> None:
        self._filter_treeview(self._ptm_tv, self._ptm_tv_all_rows, self._results_ptm_search_var.get())

    def _filter_mut_tv(self, *_args) -> None:
        self._filter_treeview(self._mut_tv, self._mut_tv_all_rows, self._results_mut_search_var.get())

    def _focus_results_search(self, event=None):
        """Ctrl+F: focus the PTM search box, but only while the Results tab is showing."""
        if self._tabview.get() == "Results":
            self._ptm_search_entry.focus_set()
            return "break"

    def _load_results(self) -> None:
        import pandas as pd

        wide_path = self._output_dir / "ptm_mutation_proximity_db.tsv"
        long_path = self._output_dir / "ptm_mutation_proximity_long.tsv"

        if not wide_path.exists():
            self._results_status.configure(
                text=f"No output found in {self._output_dir.name}/",
                text_color=_RED,
            )
            self._clear_treeview_fully(self._ptm_tv, self._ptm_tv_all_rows)
            self._clear_treeview_fully(self._mut_tv, self._mut_tv_all_rows)
            self._ptm_tv_all_rows = []
            self._mut_tv_all_rows = []
            self._results_df_wide = None
            self._results_df_long = None
            self._refresh_viz_selector(pd.DataFrame(columns=["gene", "ptm_site", "UniProt"]))
            return

        try:
            df_wide = pd.read_csv(wide_path, sep="\t", encoding="utf-16",
                                   dtype=str, keep_default_na=False)
        except Exception as exc:
            self._results_status.configure(text=f"Error loading file: {exc}",
                                            text_color=_RED)
            return

        df_long = None
        if long_path.exists():
            try:
                df_long = pd.read_csv(long_path, sep="\t", encoding="utf-16",
                                       dtype=str, keep_default_na=False)
            except Exception:
                pass

        self._results_df_wide = df_wide
        self._results_df_long = df_long

        n_sites = len(df_wide)
        n_proteins = df_wide["UniProt"].nunique() if "UniProt" in df_wide.columns else "?"
        long_note = (" · long format available" if df_long is not None
                     else " · enable long format for per-mutation detail")
        self._results_status.configure(
            text=f"{n_sites} PTM sites · {n_proteins} proteins{long_note}",
            text_color="gray60",
        )
        self._populate_ptm_tv(df_wide)
        self._clear_treeview_fully(self._mut_tv, self._mut_tv_all_rows)
        self._mut_tv_all_rows = []
        self._refresh_viz_selector(df_wide)

    def _populate_ptm_tv(self, df) -> None:
        tv = self._ptm_tv
        self._clear_treeview_fully(tv, self._ptm_tv_all_rows)
        for i, (_, row) in enumerate(df.iterrows(), 1):
            try:
                near_pts = int(float(row.get("nearby_muts_total_patient_count", "") or "0"))
                far_pts  = int(float(row.get("distant_muts_total_patient_count", "") or "0"))
                total_pts: int | str = near_pts + far_pts
            except ValueError:
                total_pts = ""
            try:
                near_unique = int(float(row.get("unique_mutation_position_count_within_5_positions", "") or "0"))
                far_unique  = int(float(row.get("unique_mutation_position_count_more_than_5_positions", "") or "0"))
                total_muts: int | str = near_unique + far_unique
            except ValueError:
                total_muts = ""
            linear_dists = [d for d in row.get("morethan5_linear_distance", "").split(",") if d.strip()]
            try:
                max_linear_dist: int | str = max(int(d) for d in linear_dists)
            except ValueError:
                max_linear_dist = ""
            tv.insert("", "end", iid=str(i), values=(
                i,
                row.get("UniProt", ""),
                row.get("gene", ""),
                row.get("ptm_site", ""),
                row.get("ptm_type", ""),
                row.get("mutation_count_within_5_positions", ""),
                row.get("mutation_count_more_than_5_positions", ""),
                total_muts,
                total_pts,
                row.get("total_cosmic_missense_patients", ""),
                row.get("mutation_at_ptm_site", ""),
                row.get("1433pred_binding_site", ""),
                row.get("ptm_is_disordered", ""),
                row.get("ptm_is_binding", ""),
                max_linear_dist,
            ), tags=("odd" if i % 2 else "even",))
        self._ptm_tv_all_rows = self._capture_tv_rows(tv)
        self._filter_ptm_tv()

    def _on_ptm_select(self, *_) -> None:
        sel = self._ptm_tv.selection()
        if not sel or self._results_df_wide is None:
            return
        row = self._results_df_wide.iloc[int(sel[0]) - 1]
        if self._results_df_long is not None:
            uid  = row.get("UniProt", "")
            site = row.get("ptm_site", "")
            mask = (
                (self._results_df_long.get("uniprot_id", "") == uid) &
                (self._results_df_long.get("ptm_position", "") == site)
            )
            self._populate_mut_tv_long(self._results_df_long[mask])
        else:
            self._populate_mut_tv_wide(row)

    def _visualize_selected_ptm(self) -> None:
        """Jump to the Visualization tab and render the lollipop plot for the selected PTM row."""
        sel = self._ptm_tv.selection()
        if not sel or self._results_df_wide is None:
            self._results_status.configure(
                text="Select a PTM site in the table first.", text_color=_YELLOW,
            )
            return
        row = self._results_df_wide.iloc[int(sel[0]) - 1]
        label = f"{row.get('gene', '?')}  {row.get('ptm_site', '?')}  ({row.get('UniProt', '?')})"

        self._tabview.set("Visualization")
        self._viz_search_var.set("")
        if label in self._viz_ptm_rows:
            self._viz_combo.set(label)
            self._generate_lollipop_plot()
        else:
            self._viz_status.configure(
                text=f"Could not find '{label}' in the Visualization selector.", text_color=_RED,
            )

    def _populate_mut_tv_long(self, df) -> None:
        tv = self._mut_tv
        self._clear_treeview_fully(tv, self._mut_tv_all_rows)
        for i, (_, r) in enumerate(df.iterrows(), 1):
            tv.insert("", "end", iid=str(i), values=(
                i,
                r.get("mutation", ""),
                r.get("sequence_distance", ""),
                r.get("distance_angstrom", ""),
                r.get("mut_is_binding", ""),
                r.get("mut_is_disordered", ""),
                r.get("polyphen_class", ""),
                r.get("mutation_plddt", ""),
                r.get("pair_pae", ""),
                r.get("patient_count", ""),
            ), tags=("odd" if i % 2 else "even",))
        self._mut_tv_all_rows = self._capture_tv_rows(tv)
        self._filter_mut_tv()

    def _populate_mut_tv_wide(self, row) -> None:
        import re as _re
        tv = self._mut_tv
        self._clear_treeview_fully(tv, self._mut_tv_all_rows)
        ptm_m = _re.search(r"(\d+)", str(row.get("ptm_site", "")))
        ptm_pos = int(ptm_m.group(1)) if ptm_m else None
        i = 0
        for col_key, _ in [
            ("mutations_within_5_positions",  "≤5 pos"),
            ("mutations_more_than_5_positions", ">5 pos"),
        ]:
            for entry in (row.get(col_key, "") or "").split(", "):
                entry = entry.strip()
                if not entry:
                    continue
                m = _MUT_ENTRY_RE.match(entry)
                if not m:
                    continue
                i += 1
                pp_code = m.group(2) or ""
                mut_m   = _re.search(r"\d+", m.group(1))
                mut_pos = int(mut_m.group()) if mut_m else None
                seq_d   = abs(mut_pos - ptm_pos) if (mut_pos is not None and ptm_pos is not None) else ""
                tv.insert("", "end", iid=str(i), values=(
                    i,
                    m.group(1),
                    seq_d,
                    m.group(4),
                    "",  # is_binding (unavailable in wide format)
                    "",  # is_disordered (unavailable in wide format)
                    _PP_LABEL.get(pp_code, ""),
                    "",  # mut pLDDT (unavailable in wide format)
                    m.group(5) or "",
                    "",  # patient count (unavailable in wide format)
                ), tags=("odd" if i % 2 else "even",))
        self._mut_tv_all_rows = self._capture_tv_rows(tv)
        self._filter_mut_tv()
