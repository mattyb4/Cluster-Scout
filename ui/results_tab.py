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

from ui.common import (
    _MUT_ENTRY_RE, _PP_LABEL, _PTM_TV_COLS, _MUT_TV_COLS, _MUT_LONG_SRC_MAP,
    _RED, _GREEN, _YELLOW, _BLUE,
    _load_column_prefs, _save_column_prefs,
)

_PTM_TV_SRC_IDS = [c[1] for c in _PTM_TV_COLS if c[1] != "#col"]
_MUT_TV_SRC_IDS = [c[1] for c in _MUT_TV_COLS if c[1] != "#col"]


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

    def _make_treeview(self, parent, col_defs: list, visible_ids: list | None = None):
        """Build a treeview with every column in *col_defs* defined, but only
        *visible_ids* (default: the default_visible=True ones) actually shown.

        All columns always exist so `values=` tuples passed to `insert()` stay
        a fixed shape; `displaycolumns` is what the Columns picker toggles at
        runtime, without needing to rebuild the widget.
        """
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

        for display, col_id, width, numeric, default in col_defs:
            anchor = "e" if numeric else "w"
            tv.heading(col_id, text=display,
                       command=lambda c=col_id: self._sort_tv(tv, c, False))
            tv.column(col_id, width=width, minwidth=30, stretch=False, anchor=anchor)

        if visible_ids is None:
            visible_ids = [c[1] for c in col_defs if c[4]]
        tv["displaycolumns"] = visible_ids
        self._apply_tv_stretch(tv)

        vsb = ttk.Scrollbar(frame, orient="vertical", command=tv.yview)
        hsb = ttk.Scrollbar(frame, orient="horizontal", command=tv.xview)
        tv.configure(yscrollcommand=vsb.set, xscrollcommand=hsb.set)

        tv.tag_configure("odd",  background="#2b2b2b")
        tv.tag_configure("even", background="#313131")

        tv.grid(row=0, column=0, sticky="nsew")
        vsb.grid(row=0, column=1, sticky="ns")
        hsb.grid(row=1, column=0, sticky="ew")

        return tv

    def _apply_tv_stretch(self, tv) -> None:
        """Make the last currently-visible column absorb extra width."""
        for c in tv["columns"]:
            tv.column(c, stretch=False)
        display = list(tv["displaycolumns"])
        if display:
            tv.column(display[-1], stretch=True)

    def _set_visible_columns(self, which: str, col_ids: list) -> None:
        tv = self._ptm_tv if which == "ptm" else self._mut_tv
        tv.configure(displaycolumns=col_ids)
        self._apply_tv_stretch(tv)
        if which == "ptm":
            self._ptm_visible_cols = col_ids
        else:
            self._mut_visible_cols = col_ids
        _save_column_prefs(which, col_ids)

    def _open_column_picker(self, which: str) -> None:
        registry = _PTM_TV_COLS if which == "ptm" else _MUT_TV_COLS
        current = set(self._ptm_visible_cols if which == "ptm" else self._mut_visible_cols)
        title = "PTM Sites Columns" if which == "ptm" else "Mutation Details Columns"

        win = ctk.CTkToplevel(self)
        win.title(title)
        win.geometry("320x480")
        win.transient(self)
        win.grab_set()

        scroll = ctk.CTkScrollableFrame(win, label_text="Show columns")
        scroll.pack(fill="both", expand=True, padx=10, pady=(10, 4))

        vars_by_id: dict = {}
        for display, col_id, _width, _numeric, _default in registry:
            if col_id == "#col":
                continue
            var = ctk.BooleanVar(value=col_id in current)
            vars_by_id[col_id] = var
            ctk.CTkCheckBox(scroll, text=display, variable=var).pack(anchor="w", pady=2, padx=4)

        btn_row = ctk.CTkFrame(win, fg_color="transparent")
        btn_row.pack(fill="x", padx=10, pady=(4, 10))

        def _reset():
            for _display, col_id, _width, _numeric, default in registry:
                if col_id in vars_by_id:
                    vars_by_id[col_id].set(default)

        def _apply():
            selected = {col_id for col_id, var in vars_by_id.items() if var.get()}
            ordered = ["#col"] + [c[1] for c in registry if c[1] in selected]
            self._set_visible_columns(which, ordered)
            win.destroy()

        ctk.CTkButton(btn_row, text="Reset to defaults", width=130,
                      command=_reset).pack(side="left")
        ctk.CTkButton(btn_row, text="Cancel", width=80,
                      command=win.destroy).pack(side="right")
        ctk.CTkButton(btn_row, text="Apply", width=80,
                      command=_apply).pack(side="right", padx=(0, 6))

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

        _ptm_known = {c[1] for c in _PTM_TV_COLS}
        _mut_known = {c[1] for c in _MUT_TV_COLS}
        _ptm_saved = _load_column_prefs("ptm")
        _mut_saved = _load_column_prefs("mut")
        self._ptm_visible_cols = (
            [c for c in _ptm_saved if c in _ptm_known] if _ptm_saved
            else [c[1] for c in _PTM_TV_COLS if c[4]]
        )
        self._mut_visible_cols = (
            [c for c in _mut_saved if c in _mut_known] if _mut_saved
            else [c[1] for c in _MUT_TV_COLS if c[4]]
        )

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
        ctk.CTkButton(top_header, text="⚙  Columns", width=95, height=28,
                       font=ctk.CTkFont(size=12),
                       command=lambda: self._open_column_picker("ptm")).pack(side=tk.RIGHT, padx=(0, 6))
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
        self._ptm_tv = self._make_treeview(top_frame, _PTM_TV_COLS, self._ptm_visible_cols)
        self._ptm_tv.bind("<<TreeviewSelect>>", self._on_ptm_select)
        self._ptm_tv.bind("<Double-Button-1>", lambda _e: self._visualize_selected_ptm())

        # ── Detail panel header ──
        bot_header = ctk.CTkFrame(outer, fg_color="transparent")
        bot_header.grid(row=2, column=0, sticky="ew", padx=8, pady=(4, 2))
        ctk.CTkLabel(bot_header, text="Mutation Details",
                     font=ctk.CTkFont(size=14, weight="bold")).pack(side=tk.LEFT, padx=(8, 0))
        ctk.CTkButton(bot_header, text="⚙  Columns", width=95, height=28,
                       font=ctk.CTkFont(size=12),
                       command=lambda: self._open_column_picker("mut")).pack(side=tk.RIGHT, padx=(0, 6))
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
        self._mut_tv = self._make_treeview(bot_frame, _MUT_TV_COLS, self._mut_visible_cols)

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
                near_pts = far_pts = ""
                total_pts = ""
            try:
                near_unique = int(float(row.get("unique_mutation_position_count_within_5_positions", "") or "0"))
                far_unique  = int(float(row.get("unique_mutation_position_count_more_than_5_positions", "") or "0"))
                total_muts: int | str = near_unique + far_unique
            except ValueError:
                near_unique = far_unique = ""
                total_muts = ""
            linear_dists = [d for d in row.get("morethan5_linear_distance", "").split(",") if d.strip()]
            try:
                max_linear_dist: int | str = max(int(d) for d in linear_dists)
            except ValueError:
                max_linear_dist = ""
            values_map = {
                "uniprot": row.get("UniProt", ""),
                "gene": row.get("gene", ""),
                "site": row.get("ptm_site", ""),
                "type": row.get("ptm_type", ""),
                "near": row.get("mutation_count_within_5_positions", ""),
                "far": row.get("mutation_count_more_than_5_positions", ""),
                "near_pts": near_pts,
                "far_pts": far_pts,
                "near_unique": near_unique,
                "far_unique": far_unique,
                "total": total_muts,
                "pts": total_pts,
                "cosmic": row.get("total_cosmic_missense_patients", ""),
                "atptm": row.get("mutation_at_ptm_site", ""),
                "confirmed_disrupt": row.get("confirmed_disrupting_mutations", ""),
                "diseases": row.get("ptm_diseases", ""),
                "pred14": row.get("1433pred_binding_site", ""),
                "pred14_consensus": row.get("1433pred_consensus", ""),
                "conf14": row.get("1433_confirmed_site", ""),
                "conf14_pmid": row.get("1433_confirmed_pmid", ""),
                "kinases": row.get("kinase_predictions", ""),
                "aiupred_gen": row.get("ptm_aiupred_general", ""),
                "aiupred_bind": row.get("ptm_aiupred_binding", ""),
                "disord": row.get("ptm_is_disordered", ""),
                "bind": row.get("ptm_is_binding", ""),
                "maxlin": max_linear_dist,
                "near_muts_raw": row.get("mutations_within_5_positions", ""),
                "far_muts_raw": row.get("mutations_more_than_5_positions", ""),
                "lin_dist_raw": row.get("morethan5_linear_distance", ""),
            }
            values = [i] + [values_map.get(c, "") for c in _PTM_TV_SRC_IDS]
            tv.insert("", "end", iid=str(i), values=values, tags=("odd" if i % 2 else "even",))
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
            values = [i] + [r.get(_MUT_LONG_SRC_MAP[c], "") for c in _MUT_TV_SRC_IDS]
            tv.insert("", "end", iid=str(i), values=values, tags=("odd" if i % 2 else "even",))
        self._mut_tv_all_rows = self._capture_tv_rows(tv)
        self._filter_mut_tv()

    def _populate_mut_tv_wide(self, row) -> None:
        """Populate the Mutation Details tv from a wide-format PTM row.

        Per-mutation fields (binding/disordered/pLDDT/patient count/etc.) aren't
        available at this granularity in the wide format, so they're left blank;
        PTM-level fields (diseases, kinases, aiupred, ...) are shared across all
        mutation rows derived from the same PTM site.
        """
        import re as _re
        tv = self._mut_tv
        self._clear_treeview_fully(tv, self._mut_tv_all_rows)
        ptm_m = _re.search(r"(\d+)", str(row.get("ptm_site", "")))
        ptm_pos = int(ptm_m.group(1)) if ptm_m else None

        shared = {
            "confirmed_disrupt": row.get("confirmed_disrupting_mutations", ""),
            "diseases": row.get("ptm_diseases", ""),
            "pred14": row.get("1433pred_binding_site", ""),
            "pred14_consensus": row.get("1433pred_consensus", ""),
            "conf14": row.get("1433_confirmed_site", ""),
            "kinases": row.get("kinase_predictions", ""),
            "ptm_aiupred_gen": row.get("ptm_aiupred_general", ""),
            "ptm_aiupred_bind": row.get("ptm_aiupred_binding", ""),
            "ptm_disord": row.get("ptm_is_disordered", ""),
            "ptm_bind": row.get("ptm_is_binding", ""),
            "cosmic": row.get("total_cosmic_missense_patients", ""),
            "gene": row.get("gene", ""),
            "uniprot": row.get("UniProt", ""),
            "ptm_position": row.get("ptm_site", ""),
            "ptm_type_l": row.get("ptm_type", ""),
        }

        i = 0
        for col_key in ("mutations_within_5_positions", "mutations_more_than_5_positions"):
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
                per_row = {
                    "mut": m.group(1),
                    "seqd": seq_d,
                    "dist": m.group(4),
                    "isbnd": "",
                    "isdis": "",
                    "ppc": _PP_LABEL.get(pp_code, ""),
                    "pps": m.group(3) or "",
                    "mpld": "",
                    "pae": m.group(5) or "",
                    "pts": "",
                    "total_near_pts": "",
                    "near_mut_count": "",
                    "ptm_plddt": "",
                    "mut_aiupred_gen": "",
                    "mut_aiupred_bind": "",
                    **shared,
                }
                values = [i] + [per_row.get(c, "") for c in _MUT_TV_SRC_IDS]
                tv.insert("", "end", iid=str(i), values=values, tags=("odd" if i % 2 else "even",))
        self._mut_tv_all_rows = self._capture_tv_rows(tv)
        self._filter_mut_tv()
