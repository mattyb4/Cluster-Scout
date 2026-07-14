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
    _PTM_COL_HELP, _MUT_COL_HELP,
    _RED, _GREEN, _YELLOW, _BLUE,
    _load_column_prefs, _save_column_prefs, help_icon,
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
                relief="groove", borderwidth=1,
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
        self._disable_tv_stretch(tv)

        vsb = ttk.Scrollbar(frame, orient="vertical", command=tv.yview)
        hsb = ttk.Scrollbar(frame, orient="horizontal", command=tv.xview)
        tv.configure(yscrollcommand=vsb.set, xscrollcommand=hsb.set)

        tv.tag_configure("odd",  background="#2b2b2b")
        tv.tag_configure("even", background="#313131")

        tv.grid(row=0, column=0, sticky="nsew")
        vsb.grid(row=0, column=1, sticky="ns")
        hsb.grid(row=1, column=0, sticky="ew")

        return tv

    def _bind_treeview_zoom_override(self, tv) -> None:
        """Intercept Ctrl+wheel on this treeview before ttk's own unmodified
        `<MouseWheel>` class-binding scrolls it (widget-level bindings run
        before class-level ones, so this reliably wins). See
        App._on_ctrl_scroll_zoom for why this override is needed at all.
        """
        tv.bind("<Control-MouseWheel>", self._on_ctrl_scroll_zoom)
        tv.bind("<Control-Button-4>", self._on_ctrl_scroll_zoom)
        tv.bind("<Control-Button-5>", self._on_ctrl_scroll_zoom)

    def _disable_tv_stretch(self, tv) -> None:
        """Keep every column at exactly its configured width, always.

        ttk.Treeview's `stretch=True` doesn't just grow a column into extra
        room when there's space to spare -- if the visible columns'
        configured widths already add up to more than the widget's actual
        width, ttk instead shrinks the stretchy column below its configured
        width (down to `minwidth`) to force everything to fit. Rather than
        try to distinguish those two cases, no column ever stretches: any
        leftover space after the last column is just left blank, and
        overflow is handled by the horizontal scrollbar. Every column then
        reliably renders at the width set for it in the registry.
        """
        for c in tv["columns"]:
            tv.column(c, stretch=False)

    def _set_visible_columns(self, which: str, col_ids: list) -> None:
        tv = self._ptm_tv if which == "ptm" else self._mut_tv
        tv.configure(displaycolumns=col_ids)
        self._disable_tv_stretch(tv)
        if which == "ptm":
            self._ptm_visible_cols = col_ids
        else:
            self._mut_visible_cols = col_ids
        _save_column_prefs(which, col_ids)

    def _open_column_picker(self, which: str) -> None:
        registry = _PTM_TV_COLS if which == "ptm" else _MUT_TV_COLS
        col_help = _PTM_COL_HELP if which == "ptm" else _MUT_COL_HELP
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
            row = ctk.CTkFrame(scroll, fg_color="transparent")
            row.pack(anchor="w", fill="x", pady=2, padx=4)
            ctk.CTkCheckBox(row, text=display, variable=var).pack(side="left")
            help_text = col_help.get(col_id)
            if help_text:
                help_icon(row, help_text).pack(side="left", padx=(6, 0))

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

    def _tv_overlay(self, tv):
        """Lazily create (and cache) a centered message label floating over *tv*,
        used to show 'Loading data…' / error / empty-state text in place of rows.
        """
        label = self._tv_placeholder_labels.get(tv)
        if label is None:
            label = ctk.CTkLabel(
                tv.master, font=ctk.CTkFont(size=13), fg_color="#242424",
                corner_radius=6, padx=16, pady=10, wraplength=380, justify="center",
            )
            self._tv_placeholder_labels[tv] = label
        return label

    def _show_tv_message(self, tv, text: str, color: str = "gray60") -> None:
        overlay = self._tv_overlay(tv)
        overlay.configure(text=text, text_color=color)
        overlay.place(relx=0.5, rely=0.5, anchor="center")
        overlay.lift()

    def _hide_tv_message(self, tv) -> None:
        label = self._tv_placeholder_labels.get(tv)
        if label is not None:
            label.place_forget()

    def _capture_tv_rows(self, tv) -> list:
        """Snapshot (iid, values, tags) for every row, to filter against later.

        Must be called immediately after a full (unfiltered) populate, while
        every row is still attached.
        """
        return [(iid, tv.item(iid, "values"), tv.item(iid, "tags")) for iid in tv.get_children("")]

    _NUMERIC_FILTER_OPS = (">", ">=", "<", "<=", "=", "≠")
    _TEXT_FILTER_OPS = ("contains", "does not contain", "equals")

    def _filter_treeview(self, tv, all_rows: list, query: str, filters: list | None = None) -> None:
        """Show only rows matching the search *query* (substring, case-insensitive)
        AND every rule in *filters* (col_id/op/value dicts — all must match, AND).

        Uses detach()/move() rather than delete(), so hidden rows keep their
        iid and can reappear — this preserves the iid-as-dataframe-position
        scheme that selection handlers rely on.
        """
        query = query.strip().lower()
        col_ids = list(tv["columns"])
        attached = set(tv.get_children(""))
        shown = 0
        for iid, values, _tags in all_rows:
            haystack = " ".join(str(v) for v in values).lower()
            text_ok = not query or query in haystack
            rules_ok = not filters or all(
                self._eval_filter_rule(values, col_ids, rule) for rule in filters
            )
            if text_ok and rules_ok:
                tv.move(iid, "", shown)
                shown += 1
            elif iid in attached:
                tv.detach(iid)

    def _eval_filter_rule(self, values, col_ids: list, rule: dict) -> bool:
        try:
            idx = col_ids.index(rule["col_id"])
        except ValueError:
            return True  # stale/unknown column (e.g. a saved rule from a removed column) — don't hide everything
        cell = values[idx] if idx < len(values) else ""
        op = rule["op"]
        target = rule["value"]

        if op in self._NUMERIC_FILTER_OPS:
            try:
                cell_num = float(cell)
                target_num = float(target)
            except (ValueError, TypeError):
                return False  # can't compare non-numeric cell numerically -> excluded
            if op == ">":
                return cell_num > target_num
            if op == ">=":
                return cell_num >= target_num
            if op == "<":
                return cell_num < target_num
            if op == "<=":
                return cell_num <= target_num
            if op == "=":
                return cell_num == target_num
            return cell_num != target_num  # "≠"

        cell_l = str(cell).lower()
        target_l = str(target).lower()
        if op == "contains":
            return target_l in cell_l
        if op == "does not contain":
            return target_l not in cell_l
        return cell_l == target_l  # "equals"

    def _ops_for_column(self, registry: list, col_id: str) -> tuple:
        for _label, cid, _width, numeric, _default in registry:
            if cid == col_id:
                return self._NUMERIC_FILTER_OPS if numeric else self._TEXT_FILTER_OPS
        return self._TEXT_FILTER_OPS

    def _update_filter_button_label(self, which: str) -> None:
        btn = self._ptm_filter_button if which == "ptm" else self._mut_filter_button
        n = len(self._ptm_filters if which == "ptm" else self._mut_filters)
        btn.configure(text=f"▽  Filter ({n})" if n else "▽  Filter")

    def _open_filter_picker(self, which: str) -> None:
        registry = [c for c in (_PTM_TV_COLS if which == "ptm" else _MUT_TV_COLS) if c[1] != "#col"]
        current_filters = self._ptm_filters if which == "ptm" else self._mut_filters
        title = "Filter PTM Sites" if which == "ptm" else "Filter Mutation Details"
        labels_by_id = {c[1]: c[0] for c in registry}
        id_by_label = {c[0]: c[1] for c in registry}
        col_labels = [c[0] for c in registry]

        win = ctk.CTkToplevel(self)
        win.title(title)
        win.geometry("560x440")
        win.transient(self)
        win.grab_set()

        scroll = ctk.CTkScrollableFrame(win, label_text="Filter rules (all must match)")
        scroll.pack(fill="both", expand=True, padx=10, pady=(10, 4))

        rows: list[dict] = []

        def _add_row(col_id: str | None = None, op: str | None = None, value: str = ""):
            col_id = col_id if col_id in labels_by_id else registry[0][1]
            ops = self._ops_for_column(registry, col_id)
            op = op if op in ops else ops[0]

            row_frame = ctk.CTkFrame(scroll, fg_color="transparent")
            row_frame.pack(fill="x", pady=3)

            col_var = ctk.StringVar(value=labels_by_id[col_id])
            op_var = ctk.StringVar(value=op)
            val_var = ctk.StringVar(value=value)

            def _on_col_change(selected_label):
                new_ops = self._ops_for_column(registry, id_by_label[selected_label])
                op_menu.configure(values=list(new_ops))
                if op_var.get() not in new_ops:
                    op_var.set(new_ops[0])

            col_menu = ctk.CTkOptionMenu(row_frame, values=col_labels, variable=col_var,
                                          width=170, command=_on_col_change)
            col_menu.pack(side="left", padx=(0, 4))

            op_menu = ctk.CTkOptionMenu(row_frame, values=list(ops), variable=op_var, width=130)
            op_menu.pack(side="left", padx=4)

            val_entry = ctk.CTkEntry(row_frame, textvariable=val_var, width=100)
            val_entry.pack(side="left", padx=4, fill="x", expand=True)

            def _remove():
                row_frame.destroy()
                rows[:] = [r for r in rows if r["frame"] is not row_frame]

            ctk.CTkButton(row_frame, text="✕", width=28, fg_color="transparent",
                          hover_color="#3a3a3a", command=_remove).pack(side="left", padx=(4, 0))

            rows.append({"frame": row_frame, "col_var": col_var, "op_var": op_var, "val_var": val_var})

        if current_filters:
            for f in current_filters:
                _add_row(f["col_id"], f["op"], f["value"])
        else:
            _add_row()

        ctk.CTkButton(win, text="+ Add filter", width=110,
                      command=lambda: _add_row()).pack(anchor="w", padx=14, pady=(0, 6))

        btn_row = ctk.CTkFrame(win, fg_color="transparent")
        btn_row.pack(fill="x", padx=10, pady=(4, 10))

        def _clear_all():
            current_filters.clear()
            self._update_filter_button_label(which)
            (self._filter_ptm_tv if which == "ptm" else self._filter_mut_tv)()
            win.destroy()

        def _apply():
            new_filters = []
            for r in rows:
                value = r["val_var"].get().strip()
                if not value:
                    continue
                new_filters.append({
                    "col_id": id_by_label[r["col_var"].get()],
                    "op": r["op_var"].get(),
                    "value": value,
                })
            if which == "ptm":
                self._ptm_filters = new_filters
            else:
                self._mut_filters = new_filters
            self._update_filter_button_label(which)
            (self._filter_ptm_tv if which == "ptm" else self._filter_mut_tv)()
            win.destroy()

        ctk.CTkButton(btn_row, text="Clear all", width=90, command=_clear_all).pack(side="left")
        ctk.CTkButton(btn_row, text="Cancel", width=80, command=win.destroy).pack(side="right")
        ctk.CTkButton(btn_row, text="Apply", width=80, command=_apply).pack(side="right", padx=(0, 6))

    def _build_results_tab(self, tab) -> None:
        import tkinter as tk

        self._results_df_wide = None
        self._results_df_long = None
        self._results_loaded_key = None
        self._ptm_tv_all_rows: list = []
        self._mut_tv_all_rows: list = []
        self._tv_placeholder_labels: dict = {}
        self._ptm_filters: list = []
        self._mut_filters: list = []
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
        self._refresh_button = ctk.CTkButton(
            top_header, text="↺  Refresh", width=90, height=28,
            font=ctk.CTkFont(size=12),
            command=lambda: self._load_results(force=True),
        )
        self._refresh_button.pack(side=tk.RIGHT)
        ctk.CTkButton(top_header, text="⚙  Columns", width=95, height=28,
                       font=ctk.CTkFont(size=12),
                       command=lambda: self._open_column_picker("ptm")).pack(side=tk.RIGHT, padx=(0, 6))
        self._ptm_filter_button = ctk.CTkButton(
            top_header, text="▽  Filter", width=90, height=28,
            font=ctk.CTkFont(size=12),
            command=lambda: self._open_filter_picker("ptm"),
        )
        self._ptm_filter_button.pack(side=tk.RIGHT, padx=(0, 6))
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
        self._bind_treeview_zoom_override(self._ptm_tv)

        # ── Detail panel header ──
        bot_header = ctk.CTkFrame(outer, fg_color="transparent")
        bot_header.grid(row=2, column=0, sticky="ew", padx=8, pady=(4, 2))
        ctk.CTkLabel(bot_header, text="Mutation Details",
                     font=ctk.CTkFont(size=14, weight="bold")).pack(side=tk.LEFT, padx=(8, 0))
        ctk.CTkButton(bot_header, text="⚙  Columns", width=95, height=28,
                       font=ctk.CTkFont(size=12),
                       command=lambda: self._open_column_picker("mut")).pack(side=tk.RIGHT, padx=(0, 6))
        self._mut_filter_button = ctk.CTkButton(
            bot_header, text="▽  Filter", width=90, height=28,
            font=ctk.CTkFont(size=12),
            command=lambda: self._open_filter_picker("mut"),
        )
        self._mut_filter_button.pack(side=tk.RIGHT, padx=(0, 6))
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
        self._bind_treeview_zoom_override(self._mut_tv)

        # Ctrl+F focuses the PTM search box whenever the Results tab is active
        self.bind_all("<Control-f>", self._focus_results_search)

    def _filter_ptm_tv(self, *_args) -> None:
        self._filter_treeview(self._ptm_tv, self._ptm_tv_all_rows,
                               self._results_ptm_search_var.get(), self._ptm_filters)

    def _filter_mut_tv(self, *_args) -> None:
        self._filter_treeview(self._mut_tv, self._mut_tv_all_rows,
                               self._results_mut_search_var.get(), self._mut_filters)

    def _focus_results_search(self, event=None):
        """Ctrl+F: focus the PTM search box, but only while the Results tab is showing."""
        if self._tabview.get() == "Results":
            self._ptm_search_entry.focus_set()
            return "break"

    def _load_results(self, force: bool = False) -> None:
        """Show a 'Loading data…' placeholder immediately, then do the actual
        file read/populate — unless *force* is False and the output file
        hasn't changed since the last load, in which case this is a no-op.

        Every switch to the Results tab calls this, so without the mtime
        check, revisiting a tab you'd already loaded would still re-read and
        re-populate a multi-thousand-row TSV from disk each time — a ~1.5s
        block that's laggy no matter how well a loading placeholder paints.
        The "↺ Refresh" button passes force=True to always bypass this.

        `after(..., callback)` only *schedules* the callback — it doesn't
        guarantee the pending tab-switch/placeholder redraw is flushed to
        screen before the callback (often blocking-fast) timer fires, so that
        approach still let the tab switch look frozen. `update_idletasks()`
        alone isn't enough either — it flushes Tk's internal geometry/idle
        queue but doesn't pump the window system's own redraw/expose events,
        so on some triggers (e.g. clicking Refresh with no tab-raise involved)
        the canvas changes were queued but never actually painted to screen
        before the blocking read started. A full `update()` processes those
        pending redraw events too, guaranteeing the placeholder is actually
        visible before we block on the pandas read here.
        """
        wide_path = self._output_dir / "ptm_mutation_proximity_db.tsv"
        if not force and self._results_df_wide is not None:
            try:
                key = (wide_path, wide_path.stat().st_mtime)
            except OSError:
                key = (wide_path, None)
            if key == self._results_loaded_key:
                return

        for tv, attr in ((self._ptm_tv, "_ptm_tv_all_rows"), (self._mut_tv, "_mut_tv_all_rows")):
            self._clear_treeview_fully(tv, getattr(self, attr))
            setattr(self, attr, [])
            self._show_tv_message(tv, "Loading data…")
        self._results_status.configure(text="Loading data…", text_color="gray60")
        self._refresh_button.configure(state="disabled")
        self.update()
        try:
            self._load_results_now()
        finally:
            self._refresh_button.configure(state="normal")

    def _load_results_now(self) -> None:
        import pandas as pd

        wide_path = self._output_dir / "ptm_mutation_proximity_db.tsv"
        long_path = self._output_dir / "ptm_mutation_proximity_long.tsv"

        if not wide_path.exists():
            msg = f"No output found in {self._output_dir.name}/"
            self._results_status.configure(text=msg, text_color=_RED)
            self._show_tv_message(self._ptm_tv, msg, _RED)
            self._show_tv_message(self._mut_tv, msg, _RED)
            self._results_df_wide = None
            self._results_df_long = None
            self._results_loaded_key = None
            self._refresh_viz_selector(pd.DataFrame(columns=["gene", "ptm_site", "UniProt"]))
            return

        try:
            df_wide = pd.read_csv(wide_path, sep="\t", encoding="utf-16",
                                   dtype=str, keep_default_na=False)

            df_long = None
            if long_path.exists():
                try:
                    df_long = pd.read_csv(long_path, sep="\t", encoding="utf-16",
                                           dtype=str, keep_default_na=False)
                except Exception:
                    pass

            self._results_df_wide = df_wide
            self._results_df_long = df_long
            self._results_loaded_key = (wide_path, wide_path.stat().st_mtime)

            n_sites = len(df_wide)
            n_proteins = df_wide["UniProt"].nunique() if "UniProt" in df_wide.columns else "?"
            long_note = (" · long format available" if df_long is not None
                         else " · enable long format for per-mutation detail")
            self._results_status.configure(
                text=f"{n_sites} PTM sites · {n_proteins} proteins{long_note}",
                text_color="gray60",
            )
            self._hide_tv_message(self._ptm_tv)
            self._hide_tv_message(self._mut_tv)
            self._populate_ptm_tv(df_wide)
            self._clear_treeview_fully(self._mut_tv, self._mut_tv_all_rows)
            self._mut_tv_all_rows = []
            self._refresh_viz_selector(df_wide)
        except Exception as exc:
            self._results_df_wide = None
            self._results_df_long = None
            self._results_loaded_key = None
            msg = f"Error loading results: {exc}"
            self._results_status.configure(text=msg, text_color=_RED)
            self._show_tv_message(self._ptm_tv, msg, _RED)
            self._show_tv_message(self._mut_tv, msg, _RED)

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
        available at this granularity in the wide format, so they're left blank.
        "PTM pLDDT" is the one PTM-level exception kept in this table (it has no
        equivalent in the PTM Sites table), so it's also left blank here since
        the wide format doesn't carry it per-row either.
        """
        import re as _re
        tv = self._mut_tv
        self._clear_treeview_fully(tv, self._mut_tv_all_rows)
        ptm_m = _re.search(r"(\d+)", str(row.get("ptm_site", "")))
        ptm_pos = int(ptm_m.group(1)) if ptm_m else None

        # Mutation names confirmed to disrupt this PTM site, so each row below
        # can report "yes"/"no" for itself rather than repeating the whole list.
        confirmed_muts = set()
        for entry in (row.get("confirmed_disrupting_mutations", "") or "").split(", "):
            cm = _MUT_ENTRY_RE.match(entry.strip())
            if cm:
                confirmed_muts.add(cm.group(1))

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
                    "confirmed_disrupt": "yes" if m.group(1) in confirmed_muts else "no",
                    "ptm_plddt": "",
                    "mut_aiupred_gen": "",
                    "mut_aiupred_bind": "",
                }
                values = [i] + [per_row.get(c, "") for c in _MUT_TV_SRC_IDS]
                tv.insert("", "end", iid=str(i), values=values, tags=("odd" if i % 2 else "even",))
        self._mut_tv_all_rows = self._capture_tv_rows(tv)
        self._filter_mut_tv()
