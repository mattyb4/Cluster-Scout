"""Visualization tab: lollipop (needle) plot of mutations near a PTM site.

Reads `self._results_df_wide`/`self._results_df_long` directly (ResultsTabMixin
state) — these two tabs are tightly bidirectionally coupled by design, since
selecting a PTM site in Results drives what's plotted here. `_style_lollipop_axis`
is also called by AnalysisToolsTabMixin's `_style_dark_figure`.
"""
from __future__ import annotations

import re

import customtkinter as ctk

from ui.common import _MUT_ENTRY_RE, _PP_LABEL, _PP_COLORS, _PTM_MARKER_COLOR, _NEEDLE_DEFAULT_COLOR, _RED, _YELLOW, _GREEN


class VisualizationTabMixin:
    def _build_viz_tab(self, tab) -> None:
        import matplotlib
        matplotlib.use("TkAgg")
        from matplotlib.figure import Figure
        from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

        tab.grid_columnconfigure(0, weight=1)
        tab.grid_rowconfigure(2, weight=1)

        self._viz_ptm_rows: dict[str, int] = {}
        self._viz_all_labels: list[str] = []

        # ── Controls row 1: PTM selection ──
        controls = ctk.CTkFrame(tab)
        controls.grid(row=0, column=0, padx=12, pady=(12, 6), sticky="ew")

        ctk.CTkLabel(
            controls, text="PTM site:", font=ctk.CTkFont(weight="bold"),
        ).pack(side="left", padx=(12, 6), pady=10)

        self._viz_search_var = ctk.StringVar(value="")
        search_entry = ctk.CTkEntry(
            controls, textvariable=self._viz_search_var, width=160,
            placeholder_text="Search gene or site…",
        )
        search_entry.pack(side="left", padx=(0, 6), pady=10)
        self._viz_search_var.trace_add("write", self._on_viz_search)

        self._viz_combo = ctk.CTkComboBox(controls, width=220, values=[])
        self._viz_combo.pack(side="left", padx=(0, 12), pady=10)

        ctk.CTkLabel(controls, text="±aa window:").pack(side="left", padx=(8, 4), pady=10)
        self._viz_window_var = ctk.StringVar(value="15")
        ctk.CTkEntry(
            controls, textvariable=self._viz_window_var, width=50,
        ).pack(side="left", padx=(0, 12), pady=10)

        # ── Controls row 2: display mode + actions ──
        controls2 = ctk.CTkFrame(tab)
        controls2.grid(row=1, column=0, padx=12, pady=(0, 6), sticky="ew")

        ctk.CTkLabel(controls2, text="Show:", font=ctk.CTkFont(weight="bold")).pack(
            side="left", padx=(12, 6), pady=10,
        )
        self._viz_mode_var = ctk.StringVar(value="All mutations")
        ctk.CTkSegmentedButton(
            controls2, values=["All mutations", "Unique per position"],
            variable=self._viz_mode_var,
            command=lambda _v: self._generate_lollipop_plot(),
        ).pack(side="left", padx=(0, 12), pady=10)

        ctk.CTkButton(
            controls2, text="Generate", width=100,
            command=self._generate_lollipop_plot,
        ).pack(side="left", padx=(0, 8), pady=10)

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

        # ── Plot canvas ──
        canvas_frame = ctk.CTkFrame(tab, fg_color="#2b2b2b")
        canvas_frame.grid(row=2, column=0, padx=12, pady=(0, 12), sticky="nsew")
        canvas_frame.grid_columnconfigure(0, weight=1)
        canvas_frame.grid_rowconfigure(0, weight=1)

        self._viz_fig = Figure(figsize=(10, 6), dpi=100, facecolor="#2b2b2b")
        self._viz_canvas = FigureCanvasTkAgg(self._viz_fig, master=canvas_frame)
        self._viz_canvas.get_tk_widget().grid(row=0, column=0, sticky="nsew")

    # ── Visualization: PTM selection ────────────────────────────────────────

    def _refresh_viz_selector(self, df) -> None:
        """(Re)populate the PTM-site combobox from the currently loaded results."""
        self._viz_ptm_rows = {}
        labels: list[str] = []
        for idx, row in df.iterrows():
            gene = row.get("gene", "?")
            site = row.get("ptm_site", "?")
            uid = row.get("UniProt", "?")
            label = f"{gene}  {site}  ({uid})"
            labels.append(label)
            self._viz_ptm_rows[label] = idx
        self._viz_all_labels = labels
        current = self._viz_combo.get()
        self._viz_combo.configure(values=labels)
        if labels and current not in labels:
            self._viz_combo.set(labels[0])
        elif not labels:
            self._viz_combo.set("")

    def _on_viz_search(self, *_args) -> None:
        query = self._viz_search_var.get().strip().lower()
        filtered = (
            [label for label in self._viz_all_labels if query in label.lower()]
            if query else self._viz_all_labels
        )
        self._viz_combo.configure(values=filtered)
        if filtered and self._viz_combo.get() not in filtered:
            self._viz_combo.set(filtered[0])

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

        for i, (_, r) in enumerate(df.iterrows()):
            x = x_plot[i]
            count = r["patient_count"]
            color = _PP_COLORS.get(str(r.get("polyphen_class", "")).strip(), _NEEDLE_DEFAULT_COLOR)
            ax.vlines(x, 0, count, color=color, linewidth=2, zorder=2)
            ax.scatter([x], [count], s=90, color=color, edgecolor="black", zorder=3)
            # The far (categorical) panel already names each mutation via its x-tick
            # label, so only annotate above the marker for the true-scale local panel.
            if jitter:
                ax.annotate(
                    str(r["mutation"]), (x, count), xytext=(4, 10),
                    textcoords="offset points", ha="left", va="bottom",
                    fontsize=8, color="#dcdcdc", rotation=45,
                )

    def _draw_lollipop(self, gene: str, ptm_site: str, ptm_pos: int, mut_df, local_window: float) -> None:
        from matplotlib.lines import Line2D
        from matplotlib.patches import Patch

        fig = self._viz_fig
        fig.clf()
        fig.patch.set_facecolor("#2b2b2b")

        local_df = mut_df[(mut_df["mutation_position"] - ptm_pos).abs() <= local_window] \
            .sort_values("mutation_position")
        far_df = mut_df[(mut_df["mutation_position"] - ptm_pos).abs() > local_window] \
            .sort_values("mutation_position").reset_index(drop=True)
        has_far = not far_df.empty

        max_count = max(float(mut_df["patient_count"].max()), 1.0)
        ptm_height = max_count * 1.15
        y_top = ptm_height * 1.45

        if has_far:
            gs = fig.add_gridspec(
                1, 3, width_ratios=[0.68, 0.05, 0.27], wspace=0.05,
                left=0.08, right=0.97, top=0.86, bottom=0.28,
            )
            ax_local = fig.add_subplot(gs[0])
            ax_gap = fig.add_subplot(gs[1])
            ax_far = fig.add_subplot(gs[2])
        else:
            gs = fig.add_gridspec(1, 1, left=0.08, right=0.97, top=0.86, bottom=0.28)
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

        legend_handles = [
            Line2D([0], [0], marker="*", color="none", markerfacecolor=_PTM_MARKER_COLOR,
                   markeredgecolor="white", markersize=14, label="PTM site"),
            Patch(facecolor=_PP_COLORS["benign"], label="Benign"),
            Patch(facecolor=_PP_COLORS["possibly_damaging"], label="Possibly damaging"),
            Patch(facecolor=_PP_COLORS["probably_damaging"], label="Probably damaging"),
            Patch(facecolor=_NEEDLE_DEFAULT_COLOR, label="Unknown / not scored"),
        ]
        ax_local.legend(
            handles=legend_handles, loc="upper left", fontsize=8,
            facecolor="#3a3a3a", edgecolor="#555555", labelcolor="#dcdcdc",
        )

        fig.suptitle(f"{gene} — {ptm_site} and nearby mutations", color="white",
                     fontsize=14, fontweight="bold")

        self._viz_canvas.draw()

    def _generate_lollipop_plot(self) -> None:
        if self._results_df_wide is None:
            self._viz_status.configure(
                text="No results loaded — run the pipeline or open the Results tab first.",
                text_color=_RED,
            )
            return

        label = self._viz_combo.get()
        idx = self._viz_ptm_rows.get(label)
        if idx is None:
            self._viz_status.configure(text="Select a PTM site to plot.", text_color=_RED)
            return

        row = self._results_df_wide.loc[idx]
        gene = row.get("gene", "?")
        ptm_site = row.get("ptm_site", "?")
        uid = row.get("UniProt", "?")

        ptm_m = re.search(r"\d+", str(ptm_site))
        if not ptm_m:
            self._viz_status.configure(
                text=f"Could not parse a residue position from '{ptm_site}'.", text_color=_RED,
            )
            return
        ptm_pos = int(ptm_m.group())

        mut_df, note = self._get_viz_mutation_df(uid, ptm_site, row)
        if mut_df.empty:
            self._viz_fig.clf()
            self._viz_canvas.draw()
            self._viz_status.configure(
                text=f"No nearby mutations found for {gene} {ptm_site}.", text_color=_YELLOW,
            )
            return

        unique_only = self._viz_mode_var.get() == "Unique per position"
        if unique_only:
            mut_df = self._collapse_to_unique_positions(mut_df)

        try:
            local_window = float(self._viz_window_var.get().strip() or "15")
        except ValueError:
            local_window = 15.0

        self._draw_lollipop(gene, ptm_site, ptm_pos, mut_df, local_window)

        kind = "unique position(s)" if unique_only else "nearby mutation(s)"
        status = f"{len(mut_df)} {kind} for {gene} {ptm_site}"
        if note:
            status += f" — {note}"
        self._viz_status.configure(text=status, text_color="gray60")

    def _collapse_to_unique_positions(self, df):
        """Merge mutations that share a residue position into one entry per position.

        Patient counts are summed across substitutions at that residue; the
        displayed color reflects the most severe PolyPhen class present.
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

    def _save_viz_plot(self) -> None:
        if not self._viz_fig.axes:
            self._viz_status.configure(text="Generate a plot before saving.", text_color=_RED)
            return
        label = self._viz_combo.get()
        safe = re.sub(r"[^\w-]+", "_", label).strip("_") or "lollipop"
        out_dir = self._output_dir
        out_dir.mkdir(parents=True, exist_ok=True)
        out_path = out_dir / f"lollipop_{safe}.png"
        self._viz_fig.savefig(out_path, dpi=200, facecolor=self._viz_fig.get_facecolor(), bbox_inches="tight")
        self._viz_status.configure(text=f"Saved to {out_path}", text_color=_GREEN)
