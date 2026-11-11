"""Structure Heatmaps tab: parameters, validation, the export run, and the list of
exported proteins shown in the results area.

Works from any UniProt-resolvable protein (or one uploaded CIF), independent of the
pipeline's outputs, so it has its own tab rather than a place among the Pipeline
modes. The parameter fields are rebuilt when a color reset needs the swatches
redrawn, so every setting lives in a persistent `self._heatmap_*` variable created on
first use. Runs through the app's single background-run slot (`self._running`) and
this tab's own status/progress/log widgets; `_poll_queue` (pipeline_runner.py) routes
queue messages here instead of to the Pipeline tab whenever `self._hm_run_active` is set.
"""
from __future__ import annotations

import queue
import sys
import threading
import time
from pathlib import Path
from tkinter import filedialog

import customtkinter as ctk

from ui import theme
from ui.common import (
    _BLUE,
    _BTN2,
    _BTN2_HOVER,
    _CHIP,
    _CHIP_HOVER,
    _GRAY,
    _GREEN,
    _PLOT_BG,
    _RED,
    _YELLOW,
    OUTPUT_DIR,
    _fmt_time,
    add_resize_grip,
    color_swatch_button,
    extract_uniprot_from_cif,
    help_icon,
    isolate_textbox_scroll,
    section_header,
)

_SIDE_W = 400  # px; fixed width of the parameters panel

_HEATMAP_DEFAULT_OUTPUT_DIR = OUTPUT_DIR / "coordinates"


class StructureHeatmapsMixin:
    # ── Tab layout ───────────────────────────────────────────────────────────

    def _build_heatmaps_tab(self, tab) -> None:
        """Parameters on the left (Run pinned at the bottom), a divider, and the
        results list with a collapsible log on the right."""
        tab.grid_columnconfigure(2, weight=1)
        tab.grid_rowconfigure(0, weight=1)

        side = ctk.CTkFrame(tab, fg_color="transparent", corner_radius=0, border_width=0, width=_SIDE_W)
        side.grid(row=0, column=0, sticky="ns", padx=(14, 10), pady=12)
        side.grid_propagate(False)
        side.grid_columnconfigure(0, weight=1)
        side.grid_rowconfigure(1, weight=1)
        section_header(side, "Parameters", row=0, pady=(0, 2))

        # Explicit fg_color (not "transparent"): children of a transparent CTkScrollableFrame
        # capture a static bg color at construction and show stale corners after a light/dark switch.
        params_scroll = ctk.CTkScrollableFrame(side, fg_color=theme.pair("bg"))
        params_scroll.grid(row=1, column=0, sticky="nsew")
        params_scroll.grid_columnconfigure(0, weight=1)
        self._hm_params_frame = ctk.CTkFrame(params_scroll, fg_color="transparent")
        self._hm_params_frame.grid(row=0, column=0, sticky="ew")
        self._hm_params_frame.grid_columnconfigure(0, weight=1)

        run_area = ctk.CTkFrame(side, fg_color="transparent")
        run_area.grid(row=2, column=0, sticky="ew", pady=(8, 0))
        run_area.grid_columnconfigure(0, weight=1)
        self._hm_progress_bar = ctk.CTkProgressBar(run_area, height=10)
        self._hm_progress_bar.set(0)
        self._hm_progress_bar.grid(row=0, column=0, sticky="ew", pady=(0, 8))
        self._hm_progress_bar.grid_remove()
        self._hm_run_btn = ctk.CTkButton(
            run_area, text="▶  Run", height=36, font=ctk.CTkFont(size=14, weight="bold"),
            command=self._start_heatmap_run,
        )
        self._hm_run_btn.grid(row=1, column=0, sticky="ew")

        ctk.CTkFrame(
            tab, width=1, corner_radius=0, border_width=0, fg_color=theme.pair("line"),
        ).grid(row=0, column=1, sticky="ns", pady=12)

        main = ctk.CTkFrame(tab, fg_color="transparent", corner_radius=0, border_width=0)
        main.grid(row=0, column=2, sticky="nsew", padx=14, pady=12)
        main.grid_columnconfigure(0, weight=1)
        main.grid_rowconfigure(1, weight=1)

        # Results heading: title, live status, hairline
        head = ctk.CTkFrame(main, fg_color="transparent")
        head.grid(row=0, column=0, sticky="ew", pady=(0, 6))
        head.grid_columnconfigure(2, weight=1)
        ctk.CTkLabel(
            head, text="RESULTS", text_color=_GRAY, font=ctk.CTkFont(size=11, weight="bold"),
        ).grid(row=0, column=0, padx=(0, 10))
        self._hm_status_label = ctk.CTkLabel(
            head, text="Set parameters on the left, then click Run.",
            text_color=_GRAY, font=ctk.CTkFont(size=11), anchor="w",
        )
        self._hm_status_label.grid(row=0, column=1, padx=(0, 10))
        ctk.CTkFrame(
            head, width=0, height=1, corner_radius=0, border_width=0, fg_color=theme.pair("line"),
        ).grid(row=0, column=2, sticky="ew")

        results_area = ctk.CTkFrame(main, fg_color=_PLOT_BG)
        results_area.grid(row=1, column=0, sticky="nsew")
        results_area.grid_columnconfigure(0, weight=1)
        results_area.grid_rowconfigure(0, weight=1)
        # Shown over the results until a run has produced some
        self._hm_empty_label = ctk.CTkLabel(
            results_area, text="Add proteins and click Run to export coordinates and ChimeraX scripts.",
            text_color=_GRAY, font=ctk.CTkFont(size=13),
        )
        self._build_heatmap_results_panel(results_area)

        # Collapsed log drawer under the results
        self._hm_log_visible = False
        self._hm_log_frame = ctk.CTkFrame(main, fg_color="transparent")
        self._hm_log = ctk.CTkTextbox(
            self._hm_log_frame, height=140, font=ctk.CTkFont(family="Courier New", size=12),
            wrap="word", state="disabled", border_width=1, border_color=theme.pair("line"),
        )
        self._hm_log.pack(fill="both", expand=True)
        isolate_textbox_scroll(self._hm_log)
        add_resize_grip(self._hm_log).pack(fill="x")
        self._hm_log_toggle = ctk.CTkButton(
            main, text="▸  Show details", height=26, anchor="w", width=0,
            fg_color="transparent", hover_color=_CHIP_HOVER,
            font=ctk.CTkFont(size=12), command=self._hm_toggle_log,
        )
        self._hm_log_toggle.grid(row=3, column=0, sticky="w", pady=(6, 0))

        self._rebuild_heatmap_fields()
        self._refresh_heatmap_results()

    def _hm_toggle_log(self) -> None:
        """Show or hide this tab's own collapsible log panel."""
        if self._hm_log_visible:
            self._hm_log_frame.grid_remove()
            self._hm_log_toggle.configure(text="▸  Show details")
            self._hm_log_visible = False
        else:
            self._hm_log_frame.grid(row=2, column=0, sticky="ew", pady=(8, 0))
            self._hm_log_toggle.configure(text="▾  Hide details")
            self._hm_log_visible = True

    def _hm_field_label(self, row: int, text: str, help_text: str | None = None, parent=None) -> None:
        """Small muted label (plus optional help icon) above a field in the params panel
        (or in *parent*, a frame nested inside it)."""
        frame = ctk.CTkFrame(parent or self._hm_params_frame, fg_color="transparent", width=0, height=0)
        frame.grid(row=row, column=0, sticky="w", pady=(12, 3))
        ctk.CTkLabel(
            frame, text=text, anchor="w", text_color=_GRAY, font=ctk.CTkFont(size=12),
        ).pack(side="left")
        if help_text:
            help_icon(frame, help_text).pack(side="left", padx=(5, 0))

    def _rebuild_heatmap_fields(self) -> None:
        """Clear and rebuild the params frame (e.g. to redraw color swatches after a reset)."""
        for w in self._hm_params_frame.winfo_children():
            w.destroy()
        self._hm_params_frame.grid_columnconfigure(0, weight=1)
        self._build_structure_heatmaps_fields()

    # ── Fields ───────────────────────────────────────────────────────────────

    def _build_structure_heatmaps_fields(self) -> None:
        """Structure Heatmaps parameter fields, built into self._hm_params_frame."""
        f = self._hm_params_frame

        from export_ca_coordinates import (
            MUTATION_DEFAULT_HIGH_COLOR,
            MUTATION_DEFAULT_LOW_COLOR,
            MUTATION_DEFAULT_MID_COLOR,
            MUTATION_MARKER_DEFAULT_COLOR,
            PLDDT_DEFAULT_HIGH_COLOR,
            PLDDT_DEFAULT_LOW_COLOR,
            PLDDT_DEFAULT_MID_COLOR,
            PTM_MARKER_DEFAULT_COLOR,
        )
        if not hasattr(self, "_heatmap_source_var"):
            self._heatmap_source_var = ctk.StringVar(value="Database")
            self._heatmap_proteins: list[str] = []
            self._heatmap_items: list[tuple[str, Path | None, str | None]] = []
            self._heatmap_custom_cif_var = ctk.StringVar(value="")
            self._heatmap_custom_cif_uniprot_var = ctk.StringVar(value="")
            self._heatmap_mutation_var = ctk.BooleanVar(value=False)
            self._heatmap_residue_var = ctk.BooleanVar(value=False)
            self._heatmap_log_scale_var = ctk.BooleanVar(value=False)
            self._heatmap_dim_confidence_var = ctk.BooleanVar(value=False)
            self._heatmap_plddt_var = ctk.BooleanVar(value=False)
            self._heatmap_mark_ptm_var = ctk.BooleanVar(value=False)
            self._heatmap_mark_mutations_var = ctk.BooleanVar(value=False)
            self._heatmap_mutation_low_var = ctk.StringVar(value=MUTATION_DEFAULT_LOW_COLOR)
            self._heatmap_mutation_mid_var = ctk.StringVar(value=MUTATION_DEFAULT_MID_COLOR)
            self._heatmap_mutation_high_var = ctk.StringVar(value=MUTATION_DEFAULT_HIGH_COLOR)
            self._heatmap_plddt_low_var = ctk.StringVar(value=PLDDT_DEFAULT_LOW_COLOR)
            self._heatmap_plddt_mid_var = ctk.StringVar(value=PLDDT_DEFAULT_MID_COLOR)
            self._heatmap_plddt_high_var = ctk.StringVar(value=PLDDT_DEFAULT_HIGH_COLOR)
            self._heatmap_ptm_marker_color_var = ctk.StringVar(value=PTM_MARKER_DEFAULT_COLOR)
            self._heatmap_mutation_marker_color_var = ctk.StringVar(value=MUTATION_MARKER_DEFAULT_COLOR)
            self._heatmap_output_dir_var = ctk.StringVar(value=str(_HEATMAP_DEFAULT_OUTPUT_DIR))

        def _reset_colors(*pairs: tuple[ctk.StringVar, str]) -> None:
            for var, default in pairs:
                var.set(default)
            self._rebuild_heatmap_fields()  # redraws the swatches with the restored colors

        def _checkbox_row(row: int, text: str, var: ctk.BooleanVar, help_text: str,
                          indent: int = 0, command=None, pady=(2, 0), font=None) -> ctk.CTkCheckBox:
            frame = ctk.CTkFrame(f, fg_color="transparent")
            frame.grid(row=row, column=0, sticky="w", padx=(indent, 0), pady=pady)
            box = ctk.CTkCheckBox(
                frame, text=text, variable=var, checkbox_width=18, checkbox_height=18, command=command, font=font,
            )
            box.pack(side="left")
            help_icon(frame, help_text).pack(side="left", padx=(4, 0))
            return box

        def _colors_row(row: int, vars_: tuple[ctk.StringVar, ...], defaults: tuple[str, ...],
                        help_text: str, pady=(2, 0)) -> None:
            frame = ctk.CTkFrame(f, fg_color="transparent")
            frame.grid(row=row, column=0, sticky="w", padx=(24, 0), pady=pady)
            ctk.CTkLabel(frame, text="Colors:", text_color=_GRAY, font=ctk.CTkFont(size=12)).pack(
                side="left", padx=(0, 6),
            )
            for i, var in enumerate(vars_):
                if i:
                    ctk.CTkLabel(frame, text="→").pack(side="left", padx=3)
                color_swatch_button(frame, var).pack(side="left")
            ctk.CTkButton(
                frame, text="↺", width=24, height=22, fg_color=_BTN2, hover_color=_BTN2_HOVER,
                command=lambda: _reset_colors(*zip(vars_, defaults)),
            ).pack(side="left", padx=(6, 0))
            help_icon(frame, help_text).pack(side="left", padx=(4, 0))

        def _marker_row(row: int, text: str, var: ctk.BooleanVar, help_text: str,
                        color_var: ctk.StringVar, default: str, pady=(2, 0)) -> None:
            frame = ctk.CTkFrame(f, fg_color="transparent")
            frame.grid(row=row, column=0, sticky="w", pady=pady)
            ctk.CTkCheckBox(
                frame, text=text, variable=var, checkbox_width=18, checkbox_height=18,
            ).pack(side="left")
            help_icon(frame, help_text).pack(side="left", padx=(4, 0))
            color_swatch_button(frame, color_var).pack(side="left", padx=(8, 0))
            ctk.CTkButton(
                frame, text="↺", width=24, height=22, fg_color=_BTN2, hover_color=_BTN2_HOVER,
                command=lambda: _reset_colors((color_var, default)),
            ).pack(side="left", padx=(4, 0))

        # Source toggle: batch of database-backed proteins, or one caller-provided CIF
        self._hm_field_label(
            0, "Source",
            "Database exports one or more proteins by gene/UniProt, using "
            "AlphaFold DB's model (downloaded automatically if needed). "
            "Upload CIF file instead uses a specific .cif you already have "
            "-- e.g. a seeded AlphaFold Server prediction from the CIF "
            "Variance tool's \"Generate AlphaFold Seeds JSON\" option -- for "
            "a single protein, with nothing downloaded or added to the "
            "shared structure cache.",
        )
        ctk.CTkSegmentedButton(
            f, values=["Database", "Upload CIF file"],
            variable=self._heatmap_source_var, command=self._on_heatmap_source_change,
        ).grid(row=1, column=0, sticky="ew")

        # ── Database source: proteins list (gene symbols and/or UniProt accessions) ──
        self._heatmap_database_frame = ctk.CTkFrame(f, fg_color="transparent")
        self._heatmap_database_frame.grid(row=2, column=0, sticky="ew")
        self._heatmap_database_frame.grid_columnconfigure(0, weight=1)
        db = self._heatmap_database_frame

        self._hm_field_label(
            0, "Proteins",
            "Add one or more gene symbols and/or UniProt accessions - each "
            "is exported in turn (its own {gene}_{UniProt}/ folder inside "
            "the output folder), with the same options below applied to "
            "every one. Each entry is validated against the UniProt API and "
            "AlphaFold DB when you click Add (or press Enter), the same as "
            "Radius Sweep's gene list, so a typo or an unmodeled protein is "
            "reported immediately instead of only failing partway through "
            "a run.",
            parent=db,
        )
        input_frame = ctk.CTkFrame(db, fg_color="transparent")
        input_frame.grid(row=1, column=0, sticky="ew")
        input_frame.grid_columnconfigure(0, weight=1)
        self._heatmap_protein_input_var = ctk.StringVar(value="")
        self._heatmap_protein_entry = ctk.CTkEntry(
            input_frame, textvariable=self._heatmap_protein_input_var,
            placeholder_text="e.g. P04637 or TP53",
        )
        self._heatmap_protein_entry.grid(row=0, column=0, sticky="ew", padx=(0, 6))
        self._heatmap_protein_entry.bind("<Return>", lambda _e: self._add_heatmap_protein())
        self._heatmap_protein_add_btn = ctk.CTkButton(
            input_frame, text="+ Add", width=64, height=28, command=self._add_heatmap_protein,
        )
        self._heatmap_protein_add_btn.grid(row=0, column=1)

        # Feedback for _add_heatmap_protein — hidden until there's something to say
        self._heatmap_protein_error_label = ctk.CTkLabel(
            db, text="", text_color=_RED, font=ctk.CTkFont(size=11), anchor="w",
            justify="left", wraplength=_SIDE_W - 50,
        )
        self._heatmap_protein_error_label.grid(row=2, column=0, sticky="w", pady=(4, 0))
        self._heatmap_protein_error_label.grid_remove()

        self._heatmap_proteins_list_frame = ctk.CTkFrame(db, fg_color="transparent")
        self._heatmap_proteins_list_frame.grid(row=3, column=0, sticky="ew", pady=(6, 0))
        self._refresh_heatmap_protein_chips()

        # ── Upload-CIF source: one specific file + its UniProt ID ──
        self._heatmap_custom_cif_frame = ctk.CTkFrame(f, fg_color="transparent")
        self._heatmap_custom_cif_frame.grid(row=2, column=0, sticky="ew")
        self._heatmap_custom_cif_frame.grid_columnconfigure(0, weight=1)
        cf = self._heatmap_custom_cif_frame

        self._hm_field_label(0, "CIF file", parent=cf)
        cif_row = ctk.CTkFrame(cf, fg_color="transparent")
        cif_row.grid(row=1, column=0, sticky="ew")
        cif_row.grid_columnconfigure(0, weight=1)
        ctk.CTkEntry(cif_row, textvariable=self._heatmap_custom_cif_var).grid(
            row=0, column=0, sticky="ew", padx=(0, 6),
        )
        ctk.CTkButton(
            cif_row, text="Browse", width=70, height=26, font=ctk.CTkFont(size=12),
            command=self._browse_heatmap_custom_cif,
        ).grid(row=0, column=1)
        self._hm_field_label(2, "UniProt ID", parent=cf)
        ctk.CTkEntry(
            cf, textvariable=self._heatmap_custom_cif_uniprot_var,
            placeholder_text="auto-detected, or enter if not found",
        ).grid(row=3, column=0, sticky="ew")

        # Both frames default to gridded -- show only the selected source's
        self._on_heatmap_source_change()

        # ── Heatmaps and markers ──
        self._hm_field_label(3, "Heatmaps and markers")

        self._heatmap_mutation_checkbox = _checkbox_row(
            4, "Mutation heatmap (patients within 10 Å)", self._heatmap_mutation_var,
            "Colors the ChimeraX cartoon by COSMIC patient count summed "
            "within 10 Å of each residue - a mutation-hotspot heatmap, "
            "using a sequential red palette. Single-fragment proteins only.",
            indent=0, command=self._on_heatmap_mutation_toggle,
            # CustomTkinter's bundled Roboto clips the ring of the Å on Windows; Segoe UI doesn't
            font=ctk.CTkFont(family="Segoe UI", size=13) if sys.platform == "win32" else None,
        )
        _checkbox_row(
            5, "Mutation heatmap (patients at residue only)", self._heatmap_residue_var,
            "Colors the ChimeraX cartoon by the COSMIC patient count at each "
            "residue itself - the patients whose mutation is at that exact "
            "position, with no neighbors counted - so only mutated residues "
            "are colored and every other residue stays at the low color. "
            "Written to its own script, separate from the 10 Å heatmap above, "
            "so either or both can be on. Uses the same colors, Log-scale and "
            "Dim options as that heatmap. Single-fragment proteins only.",
            command=self._on_heatmap_mutation_toggle,
        )
        _colors_row(
            6, (self._heatmap_mutation_low_var, self._heatmap_mutation_mid_var, self._heatmap_mutation_high_var),
            (MUTATION_DEFAULT_LOW_COLOR, MUTATION_DEFAULT_MID_COLOR, MUTATION_DEFAULT_HIGH_COLOR),
            "Low, middle and high colors of the mutation heatmap's scale (no "
            "nearby patients, half the protein's highest count, and the "
            "highest count). Click a swatch to choose a color; ↺ resets all "
            "three back to the default red scale. Applies to both mutation "
            "heatmaps.",
        )
        self._heatmap_log_scale_checkbox = _checkbox_row(
            7, "Log-scale", self._heatmap_log_scale_var,
            "Only applies to the mutation heatmaps above. COSMIC patient "
            "counts are often heavily right-skewed - a few hotspot residues "
            "can have counts far higher than the rest, so a linear color "
            "scale crushes nearly the whole protein into one flat color. "
            "Log-scaling compresses that range so variation among "
            "lower-count residues becomes visible too.",
            indent=24,
        )
        self._heatmap_dim_confidence_checkbox = _checkbox_row(
            8, "Dim low-confidence residues", self._heatmap_dim_confidence_var,
            "Only applies to the mutation heatmaps above. Fades each "
            "residue's heatmap color in proportion to how low its AlphaFold "
            "confidence (pLDDT) is, so a mutation hotspot in a poorly-"
            "modeled region reads as less trustworthy than an equally hot "
            "one in a well-modeled region, instead of both looking equally "
            "certain. Switches the script's lighting from \"soft\" to "
            "\"simple\", since ChimeraX's soft ambient shadows don't render "
            "correctly once part of the model is transparent.",
            indent=24, pady=(2, 8),
        )
        self._on_heatmap_mutation_toggle()  # sync initial enabled/disabled state

        _checkbox_row(
            9, "pLDDT heatmap", self._heatmap_plddt_var,
            "Colors the ChimeraX cartoon by AlphaFold's per-residue "
            "confidence score (pLDDT), using ChimeraX's built-in AlphaFold "
            "color scheme (blue = very high confidence, down to orange = "
            "very low). Useful for judging how reliable the modeled "
            "structure is in a region, independent of mutation data. "
            "Single-fragment proteins only.",
        )
        _colors_row(
            10, (self._heatmap_plddt_low_var, self._heatmap_plddt_mid_var, self._heatmap_plddt_high_var),
            (PLDDT_DEFAULT_LOW_COLOR, PLDDT_DEFAULT_MID_COLOR, PLDDT_DEFAULT_HIGH_COLOR),
            "Low, middle and high colors of the pLDDT heatmap's scale "
            "(pLDDT 0, 50 and 100). Click a swatch to choose a color; ↺ "
            "resets all three back to AlphaFold's own color scheme.",
            pady=(2, 8),
        )

        # Independent markers, layered on top of the heatmaps rather than recoloring
        _marker_row(
            11, "Mark PTM sites", self._heatmap_mark_ptm_var,
            "Marks each known PTM site (from the pipeline's PTM/mutation "
            "data) with a small sphere at its CA coordinate, layered on "
            "top of whichever heatmap(s) above are selected (or a plain "
            "cartoon if neither is). This is a separate marker, not a "
            "recoloring, so it never overwrites a heatmap's own color at "
            "that residue. Single-fragment proteins only.",
            self._heatmap_ptm_marker_color_var, PTM_MARKER_DEFAULT_COLOR,
        )
        _marker_row(
            12, "Show mutation markers", self._heatmap_mark_mutations_var,
            "Shows each COSMIC mutation position's side chain as a colored "
            "stick, layered on top of whichever heatmap(s) above are "
            "selected (or a plain cartoon if neither is). Like Mark PTM "
            "sites, this reveals real side-chain atoms rather than "
            "recoloring the cartoon, so it never overwrites a heatmap's own "
            "color at that residue. Single-fragment proteins only.",
            self._heatmap_mutation_marker_color_var, MUTATION_MARKER_DEFAULT_COLOR,
        )

        # Output folder
        self._hm_field_label(
            13, "Output folder",
            "Where each protein's {gene}_{UniProt}/ folder (coordinate tables, "
            "ChimeraX scripts) is written. Defaults to Output/coordinates.",
        )
        ctk.CTkEntry(
            f, textvariable=self._heatmap_output_dir_var, state="readonly",
        ).grid(row=14, column=0, sticky="ew")
        buttons = ctk.CTkFrame(f, fg_color="transparent")
        buttons.grid(row=15, column=0, sticky="w", pady=(6, 0))
        ctk.CTkButton(
            buttons, text="Browse", width=70, height=26, font=ctk.CTkFont(size=12),
            command=self._browse_heatmap_output_dir,
        ).pack(side="left", padx=(0, 6))
        ctk.CTkButton(
            buttons, text="Reset", width=60, height=26, font=ctk.CTkFont(size=12),
            fg_color=_BTN2, hover_color=_BTN2_HOVER,
            command=lambda: self._heatmap_output_dir_var.set(str(_HEATMAP_DEFAULT_OUTPUT_DIR)),
        ).pack(side="left", padx=(0, 6))
        ctk.CTkButton(
            buttons, text="Open Folder", width=90, height=26, font=ctk.CTkFont(size=12),
            fg_color=_BTN2, hover_color=_BTN2_HOVER,
            command=lambda: self._open_folder(Path(self._heatmap_output_dir_var.get())),
        ).pack(side="left")

    def _browse_heatmap_output_dir(self) -> None:
        path = filedialog.askdirectory(
            title="Select Structure Heatmaps output folder",
            initialdir=self._heatmap_output_dir_var.get(),
        )
        if path:
            self._heatmap_output_dir_var.set(path)

    def _on_heatmap_source_change(self, _value: str = "") -> None:
        if self._heatmap_source_var.get() == "Upload CIF file":
            self._heatmap_database_frame.grid_remove()
            self._heatmap_custom_cif_frame.grid()
        else:
            self._heatmap_custom_cif_frame.grid_remove()
            self._heatmap_database_frame.grid()

    def _browse_heatmap_custom_cif(self) -> None:
        """Open a file dialog for a user-supplied CIF and try to auto-fill its
        UniProt ID from embedded metadata.

        Unlike Single Protein mode's _browse_cif, this does NOT fall back to
        the file's parent folder name when detection fails -- that heuristic
        only makes sense for cif_models/{UniProt}/... paths, but a custom
        upload (e.g. fresh out of an AlphaFold Server download or a
        data/cif_comparison/ folder) can live anywhere, so a folder-name
        guess here would usually just be wrong rather than merely unhelpful.
        """
        path = filedialog.askopenfilename(
            title="Select CIF structure file",
            filetypes=[("CIF files", "*.cif"), ("All files", "*.*")],
        )
        if not path:
            return
        self._heatmap_custom_cif_var.set(path)
        uid = extract_uniprot_from_cif(Path(path))
        if uid:
            self._heatmap_custom_cif_uniprot_var.set(uid)

    def _on_heatmap_mutation_toggle(self) -> None:
        """Keep the log-scale and dim-confidence checkboxes in sync with the
        mutation-heatmap toggles they modify -- disabled (and visually
        grayed) when there's no mutation heatmap left for them to affect.
        """
        state = "normal" if self._heatmap_mutation_var.get() or self._heatmap_residue_var.get() else "disabled"
        self._heatmap_log_scale_checkbox.configure(state=state)
        self._heatmap_dim_confidence_checkbox.configure(state=state)

    # ── Protein list ─────────────────────────────────────────────────────────

    def _add_heatmap_protein(self) -> None:
        """Add a protein token (gene symbol or UniProt accession) to the
        batch list, after validating it live against the UniProt API and
        AlphaFold DB -- the same "catch bad input at Add time, not mid-run"
        philosophy as Radius Sweep's own gene picker (_add_radius_gene),
        adapted to this tool's own data dependency: Radius Sweep
        validates against the pipeline's local PTM_TSV (cheap, synchronous),
        but Structure Heatmaps works from any UniProt-resolvable protein via
        live lookups, so the check runs off the main thread (same
        queue.Queue + self.after() polling pattern as
        AnalysisToolsTabMixin._generate_alphafold_seed_json) to avoid
        freezing the UI on every click, with the input/button disabled and
        a "Checking…" status shown meanwhile.
        """
        token = self._heatmap_protein_input_var.get().strip()
        if not token:
            return

        if token.upper() in self._heatmap_proteins:
            self._heatmap_protein_error_label.configure(
                text=f"⚠  {token.upper()} is already in the list.", text_color=_YELLOW,
            )
            self._heatmap_protein_error_label.grid()
            self._heatmap_protein_input_var.set("")
            return

        self._heatmap_protein_error_label.configure(text=f"⏳  Checking {token}…", text_color=_GRAY)
        self._heatmap_protein_error_label.grid()
        self._heatmap_protein_add_btn.configure(state="disabled")
        self._heatmap_protein_entry.configure(state="disabled")

        result_q: queue.Queue = queue.Queue()

        def _worker() -> None:
            from export_ca_coordinates import resolve_and_check_protein
            result_q.put(resolve_and_check_protein(token, log_cb=lambda *_: None))

        threading.Thread(target=_worker, daemon=True).start()
        self.after(100, lambda: self._poll_heatmap_protein_check(result_q))

    def _poll_heatmap_protein_check(self, result_q: "queue.Queue") -> None:
        try:
            result = result_q.get_nowait()
        except queue.Empty:
            self.after(100, lambda: self._poll_heatmap_protein_check(result_q))
            return

        # The fields are torn down and rebuilt on every sub-tool switch -- bail
        # out quietly if the user switched away before this check finished.
        if not self._heatmap_protein_entry.winfo_exists():
            return

        self._heatmap_protein_add_btn.configure(state="normal")
        self._heatmap_protein_entry.configure(state="normal")

        if not result["ok"]:
            self._heatmap_protein_error_label.configure(text=f"⚠  {result['error']}", text_color=_RED)
            self._heatmap_protein_error_label.grid()
            return

        gene = result["gene"]
        if gene in self._heatmap_proteins:
            self._heatmap_protein_error_label.configure(
                text=f"⚠  {gene} is already in the list.", text_color=_YELLOW,
            )
            self._heatmap_protein_error_label.grid()
            self._heatmap_protein_input_var.set("")
            return

        self._heatmap_proteins.append(gene)
        self._refresh_heatmap_protein_chips()
        self._heatmap_protein_input_var.set("")

        # Non-blocking warning: coordinate export still covers every fragment,
        # only the ChimeraX heatmap/marker scripts are skipped for it.
        if result["fragment_count"] > 1:
            self._heatmap_protein_error_label.configure(
                text=f"⚠  {gene} ({result['uniprot']}) spans multiple AlphaFold "
                     f"fragments — ChimeraX heatmaps/markers will be skipped for it "
                     f"(coordinate export still covers all fragments).",
                text_color=_YELLOW,
            )
            self._heatmap_protein_error_label.grid()
        else:
            self._heatmap_protein_error_label.grid_remove()

    def _remove_heatmap_protein(self, token: str) -> None:
        if token in self._heatmap_proteins:
            self._heatmap_proteins.remove(token)
            self._refresh_heatmap_protein_chips()

    def _refresh_heatmap_protein_chips(self) -> None:
        """Redraw the added-proteins chip list, each removable via its own ✕ button."""
        for w in self._heatmap_proteins_list_frame.winfo_children():
            w.destroy()

        if not self._heatmap_proteins:
            ctk.CTkLabel(
                self._heatmap_proteins_list_frame, text="None added yet.",
                text_color=_GRAY, font=ctk.CTkFont(size=11),
            ).pack(anchor="w")
            return

        row = None
        for i, token in enumerate(self._heatmap_proteins):
            if i % 3 == 0:
                row = ctk.CTkFrame(self._heatmap_proteins_list_frame, fg_color="transparent")
                row.pack(anchor="w", pady=(0, 4))
            chip = ctk.CTkFrame(row, fg_color=_CHIP, corner_radius=6)
            chip.pack(side="left", padx=(0, 6))
            ctk.CTkLabel(chip, text=token, font=ctk.CTkFont(size=12)).pack(
                side="left", padx=(8, 4), pady=4,
            )
            ctk.CTkButton(
                chip, text="✕", width=20, height=20, fg_color="transparent",
                hover_color=_CHIP_HOVER, font=ctk.CTkFont(size=11),
                command=lambda t=token: self._remove_heatmap_protein(t),
            ).pack(side="left", padx=(0, 6), pady=4)

    # ── Results area ─────────────────────────────────────────────────────────

    def _build_heatmap_results_panel(self, parent) -> None:
        """The scrollable list of exported proteins, shown over the empty-state hint."""
        self._heatmap_results_frame = ctk.CTkScrollableFrame(parent, fg_color=_PLOT_BG)
        self._heatmap_results_frame.grid(row=0, column=0, sticky="nsew")
        self._heatmap_results_frame.grid_columnconfigure(1, weight=1)
        self._heatmap_results_frame.grid_remove()
        if not hasattr(self, "_heatmap_items"):
            self._heatmap_items = []

    def _refresh_heatmap_results(self) -> None:
        """Redraw the exported-proteins list from self._heatmap_items."""
        frame = self._heatmap_results_frame
        for w in frame.winfo_children():
            w.destroy()
        for i, (label, folder, error) in enumerate(self._heatmap_items):
            ok = error is None
            ctk.CTkLabel(
                frame, text="✓" if ok else "✗", text_color=_GREEN if ok else _RED,
                font=ctk.CTkFont(size=14, weight="bold"), width=24,
            ).grid(row=i, column=0, padx=(10, 4), pady=6, sticky="nw")
            ctk.CTkLabel(
                frame, text=label if ok else f"{label}\n{error}", anchor="w", justify="left",
                text_color=None if ok else _RED, wraplength=520,
            ).grid(row=i, column=1, padx=4, pady=6, sticky="w")
            if folder is not None:
                ctk.CTkButton(
                    frame, text="Open folder", width=90, height=26, font=ctk.CTkFont(size=12),
                    fg_color=_BTN2, hover_color=_BTN2_HOVER,
                    command=lambda p=folder: self._open_folder(p),
                ).grid(row=i, column=2, padx=10, pady=6, sticky="ne")
        if self._heatmap_items:
            self._hm_empty_label.place_forget()
            frame.grid()
        else:
            frame.grid_remove()
            self._hm_empty_label.place(relx=0.5, rely=0.5, anchor="center")
            self._hm_empty_label.lift()

    # ── Validation and run ───────────────────────────────────────────────────

    def _validate_heatmap_run(self) -> bool:
        """Check the inputs a Structure Heatmaps run needs, showing a dialog and
        returning False for the first problem found."""
        from tkinter import messagebox

        input_problems = self._validate_input_files("structure-heatmaps")
        if input_problems:
            messagebox.showerror(
                "Input File Problem",
                "Cannot start the export — one or more input files look invalid:\n\n"
                + "\n".join(f"  • {p}" for p in input_problems)
                + "\n\nCheck that the correct file is in the input folder and try again.",
            )
            return False

        if self._heatmap_source_var.get() == "Upload CIF file":
            if not self._heatmap_custom_cif_var.get().strip():
                messagebox.showwarning("Missing input", "Select a CIF file to upload.")
                return False
            if not self._heatmap_custom_cif_uniprot_var.get().strip():
                messagebox.showwarning(
                    "Missing input",
                    "Enter the UniProt ID for this CIF file (could not auto-detect it).",
                )
                return False
        elif not self._heatmap_proteins:
            messagebox.showwarning("Missing input", "Add at least one gene symbol or UniProt accession.")
            return False
        return True

    def _start_heatmap_run(self) -> None:
        """Validate and launch a Structure Heatmaps export."""
        if self._running or not self._validate_heatmap_run():
            return

        self._running = True
        self._hm_run_active = True
        self._stop_requested = False
        self._suspended = False
        self._set_run_controls(True)

        self._heatmap_items = []  # the previous export's list no longer applies
        self._refresh_heatmap_results()
        self._hm_progress_bar.set(0)
        self._hm_progress_bar.grid_remove()
        self._hm_status_label.configure(text="●  Starting…", text_color=_GRAY)
        self._hm_log.configure(state="normal")
        self._hm_log.delete("1.0", "end")
        self._hm_log.configure(state="disabled")

        threading.Thread(target=self._run_structure_heatmaps, daemon=True).start()

    def _heatmap_kwargs(self) -> dict:
        """Heatmap/marker options shared by both the batch and custom-CIF
        export paths, read fresh from the fields' variables.
        """
        return dict(
            mutation_heatmap=self._heatmap_mutation_var.get(),
            residue_heatmap=self._heatmap_residue_var.get(),
            plddt_heatmap=self._heatmap_plddt_var.get(),
            mark_ptm_sites=self._heatmap_mark_ptm_var.get(),
            mark_mutations=self._heatmap_mark_mutations_var.get(),
            log_scale=self._heatmap_log_scale_var.get(),
            dim_low_confidence=self._heatmap_dim_confidence_var.get(),
            mutation_low_color=self._heatmap_mutation_low_var.get(),
            mutation_mid_color=self._heatmap_mutation_mid_var.get(),
            mutation_high_color=self._heatmap_mutation_high_var.get(),
            plddt_low_color=self._heatmap_plddt_low_var.get(),
            plddt_mid_color=self._heatmap_plddt_mid_var.get(),
            plddt_high_color=self._heatmap_plddt_high_var.get(),
            ptm_marker_color=self._heatmap_ptm_marker_color_var.get(),
            mutation_marker_color=self._heatmap_mutation_marker_color_var.get(),
        )

    def _run_structure_heatmaps(self) -> None:
        """Dispatch to the batch (database-backed) or single custom-CIF
        export, based on the Source toggle. Runs in the background thread.
        """
        if self._heatmap_source_var.get() == "Upload CIF file":
            self._run_heatmaps_custom_cif()
        else:
            self._run_heatmaps_batch()

    def _run_heatmaps_batch(self) -> None:
        """Run the CA-coordinate export in-process, in the background thread.

        Batches over every protein added to the list, sharing one progress
        bar/status label across the whole batch (this tool has no fixed
        multi-step layout to give each protein its own row, unlike the main
        pipeline) -- see export_ca_coordinates.run_batch_export for how a
        single protein's failure is contained rather than aborting the rest.
        """
        tokens = list(self._heatmap_proteins)
        output_dir = Path(self._heatmap_output_dir_var.get())
        kwargs = self._heatmap_kwargs()

        self._q("pipeline_start", 1, "structure-heatmaps", "warm")
        self._q("show_log")
        self._q("status", 0, "▶  Exporting…", _BLUE)
        self._q("show_progress", 0)
        self._q("log", f"Exporting CA coordinates for {len(tokens)} protein(s): {', '.join(tokens)}")
        self._q("log", "")

        def _progress(i: int, total: int, token: str) -> None:
            self._q("progress", 0, i / total, f"▶  Exporting {i}/{total}: {token}…")

        t0 = time.time()
        try:
            from export_ca_coordinates import run_batch_export
            items = run_batch_export(
                tokens, output_dir=output_dir, progress_cb=_progress,
                log_cb=lambda line: self._q("log", line), **kwargs,
            )
        except ImportError as exc:
            self._fail_heatmap_run(f"Missing dependency: {exc}. Run: uv sync")
            return
        except (FileNotFoundError, ValueError, RuntimeError) as exc:
            self._fail_heatmap_run(f"Error: {exc}")
            return
        except Exception as exc:
            self._fail_heatmap_run(f"Unexpected error: {exc}")
            return

        elapsed = time.time() - t0
        n_ok = sum(1 for item in items if item.result is not None)
        n_total = len(items)
        self._q("hide_progress", 0)
        if n_ok == n_total:
            self._q("status", 0, f"✓  {n_ok}/{n_total} in {_fmt_time(elapsed)}", _GREEN)
        elif n_ok > 0:
            self._q("status", 0, f"⚠  {n_ok}/{n_total} succeeded", _YELLOW)
        else:
            self._q("status", 0, "✗  Failed (0 succeeded)", _RED)
        self._q("heatmap_results", [
            (f"{item.result.gene} ({item.result.uid})", item.result.all_out.parent, None)
            if item.result is not None else (item.token, None, item.error or "Export failed")
            for item in items
        ])
        self._q("finished")

    def _run_heatmaps_custom_cif(self) -> None:
        """Run a single export against a caller-provided CIF file (e.g. a
        seeded AlphaFold Server prediction), in-process, in the background
        thread -- see export_ca_coordinates.run_export's custom_cif_path
        docstring for why this bypasses the AlphaFold DB download/cache
        entirely.
        """
        cif_path = Path(self._heatmap_custom_cif_var.get().strip())
        uniprot = self._heatmap_custom_cif_uniprot_var.get().strip()
        output_dir = Path(self._heatmap_output_dir_var.get())
        kwargs = self._heatmap_kwargs()

        self._q("pipeline_start", 1, "structure-heatmaps", "warm")
        self._q("show_log")
        self._q("status", 0, "▶  Exporting…", _BLUE)
        self._q("show_progress", 0)
        self._q("log", f"Exporting Structure Heatmaps for custom CIF: {cif_path}")
        self._q("log", "")

        t0 = time.time()
        try:
            from export_ca_coordinates import run_export
            result = run_export(
                uniprot=uniprot, custom_cif_path=cif_path, output_dir=output_dir,
                log_cb=lambda line: self._q("log", line), **kwargs,
            )
        except ImportError as exc:
            self._fail_heatmap_run(f"Missing dependency: {exc}. Run: uv sync")
            return
        except (FileNotFoundError, ValueError, RuntimeError) as exc:
            self._fail_heatmap_run(f"Error: {exc}")
            return
        except Exception as exc:
            self._fail_heatmap_run(f"Unexpected error: {exc}")
            return

        elapsed = time.time() - t0
        self._q("hide_progress", 0)
        self._q("status", 0, f"✓  {result.gene} ({result.uid}) in {_fmt_time(elapsed)}", _GREEN)
        self._q("heatmap_results", [(f"{result.gene} ({result.uid})", result.all_out.parent, None)])
        self._q("finished")

    def _fail_heatmap_run(self, message: str) -> None:
        """Report a run that stopped before producing any results."""
        self._q("status", 0, "✗  Failed", _RED)
        self._q("log", message)
        self._q("hide_progress", 0)
        self._q("finished")
