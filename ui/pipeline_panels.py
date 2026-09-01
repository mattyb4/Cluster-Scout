"""Pipeline tab widget building: mode panels, settings, browse dialogs.

The execution engine (run/stop/resume, subprocess streaming, the progress
queue) lives in ui/pipeline_runner.py — this file only builds widgets and
handles simple input events (browsing for a file, toggling the log).
"""
from __future__ import annotations

import shutil
from pathlib import Path
from tkinter import filedialog

import customtkinter as ctk
from PIL import Image

from ui.common import (
    _GRAY,
    _GREEN,
    _INPUT_FOLDERS,
    _MODE_HELP,
    _RED,
    _YELLOW,
    MUTATION_CLUSTERING_STEPS,
    OUTPUT_DIR,
    PROJECT_ROOT,
    PTM_PROXIMITY_STEPS,
    add_resize_grip,
    color_swatch_button,
    extract_uniprot_from_cif,
    help_icon,
    isolate_textbox_scroll,
    resolve_input_file,
)


class PipelineTabMixin:
    def _build_pipeline_tab(self, tab) -> None:
        tab.grid_columnconfigure(0, weight=1)
        tab.grid_rowconfigure(0, weight=1)

        scroll = ctk.CTkScrollableFrame(tab, fg_color="transparent")
        scroll.grid(row=0, column=0, sticky="nsew")
        scroll.grid_columnconfigure(0, weight=1)

        p = scroll  # all pipeline widgets go in the scrollable frame

        # Title
        _logo_path_dark = PROJECT_ROOT / "cluster_scout_logo_dark.png"
        if _logo_path_dark.exists():
            _pil_dark = Image.open(_logo_path_dark)
            _h = 200
            _w = int(_pil_dark.width * _h / _pil_dark.height)
            _logo_img = ctk.CTkImage(
                light_image=_pil_dark,
                dark_image=_pil_dark,
                size=(_w, _h),
            )
            ctk.CTkLabel(p, image=_logo_img, text="").grid(
                row=0, column=0, padx=24, pady=(12, 4), sticky="w"
            )
        else:
            ctk.CTkLabel(
                p,
                text="Cluster-Scout",
                font=ctk.CTkFont(size=22, weight="bold"),
            ).grid(row=0, column=0, padx=24, pady=(12, 4), sticky="w")

        # Data-file status bar with Browse buttons
        self._file_frame = ctk.CTkFrame(p)
        self._file_frame.grid(row=1, column=0, padx=24, pady=4, sticky="ew")
        self._file_frame.grid_columnconfigure(1, weight=1)
        self._file_indicators: dict[str, ctk.CTkLabel] = {}
        self._file_buttons: dict[str, ctk.CTkButton] = {}

        ctk.CTkLabel(
            self._file_frame,
            text="Input files:",
            font=ctk.CTkFont(weight="bold"),
        ).grid(row=0, column=0, columnspan=3, padx=12, pady=(8, 2), sticky="w")

        for i, (name, (folder, exts, desc, validator)) in enumerate(_INPUT_FOLDERS.items(), 1):
            lbl = ctk.CTkLabel(self._file_frame, text=f"{name} …", anchor="w")
            lbl.grid(row=i, column=0, columnspan=2, padx=(12, 6), pady=3, sticky="ew")
            self._file_indicators[name] = lbl

            filetypes = [(desc, " ".join(f"*{e}" for e in exts))]
            btn = ctk.CTkButton(
                self._file_frame,
                text="Browse",
                width=70,
                height=26,
                font=ctk.CTkFont(size=12),
                command=lambda n=name, f=folder, ft=filetypes, v=validator: self._browse_file(n, f, ft, v),
            )
            btn.grid(row=i, column=2, padx=12, pady=3, sticky="e")
            self._file_buttons[name] = btn

        # Mode selection
        mode_frame = ctk.CTkFrame(p)
        mode_frame.grid(row=2, column=0, padx=24, pady=4, sticky="ew")

        ctk.CTkLabel(
            mode_frame, text="Mode:", font=ctk.CTkFont(weight="bold")
        ).pack(side="left", padx=(12, 8), pady=10)

        self._mode = ctk.StringVar(value="ptm-proximity")
        for label, value in [
            ("PTM Proximity", "ptm-proximity"),
            ("Mutation Clustering", "mutation-clustering"),
            ("Single Protein", "single-protein"),
            ("Structure Heatmaps", "ca-coordinates"),
        ]:
            ctk.CTkRadioButton(
                mode_frame,
                text=label,
                variable=self._mode,
                value=value,
                command=self._rebuild_step_rows,
            ).pack(side="left", padx=(8, 0), pady=10)
            help_icon(mode_frame, _MODE_HELP[value]).pack(side="left", padx=(4, 8), pady=10)

        # Output folder selector
        out_frame = ctk.CTkFrame(p)
        out_frame.grid(row=3, column=0, padx=24, pady=4, sticky="ew")
        out_frame.grid_columnconfigure(1, weight=1)

        ctk.CTkLabel(
            out_frame, text="Output folder:", font=ctk.CTkFont(weight="bold"),
        ).grid(row=0, column=0, padx=(12, 6), pady=8, sticky="w")

        self._output_dir_var = ctk.StringVar(value=str(OUTPUT_DIR))
        self._output_dir_entry = ctk.CTkEntry(
            out_frame, textvariable=self._output_dir_var, state="readonly",
        )
        self._output_dir_entry.grid(row=0, column=1, padx=6, pady=8, sticky="ew")

        ctk.CTkButton(
            out_frame, text="Change", width=70, height=26,
            font=ctk.CTkFont(size=12),
            command=self._browse_output_dir,
        ).grid(row=0, column=2, padx=(6, 4), pady=8, sticky="e")

        ctk.CTkButton(
            out_frame, text="Reset", width=60, height=26,
            font=ctk.CTkFont(size=12), fg_color="gray30", hover_color="gray40",
            command=lambda: self._output_dir_var.set(str(OUTPUT_DIR)),
        ).grid(row=0, column=3, padx=(0, 12), pady=8, sticky="e")

        # Pipeline settings (ptm-proximity, mutation-clustering, AND single-protein
        # -- analyze_single_cif_nearby_mutations.py accepts all four of these too)
        settings_frame = self._settings_frame = ctk.CTkFrame(p)
        settings_frame.grid(row=4, column=0, padx=24, pady=4, sticky="ew")

        ctk.CTkLabel(
            settings_frame, text="Settings:", font=ctk.CTkFont(weight="bold"),
        ).pack(side="left", padx=(12, 8), pady=8)

        ctk.CTkLabel(settings_frame, text="Cutoff (Å):").pack(side="left", padx=(8, 4), pady=8)
        help_icon(
            settings_frame,
            "Maximum 3D distance, in Ångströms, between a PTM site and a "
            "mutation for them to be counted as \"nearby.\" Default: 10.0 Å.",
        ).pack(side="left", padx=(0, 4), pady=8)
        self._cutoff_var = ctk.StringVar(value="10.0")
        ctk.CTkEntry(
            settings_frame, textvariable=self._cutoff_var, width=60,
        ).pack(side="left", padx=(0, 16), pady=8)

        ctk.CTkLabel(settings_frame, text="Min samples:").pack(side="left", padx=(8, 4), pady=8)
        help_icon(
            settings_frame,
            "Minimum number of distinct COSMIC patient samples a mutation "
            "must appear in to count as a recurrent hotspot. Lower values "
            "include rarer mutations; higher values restrict to mutations "
            "seen more often. In Single Protein mode, this can only tighten "
            "the threshold already applied when the input TSV was built -- "
            "mutations below that original threshold aren't in the data at all.",
        ).pack(side="left", padx=(0, 4), pady=8)
        self._min_samples_var = ctk.StringVar(value="3")
        ctk.CTkEntry(
            settings_frame, textvariable=self._min_samples_var, width=60,
        ).pack(side="left", padx=(0, 16), pady=8)

        ctk.CTkLabel(settings_frame, text="Min pLDDT:").pack(side="left", padx=(8, 4), pady=8)
        help_icon(
            settings_frame,
            "Minimum AlphaFold per-residue confidence score (0-100). "
            "Residues below this threshold are excluded, since the model "
            "is less certain about their position. Leave blank to disable.",
        ).pack(side="left", padx=(0, 4), pady=8)
        self._min_plddt_var = ctk.StringVar(value="")
        ctk.CTkEntry(
            settings_frame, textvariable=self._min_plddt_var, width=60,
            placeholder_text="off",
        ).pack(side="left", padx=(0, 16), pady=8)

        ctk.CTkLabel(settings_frame, text="Max PAE:").pack(side="left", padx=(8, 4), pady=8)
        help_icon(
            settings_frame,
            "Maximum Predicted Aligned Error, in Ångströms, allowed between "
            "a PTM site and mutation pair. Filters out pairs where AlphaFold "
            "isn't confident about their relative position, even if the raw "
            "distance looks close. Leave blank to disable.",
        ).pack(side="left", padx=(0, 4), pady=8)
        self._max_pae_var = ctk.StringVar(value="")
        ctk.CTkEntry(
            settings_frame, textvariable=self._max_pae_var, width=60,
            placeholder_text="off",
        ).pack(side="left", padx=(0, 12), pady=8)

        # PolyPhen filter
        pp_frame = self._pp_frame = ctk.CTkFrame(p)
        pp_frame.grid(row=5, column=0, padx=24, pady=4, sticky="ew")

        ctk.CTkLabel(
            pp_frame, text="PolyPhen filter:", font=ctk.CTkFont(weight="bold"),
        ).pack(side="left", padx=(12, 4), pady=8)

        self._pp_benign_var = ctk.BooleanVar(value=True)
        ctk.CTkCheckBox(
            pp_frame, text="Benign",
            variable=self._pp_benign_var,
            checkbox_width=18, checkbox_height=18,
        ).pack(side="left", padx=8, pady=8)

        self._pp_possibly_var = ctk.BooleanVar(value=True)
        ctk.CTkCheckBox(
            pp_frame, text="Possibly damaging",
            variable=self._pp_possibly_var,
            checkbox_width=18, checkbox_height=18,
        ).pack(side="left", padx=8, pady=8)

        self._pp_probably_var = ctk.BooleanVar(value=True)
        ctk.CTkCheckBox(
            pp_frame, text="Probably damaging",
            variable=self._pp_probably_var,
            checkbox_width=18, checkbox_height=18,
        ).pack(side="left", padx=(8, 12), pady=8)

        # Steps panel
        self._steps_outer = ctk.CTkFrame(p)
        self._steps_outer.grid(row=7, column=0, padx=24, pady=4, sticky="ew")
        self._steps_outer.grid_columnconfigure(1, weight=1)
        self._rebuild_step_rows()

        # Buttons
        btn_frame = ctk.CTkFrame(p, fg_color="transparent")
        btn_frame.grid(row=8, column=0, padx=24, pady=8, sticky="ew")

        self._run_btn = ctk.CTkButton(
            btn_frame,
            text="▶  Run Pipeline",
            command=self._start_pipeline,
            width=160,
            height=44,
            font=ctk.CTkFont(size=15, weight="bold"),
        )
        self._run_btn.pack(side="left", padx=(0, 8))

        self._stop_btn = ctk.CTkButton(
            btn_frame,
            text="■  Stop",
            command=self._stop_pipeline,
            width=110,
            height=44,
            font=ctk.CTkFont(size=15, weight="bold"),
            state="disabled",
            fg_color="gray30",
            hover_color=_RED,
        )
        self._stop_btn.pack(side="left", padx=(0, 12))

        self._open_btn = ctk.CTkButton(
            btn_frame,
            text="Open Output Folder",
            command=self._open_output_folder,
            width=180,
            height=44,
        )
        self._open_btn.pack(side="left")

        ctk.CTkButton(
            btn_frame,
            text="Manage Cache",
            command=self._manage_cache,
            width=140,
            height=44,
            fg_color="gray30",
            hover_color="gray40",
        ).pack(side="left", padx=(8, 0))

        self._timer_label = ctk.CTkLabel(
            btn_frame,
            text="",
            text_color=_GRAY,
            font=ctk.CTkFont(size=13),
        )
        self._timer_label.pack(side="right", padx=12)

        import tkinter as tk
        self._activity_width = 40
        self._activity_height = 6
        self._activity_chunk = 12
        self._activity_canvas = tk.Canvas(
            btn_frame, width=self._activity_width, height=self._activity_height,
            highlightthickness=0, bg="#333333", bd=0,
        )
        self._activity_canvas.pack(side="right", padx=(0, 4))
        self._activity_canvas.pack_forget()
        self._activity_animating = False
        self._activity_pos = -self._activity_chunk

        # Log (collapsible)
        self._log_visible = False
        self._log_toggle = ctk.CTkButton(
            p,
            text="Show Details",
            width=120,
            height=28,
            font=ctk.CTkFont(size=12),
            fg_color="gray30",
            hover_color="gray40",
            command=self._toggle_log,
        )
        self._log_toggle.grid(row=9, column=0, padx=24, pady=(8, 0), sticky="w")

        self._log_frame = ctk.CTkFrame(p, fg_color="transparent")
        self._log = ctk.CTkTextbox(
            self._log_frame,
            height=140,
            font=ctk.CTkFont(family="Courier New", size=12),
            wrap="word",
            state="disabled",
        )
        self._log.pack(fill="both", expand=True)
        isolate_textbox_scroll(self._log)
        add_resize_grip(self._log).pack(fill="x")
        self._toggle_log()  # visible by default

        # Hide the scrollbar when content fits; show it only when scrolling is needed
        def _update_scrollbar(*_):
            try:
                canvas = scroll._parent_canvas
                sr = canvas.cget("scrollregion")
                if not sr:
                    scroll._scrollbar.grid_remove()
                    return
                content_h = int(float(sr.split()[3]))
                if content_h > canvas.winfo_height():
                    scroll._scrollbar.grid()
                else:
                    scroll._scrollbar.grid_remove()
            except Exception:
                pass

        # add="+" is required here: CTkScrollableFrame already binds its own
        # <Configure> handler (to keep the canvas scrollregion in sync with
        # content size) — a plain .bind() would replace it instead of adding
        # to it, freezing the scrollregion at whatever it was when this ran.
        scroll.bind("<Configure>", _update_scrollbar, add="+")
        tab.bind("<Configure>", _update_scrollbar, add="+")
        self.after(200, _update_scrollbar)

    def _rebuild_step_rows(self):
        for w in self._steps_outer.winfo_children():
            w.destroy()
        self._step_status_labels = []
        self._step_progress_bars: list[ctk.CTkProgressBar] = []
        self._steps_outer.grid_columnconfigure(1, weight=1)

        mode = self._mode.get()

        # Settings frame also applies to single-protein (accepted by
        # analyze_single_cif_nearby_mutations.py); PolyPhen filter feeds step 4 of
        # both ptm-proximity and mutation-clustering. Neither applies to ca-coordinates.
        if mode in ("ptm-proximity", "mutation-clustering", "single-protein"):
            self._settings_frame.grid()
        else:
            self._settings_frame.grid_remove()
        if mode in ("ptm-proximity", "mutation-clustering"):
            self._pp_frame.grid()
        else:
            self._pp_frame.grid_remove()

        if mode == "single-protein":
            self._build_single_protein_panel()
            return
        if mode == "ca-coordinates":
            self._build_ca_coordinates_panel()
            return

        steps = PTM_PROXIMITY_STEPS if mode == "ptm-proximity" else MUTATION_CLUSTERING_STEPS

        ctk.CTkLabel(
            self._steps_outer,
            text="Steps",
            font=ctk.CTkFont(weight="bold"),
        ).grid(row=0, column=0, columnspan=4, padx=12, pady=(8, 2), sticky="w")

        for i, (panel_label, _log_label) in enumerate(steps, 1):
            ctk.CTkLabel(self._steps_outer, text=f"  {i}.", width=28).grid(
                row=i, column=0, padx=(12, 0), pady=5, sticky="w"
            )
            ctk.CTkLabel(self._steps_outer, text=panel_label, anchor="w").grid(
                row=i, column=1, padx=6, pady=5, sticky="ew"
            )

            bar = ctk.CTkProgressBar(self._steps_outer, width=120, height=14)
            bar.set(0)
            bar.grid(row=i, column=2, padx=6, pady=5, sticky="e")
            bar.grid_remove()
            self._step_progress_bars.append(bar)

            status = ctk.CTkLabel(
                self._steps_outer,
                text="●  Waiting",
                width=100,
                anchor="e",
                text_color=_GRAY,
            )
            status.grid(row=i, column=3, padx=12, pady=5, sticky="e")
            self._step_status_labels.append(status)

    def _build_single_protein_panel(self):
        """Build the input fields for single-protein analysis mode."""
        ctk.CTkLabel(
            self._steps_outer,
            text="Single Protein Analysis",
            font=ctk.CTkFont(weight="bold"),
        ).grid(row=0, column=0, columnspan=3, padx=12, pady=(8, 2), sticky="w")

        # CIF file picker
        ctk.CTkLabel(self._steps_outer, text="CIF file:", anchor="w").grid(
            row=1, column=0, padx=(12, 6), pady=6, sticky="w"
        )
        self._single_cif_var = ctk.StringVar(value="")
        self._single_cif_entry = ctk.CTkEntry(
            self._steps_outer, textvariable=self._single_cif_var, width=400,
        )
        self._single_cif_entry.grid(row=1, column=1, padx=6, pady=6, sticky="ew")
        ctk.CTkButton(
            self._steps_outer, text="Browse", width=70, height=26,
            font=ctk.CTkFont(size=12),
            command=self._browse_cif,
        ).grid(row=1, column=2, padx=12, pady=6, sticky="e")

        # UniProt ID
        ctk.CTkLabel(self._steps_outer, text="UniProt ID:", anchor="w").grid(
            row=2, column=0, padx=(12, 6), pady=6, sticky="w"
        )
        self._single_uniprot_var = ctk.StringVar(value="")
        ctk.CTkEntry(
            self._steps_outer, textvariable=self._single_uniprot_var, width=200,
            placeholder_text="Edit if CIF is not in a UniProt-named folder",
        ).grid(row=2, column=1, padx=6, pady=6, sticky="w")

        # Status label (reuse the step status pattern so the pipeline runner can update it)
        status = ctk.CTkLabel(
            self._steps_outer, text="●  Ready", width=100,
            anchor="e", text_color=_GRAY,
        )
        status.grid(row=3, column=1, columnspan=2, padx=12, pady=6, sticky="e")
        self._step_status_labels.append(status)

        bar = ctk.CTkProgressBar(self._steps_outer, width=120, height=14)
        bar.set(0)
        bar.grid(row=3, column=0, padx=12, pady=6, sticky="w")
        bar.grid_remove()
        self._step_progress_bars.append(bar)

    def _browse_cif(self):
        """Open a file dialog for selecting a CIF file and auto-fill the UniProt ID."""
        path = filedialog.askopenfilename(
            title="Select CIF structure file",
            filetypes=[("CIF files", "*.cif"), ("All files", "*.*")],
        )
        if not path:
            return
        self._single_cif_var.set(path)
        uid = extract_uniprot_from_cif(Path(path))
        if uid:
            self._single_uniprot_var.set(uid)
        else:
            parent_name = Path(path).parent.name
            self._single_uniprot_var.set(parent_name)

    def _build_ca_coordinates_panel(self):
        """Build the input fields for Structure Heatmaps mode (CA-coordinate export)."""
        ctk.CTkLabel(
            self._steps_outer,
            text="Structure Heatmaps",
            font=ctk.CTkFont(weight="bold"),
        ).grid(row=0, column=0, columnspan=3, padx=12, pady=(8, 2), sticky="w")

        # Source toggle: batch of database-backed proteins, or one caller-provided CIF
        source_frame = ctk.CTkFrame(self._steps_outer, fg_color="transparent")
        source_frame.grid(row=1, column=0, columnspan=3, padx=12, pady=(0, 4), sticky="w")
        ctk.CTkLabel(source_frame, text="Source:", font=ctk.CTkFont(weight="bold")).pack(side="left", padx=(0, 8))
        if not hasattr(self, "_ca_source_var"):
            self._ca_source_var = ctk.StringVar(value="Database")
        ctk.CTkSegmentedButton(
            source_frame, values=["Database", "Upload CIF file"],
            variable=self._ca_source_var, command=self._on_ca_source_change,
        ).pack(side="left")
        help_icon(
            source_frame,
            "Database exports one or more proteins by gene/UniProt, using "
            "AlphaFold DB's model (downloaded automatically if needed). "
            "Upload CIF file instead uses a specific .cif you already have "
            "-- e.g. a seeded AlphaFold Server prediction from the CIF "
            "Variance tool's \"Generate AlphaFold Seeds JSON\" option -- for "
            "a single protein, with nothing downloaded or added to the "
            "shared structure cache.",
        ).pack(side="left", padx=(6, 0))

        # ── Database source: proteins list (gene symbols and/or UniProt accessions) ──
        self._ca_database_frame = ctk.CTkFrame(self._steps_outer, fg_color="transparent")
        self._ca_database_frame.grid(row=2, column=0, columnspan=3, sticky="ew")
        self._ca_database_frame.grid_columnconfigure(1, weight=1)

        proteins_label_frame = ctk.CTkFrame(self._ca_database_frame, fg_color="transparent")
        proteins_label_frame.grid(row=0, column=0, padx=(12, 6), pady=6, sticky="nw")
        ctk.CTkLabel(proteins_label_frame, text="Proteins:", anchor="w").pack(side="left")
        help_icon(
            proteins_label_frame,
            "Add one or more gene symbols and/or UniProt accessions - each "
            "is exported in turn (its own Output/coordinates/{gene}_{UniProt}/ "
            "folder), with the same options below applied to every one. A "
            "gene's UniProt accession, or a UniProt accession's gene "
            "symbol, is resolved automatically as needed. One protein "
            "failing (e.g. no AlphaFold model) doesn't stop the rest.",
        ).pack(side="left", padx=(4, 0))
        if not hasattr(self, "_ca_proteins"):
            self._ca_proteins: list[str] = []

        protein_input_frame = ctk.CTkFrame(self._ca_database_frame, fg_color="transparent")
        protein_input_frame.grid(row=0, column=1, columnspan=2, padx=6, pady=6, sticky="w")

        self._ca_protein_input_var = ctk.StringVar(value="")
        protein_entry = ctk.CTkEntry(
            protein_input_frame, textvariable=self._ca_protein_input_var, width=180,
            placeholder_text="e.g. P04637 or TP53",
        )
        protein_entry.pack(side="left", padx=(0, 6))
        protein_entry.bind("<Return>", lambda _e: self._add_ca_protein())

        ctk.CTkButton(
            protein_input_frame, text="+ Add", width=70, height=28,
            command=self._add_ca_protein,
        ).pack(side="left")

        # Feedback for _add_ca_protein — hidden until there's something to say
        self._ca_protein_error_label = ctk.CTkLabel(
            self._ca_database_frame, text="", text_color=_RED,
            font=ctk.CTkFont(size=11), anchor="w",
        )
        self._ca_protein_error_label.grid(row=1, column=1, columnspan=2, padx=6, pady=(0, 4), sticky="w")
        self._ca_protein_error_label.grid_remove()

        self._ca_proteins_list_frame = ctk.CTkFrame(self._ca_database_frame, fg_color="transparent")
        self._ca_proteins_list_frame.grid(row=2, column=0, columnspan=3, padx=12, pady=(0, 6), sticky="ew")
        self._refresh_ca_protein_chips()

        # ── Upload-CIF source: one specific file + its UniProt ID ──
        self._ca_custom_cif_frame = ctk.CTkFrame(self._steps_outer, fg_color="transparent")
        self._ca_custom_cif_frame.grid(row=2, column=0, columnspan=3, sticky="ew")
        self._ca_custom_cif_frame.grid_columnconfigure(1, weight=1)
        self._ca_custom_cif_frame.grid_remove()  # Database is the default source

        ctk.CTkLabel(self._ca_custom_cif_frame, text="CIF file:", anchor="w").grid(
            row=0, column=0, padx=(12, 6), pady=6, sticky="w"
        )
        if not hasattr(self, "_ca_custom_cif_var"):
            self._ca_custom_cif_var = ctk.StringVar(value="")
        ctk.CTkEntry(
            self._ca_custom_cif_frame, textvariable=self._ca_custom_cif_var,
        ).grid(row=0, column=1, padx=6, pady=6, sticky="ew")
        ctk.CTkButton(
            self._ca_custom_cif_frame, text="Browse", width=70, height=26,
            font=ctk.CTkFont(size=12),
            command=self._browse_ca_custom_cif,
        ).grid(row=0, column=2, padx=12, pady=6, sticky="e")

        ctk.CTkLabel(self._ca_custom_cif_frame, text="UniProt ID:", anchor="w").grid(
            row=1, column=0, padx=(12, 6), pady=(0, 6), sticky="w"
        )
        if not hasattr(self, "_ca_custom_cif_uniprot_var"):
            self._ca_custom_cif_uniprot_var = ctk.StringVar(value="")
        ctk.CTkEntry(
            self._ca_custom_cif_frame, textvariable=self._ca_custom_cif_uniprot_var, width=200,
            placeholder_text="auto-detected, or enter if not found",
        ).grid(row=1, column=1, padx=6, pady=(0, 6), sticky="w")

        # Both frames above default to gridded/removed for a fresh "Database"
        # source -- resync to whichever source was actually last selected,
        # since this whole panel (and both frames) is torn down and rebuilt
        # on every mode switch, but self._ca_source_var itself persists.
        self._on_ca_source_change()

        from export_ca_coordinates import (
            MUTATION_DEFAULT_HIGH_COLOR,
            MUTATION_DEFAULT_LOW_COLOR,
            MUTATION_MARKER_DEFAULT_COLOR,
            PLDDT_DEFAULT_HIGH_COLOR,
            PLDDT_DEFAULT_LOW_COLOR,
            PTM_MARKER_DEFAULT_COLOR,
        )
        if not hasattr(self, "_ca_mutation_low_var"):
            self._ca_mutation_low_var = ctk.StringVar(value=MUTATION_DEFAULT_LOW_COLOR)
            self._ca_mutation_high_var = ctk.StringVar(value=MUTATION_DEFAULT_HIGH_COLOR)
            self._ca_plddt_low_var = ctk.StringVar(value=PLDDT_DEFAULT_LOW_COLOR)
            self._ca_plddt_high_var = ctk.StringVar(value=PLDDT_DEFAULT_HIGH_COLOR)
            self._ca_ptm_marker_color_var = ctk.StringVar(value=PTM_MARKER_DEFAULT_COLOR)
            self._ca_mutation_marker_color_var = ctk.StringVar(value=MUTATION_MARKER_DEFAULT_COLOR)

        def _reset_colors(*pairs: tuple[ctk.StringVar, str]) -> None:
            for var, default in pairs:
                var.set(default)
            self._rebuild_step_rows()

        # Heatmaps section
        ctk.CTkLabel(
            self._steps_outer, text="Heatmaps:", font=ctk.CTkFont(weight="bold"),
        ).grid(row=3, column=0, columnspan=3, padx=12, pady=(10, 2), sticky="w")

        # Mutation heatmap
        if not hasattr(self, "_ca_mutation_heatmap_var"):
            self._ca_mutation_heatmap_var = ctk.BooleanVar(value=True)
        mut_heatmap_frame = ctk.CTkFrame(self._steps_outer, fg_color="transparent")
        mut_heatmap_frame.grid(row=4, column=0, columnspan=3, padx=24, pady=2, sticky="w")
        ctk.CTkCheckBox(
            mut_heatmap_frame, text="Mutation heatmap (patients within 10 Å)",
            variable=self._ca_mutation_heatmap_var,
            checkbox_width=18, checkbox_height=18,
            command=self._on_ca_mutation_heatmap_toggle,
        ).pack(side="left")
        help_icon(
            mut_heatmap_frame,
            "Colors the ChimeraX cartoon by COSMIC patient count summed "
            "within 10 Å of each residue - a mutation-hotspot heatmap, "
            "using a sequential red palette. Single-fragment proteins only.",
        ).pack(side="left", padx=(4, 0))

        # Mutation heatmap colors (low/high swatches)
        mut_colors_frame = ctk.CTkFrame(self._steps_outer, fg_color="transparent")
        mut_colors_frame.grid(row=5, column=0, columnspan=3, padx=(48, 12), pady=2, sticky="w")
        ctk.CTkLabel(mut_colors_frame, text="Colors:").pack(side="left", padx=(0, 6))
        color_swatch_button(mut_colors_frame, self._ca_mutation_low_var).pack(side="left")
        ctk.CTkLabel(mut_colors_frame, text="→").pack(side="left", padx=4)
        color_swatch_button(mut_colors_frame, self._ca_mutation_high_var).pack(side="left")
        ctk.CTkButton(
            mut_colors_frame, text="↺", width=24, height=22,
            fg_color="gray30", hover_color="gray40",
            command=lambda: _reset_colors(
                (self._ca_mutation_low_var, MUTATION_DEFAULT_LOW_COLOR),
                (self._ca_mutation_high_var, MUTATION_DEFAULT_HIGH_COLOR),
            ),
        ).pack(side="left", padx=(6, 0))
        help_icon(
            mut_colors_frame,
            "Low/high ends of the mutation heatmap's color scale (few nearby "
            "patients to many). Click a swatch to choose a color; ↺ resets "
            "both back to the default red scale shown above.",
        ).pack(side="left", padx=(4, 0))

        # Log-scale (mutation heatmap only)
        if not hasattr(self, "_ca_log_scale_var"):
            self._ca_log_scale_var = ctk.BooleanVar(value=False)
        log_scale_frame = ctk.CTkFrame(self._steps_outer, fg_color="transparent")
        log_scale_frame.grid(row=6, column=0, columnspan=3, padx=(48, 12), pady=2, sticky="w")
        self._ca_log_scale_checkbox = ctk.CTkCheckBox(
            log_scale_frame, text="Log-scale",
            variable=self._ca_log_scale_var,
            checkbox_width=18, checkbox_height=18,
        )
        self._ca_log_scale_checkbox.pack(side="left")
        help_icon(
            log_scale_frame,
            "Only applies to the mutation heatmap above. COSMIC patient "
            "counts are often heavily right-skewed - a few hotspot residues "
            "can have counts far higher than the rest, so a linear color "
            "scale crushes nearly the whole protein into one flat color. "
            "Log-scaling compresses that range so variation among "
            "lower-count residues becomes visible too.",
        ).pack(side="left", padx=(4, 0))

        # Dim low-confidence residues (mutation heatmap only)
        if not hasattr(self, "_ca_dim_confidence_var"):
            self._ca_dim_confidence_var = ctk.BooleanVar(value=False)
        dim_frame = ctk.CTkFrame(self._steps_outer, fg_color="transparent")
        dim_frame.grid(row=7, column=0, columnspan=3, padx=(48, 12), pady=(2, 8), sticky="w")
        self._ca_dim_confidence_checkbox = ctk.CTkCheckBox(
            dim_frame, text="Dim low-confidence residues",
            variable=self._ca_dim_confidence_var,
            checkbox_width=18, checkbox_height=18,
        )
        self._ca_dim_confidence_checkbox.pack(side="left")
        help_icon(
            dim_frame,
            "Only applies to the mutation heatmap above. Fades each "
            "residue's heatmap color in proportion to how low its AlphaFold "
            "confidence (pLDDT) is, so a mutation hotspot in a poorly-"
            "modeled region reads as less trustworthy than an equally hot "
            "one in a well-modeled region, instead of both looking equally "
            "certain. Switches the script's lighting from \"soft\" to "
            "\"simple\", since ChimeraX's soft ambient shadows don't render "
            "correctly once part of the model is transparent.",
        ).pack(side="left", padx=(4, 0))
        self._on_ca_mutation_heatmap_toggle()  # sync initial enabled/disabled state

        # pLDDT heatmap
        if not hasattr(self, "_ca_plddt_heatmap_var"):
            self._ca_plddt_heatmap_var = ctk.BooleanVar(value=False)
        plddt_frame = ctk.CTkFrame(self._steps_outer, fg_color="transparent")
        plddt_frame.grid(row=8, column=0, columnspan=3, padx=24, pady=2, sticky="w")
        ctk.CTkCheckBox(
            plddt_frame, text="pLDDT heatmap",
            variable=self._ca_plddt_heatmap_var,
            checkbox_width=18, checkbox_height=18,
        ).pack(side="left")
        help_icon(
            plddt_frame,
            "Colors the ChimeraX cartoon by AlphaFold's per-residue "
            "confidence score (pLDDT), using ChimeraX's built-in AlphaFold "
            "color scheme (blue = very high confidence, down to orange = "
            "very low). Useful for judging how reliable the modeled "
            "structure is in a region, independent of mutation data. "
            "Single-fragment proteins only.",
        ).pack(side="left", padx=(4, 0))

        # pLDDT heatmap colors (low/high swatches)
        plddt_colors_frame = ctk.CTkFrame(self._steps_outer, fg_color="transparent")
        plddt_colors_frame.grid(row=9, column=0, columnspan=3, padx=(48, 12), pady=(2, 8), sticky="w")
        ctk.CTkLabel(plddt_colors_frame, text="Colors:").pack(side="left", padx=(0, 6))
        color_swatch_button(plddt_colors_frame, self._ca_plddt_low_var).pack(side="left")
        ctk.CTkLabel(plddt_colors_frame, text="→").pack(side="left", padx=4)
        color_swatch_button(plddt_colors_frame, self._ca_plddt_high_var).pack(side="left")
        ctk.CTkButton(
            plddt_colors_frame, text="↺", width=24, height=22,
            fg_color="gray30", hover_color="gray40",
            command=lambda: _reset_colors(
                (self._ca_plddt_low_var, PLDDT_DEFAULT_LOW_COLOR),
                (self._ca_plddt_high_var, PLDDT_DEFAULT_HIGH_COLOR),
            ),
        ).pack(side="left", padx=(6, 0))
        help_icon(
            plddt_colors_frame,
            "Low/high ends of the pLDDT heatmap's color scale (low "
            "confidence to high). Click a swatch to choose a color; ↺ "
            "resets both back to AlphaFold's own color scheme shown above.",
        ).pack(side="left", padx=(4, 0))

        # Mark PTM sites (independent marker, not a heatmap)
        if not hasattr(self, "_ca_mark_ptm_var"):
            self._ca_mark_ptm_var = ctk.BooleanVar(value=False)
        ptm_frame = ctk.CTkFrame(self._steps_outer, fg_color="transparent")
        ptm_frame.grid(row=10, column=0, columnspan=3, padx=24, pady=2, sticky="w")
        ctk.CTkCheckBox(
            ptm_frame, text="Mark PTM sites",
            variable=self._ca_mark_ptm_var,
            checkbox_width=18, checkbox_height=18,
        ).pack(side="left")
        help_icon(
            ptm_frame,
            "Marks each known PTM site (from the pipeline's PTM/mutation "
            "data) with a small sphere at its CA coordinate, layered on "
            "top of whichever heatmap(s) above are selected (or a plain "
            "cartoon if neither is). This is a separate marker, not a "
            "recoloring, so it never overwrites a heatmap's own color at "
            "that residue. Single-fragment proteins only.",
        ).pack(side="left", padx=(4, 0))
        color_swatch_button(ptm_frame, self._ca_ptm_marker_color_var).pack(side="left", padx=(8, 0))
        ctk.CTkButton(
            ptm_frame, text="↺", width=24, height=22,
            fg_color="gray30", hover_color="gray40",
            command=lambda: _reset_colors((self._ca_ptm_marker_color_var, PTM_MARKER_DEFAULT_COLOR)),
        ).pack(side="left", padx=(4, 0))

        # Mark mutations (independent marker, not a heatmap)
        if not hasattr(self, "_ca_mark_mutations_var"):
            self._ca_mark_mutations_var = ctk.BooleanVar(value=False)
        mut_marker_frame = ctk.CTkFrame(self._steps_outer, fg_color="transparent")
        mut_marker_frame.grid(row=11, column=0, columnspan=3, padx=24, pady=(2, 8), sticky="w")
        ctk.CTkCheckBox(
            mut_marker_frame, text="Show mutation markers",
            variable=self._ca_mark_mutations_var,
            checkbox_width=18, checkbox_height=18,
        ).pack(side="left")
        help_icon(
            mut_marker_frame,
            "Shows each COSMIC mutation position's side chain as a colored "
            "stick, layered on top of whichever heatmap(s) above are "
            "selected (or a plain cartoon if neither is). Like Mark PTM "
            "sites, this reveals real side-chain atoms rather than "
            "recoloring the cartoon, so it never overwrites a heatmap's own "
            "color at that residue. Single-fragment proteins only.",
        ).pack(side="left", padx=(4, 0))
        color_swatch_button(mut_marker_frame, self._ca_mutation_marker_color_var).pack(side="left", padx=(8, 0))
        ctk.CTkButton(
            mut_marker_frame, text="↺", width=24, height=22,
            fg_color="gray30", hover_color="gray40",
            command=lambda: _reset_colors(
                (self._ca_mutation_marker_color_var, MUTATION_MARKER_DEFAULT_COLOR),
            ),
        ).pack(side="left", padx=(4, 0))

        # Status label + progress bar (reuse the step status pattern)
        status = ctk.CTkLabel(
            self._steps_outer, text="●  Ready", width=100,
            anchor="e", text_color=_GRAY,
        )
        status.grid(row=12, column=1, columnspan=2, padx=12, pady=6, sticky="e")
        self._step_status_labels.append(status)

        bar = ctk.CTkProgressBar(self._steps_outer, width=120, height=14)
        bar.set(0)
        bar.grid(row=12, column=0, padx=12, pady=6, sticky="w")
        bar.grid_remove()
        self._step_progress_bars.append(bar)

    def _on_ca_source_change(self, _value: str = "") -> None:
        if self._ca_source_var.get() == "Upload CIF file":
            self._ca_database_frame.grid_remove()
            self._ca_custom_cif_frame.grid()
        else:
            self._ca_custom_cif_frame.grid_remove()
            self._ca_database_frame.grid()

    def _browse_ca_custom_cif(self) -> None:
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
        self._ca_custom_cif_var.set(path)
        uid = extract_uniprot_from_cif(Path(path))
        if uid:
            self._ca_custom_cif_uniprot_var.set(uid)

    def _on_ca_mutation_heatmap_toggle(self) -> None:
        """Keep the log-scale and dim-confidence checkboxes in sync with the
        mutation-heatmap toggle they modify -- disabled (and visually
        grayed) when there's no mutation heatmap left for them to affect.
        """
        state = "normal" if self._ca_mutation_heatmap_var.get() else "disabled"
        self._ca_log_scale_checkbox.configure(state=state)
        self._ca_dim_confidence_checkbox.configure(state=state)

    def _add_ca_protein(self) -> None:
        """Add a protein token (gene symbol or UniProt accession) to the
        batch list. Unlike Radius Sweep's gene picker, this doesn't validate
        against local pipeline data or resolve it up front -- Structure
        Heatmaps works from raw COSMIC + a live UniProt lookup, not the pipeline's own
        intermediate TSVs, so there's nothing to check locally, and a live
        network call on every "Add" click would make the button feel slow.
        Any token that fails to resolve is instead reported per-protein when
        the batch actually runs.
        """
        token = self._ca_protein_input_var.get().strip().upper()
        if not token:
            return

        if token in self._ca_proteins:
            self._ca_protein_error_label.configure(
                text=f"⚠  {token} is already in the list.", text_color=_YELLOW,
            )
            self._ca_protein_error_label.grid()
            self._ca_protein_input_var.set("")
            return

        self._ca_protein_error_label.grid_remove()
        self._ca_proteins.append(token)
        self._refresh_ca_protein_chips()
        self._ca_protein_input_var.set("")

    def _remove_ca_protein(self, token: str) -> None:
        if token in self._ca_proteins:
            self._ca_proteins.remove(token)
            self._refresh_ca_protein_chips()

    def _refresh_ca_protein_chips(self) -> None:
        """Redraw the added-proteins chip list, each removable via its own ✕ button."""
        for w in self._ca_proteins_list_frame.winfo_children():
            w.destroy()

        if not self._ca_proteins:
            ctk.CTkLabel(
                self._ca_proteins_list_frame, text="None added yet.",
                text_color=_GRAY, font=ctk.CTkFont(size=11),
            ).pack(anchor="w")
            return

        row = None
        for i, token in enumerate(self._ca_proteins):
            if i % 6 == 0:
                row = ctk.CTkFrame(self._ca_proteins_list_frame, fg_color="transparent")
                row.pack(anchor="w", pady=(0, 4))
            chip = ctk.CTkFrame(row, fg_color="#3a3a3a", corner_radius=6)
            chip.pack(side="left", padx=(0, 6))
            ctk.CTkLabel(chip, text=token, font=ctk.CTkFont(size=12)).pack(
                side="left", padx=(8, 4), pady=4,
            )
            ctk.CTkButton(
                chip, text="✕", width=20, height=20, fg_color="transparent",
                hover_color="#4a4a4a", font=ctk.CTkFont(size=11),
                command=lambda t=token: self._remove_ca_protein(t),
            ).pack(side="left", padx=(0, 6), pady=4)

    # ── File-status bar ──────────────────────────────────────────────────────

    def _refresh_file_status(self):
        """Update the status indicator for each input folder."""
        for name, (folder, exts, _desc, _validator) in _INPUT_FOLDERS.items():
            lbl = self._file_indicators[name]
            try:
                f = resolve_input_file(folder, exts)
                lbl.configure(text=f"✓  {name}: {f.name}", text_color=_GREEN)
            except FileNotFoundError:
                lbl.configure(text=f"✗  {name}: no file", text_color=_RED)
            except RuntimeError:
                lbl.configure(text=f"⚠  {name}: multiple files", text_color=_YELLOW)

    def _browse_file(self, name: str, folder: Path, filetypes: list, validator) -> None:
        """Open a file dialog, validate the selected file's content, then swap it
        into the input folder and refresh status.

        Copies to a hidden staging name in `folder` first and only clears the
        existing file(s) once the new file is confirmed copied and valid --
        never deletes-then-copies, which would empty the folder on a failed
        copy or if the source file already lives inside `folder`.
        """
        path = filedialog.askopenfilename(
            title=f"Select {name} input file",
            filetypes=filetypes + [("All files", "*.*")],
        )
        if not path:
            return

        src = Path(path).resolve()
        folder = folder.resolve()
        folder.mkdir(parents=True, exist_ok=True)

        staged = folder / f".browsing_{src.name}"
        try:
            shutil.copy2(src, staged)
        except OSError as exc:
            from tkinter import messagebox
            messagebox.showerror("Copy Failed", f"Could not copy {src.name}:\n\n{exc}")
            staged.unlink(missing_ok=True)
            return

        problems = [p.replace(staged.name, src.name) for p in validator(staged)]
        if problems:
            staged.unlink(missing_ok=True)
            from tkinter import messagebox
            messagebox.showerror(
                "Invalid File",
                f"{src.name} doesn't look like a valid {name} file — it was NOT "
                f"copied in (your existing file, if any, is untouched):\n\n"
                + "\n".join(f"  • {p}" for p in problems),
            )
            return

        for existing in folder.iterdir():
            if existing.is_file() and existing != staged:
                existing.unlink()
        staged.replace(folder / src.name)

        self._refresh_file_status()

    def _toggle_log(self):
        """Show or hide the raw log output panel (a normal part of the scrollable content)."""
        if self._log_visible:
            self._log_frame.grid_remove()
            self._log_toggle.configure(text="Show Details")
            self._log_visible = False
        else:
            self._log_frame.grid(row=10, column=0, padx=24, pady=(6, 12), sticky="ew")
            self._log_toggle.configure(text="Hide Details")
            self._log_visible = True

    def _browse_output_dir(self):
        """Let the user pick a custom output folder."""
        path = filedialog.askdirectory(
            title="Select output folder",
            initialdir=self._output_dir_var.get(),
        )
        if path:
            self._output_dir_var.set(path)
