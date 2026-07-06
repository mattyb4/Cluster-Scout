"""Pipeline tab widget building: mode panels, settings, browse dialogs.

The execution engine (run/stop/resume, subprocess streaming, the progress
queue) lives in ui/pipeline_runner.py — this file only builds widgets and
handles simple input events (browsing for a file, toggling the log).
"""
from __future__ import annotations

from pathlib import Path
from tkinter import filedialog
import shutil

import customtkinter as ctk
from PIL import Image

from ui.common import (
    PROJECT_ROOT, OUTPUT_DIR, _INPUT_FOLDERS,
    _GRAY, _RED, _GREEN, _YELLOW,
    PTM_PROXIMITY_STEPS, MUTATION_CLUSTERING_STEPS,
    resolve_input_file, extract_uniprot_from_cif,
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

        for i, (name, (folder, exts, desc)) in enumerate(_INPUT_FOLDERS.items(), 1):
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
                command=lambda n=name, f=folder, ft=filetypes: self._browse_file(n, f, ft),
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
            ("Analysis Tools", "analysis-tools"),
            ("CA Coordinates", "ca-coordinates"),
        ]:
            ctk.CTkRadioButton(
                mode_frame,
                text=label,
                variable=self._mode,
                value=value,
                command=self._rebuild_step_rows,
            ).pack(side="left", padx=8, pady=10)

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

        # Pipeline settings
        settings_frame = ctk.CTkFrame(p)
        settings_frame.grid(row=4, column=0, padx=24, pady=4, sticky="ew")

        ctk.CTkLabel(
            settings_frame, text="Settings:", font=ctk.CTkFont(weight="bold"),
        ).pack(side="left", padx=(12, 8), pady=8)

        ctk.CTkLabel(settings_frame, text="Cutoff (Å):").pack(side="left", padx=(8, 4), pady=8)
        self._cutoff_var = ctk.StringVar(value="10.0")
        ctk.CTkEntry(
            settings_frame, textvariable=self._cutoff_var, width=60,
        ).pack(side="left", padx=(0, 16), pady=8)

        ctk.CTkLabel(settings_frame, text="Min samples:").pack(side="left", padx=(8, 4), pady=8)
        self._min_samples_var = ctk.StringVar(value="3")
        ctk.CTkEntry(
            settings_frame, textvariable=self._min_samples_var, width=60,
        ).pack(side="left", padx=(0, 16), pady=8)

        ctk.CTkLabel(settings_frame, text="Min pLDDT:").pack(side="left", padx=(8, 4), pady=8)
        self._min_plddt_var = ctk.StringVar(value="")
        ctk.CTkEntry(
            settings_frame, textvariable=self._min_plddt_var, width=60,
            placeholder_text="off",
        ).pack(side="left", padx=(0, 16), pady=8)

        ctk.CTkLabel(settings_frame, text="Max PAE:").pack(side="left", padx=(8, 4), pady=8)
        self._max_pae_var = ctk.StringVar(value="")
        ctk.CTkEntry(
            settings_frame, textvariable=self._max_pae_var, width=60,
            placeholder_text="off",
        ).pack(side="left", padx=(0, 12), pady=8)

        # PolyPhen filter
        pp_frame = ctk.CTkFrame(p)
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

        self._log = ctk.CTkTextbox(
            p,
            height=140,
            font=ctk.CTkFont(family="Courier New", size=12),
            wrap="word",
            state="disabled",
        )
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

        if mode == "single-protein":
            self._build_single_protein_panel()
            return
        if mode == "analysis-tools":
            self._build_analysis_tools_panel()
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

    def _build_analysis_tools_panel(self):
        """Build the combined Radius Sweep / CIF Variance analysis panel.

        Both tools share one mode and one segmented-button toggle (also used by
        the Analysis Tools results tab) since they're both auxiliary structural
        analyses rather than pipeline steps — only the parameter fields below the
        toggle, and which script actually runs, differ per sub-tool.
        """
        ctk.CTkLabel(
            self._steps_outer,
            text="Analysis Tools",
            font=ctk.CTkFont(weight="bold"),
        ).grid(row=0, column=0, columnspan=4, padx=12, pady=(8, 2), sticky="w")

        if not hasattr(self, "_analysis_subtool_var"):
            self._analysis_subtool_var = ctk.StringVar(value="Radius Sweep")
        ctk.CTkSegmentedButton(
            self._steps_outer, values=["Radius Sweep", "CIF Variance"],
            variable=self._analysis_subtool_var,
            command=lambda _v: self._rebuild_step_rows(),
        ).grid(row=1, column=0, columnspan=4, padx=12, pady=(0, 8), sticky="w")

        if self._analysis_subtool_var.get() == "Radius Sweep":
            self._build_radius_sweep_fields()
        else:
            self._build_cif_variance_fields()

    def _build_radius_sweep_fields(self):
        """Radius Sweep parameter fields — rows 2-5 of the Analysis Tools panel."""
        # Genes
        ctk.CTkLabel(self._steps_outer, text="Genes:", anchor="w").grid(
            row=2, column=0, padx=(12, 6), pady=6, sticky="w"
        )
        if not hasattr(self, "_radius_genes_var"):
            from radius_sweep import DEFAULT_GENES
            self._radius_genes_var = ctk.StringVar(value=" ".join(DEFAULT_GENES))
        ctk.CTkEntry(
            self._steps_outer, textvariable=self._radius_genes_var, width=400,
            placeholder_text="Space-separated gene symbols",
        ).grid(row=2, column=1, columnspan=3, padx=6, pady=6, sticky="ew")

        # Radius range
        ctk.CTkLabel(self._steps_outer, text="Radius range (Å):", anchor="w").grid(
            row=3, column=0, padx=(12, 6), pady=6, sticky="w"
        )
        range_frame = ctk.CTkFrame(self._steps_outer, fg_color="transparent")
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

        # Unfiltered comparison
        if not hasattr(self, "_radius_unfiltered_var"):
            self._radius_unfiltered_var = ctk.BooleanVar(value=False)
        ctk.CTkCheckBox(
            self._steps_outer, text="Include unfiltered COSMIC comparison",
            variable=self._radius_unfiltered_var,
            checkbox_width=18, checkbox_height=18,
        ).grid(row=4, column=0, columnspan=4, padx=12, pady=6, sticky="w")

        # Status label + progress bar (reuse the step status pattern)
        status = ctk.CTkLabel(
            self._steps_outer, text="●  Ready", width=100,
            anchor="e", text_color=_GRAY,
        )
        status.grid(row=5, column=1, columnspan=3, padx=12, pady=6, sticky="e")
        self._step_status_labels.append(status)

        bar = ctk.CTkProgressBar(self._steps_outer, width=120, height=14)
        bar.set(0)
        bar.grid(row=5, column=0, padx=12, pady=6, sticky="w")
        bar.grid_remove()
        self._step_progress_bars.append(bar)

    def _build_cif_variance_fields(self):
        """CIF Variance parameter fields — rows 2-9 of the Analysis Tools panel."""
        # Input folder
        ctk.CTkLabel(self._steps_outer, text="Input folder:", anchor="w").grid(
            row=2, column=0, padx=(12, 6), pady=6, sticky="w"
        )
        if not hasattr(self, "_variance_input_dir_var"):
            from cif_variance import DEFAULT_INPUT_DIR
            self._variance_input_dir_var = ctk.StringVar(value=str(DEFAULT_INPUT_DIR))
            self._variance_input_dir_var.trace_add("write", self._update_variance_cif_count)
        ctk.CTkEntry(
            self._steps_outer, textvariable=self._variance_input_dir_var, width=340,
        ).grid(row=2, column=1, padx=6, pady=6, sticky="ew")
        ctk.CTkButton(
            self._steps_outer, text="Browse", width=70, height=26,
            font=ctk.CTkFont(size=12),
            command=self._browse_variance_input_dir,
        ).grid(row=2, column=2, padx=12, pady=6, sticky="e")

        self._variance_cif_count_label = ctk.CTkLabel(
            self._steps_outer, text="", anchor="w", font=ctk.CTkFont(size=11),
        )
        self._variance_cif_count_label.grid(row=3, column=1, columnspan=2, padx=6, pady=(0, 6), sticky="w")
        self._update_variance_cif_count()

        # Top N
        ctk.CTkLabel(self._steps_outer, text="Top N residues:", anchor="w").grid(
            row=4, column=0, padx=(12, 6), pady=6, sticky="w"
        )
        if not hasattr(self, "_variance_top_var"):
            self._variance_top_var = ctk.StringVar(value="10")
        ctk.CTkEntry(
            self._steps_outer, textvariable=self._variance_top_var, width=60,
        ).grid(row=4, column=1, padx=6, pady=6, sticky="w")

        # Report range
        ctk.CTkLabel(self._steps_outer, text="Report range:", anchor="w").grid(
            row=5, column=0, padx=(12, 6), pady=6, sticky="w"
        )
        report_frame = ctk.CTkFrame(self._steps_outer, fg_color="transparent")
        report_frame.grid(row=5, column=1, columnspan=2, padx=6, pady=6, sticky="w")
        if not hasattr(self, "_variance_range_start_var"):
            self._variance_range_start_var = ctk.StringVar(value="")
            self._variance_range_end_var = ctk.StringVar(value="")
        ctk.CTkEntry(report_frame, textvariable=self._variance_range_start_var, width=70,
                     placeholder_text="start").pack(side="left", padx=(0, 6))
        ctk.CTkEntry(report_frame, textvariable=self._variance_range_end_var, width=70,
                     placeholder_text="end (blank = all)").pack(side="left")

        # Align range
        ctk.CTkLabel(self._steps_outer, text="Align range:", anchor="w").grid(
            row=6, column=0, padx=(12, 6), pady=6, sticky="w"
        )
        align_frame = ctk.CTkFrame(self._steps_outer, fg_color="transparent")
        align_frame.grid(row=6, column=1, columnspan=2, padx=6, pady=6, sticky="w")
        if not hasattr(self, "_variance_align_start_var"):
            self._variance_align_start_var = ctk.StringVar(value="")
            self._variance_align_end_var = ctk.StringVar(value="")
        ctk.CTkEntry(align_frame, textvariable=self._variance_align_start_var, width=70,
                     placeholder_text="start").pack(side="left", padx=(0, 6))
        ctk.CTkEntry(align_frame, textvariable=self._variance_align_end_var, width=70,
                     placeholder_text="end (blank = same as report range)").pack(side="left")

        # UniProt / gene overrides
        ctk.CTkLabel(self._steps_outer, text="UniProt override:", anchor="w").grid(
            row=7, column=0, padx=(12, 6), pady=6, sticky="w"
        )
        if not hasattr(self, "_variance_uniprot_var"):
            self._variance_uniprot_var = ctk.StringVar(value="")
        ctk.CTkEntry(
            self._steps_outer, textvariable=self._variance_uniprot_var, width=150,
            placeholder_text="auto-detected from CIF",
        ).grid(row=7, column=1, padx=6, pady=6, sticky="w")

        ctk.CTkLabel(self._steps_outer, text="Gene override:", anchor="w").grid(
            row=8, column=0, padx=(12, 6), pady=6, sticky="w"
        )
        if not hasattr(self, "_variance_gene_var"):
            self._variance_gene_var = ctk.StringVar(value="")
        ctk.CTkEntry(
            self._steps_outer, textvariable=self._variance_gene_var, width=150,
            placeholder_text="optional, for UniProt lookup",
        ).grid(row=8, column=1, padx=6, pady=6, sticky="w")

        # Status label + progress bar (reuse the step status pattern)
        status = ctk.CTkLabel(
            self._steps_outer, text="●  Ready", width=100,
            anchor="e", text_color=_GRAY,
        )
        status.grid(row=9, column=1, columnspan=2, padx=12, pady=6, sticky="e")
        self._step_status_labels.append(status)

        bar = ctk.CTkProgressBar(self._steps_outer, width=120, height=14)
        bar.set(0)
        bar.grid(row=9, column=0, padx=12, pady=6, sticky="w")
        bar.grid_remove()
        self._step_progress_bars.append(bar)

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

    def _build_ca_coordinates_panel(self):
        """Build the input fields for CA-coordinate export mode."""
        ctk.CTkLabel(
            self._steps_outer,
            text="Export Alpha-Carbon Coordinates",
            font=ctk.CTkFont(weight="bold"),
        ).grid(row=0, column=0, columnspan=3, padx=12, pady=(8, 2), sticky="w")

        # UniProt ID
        ctk.CTkLabel(self._steps_outer, text="UniProt ID:", anchor="w").grid(
            row=1, column=0, padx=(12, 6), pady=6, sticky="w"
        )
        if not hasattr(self, "_ca_uniprot_var"):
            self._ca_uniprot_var = ctk.StringVar(value="")
        ctk.CTkEntry(
            self._steps_outer, textvariable=self._ca_uniprot_var, width=200,
            placeholder_text="e.g. P04637",
        ).grid(row=1, column=1, padx=6, pady=6, sticky="w")

        # Gene override
        ctk.CTkLabel(self._steps_outer, text="Gene override:", anchor="w").grid(
            row=2, column=0, padx=(12, 6), pady=6, sticky="w"
        )
        if not hasattr(self, "_ca_gene_var"):
            self._ca_gene_var = ctk.StringVar(value="")
        ctk.CTkEntry(
            self._steps_outer, textvariable=self._ca_gene_var, width=200,
            placeholder_text="optional, skips UniProt API lookup",
        ).grid(row=2, column=1, padx=6, pady=6, sticky="w")

        # Status label + progress bar (reuse the step status pattern)
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

    # ── File-status bar ──────────────────────────────────────────────────────

    def _refresh_file_status(self):
        """Update the status indicator for each input folder."""
        for name, (folder, exts, _desc) in _INPUT_FOLDERS.items():
            lbl = self._file_indicators[name]
            try:
                f = resolve_input_file(folder, exts)
                lbl.configure(text=f"✓  {name}: {f.name}", text_color=_GREEN)
            except FileNotFoundError:
                lbl.configure(text=f"✗  {name}: no file", text_color=_RED)
            except RuntimeError:
                lbl.configure(text=f"⚠  {name}: multiple files", text_color=_YELLOW)

    def _browse_file(self, name: str, folder: Path, filetypes: list) -> None:
        """Open a file dialog, copy the selected file into the input folder, and refresh status."""
        path = filedialog.askopenfilename(
            title=f"Select {name} input file",
            filetypes=filetypes + [("All files", "*.*")],
        )
        if not path:
            return

        src = Path(path)
        folder.mkdir(parents=True, exist_ok=True)

        for existing in folder.iterdir():
            if existing.is_file():
                existing.unlink()

        shutil.copy2(src, folder / src.name)
        self._refresh_file_status()

    def _toggle_log(self):
        """Show or hide the raw log output panel (a normal part of the scrollable content)."""
        if self._log_visible:
            self._log.grid_remove()
            self._log_toggle.configure(text="Show Details")
            self._log_visible = False
        else:
            self._log.grid(row=10, column=0, padx=24, pady=(6, 12), sticky="ew")
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
