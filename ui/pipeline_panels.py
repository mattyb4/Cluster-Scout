"""Pipeline tab widget building: mode panels, settings, browse dialogs.

The execution engine (run/stop/resume, subprocess streaming, the progress
queue) lives in ui/pipeline_runner.py — this file only builds widgets and
handles simple input events (browsing for a file, toggling the log).
"""
from __future__ import annotations

import shutil
import sys
from pathlib import Path
from tkinter import filedialog

import customtkinter as ctk

from ui import theme
from ui.common import (
    _BTN2,
    _BTN2_HOVER,
    _CLINVAR_CLASS_OPTIONS,
    _DANGER,
    _GRAY,
    _GREEN,
    _INPUT_FOLDERS,
    _MODE_HELP,
    _MUTATION_SOURCE_BY_LABEL,
    _MUTATION_SOURCE_OPTIONS,
    _PSP_INPUT_FOLDER,
    _PTM_SOURCE_BY_LABEL,
    _PTM_SOURCE_OPTIONS,
    _RED,
    _YELLOW,
    OUTPUT_DIR,
    add_resize_grip,
    discover_psp_files,
    extract_uniprot_from_cif,
    help_icon,
    isolate_textbox_scroll,
    load_logo_image,
    pipeline_steps,
    psp_folder_summary,
    resolve_input_file,
    section_header,
    validate_psp_folder,
)

_LABEL_W = 120  # shared label column width, so every control starts on the same line


class PipelineTabMixin:
    def _build_pipeline_tab(self, tab) -> None:
        tab.grid_columnconfigure(0, weight=1)
        tab.grid_rowconfigure(0, weight=1)  # scroll area; the run bar below it stays pinned

        # Explicit fg_color (not "transparent"): children of a transparent CTkScrollableFrame
        # capture a static bg color at construction and show stale corners after a light/dark switch.
        scroll = ctk.CTkScrollableFrame(tab, fg_color=theme.pair("bg"))
        scroll.grid(row=0, column=0, sticky="nsew")
        scroll.grid_columnconfigure(0, weight=1)

        # One content column spanning the window
        content = ctk.CTkFrame(scroll, fg_color="transparent", border_width=0)
        content.grid(row=0, column=0, sticky="ew", padx=24)
        content.grid_columnconfigure(0, weight=1)
        p = content  # all pipeline widgets go in the content column

        # Header: the logo
        header = ctk.CTkFrame(p, fg_color="transparent")
        header.grid(row=0, column=0, pady=(10, 0), sticky="ew")
        header.grid_columnconfigure(1, weight=1)
        logo = load_logo_image(104)
        if logo is not None:
            ctk.CTkLabel(header, image=logo, text="").grid(row=0, column=0, sticky="w")
        else:
            ctk.CTkLabel(
                header, text="Cluster-Scout", font=ctk.CTkFont(size=22, weight="bold"),
            ).grid(row=0, column=0, sticky="w")

        # Data-file status bar with Browse buttons
        # Each of the three main sections (input files, run settings, steps) is a bordered card
        file_card = ctk.CTkFrame(p)
        file_card.grid(row=1, column=0, pady=(14, 0), sticky="ew")
        file_card.grid_columnconfigure(0, weight=1)
        section_header(file_card, "Input files", row=0, pady=(12, 4), padx=14)
        self._file_frame = ctk.CTkFrame(file_card, fg_color="transparent")
        self._file_frame.grid(row=1, column=0, padx=14, pady=(0, 10), sticky="ew")
        self._file_frame.grid_columnconfigure(1, weight=1)
        self._file_indicators: dict[str, ctk.CTkLabel] = {}
        self._file_buttons: dict[str, ctk.CTkButton] = {}

        for i, (name, (folder, exts, desc, validator)) in enumerate(_INPUT_FOLDERS.items(), 1):
            lbl = ctk.CTkLabel(self._file_frame, text=f"{name} …", anchor="w")
            lbl.grid(row=i, column=0, columnspan=2, padx=(0, 6), pady=3, sticky="ew")
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
            btn.grid(row=i, column=2, padx=(12, 0), pady=3, sticky="e")
            self._file_buttons[name] = btn

        # PhosphoSitePlus is a folder of files, so it gets its own row and a
        # folder picker instead of the single-file Browse above
        psp_row = len(_INPUT_FOLDERS) + 1
        self._psp_file_indicator = ctk.CTkLabel(self._file_frame, text="PhosphoSitePlus …", anchor="w")
        self._psp_file_indicator.grid(row=psp_row, column=0, columnspan=2, padx=(0, 6), pady=3, sticky="ew")
        ctk.CTkButton(
            self._file_frame, text="Browse", width=70, height=26,
            font=ctk.CTkFont(size=12), command=self._browse_psp_folder,
        ).grid(row=psp_row, column=2, padx=(12, 0), pady=3, sticky="e")

        # Mode selection, with PTM Proximity's PTM-source selector underneath
        # (packed/unpacked in _rebuild_step_rows as the mode changes)
        settings_card = ctk.CTkFrame(p)
        settings_card.grid(row=2, column=0, pady=(14, 0), sticky="ew")
        settings_card.grid_columnconfigure(0, weight=1)
        section_header(settings_card, "Run settings", row=0, pady=(12, 4), padx=14)
        mode_area = ctk.CTkFrame(settings_card, fg_color="transparent")
        mode_area.grid(row=1, column=0, padx=14, pady=2, sticky="ew")
        ctk.CTkFrame(settings_card, fg_color="transparent", width=0, height=8).grid(row=6, column=0)  # bottom padding
        mode_frame = ctk.CTkFrame(mode_area, fg_color="transparent")
        mode_frame.pack(fill="x")

        ctk.CTkLabel(
            mode_frame, text="Mode", width=_LABEL_W, anchor="w",
        ).pack(side="left", padx=(0, 8), pady=8)

        self._mode = ctk.StringVar(value="ptm-proximity")
        for label, value in [
            ("PTM Proximity", "ptm-proximity"),
            ("Mutation Clustering", "mutation-clustering"),
            ("Single Protein", "single-protein"),
        ]:
            ctk.CTkRadioButton(
                mode_frame,
                text=label,
                variable=self._mode,
                value=value,
                command=self._rebuild_step_rows,
            ).pack(side="left", padx=(8, 0), pady=8)
            help_icon(mode_frame, _MODE_HELP[value]).pack(side="left", padx=(4, 8), pady=8)

        self._ptm_source_frame = ctk.CTkFrame(mode_area, fg_color="transparent")
        ctk.CTkLabel(
            self._ptm_source_frame, text="PTM source", width=_LABEL_W, anchor="w",
        ).pack(side="left", padx=(0, 8), pady=6)
        self._ptm_source_var = ctk.StringVar(value=_PTM_SOURCE_OPTIONS[0][0])
        ctk.CTkSegmentedButton(
            self._ptm_source_frame, values=[label for label, _ in _PTM_SOURCE_OPTIONS],
            variable=self._ptm_source_var, command=lambda _v: self._rebuild_step_rows(),
        ).pack(side="left", pady=8)
        help_icon(
            self._ptm_source_frame,
            "Where PTM Proximity gets its PTM sites. PTMD: PTMD 2.0's "
            "disease-associated PTMs. PhosphoSitePlus: every human site in a "
            "PhosphoSitePlus download, with its LTP/HTP evidence counts and any "
            "disease associations. Each source writes its own output files "
            "(PhosphoSitePlus's start with psp_), so running one never "
            "overwrites the other.",
        ).pack(side="left", padx=(6, 8), pady=8)

        # Mutation source (PTM Proximity and Mutation Clustering), and for ClinVar
        # the classifications to include -- both shown/hidden in _rebuild_step_rows
        self._mutation_source_frame = ctk.CTkFrame(mode_area, fg_color="transparent")
        ctk.CTkLabel(
            self._mutation_source_frame, text="Mutation source", width=_LABEL_W, anchor="w",
        ).pack(side="left", padx=(0, 8), pady=6)
        self._mutation_source_var = ctk.StringVar(value=_MUTATION_SOURCE_OPTIONS[0][0])
        ctk.CTkSegmentedButton(
            self._mutation_source_frame, values=[label for label, _ in _MUTATION_SOURCE_OPTIONS],
            variable=self._mutation_source_var, command=lambda _v: self._rebuild_step_rows(),
        ).pack(side="left", pady=8)
        help_icon(
            self._mutation_source_frame,
            "Where the mutations come from. COSMIC: recurrent somatic cancer "
            "mutations from the COSMIC Mutant Census. ClinVar: missense variants "
            "from ClinVar's variant_summary file (any disease), counted by the "
            "number of ClinVar submitters. Each source writes its own files "
            "(ClinVar's start with clinvar_), so running one never overwrites "
            "the other.",
        ).pack(side="left", padx=(6, 8), pady=8)

        self._clinvar_classes_frame = ctk.CTkFrame(mode_area, fg_color="transparent")
        ctk.CTkLabel(
            self._clinvar_classes_frame, text="ClinVar classes", width=_LABEL_W, anchor="w",
        ).pack(side="left", padx=(0, 8), pady=6)
        self._clinvar_class_vars: dict[str, ctk.BooleanVar] = {}
        for label, key in _CLINVAR_CLASS_OPTIONS:
            var = ctk.BooleanVar(value=key in ("P", "LP"))
            self._clinvar_class_vars[key] = var
            ctk.CTkCheckBox(
                self._clinvar_classes_frame, text=label, variable=var,
                checkbox_width=18, checkbox_height=18,
            ).pack(side="left", padx=(0, 10), pady=8)
        help_icon(
            self._clinvar_classes_frame,
            "Which ClinVar classifications count. Pathogenic + Likely pathogenic "
            "is about 4,400 genes; adding Uncertain significance makes the run "
            "far larger. A combined classification such as "
            "'Pathogenic/Likely pathogenic' counts for each part.",
        ).pack(side="left", padx=(0, 8), pady=8)

        # Output folder selector
        out_frame = ctk.CTkFrame(settings_card, fg_color="transparent")
        out_frame.grid(row=2, column=0, padx=14, pady=2, sticky="ew")
        out_frame.grid_columnconfigure(1, weight=1)

        ctk.CTkLabel(
            out_frame, text="Output folder", width=_LABEL_W, anchor="w",
        ).grid(row=0, column=0, padx=(0, 8), pady=6, sticky="w")

        # Shared by every mode (they intentionally write into the same output DB files)
        self._output_dir_var = ctk.StringVar(value=str(OUTPUT_DIR))
        self._output_dir_entry = ctk.CTkEntry(
            out_frame, textvariable=self._output_dir_var, state="readonly",
        )
        self._output_dir_entry.grid(row=0, column=1, padx=6, pady=6, sticky="ew")

        ctk.CTkButton(
            out_frame, text="Change", width=70, height=26,
            font=ctk.CTkFont(size=12),
            command=self._browse_output_dir,
        ).grid(row=0, column=2, padx=(6, 4), pady=8, sticky="e")

        ctk.CTkButton(
            out_frame, text="Reset", width=60, height=26,
            font=ctk.CTkFont(size=12), fg_color=_BTN2, hover_color=_BTN2_HOVER,
            command=lambda: self._output_dir_var.set(str(OUTPUT_DIR)),
        ).grid(row=0, column=3, padx=(0, 0), pady=6, sticky="e")

        # Pipeline settings (ptm-proximity, mutation-clustering, AND single-protein
        # -- analyze_single_cif_nearby_mutations.py accepts all four of these too)
        settings_frame = self._settings_frame = ctk.CTkFrame(settings_card, fg_color="transparent")
        settings_frame.grid(row=3, column=0, padx=14, pady=2, sticky="ew")

        ctk.CTkLabel(
            settings_frame, text="Thresholds", width=_LABEL_W, anchor="w",
        ).pack(side="left", padx=(0, 8), pady=6)

        # CustomTkinter's bundled Roboto reports an ascent too small for the
        # ring of the Å, so Windows clips it; Segoe UI's metrics include it
        cutoff_font = ctk.CTkFont(family="Segoe UI", size=13) if sys.platform == "win32" else None
        ctk.CTkLabel(settings_frame, text="Cutoff (Å):", font=cutoff_font).pack(side="left", padx=(0, 4), pady=6)
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
            "Minimum number of distinct COSMIC patient samples (or ClinVar "
            "submitters, with the ClinVar source) a mutation "
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
        pp_frame = self._pp_frame = ctk.CTkFrame(settings_card, fg_color="transparent")
        pp_frame.grid(row=4, column=0, padx=14, pady=2, sticky="ew")

        ctk.CTkLabel(
            pp_frame, text="PolyPhen filter", width=_LABEL_W, anchor="w",
        ).pack(side="left", padx=(0, 4), pady=6)

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

        # PhosphoSitePlus site filters -- only shown for PTM Proximity with the
        # PhosphoSitePlus source (see _rebuild_step_rows); applied in step 1
        psp_frame = self._psp_filter_frame = ctk.CTkFrame(settings_card, fg_color="transparent")
        psp_frame.grid(row=5, column=0, padx=14, pady=2, sticky="ew")

        ctk.CTkLabel(
            psp_frame, text="PSP site filters", width=_LABEL_W, anchor="w",
        ).pack(side="left", padx=(0, 4), pady=6)

        for kind, label, help_text in (
            ("ltp", "LTP", "Keep only sites whose LTP count - literature records from "
                           "low-throughput experiments aimed at the site itself - falls in "
                           "this range. Leave the maximum blank for no upper limit."),
            ("htp", "HTP", "Keep only sites whose HTP count - high-throughput mass-spec "
                           "papers that detected the site - falls in this range. Leave the "
                           "maximum blank for no upper limit."),
        ):
            ctk.CTkLabel(psp_frame, text=f"{label}:").pack(side="left", padx=(12, 4), pady=8)
            help_icon(psp_frame, help_text).pack(side="left", padx=(0, 4), pady=8)
            min_var = ctk.StringVar(value="0")
            max_var = ctk.StringVar(value="")
            setattr(self, f"_psp_min_{kind}_var", min_var)
            setattr(self, f"_psp_max_{kind}_var", max_var)
            ctk.CTkEntry(psp_frame, textvariable=min_var, width=50).pack(side="left", pady=8)
            ctk.CTkLabel(psp_frame, text="to").pack(side="left", padx=4, pady=8)
            ctk.CTkEntry(
                psp_frame, textvariable=max_var, width=70, placeholder_text="no limit",
            ).pack(side="left", padx=(0, 8), pady=8)

        self._psp_disease_only_var = ctk.BooleanVar(value=False)
        ctk.CTkCheckBox(
            psp_frame, text="Disease-associated sites only",
            variable=self._psp_disease_only_var,
            checkbox_width=18, checkbox_height=18,
        ).pack(side="left", padx=(12, 4), pady=8)
        help_icon(
            psp_frame,
            "Keep only sites PhosphoSitePlus lists a disease association for "
            "(needs the Disease-associated_sites file in the input folder).",
        ).pack(side="left", padx=(0, 12), pady=8)

        # Steps panel
        steps_card = ctk.CTkFrame(p)
        steps_card.grid(row=3, column=0, pady=(14, 0), sticky="ew")
        steps_card.grid_columnconfigure(0, weight=1)
        self._steps_outer = ctk.CTkFrame(steps_card, fg_color="transparent")
        self._steps_outer.grid(row=0, column=0, padx=14, pady=(0, 10), sticky="ew")
        self._steps_outer.grid_columnconfigure(1, weight=1)
        self._rebuild_step_rows()

        # Run bar, pinned under the scroll area so Run/Stop are always visible
        ctk.CTkFrame(
            tab, width=0, height=1, corner_radius=0, border_width=0, fg_color=theme.pair("line"),
        ).grid(row=1, column=0, sticky="ew")
        btn_frame = ctk.CTkFrame(tab, fg_color=theme.pair("surface"), corner_radius=0, border_width=0)
        btn_frame.grid(row=2, column=0, sticky="ew")

        self._run_btn = ctk.CTkButton(
            btn_frame,
            text="▶  Run Pipeline",
            command=self._start_pipeline,
            width=150,
            height=36,
            font=ctk.CTkFont(size=14, weight="bold"),
        )
        self._run_btn.pack(side="left", padx=(24, 8), pady=10)

        self._stop_btn = ctk.CTkButton(
            btn_frame,
            text="■  Stop",
            command=self._stop_pipeline,
            width=110,
            height=36,
            font=ctk.CTkFont(size=14, weight="bold"),
            state="disabled",
            fg_color=_BTN2,
            hover_color=_DANGER,
        )
        self._stop_btn.pack(side="left", padx=(0, 8), pady=10)

        self._open_btn = ctk.CTkButton(
            btn_frame,
            text="Open Output Folder",
            command=self._open_output_folder,
            width=160,
            height=36,
        )
        self._open_btn.pack(side="left", pady=10)

        ctk.CTkButton(
            btn_frame,
            text="Manage Cache",
            command=self._manage_cache,
            width=140,
            height=36,
            fg_color=_BTN2,
            hover_color=_BTN2_HOVER,
        ).pack(side="left", padx=(8, 0), pady=10)

        self._timer_label = ctk.CTkLabel(
            btn_frame,
            text="",
            text_color=_GRAY,
            font=ctk.CTkFont(size=13),
        )
        self._timer_label.pack(side="right", padx=(12, 24))

        import tkinter as tk
        self._activity_width = 40
        self._activity_height = 6
        self._activity_chunk = 12
        self._activity_canvas = tk.Canvas(
            btn_frame, width=self._activity_width, height=self._activity_height,
            highlightthickness=0, bg=theme.c("line"), bd=0,
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
            fg_color=_BTN2,
            hover_color=_BTN2_HOVER,
            command=self._toggle_log,
        )
        self._log_toggle.grid(row=4, column=0, pady=(14, 0), sticky="w")

        self._log_frame = ctk.CTkFrame(p, fg_color="transparent")
        self._log = ctk.CTkTextbox(
            self._log_frame,
            height=140,
            font=ctk.CTkFont(family="Courier New", size=12),
            wrap="word",
            state="disabled",
            border_width=1,
            border_color=theme.pair("line"),
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

    def _ptm_source(self) -> str:
        """The selected PTM source as a --ptm-source value ("ptmd" or "psp")."""
        return _PTM_SOURCE_BY_LABEL[self._ptm_source_var.get()]

    def _mutation_source(self) -> str:
        """The selected mutation source as a --mutation-source value ("cosmic" or "clinvar")."""
        return _MUTATION_SOURCE_BY_LABEL[self._mutation_source_var.get()]

    def _clinvar_classes(self) -> list[str]:
        """The ticked ClinVar classes as --clinvar-classes values."""
        return [key for _label, key in _CLINVAR_CLASS_OPTIONS if self._clinvar_class_vars[key].get()]

    def _psp_filter_args(self) -> tuple[list[str], list[str]]:
        """Step-1 arguments for the PhosphoSitePlus filters, plus any problems
        with the entered values (blank minimum = 0, blank maximum = no limit)."""
        args: list[str] = []
        problems: list[str] = []
        for kind, label in (("ltp", "LTP"), ("htp", "HTP")):
            lo_text = getattr(self, f"_psp_min_{kind}_var").get().strip() or "0"
            hi_text = getattr(self, f"_psp_max_{kind}_var").get().strip()
            try:
                lo = int(lo_text)
                hi = int(hi_text) if hi_text else None
            except ValueError:
                problems.append(f"{label} range must be whole numbers")
                continue
            if lo < 0 or (hi is not None and hi < 0):
                problems.append(f"{label} range can't be negative")
            elif hi is not None and hi < lo:
                problems.append(f"{label} maximum ({hi}) is below its minimum ({lo})")
            else:
                if lo:
                    args += [f"--psp-min-{kind}", str(lo)]
                if hi is not None:
                    args += [f"--psp-max-{kind}", str(hi)]
        if self._psp_disease_only_var.get():
            args.append("--psp-disease-only")
        return args, problems

    def _rebuild_step_rows(self):
        for w in self._steps_outer.winfo_children():
            w.destroy()
        self._step_status_labels = []
        self._step_progress_bars: list[ctk.CTkProgressBar] = []
        self._steps_outer.grid_columnconfigure(1, weight=1)

        mode = self._mode.get()
        # Re-pack in a fixed order, since pack() appends after whatever is showing
        for frame in (self._ptm_source_frame, self._mutation_source_frame, self._clinvar_classes_frame):
            frame.pack_forget()
        if mode == "ptm-proximity":
            self._ptm_source_frame.pack(fill="x", pady=(4, 0))
        if mode in ("ptm-proximity", "mutation-clustering"):
            self._mutation_source_frame.pack(fill="x", pady=(4, 0))
            if self._mutation_source() == "clinvar":
                self._clinvar_classes_frame.pack(fill="x", pady=(0, 0))
        if mode == "ptm-proximity" and self._ptm_source() == "psp":
            self._psp_filter_frame.grid()
        else:
            self._psp_filter_frame.grid_remove()

        # Settings frame also applies to single-protein (accepted by
        # analyze_single_cif_nearby_mutations.py); PolyPhen filter feeds step 4 of
        # both ptm-proximity and mutation-clustering.
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

        steps = pipeline_steps(mode, self._ptm_source(), self._mutation_source())

        section_header(self._steps_outer, "Steps", row=0, columnspan=4)

        for i, (panel_label, _log_label) in enumerate(steps, 1):
            ctk.CTkLabel(self._steps_outer, text=f"{i}.", width=28, anchor="w").grid(
                row=i, column=0, padx=(0, 0), pady=5, sticky="w"
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
            status.grid(row=i, column=3, padx=(12, 0), pady=5, sticky="e")
            self._step_status_labels.append(status)

    def _build_single_protein_panel(self):
        """Build the input fields for single-protein analysis mode."""
        section_header(self._steps_outer, "Single protein analysis", row=0, columnspan=3)

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

        psp_lbl = self._psp_file_indicator
        if not discover_psp_files(_PSP_INPUT_FOLDER).site:
            psp_lbl.configure(text="✗  PhosphoSitePlus: no site dataset files", text_color=_RED)
        elif validate_psp_folder(_PSP_INPUT_FOLDER):
            psp_lbl.configure(text="⚠  PhosphoSitePlus: problem with files (see Run check)",
                              text_color=_YELLOW)
        else:
            psp_lbl.configure(text=f"✓  PhosphoSitePlus: {psp_folder_summary(_PSP_INPUT_FOLDER)}",
                              text_color=_GREEN)

    def _browse_psp_folder(self) -> None:
        """Pick a PhosphoSitePlus download folder, validate it, then copy the
        files the pipeline uses (site datasets + disease file, recognized by
        their columns) into data/input/phosphositeplus/, replacing what's there.

        Copies into a staging folder first and only swaps once every copy
        succeeded, so a failed copy never leaves the input folder half-empty.
        """
        from tkinter import messagebox

        path = filedialog.askdirectory(title="Select a PhosphoSitePlus download folder")
        if not path:
            return
        src = Path(path).resolve()
        dest = _PSP_INPUT_FOLDER.resolve()
        if src == dest:
            self._refresh_file_status()
            return

        problems = validate_psp_folder(src)
        if problems:
            messagebox.showerror(
                "Invalid Folder",
                f"{src.name} can't be used as PhosphoSitePlus input — nothing was "
                f"copied (your existing files, if any, are untouched):\n\n"
                + "\n".join(f"  • {p}" for p in problems),
            )
            return

        found = discover_psp_files(src)
        staging = dest.parent / f".browsing_{dest.name}"
        shutil.rmtree(staging, ignore_errors=True)
        try:
            staging.mkdir(parents=True)
            for f in found.site + found.disease:
                shutil.copy2(f, staging / f.name)
        except OSError as exc:
            shutil.rmtree(staging, ignore_errors=True)
            messagebox.showerror("Copy Failed", f"Could not copy PhosphoSitePlus files:\n\n{exc}")
            return

        shutil.rmtree(dest, ignore_errors=True)
        staging.replace(dest)
        self._refresh_file_status()

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
            self._log_frame.grid(row=5, column=0, pady=(6, 16), sticky="ew")
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
