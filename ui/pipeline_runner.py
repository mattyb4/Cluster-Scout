"""Pipeline execution engine: run/stop/resume/cancel, subprocess streaming,
the progress queue, output backups, and cache management.

The `"viz_data"` branch in `_poll_queue` is how a background worker (radius
sweep / CIF variance) hands its result to AnalysisToolsTabMixin — the only
cross-mixin handoff in the app that goes through the queue instead of a
direct method call.
"""
from __future__ import annotations

import platform
import queue
import re
import shutil
import subprocess
import sys
import threading
import time
from pathlib import Path

import psutil

import customtkinter as ctk

from ui.common import (
    PROJECT_ROOT, SCRIPTS_DIR,
    _GRAY, _BLUE, _RED, _GREEN, _YELLOW,
    _fmt_time, _load_runtimes, _save_runtimes, _detect_run_type,
    _CACHE_ITEMS, _cache_entry_count,
    PTM_PROXIMITY_STEPS, MUTATION_CLUSTERING_STEPS,
    input_dir, resolve_input_file, COSMIC_INPUT_DIR, COSMIC_SOMATIC_STATUSES,
)


class PipelineRunnerMixin:
    def _tick(self) -> None:
        if not self._running or self._pipeline_start is None:
            return
        if self._suspended:
            return
        elapsed = time.time() - self._pipeline_start

        text = f"Elapsed: {_fmt_time(elapsed)}"

        est = getattr(self, "_precheck_estimate", None)
        hist = self._historical_times
        if est is not None:
            text += f"  |  Est. total: ~{_fmt_time(est)}"
        elif hist and len(hist) == self._total_steps:
            text += f"  |  Est. total: ~{_fmt_time(sum(hist))}"
        elif self._step_times:
            avg = sum(self._step_times) / len(self._step_times)
            text += f"  |  Est. total: ~{_fmt_time(avg * self._total_steps)}"

        self._timer_label.configure(text=text)
        self.after(1000, self._tick)

    # ── Activity animation ─────────────────────────────────────────────────

    def _start_activity_animation(self):
        self._activity_animating = True
        self._activity_pos = -self._activity_chunk
        self._animate_activity()

    def _stop_activity_animation(self):
        self._activity_animating = False
        self._activity_canvas.delete("chunk")

    def _animate_activity(self):
        if not self._activity_animating:
            return
        self._activity_canvas.delete("chunk")
        x1 = max(0, self._activity_pos)
        x2 = min(self._activity_width, self._activity_pos + self._activity_chunk)
        if x2 > x1:
            self._activity_canvas.create_rectangle(
                x1, 0, x2, self._activity_height,
                fill=_BLUE, outline="", tags="chunk",
            )
        self._activity_pos += 1
        if self._activity_pos > self._activity_width:
            self._activity_pos = -self._activity_chunk
        self.after(25, self._animate_activity)

    @property
    def _output_dir(self) -> Path:
        return Path(self._output_dir_var.get())

    # ── Pipeline execution ───────────────────────────────────────────────────

    _DEFAULT_OUTPUT_FILES = [
        "ptm_mutation_proximity_db.tsv",
        "ptm_mutation_proximity_long.tsv",
        "mutation_cluster_db.tsv",
    ]

    @property
    def _output_files(self) -> list[Path]:
        return [self._output_dir / name for name in self._DEFAULT_OUTPUT_FILES]

    def _find_locked_output_file(self) -> Path | None:
        """Return the first existing output file that cannot be opened for writing, or None."""
        for path in self._output_files:
            if not path.exists():
                continue
            try:
                with path.open("r+b"):
                    pass
            except PermissionError:
                return path
        return None

    def _backup_outputs(self) -> list[Path]:
        """Copy existing output files to .bak before the pipeline overwrites them."""
        backups = []
        for path in self._output_files:
            if path.exists():
                bak = path.with_suffix(path.suffix + ".bak")
                shutil.copy2(path, bak)
                backups.append(bak)
        return backups

    @staticmethod
    def _restore_backups(backups: list[Path]) -> None:
        """Restore .bak files over their originals (undo a partial pipeline run)."""
        for bak in backups:
            if bak.exists():
                original = bak.with_suffix("")  # strip .bak
                shutil.copy2(bak, original)

    def _remove_backups(self, backups: list[Path]) -> None:
        for bak in backups:
            try:
                if bak.exists():
                    bak.unlink()
            except Exception as exc:
                self._q("log", f"Warning: could not delete {bak.name}: {exc}")

    def _start_pipeline(self):
        if self._running:
            return

        locked = self._find_locked_output_file()
        if locked:
            from tkinter import messagebox
            messagebox.showerror(
                "Output File Not Writable",
                f"Cannot start the pipeline — the output file\n\n"
                f"    {locked.name}\n\n"
                f"cannot be written to. Common causes:\n"
                f"  • Open in another program (e.g. Excel)\n"
                f"  • File is marked read-only\n"
                f"  • A cloud sync service (OneDrive, Dropbox) has it locked\n"
                f"  • You do not have write permission on the output folder\n\n"
                f"Close or resolve the issue above and try again.",
            )
            return

        mode = self._mode.get()
        if mode == "single-protein":
            cif = getattr(self, "_single_cif_var", None)
            if not cif or not cif.get().strip():
                from tkinter import messagebox
                messagebox.showwarning("Missing input", "Please select a CIF file first.")
                return
        elif mode == "ca-coordinates":
            from tkinter import messagebox
            if not self._ca_uniprot_var.get().strip():
                messagebox.showwarning("Missing input", "Enter a UniProt accession.")
                return

        self._running = True
        self._stop_requested = False
        self._suspended = False
        self._set_run_controls(True)
        self._stop_btn.configure(state="normal", fg_color=_RED)
        self._timer_label.configure(text="")

        self._log.configure(state="normal")
        self._log.delete("1.0", "end")
        self._log.configure(state="disabled")

        for lbl in self._step_status_labels:
            lbl.configure(text="●  Waiting", text_color=_GRAY)
        for bar in self._step_progress_bars:
            bar.set(0)
            bar.grid_remove()

        mode = self._mode.get()
        threading.Thread(target=self._run_pipeline, args=(mode,), daemon=True).start()

    def _set_run_controls(self, running: bool) -> None:
        """Keep the Pipeline tab's and Analysis Tools tab's Run buttons in sync —
        only one run (of any kind) can be active at a time, app-wide, since both
        share the same `self._running` flag and background-thread execution engine.
        """
        state = "disabled" if running else "normal"
        self._run_btn.configure(state=state)
        self._at_run_btn.configure(state=state)

    def _start_analysis_tool_run(self, kind: str) -> None:
        """Validate and launch a Radius Sweep / CIF Variance run triggered from
        the Analysis Tools tab. Shares `self._running` with `_start_pipeline` —
        only one run of any kind can be active at a time.

        Unlike `_start_pipeline`, this skips the locked-output-file check and
        `_backup_outputs()`: those only guard the 3 main-pipeline TSVs
        (`_DEFAULT_OUTPUT_FILES`), which radius_sweep.tsv/cif_variance/ aren't
        part of.
        """
        if self._running:
            return

        from tkinter import messagebox
        if kind == "radius-sweep":
            genes = list(self._radius_genes)
            if not genes:
                messagebox.showwarning("Missing input", "Add at least one gene symbol.")
                return
            try:
                start = float(self._radius_start_var.get())
                stop = float(self._radius_stop_var.get())
                step = float(self._radius_step_var.get())
                if step <= 0 or stop <= start:
                    raise ValueError
            except ValueError:
                messagebox.showwarning(
                    "Invalid input",
                    "Radius start/stop/step must be numeric, with stop > start and step > 0.",
                )
                return
            try:
                min_cases = int(self._radius_min_cases_var.get())
                if min_cases < 1:
                    raise ValueError
            except ValueError:
                messagebox.showwarning(
                    "Invalid input", "Min samples must be a whole number of 1 or more.",
                )
                return
            ptm_tsv = PROJECT_ROOT / "data" / "steps" / "PTMD_COSMIC_hotspots_by_protein.tsv"
            if not ptm_tsv.exists():
                messagebox.showerror(
                    "Missing data",
                    f"Required file not found:\n{ptm_tsv}\n\n"
                    f"Run the PTM Proximity or Mutation Clustering pipeline (step 1) first.",
                )
                return
        else:  # "cif-variance"
            input_dir = Path(self._variance_input_dir_var.get().strip())
            n_cifs = len(list(input_dir.glob("*.cif"))) if input_dir.is_dir() else 0
            if n_cifs < 2:
                messagebox.showwarning(
                    "Missing input", f"Need at least 2 .cif files in\n{input_dir}\n\nFound {n_cifs}.",
                )
                return
            if self._variance_top_var.get().strip() and not self._variance_top_var.get().strip().isdigit():
                messagebox.showwarning("Invalid input", "Top N must be a whole number.")
                return
            for start_var, end_var, label in [
                (self._variance_range_start_var, self._variance_range_end_var, "Report range"),
                (self._variance_align_start_var, self._variance_align_end_var, "Align range"),
            ]:
                s, e = start_var.get().strip(), end_var.get().strip()
                if s or e:
                    if not (s.isdigit() and e.isdigit()) or int(s) >= int(e):
                        messagebox.showwarning(
                            "Invalid input", f"{label} must be two whole numbers with start < end.",
                        )
                        return

        self._running = True
        self._at_run_active = True
        self._stop_requested = False
        self._suspended = False
        self._set_run_controls(True)

        self._at_progress_bar.set(0)
        self._at_progress_bar.grid_remove()
        self._at_status_label.configure(text="●  Starting…", text_color=_GRAY)
        self._at_log.configure(state="normal")
        self._at_log.delete("1.0", "end")
        self._at_log.configure(state="disabled")

        target = self._run_radius_sweep if kind == "radius-sweep" else self._run_cif_variance
        threading.Thread(target=target, daemon=True).start()

    def _stop_pipeline(self):
        """Suspend the running subprocess immediately and show Resume/Cancel options."""
        if not self._running:
            return
        proc = self._current_proc
        actually_suspended = False
        if proc and proc.poll() is None:
            try:
                psutil.Process(proc.pid).suspend()
                self._suspended = True
                actually_suspended = True
            except psutil.NoSuchProcess:
                pass
        self._stop_btn.configure(state="disabled", fg_color="gray30")
        if actually_suspended:
            self._paused_elapsed = time.time() - self._pipeline_start
            self._stop_activity_animation()
            self._activity_canvas.pack_forget()
            self._run_btn.configure(
                state="normal", text="▶  Resume", command=self._resume_pipeline,
            )
            self._cancel_btn_ref = ctk.CTkButton(
                self._run_btn.master,
                text="✗  Cancel",
                command=self._cancel_pipeline,
                width=110,
                height=44,
                font=ctk.CTkFont(size=15, weight="bold"),
                fg_color=_RED,
                hover_color="#b82020",
            )
            self._cancel_btn_ref.pack(side="left", padx=(8, 0))
            self._q("log", "Pipeline paused.")

    def _resume_pipeline(self):
        """Resume the suspended subprocess."""
        proc = self._current_proc
        if proc and self._suspended:
            # Adjust start time so elapsed doesn't include the pause duration
            pause_duration = time.time() - (self._pipeline_start + self._paused_elapsed)
            self._pipeline_start += pause_duration
            try:
                psutil.Process(proc.pid).resume()
            except psutil.NoSuchProcess:
                pass
            self._suspended = False
        if hasattr(self, "_cancel_btn_ref"):
            self._cancel_btn_ref.destroy()
            del self._cancel_btn_ref
        self._run_btn.configure(
            state="disabled", text="▶  Run Pipeline", command=self._start_pipeline,
        )
        self._stop_btn.configure(state="normal", fg_color=_RED)
        self._activity_canvas.pack(side="right", padx=(0, 4))
        self._start_activity_animation()
        self._tick()
        self._q("log", "Pipeline resumed.")

    def _cancel_pipeline(self):
        """Kill the suspended subprocess and restore previous output."""
        self._stop_requested = True
        proc = self._current_proc
        if proc and proc.poll() is None:
            if self._suspended:
                try:
                    psutil.Process(proc.pid).resume()
                except psutil.NoSuchProcess:
                    pass
                self._suspended = False
            proc.terminate()
        if hasattr(self, "_cancel_btn_ref"):
            self._cancel_btn_ref.destroy()
            del self._cancel_btn_ref
        self._run_btn.configure(
            state="normal", text="▶  Run Pipeline", command=self._start_pipeline,
        )
        self._stop_btn.configure(state="disabled", fg_color="gray30")
        self._open_btn.configure(state="normal", fg_color=_BLUE, hover_color="#2563eb")
        for lbl in self._step_status_labels:
            lbl.configure(text="●  Waiting", text_color=_GRAY)
        for bar in self._step_progress_bars:
            bar.set(0)
            bar.grid_remove()

    # Per-item time estimates (seconds) for runtime calculation
    _TIME_PER_CIF_DOWNLOAD = 2.5
    _TIME_PER_UNIPROT_BATCH = 1.5       # ~100 IDs per batch
    _TIME_PER_1433_FETCH = 0.2          # 5 concurrent workers
    _TIME_PER_PP_FETCH = 0.05           # 10 concurrent workers
    _TIME_PER_KINASE_PREDICT = 0.5
    _TIME_PER_AIUPRED_FETCH = 2.0       # per-protein REST call, sequential
    _TIME_PER_PROTEIN_STEP3 = 0.35
    _TIME_STEP1_BASE = 40               # base time for filtering/merging
    _TIME_STEP4_BASE = 40               # base time for reading/writing the proximity DB

    def _run_precheck(self, mode: str) -> float:
        """Analyze caches and data to estimate total pipeline runtime. Returns seconds."""
        import pandas as pd

        self._q("log", "Initializing pipeline...")
        self._q("log", "")

        input_tsv = PROJECT_ROOT / "data" / "steps" / "PTMD_COSMIC_hotspots_by_protein.tsv"
        cache_dir = PROJECT_ROOT / "data" / "cache"
        models_dir = PROJECT_ROOT / "cif_models"

        # Count proteins from the intermediate TSV if a previous run produced one;
        # otherwise step 1 hasn't run yet, so estimate from the raw input files.
        n_proteins = 0
        if input_tsv.exists():
            try:
                df = pd.read_csv(input_tsv, sep="\t", usecols=["uniprot_id"], dtype=str)
                n_proteins = df["uniprot_id"].nunique()
            except Exception:
                pass

        if n_proteins == 0:
            # Estimate from COSMIC directly, applying the same somatic-status and
            # hotspot-recurrence filtering step 1 will apply, so the count approximates
            # genes that will actually survive filtering (not just genes mentioned
            # anywhere in COSMIC). This works for both modes: in ptm-proximity mode,
            # PTMD's PTM-site coverage is broad enough that COSMIC's hotspot threshold
            # — not the PTMD intersection — is the dominant bottleneck in practice, and
            # PTMD's own "Gene name" column is too sparse to use directly (most rows
            # only resolve to a gene via a UniProt API lookup, which isn't available
            # before step 1 has run).
            try:
                cosmic_file = resolve_input_file(input_dir(PROJECT_ROOT, COSMIC_INPUT_DIR), (".tsv",))
                cosmic_cols = ["GENE_SYMBOL", "MUTATION_AA", "COSMIC_SAMPLE_ID", "MUTATION_SOMATIC_STATUS"]
                cosmic_df = pd.read_csv(cosmic_file, sep="\t", usecols=cosmic_cols,
                                        dtype=str, low_memory=False)
                cosmic_df = cosmic_df[cosmic_df["MUTATION_SOMATIC_STATUS"].isin(COSMIC_SOMATIC_STATUSES)]
                try:
                    min_samples = int(self._min_samples_var.get().strip() or 3)
                except ValueError:
                    min_samples = 3
                affected_cases = cosmic_df.groupby(["GENE_SYMBOL", "MUTATION_AA"])["COSMIC_SAMPLE_ID"].nunique()
                hotspot_genes = affected_cases[affected_cases >= min_samples].index.get_level_values("GENE_SYMBOL")
                n_proteins = hotspot_genes.nunique()
                self._q("log", f"Step 1 hasn't run yet — estimated {n_proteins} proteins from COSMIC hotspots")
            except Exception:
                pass

        # Step 1: UniProt API cache
        gene_cache = cache_dir / "uniprot_gene_mapping.tsv"
        cached_genes = 0
        if gene_cache.exists():
            try:
                cached_genes = len(pd.read_csv(gene_cache, sep="\t", dtype=str))
            except Exception:
                pass
        uncached_batches = max(0, (n_proteins - cached_genes)) // 100 + 1 if n_proteins > cached_genes else 0
        step1_est = self._TIME_STEP1_BASE + uncached_batches * self._TIME_PER_UNIPROT_BATCH
        self._q("log", f"Step 1: {cached_genes} UniProt gene mappings cached")

        # Step 2: CIF downloads
        cifs_present = 0
        if models_dir.exists():
            cifs_present = sum(1 for d in models_dir.iterdir()
                               if d.is_dir() and any(d.glob("*.cif")))
        cifs_needed = max(0, n_proteins - cifs_present)
        step2_est = cifs_needed * self._TIME_PER_CIF_DOWNLOAD + 5
        self._q("log", f"Step 2: {cifs_present} CIF files present"
                + (f", ~{cifs_needed} to download" if cifs_needed else " (all cached)"))

        # Step 3: local computation
        step3_est = n_proteins * self._TIME_PER_PROTEIN_STEP3 if n_proteins else 60
        self._q("log", f"Step 3: {n_proteins} proteins to process")

        # Step 4: annotation caches (only for ptm-proximity)
        step4_est = 0
        if mode == "ptm-proximity":
            # 14-3-3
            cache_1433 = cache_dir / "1433pred"
            cached_1433 = len(list(cache_1433.glob("*.json"))) if cache_1433.exists() else 0
            uncached_1433 = max(0, n_proteins - cached_1433)

            # PolyPhen
            pp_cache = cache_dir / "polyphen.tsv"
            cached_pp = 0
            if pp_cache.exists():
                try:
                    cached_pp = len(pd.read_csv(pp_cache, sep="\t", dtype=str))
                except Exception:
                    pass

            # Kinase
            kin_cache = cache_dir / "kinase_predictions.tsv"
            cached_kin = 0
            if kin_cache.exists():
                try:
                    cached_kin = len(pd.read_csv(kin_cache, sep="\t", dtype=str))
                except Exception:
                    pass

            # AIUPred — one binding call per protein (yields general + binding scores)
            cached_aiupred = 0
            aiupred_cache = cache_dir / "aiupred_disorder.tsv"
            if aiupred_cache.exists():
                try:
                    df_ac = pd.read_csv(aiupred_cache, sep="\t", dtype=str,
                                        keep_default_na=False)
                    cached_aiupred = int((df_ac["analysis_type"] == "binding").sum())
                except Exception:
                    pass
            uncached_aiupred = max(0, n_proteins - cached_aiupred)

            self._q("log", f"Step 4: {cached_1433}/{n_proteins} 14-3-3 predictions cached, "
                    f"{cached_pp} PolyPhen pairs cached, {cached_kin} kinase windows cached, "
                    f"{cached_aiupred}/{n_proteins} AIUPred cached")

            step4_est = (uncached_1433 * self._TIME_PER_1433_FETCH
                         + max(0, n_proteins * 2 - cached_pp) * self._TIME_PER_PP_FETCH
                         + max(0, n_proteins - cached_kin) * self._TIME_PER_KINASE_PREDICT
                         + uncached_aiupred * self._TIME_PER_AIUPRED_FETCH
                         + self._TIME_STEP4_BASE)

        total_est = step1_est + step2_est + step3_est + step4_est
        self._q("log", "")
        self._q("log", f"Estimated runtime: ~{_fmt_time(total_est)}")
        self._q("log", "")

        return total_est

    def _run_pipeline(self, mode: str):
        if mode == "single-protein":
            self._run_single_protein()
            return
        if mode == "ca-coordinates":
            self._run_ca_coordinates()
            return

        python = [sys.executable, "-u"]
        input_tsv = PROJECT_ROOT / "data" / "steps" / "PTMD_COSMIC_hotspots_by_protein.tsv"
        models_dir = PROJECT_ROOT / "cif_models"

        cutoff = self._cutoff_var.get().strip() or "10.0"
        min_samples = self._min_samples_var.get().strip() or "3"
        min_plddt = self._min_plddt_var.get().strip()
        max_pae = self._max_pae_var.get().strip()
        pp_exclude = []
        if not self._pp_benign_var.get():
            pp_exclude.append("benign")
        if not self._pp_possibly_var.get():
            pp_exclude.append("possibly_damaging")
        if not self._pp_probably_var.get():
            pp_exclude.append("probably_damaging")

        cmds = [
            [*python, str(SCRIPTS_DIR / "1_filter.py"), "--mode", mode,
             "--min-samples", min_samples],
            [
                *python, str(SCRIPTS_DIR / "2_download_structures.py"),
                str(input_tsv),
                "--id_column", "uniprot_id",
                "--out_dir", str(models_dir),
                "--prefer", "cif",
                "--delay", "0.1",
                "--also_pae",
                "--logs_dir", str(self._output_dir / "logs"),
            ],
            [*python, str(SCRIPTS_DIR / "3_find_nearby_mutations.py"), "--mode", mode,
             "--output-dir", str(self._output_dir), "--cutoff", cutoff,
             *(["--min-plddt", min_plddt] if min_plddt else []),
             *(["--max-pae", max_pae] if max_pae else [])],
        ]
        if mode == "ptm-proximity":
            cmds.append([*python, str(SCRIPTS_DIR / "4_annotate.py"),
                         "--output-dir", str(self._output_dir),
                         *(["--pp-exclude"] + pp_exclude if pp_exclude else [])])

        steps = PTM_PROXIMITY_STEPS if mode == "ptm-proximity" else MUTATION_CLUSTERING_STEPS
        run_type = _detect_run_type()
        self._q("pipeline_start", len(steps), mode, run_type)

        estimated_total = self._run_precheck(mode)
        self._q("set_estimate", estimated_total)

        backups = self._backup_outputs()

        try:
            all_ok = True
            cancelled = False
            step_times: list[float] = []
            for i, (step_def, cmd) in enumerate(zip(steps, cmds)):
                _panel_label, log_label = step_def
                if self._stop_requested:
                    self._q("status", i, "●  Cancelled", _YELLOW)
                    cancelled = True
                    break

                self._q("status", i, "▶  Running…", _BLUE)
                self._q("show_progress", i)
                self._q("log", f"Step {i + 1}/{len(steps)}: {log_label}")
                self._q("step_start")

                t0 = time.time()
                ok = self._stream_cmd(cmd, i)
                elapsed = time.time() - t0

                if self._stop_requested:
                    self._q("hide_progress", i)
                    self._q("status", i, "■  Stopped", _YELLOW)
                    self._q("log", f"Step {i + 1} stopped by user")
                    cancelled = True
                    break

                if ok:
                    step_times.append(elapsed)
                    self._q("progress", i, 1.0, "Done")
                    self._q("hide_progress", i)
                    self._q("status", i, f"✓  {_fmt_time(elapsed)}", _GREEN)
                    self._q("step_complete", elapsed)
                else:
                    self._q("hide_progress", i)
                    self._q("status", i, "✗  Failed", _RED)
                    self._q("log", f"Step {i + 1} FAILED after {_fmt_time(elapsed)}")
                    all_ok = False
                    break

            restore_ok = False
            if cancelled:
                try:
                    self._restore_backups(backups)
                    self._q("log", "Pipeline cancelled — previous output restored.")
                    restore_ok = True
                except Exception:
                    self._q("log", "Pipeline cancelled. Output file may be locked — "
                            "close it and rename the .bak file to restore your previous output.")
            elif not all_ok:
                try:
                    self._restore_backups(backups)
                    self._q("log", "Pipeline failed — previous output restored.")
                    restore_ok = True
                except Exception:
                    self._q("log", "Pipeline failed. Output file may be locked — "
                            "close it and rename the .bak file to restore your previous output.")
            else:
                restore_ok = True
                self._q("save_runtimes", mode, run_type, step_times)
                self._q("log", "Pipeline complete! Output saved to the Output/ folder.")
        except Exception as exc:
            restore_ok = False
            self._q("log", f"Pipeline error: {exc}")
        finally:
            if restore_ok:
                self._remove_backups(backups)
            self._q("finished")

    def _run_single_protein(self):
        """Run the single-protein CIF analysis in the background thread."""
        cif_path = self._single_cif_var.get().strip()
        uniprot = self._single_uniprot_var.get().strip()

        if not cif_path:
            self._q("log", "Error: no CIF file selected.")
            self._q("finished")
            return

        cmd = [sys.executable, "-u", str(SCRIPTS_DIR / "analyze_single_cif_nearby_mutations.py"), cif_path]
        if uniprot:
            cmd.extend(["--uniprot", uniprot])

        self._q("pipeline_start", 1, "single-protein", "warm")
        self._q("show_log")
        self._q("status", 0, "▶  Analyzing…", _BLUE)
        self._q("show_progress", 0)
        self._q("log", f"Analyzing CIF: {Path(cif_path).name}")
        if uniprot:
            self._q("log", f"UniProt ID: {uniprot}")
        self._q("log", "")

        t0 = time.time()
        ok = self._stream_cmd(cmd, 0)
        elapsed = time.time() - t0

        if self._stop_requested:
            self._q("hide_progress", 0)
            self._q("status", 0, "■  Stopped", _YELLOW)
            self._q("finished")
            return

        if not ok:
            self._q("hide_progress", 0)
            self._q("status", 0, "✗  Failed", _RED)
            self._q("log", f"Analysis failed after {_fmt_time(elapsed)}")
            self._q("finished")
            return

        self._q("progress", 0, 1.0, "Done")
        self._q("hide_progress", 0)
        self._q("status", 0, f"✓  {_fmt_time(elapsed)}", _GREEN)
        self._q("log", "")
        self._q("log", f"Analysis complete in {_fmt_time(elapsed)}.")

        output_db = self._output_dir / "ptm_mutation_proximity_db.tsv"
        if output_db.exists():
            self._q("log", "Review the results above, then choose whether to append to the existing database.")
            self._q("offer_append", cif_path, uniprot)
        else:
            self._q("log", "No existing proximity database found — results printed above.")

        self._q("enable_open")
        self._q("finished")

    def _run_radius_sweep(self):
        """Run the radius-sweep analysis in-process, in the background thread."""
        import numpy as np

        genes = list(self._radius_genes)
        try:
            start = float(self._radius_start_var.get())
            stop = float(self._radius_stop_var.get())
            step = float(self._radius_step_var.get())
        except ValueError:
            self._q("log", "Error: invalid radius range.")
            self._q("finished")
            return
        radii = list(np.arange(start, stop + step / 2, step))
        unfiltered = self._radius_unfiltered_var.get()
        min_cases = int(self._radius_min_cases_var.get())

        self._q("pipeline_start", 1, "radius-sweep", "warm")
        self._q("show_log")
        self._q("status", 0, "▶  Running sweep…", _BLUE)
        self._q("show_progress", 0)
        self._q("log", f"Testing radii {radii[0]:.0f}-{radii[-1]:.0f} Å in {step:.0f} Å steps "
                        f"for {len(genes)} gene(s) (hotspot threshold: >= {min_cases} samples)")
        self._q("log", "")

        t0 = time.time()
        try:
            from radius_sweep import run_sweep
            result = run_sweep(
                genes, radii, min_cases=min_cases, unfiltered=unfiltered,
                output_tsv_path=self._output_dir / "radius_sweep.tsv",
                log_cb=lambda line: self._q("log", line),
            )
        except ImportError as exc:
            self._q("status", 0, "✗  Failed", _RED)
            self._q("log", f"Missing dependency: {exc}. Run: uv sync")
            self._q("hide_progress", 0)
            self._q("finished")
            return
        except (FileNotFoundError, ValueError) as exc:
            self._q("status", 0, "✗  Failed", _RED)
            self._q("log", f"Error: {exc}")
            self._q("hide_progress", 0)
            self._q("finished")
            return
        except Exception as exc:
            self._q("status", 0, "✗  Failed", _RED)
            self._q("log", f"Unexpected error: {exc}")
            self._q("hide_progress", 0)
            self._q("finished")
            return

        elapsed = time.time() - t0
        self._q("hide_progress", 0)
        self._q("status", 0, f"✓  {_fmt_time(elapsed)}", _GREEN)
        self._q("log", "")
        self._q("log", f"Sweep complete in {_fmt_time(elapsed)}.")
        self._q("viz_data", "radius-sweep", result)
        self._q("enable_open")
        self._q("finished")

    def _run_cif_variance(self):
        """Run the CIF variance analysis in-process, in the background thread."""
        input_dir = Path(self._variance_input_dir_var.get().strip())
        top = int(self._variance_top_var.get().strip() or 10)

        def _range_or_none(start_var, end_var):
            s, e = start_var.get().strip(), end_var.get().strip()
            return (int(s), int(e)) if s and e else None

        range_ = _range_or_none(self._variance_range_start_var, self._variance_range_end_var)
        align_range = _range_or_none(self._variance_align_start_var, self._variance_align_end_var)
        uniprot = self._variance_uniprot_var.get().strip() or None
        gene = self._variance_gene_var.get().strip() or None

        self._q("pipeline_start", 1, "cif-variance", "warm")
        self._q("show_log")
        self._q("status", 0, "▶  Running analysis…", _BLUE)
        self._q("show_progress", 0)
        self._q("log", f"Comparing CIF files in {input_dir}")
        self._q("log", "")

        t0 = time.time()
        try:
            from cif_variance import run_variance_analysis
            result = run_variance_analysis(
                input_dir=input_dir,
                output_dir=self._output_dir / "cif_variance",
                top=top,
                range_=range_,
                align_range=align_range,
                uniprot=uniprot,
                gene=gene,
                log_cb=lambda line: self._q("log", line),
            )
        except ImportError as exc:
            self._q("status", 0, "✗  Failed", _RED)
            self._q("log", f"Missing dependency: {exc}. Run: uv sync")
            self._q("hide_progress", 0)
            self._q("finished")
            return
        except (FileNotFoundError, ValueError) as exc:
            self._q("status", 0, "✗  Failed", _RED)
            self._q("log", f"Error: {exc}")
            self._q("hide_progress", 0)
            self._q("finished")
            return
        except Exception as exc:
            self._q("status", 0, "✗  Failed", _RED)
            self._q("log", f"Unexpected error: {exc}")
            self._q("hide_progress", 0)
            self._q("finished")
            return

        elapsed = time.time() - t0
        self._q("hide_progress", 0)
        self._q("status", 0, f"✓  {_fmt_time(elapsed)}", _GREEN)
        self._q("log", "")
        self._q("log", f"Analysis complete in {_fmt_time(elapsed)}.")
        self._q("viz_data", "cif-variance", result)
        self._q("enable_open")
        self._q("finished")

    def _run_ca_coordinates(self):
        """Run the CA-coordinate export in-process, in the background thread."""
        uniprot = self._ca_uniprot_var.get().strip()
        gene = self._ca_gene_var.get().strip() or None

        self._q("pipeline_start", 1, "ca-coordinates", "warm")
        self._q("show_log")
        self._q("status", 0, "▶  Exporting…", _BLUE)
        self._q("show_progress", 0)
        self._q("log", f"Exporting CA coordinates for {uniprot}")
        self._q("log", "")

        t0 = time.time()
        try:
            from export_ca_coordinates import run_export
            run_export(
                uniprot, gene=gene,
                output_dir=self._output_dir / "coordinates",
                log_cb=lambda line: self._q("log", line),
            )
        except ImportError as exc:
            self._q("status", 0, "✗  Failed", _RED)
            self._q("log", f"Missing dependency: {exc}. Run: uv sync")
            self._q("hide_progress", 0)
            self._q("finished")
            return
        except (FileNotFoundError, ValueError, RuntimeError) as exc:
            self._q("status", 0, "✗  Failed", _RED)
            self._q("log", f"Error: {exc}")
            self._q("hide_progress", 0)
            self._q("finished")
            return
        except Exception as exc:
            self._q("status", 0, "✗  Failed", _RED)
            self._q("log", f"Unexpected error: {exc}")
            self._q("hide_progress", 0)
            self._q("finished")
            return

        elapsed = time.time() - t0
        self._q("hide_progress", 0)
        self._q("status", 0, f"✓  {_fmt_time(elapsed)}", _GREEN)
        self._q("enable_open")
        self._q("finished")

    _TQDM_RE = re.compile(r"^(.*?):\s+(\d+)%\|[^|]*\|\s+(\d+)/(\d+)")
    _OVERALL_RE = re.compile(r"^##PROGRESS##\s+(\d+)\s+(.+)$")
    _using_overall_progress = False

    def _parse_tqdm(self, text: str, step_idx: int) -> bool:
        """Try to parse a progress line (custom ##PROGRESS## or tqdm). Returns True if matched."""
        m_overall = self._OVERALL_RE.search(text)
        if m_overall:
            self._using_overall_progress = True
            pct = int(m_overall.group(1))
            desc = m_overall.group(2)
            self._q("progress", step_idx, pct / 100, f"{pct}%")
            self._q("progress_log", desc)
            return True
        m = self._TQDM_RE.search(text)
        if not m:
            return False
        if self._using_overall_progress:
            return True
        desc, pct, current, total = m.group(1), int(m.group(2)), m.group(3), m.group(4)
        self._q("progress", step_idx, pct / 100, f"{pct}% ({current}/{total})")
        self._q("progress_log", f"{desc}: {current}/{total}")
        return True

    def _stream_cmd(self, cmd: list[str], step_idx: int) -> bool:
        """Run a subprocess, parsing tqdm output for progress bar updates."""
        env = {**__import__("os").environ, "PYTHONIOENCODING": "utf-8"}
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, env=env)
        self._current_proc = proc
        buf = b""
        while True:
            byte = proc.stdout.read(1)
            if not byte:
                break
            if self._suspended:
                buf = b""
                continue
            if byte == b"\r":
                text = buf.decode("utf-8", errors="replace")
                if self._parse_tqdm(text, step_idx):
                    buf = b""
            elif byte == b"\n":
                text = buf.decode("utf-8", errors="replace").strip()
                if text and not self._parse_tqdm(text, step_idx):
                    self._q("log", text)
                buf = b""
            else:
                buf += byte

        if buf:
            text = buf.decode("utf-8", errors="replace").strip()
            if text:
                self._q("log", text)

        proc.wait()
        return proc.returncode == 0

    # ── Queue helpers ────────────────────────────────────────────────────────

    def _q(self, *args):
        self._queue.put(args)

    def _status_target(self, idx: int):
        """Resolve which status label a "status"/"progress" queue message should
        update: the Analysis Tools tab's single persistent label when a run was
        triggered from there, otherwise the Pipeline tab's step-indexed one.

        Not a cached list reference — `_rebuild_step_rows` rebinds
        `self._step_status_labels` to a brand-new list on every Pipeline-tab
        mode change, so anything resolved once up front would go stale.
        """
        if self._at_run_active:
            return self._at_status_label
        return self._step_status_labels[idx] if 0 <= idx < len(self._step_status_labels) else None

    def _progress_target(self, idx: int):
        if self._at_run_active:
            return self._at_progress_bar
        return self._step_progress_bars[idx] if 0 <= idx < len(self._step_progress_bars) else None

    def _poll_queue(self):
        try:
            while True:
                msg = self._queue.get_nowait()
                kind = msg[0]
                if kind == "log":
                    self._append_log(msg[1])
                elif kind == "show_log":
                    if self._at_run_active:
                        if not self._at_log_visible:
                            self._at_toggle_log()
                    elif not self._log_visible:
                        self._toggle_log()
                elif kind == "status":
                    _, idx, text, color = msg
                    target = self._status_target(idx)
                    if target is not None:
                        target.configure(text=text, text_color=color)
                elif kind == "progress":
                    _, idx, pct, status_text = msg
                    bar = self._progress_target(idx)
                    status = self._status_target(idx)
                    if bar is not None:
                        bar.set(pct)
                    if status is not None:
                        status.configure(text=status_text, text_color=_BLUE)
                elif kind == "progress_log":
                    self._update_progress_line(msg[1])
                elif kind == "show_progress":
                    _, idx = msg
                    self._using_overall_progress = False
                    bar = self._progress_target(idx)
                    if bar is not None:
                        bar.set(0)
                        bar.grid()
                elif kind == "hide_progress":
                    _, idx = msg
                    bar = self._progress_target(idx)
                    if bar is not None:
                        bar.grid_remove()
                elif kind == "pipeline_start":
                    _, total, mode, run_type = msg
                    self._pipeline_start = time.time()
                    self._step_start = None
                    self._total_steps = total
                    self._steps_done = 0
                    self._step_times = []
                    self._historical_times = _load_runtimes(mode, run_type)
                    self._precheck_estimate = None
                    self._activity_canvas.pack(side="right", padx=(0, 8))
                    self._start_activity_animation()
                    self._tick()
                elif kind == "set_estimate":
                    self._precheck_estimate = msg[1]
                elif kind == "step_start":
                    self._step_start = time.time()
                elif kind == "step_complete":
                    _, elapsed = msg
                    self._steps_done += 1
                    self._step_times.append(elapsed)
                    self._step_start = None
                elif kind == "save_runtimes":
                    _, mode, run_type, times = msg
                    _save_runtimes(mode, run_type, times)
                elif kind == "offer_append":
                    _, cif_path, uniprot = msg
                    self._show_append_dialog(cif_path, uniprot)
                elif kind == "viz_data":
                    _, which, result = msg
                    self._tabview.set("Analysis Tools")
                    if which == "radius-sweep":
                        self._radius_sweep_result = result
                        self._analysis_subtool_var.set("Radius Sweep")
                        self._draw_radius_sweep_plot()
                    elif which == "cif-variance":
                        self._cif_variance_result = result
                        self._analysis_subtool_var.set("CIF Variance")
                        self._draw_cif_variance_plot()
                elif kind == "enable_open":
                    self._open_btn.configure(state="normal", fg_color=_BLUE, hover_color="#2563eb")
                elif kind == "finished":
                    was_cancelled = self._stop_requested
                    self._running = False
                    self._at_run_active = False
                    self._stop_requested = False
                    self._suspended = False
                    self._current_proc = None
                    self._stop_activity_animation()
                    self._activity_canvas.pack_forget()
                    if hasattr(self, "_cancel_btn_ref"):
                        self._cancel_btn_ref.destroy()
                        del self._cancel_btn_ref
                    self._run_btn.configure(
                        state="normal", text="▶  Run Pipeline",
                        command=self._start_pipeline,
                    )
                    self._at_run_btn.configure(state="normal")
                    self._stop_btn.configure(state="disabled", fg_color="gray30")
                    self._open_btn.configure(state="normal", fg_color=_BLUE, hover_color="#2563eb")
                    if was_cancelled:
                        for lbl in self._step_status_labels:
                            lbl.configure(text="●  Waiting", text_color=_GRAY)
                        for bar in self._step_progress_bars:
                            bar.set(0)
                            bar.grid_remove()
                    if self._pipeline_start is not None:
                        total = time.time() - self._pipeline_start
                        self._timer_label.configure(
                            text=f"Total: {_fmt_time(total)}", text_color=_GREEN
                        )
        except queue.Empty:
            pass
        self.after(100, self._poll_queue)

    def _append_log(self, line: str):
        log = self._at_log if self._at_run_active else self._log
        log.configure(state="normal")
        log.insert("end", line + "\n")
        log.see("end")
        log.configure(state="disabled")
        self._last_log_was_progress = False

    def _update_progress_line(self, text: str):
        """Replace the last line in the log if it was a progress update, otherwise append."""
        self._log.configure(state="normal")
        if getattr(self, "_last_log_was_progress", False):
            self._log.delete("end-2l", "end-1l")
        self._log.insert("end", text + "\n")
        self._log.see("end")
        self._log.configure(state="disabled")
        self._last_log_was_progress = True

    # ── Single-protein append dialog ────────────────────────────────────────

    def _check_existing_uniprot(self, uniprot: str) -> int:
        """Return the number of rows in the proximity DB that match this UniProt ID."""
        output_db = self._output_dir / "ptm_mutation_proximity_db.tsv"
        if not output_db.exists() or not uniprot:
            return 0
        try:
            import pandas as pd
            df = pd.read_csv(output_db, sep="\t", encoding="utf-16", usecols=["UniProt"],
                             dtype=str, keep_default_na=False)
            return int((df["UniProt"] == uniprot).sum())
        except Exception:
            return 0

    def _remove_uniprot_rows(self, uniprot: str) -> None:
        """Remove all rows for a UniProt ID from the proximity DB."""
        output_db = self._output_dir / "ptm_mutation_proximity_db.tsv"
        if not output_db.exists():
            return
        try:
            import pandas as pd
            df = pd.read_csv(output_db, sep="\t", encoding="utf-16", dtype=str,
                             keep_default_na=False)
            before = len(df)
            df = df[df["UniProt"] != uniprot]
            df.to_csv(output_db, sep="\t", index=False, encoding="utf-16")
            self._q("log", f"Removed {before - len(df)} existing row(s) for {uniprot}.")
        except Exception as exc:
            self._q("log", f"Warning: could not remove existing rows: {exc}")

    def _show_append_dialog(self, cif_path: str, uniprot: str):
        """Ask the user whether to append, replace, or skip based on existing data."""
        existing_count = self._check_existing_uniprot(uniprot)

        dialog = ctk.CTkToplevel(self)
        dialog.title("Save to Database?")
        dialog.resizable(False, False)
        dialog.grab_set()

        dialog.geometry("500x190")

        if existing_count > 0:
            msg = (f"The output database already contains {existing_count} row(s) "
                   f"for {uniprot}.\nHow would you like to proceed?")
        else:
            msg = (f"The output database has no existing data for {uniprot}.\n"
                   f"How would you like to proceed?")

        ctk.CTkLabel(
            dialog, text=msg, font=ctk.CTkFont(size=14), justify="center",
        ).pack(pady=(20, 16))

        btn_frame = ctk.CTkFrame(dialog, fg_color="transparent")
        btn_frame.pack()

        def do_replace():
            dialog.destroy()
            threading.Thread(
                target=self._replace_and_append, args=(cif_path, uniprot),
                daemon=True,
            ).start()

        def do_append():
            dialog.destroy()
            threading.Thread(
                target=self._append_single_protein, args=(cif_path, uniprot),
                daemon=True,
            ).start()

        if existing_count > 0:
            ctk.CTkButton(
                btn_frame, text="Replace existing", width=140,
                command=do_replace,
            ).pack(side="left", padx=6)

        ctk.CTkButton(
            btn_frame, text="Append to database", width=140,
            command=do_append,
        ).pack(side="left", padx=6)

        ctk.CTkButton(
            btn_frame, text="Skip", width=100,
            fg_color="gray30", hover_color="gray40",
            command=dialog.destroy,
        ).pack(side="left", padx=6)

    def _replace_and_append(self, cif_path: str, uniprot: str):
        """Remove existing rows for this UniProt ID, then append new results."""
        self._remove_uniprot_rows(uniprot)
        self._append_single_protein(cif_path, uniprot)

    def _append_single_protein(self, cif_path: str, uniprot: str):
        """Re-run the single-protein script with --append-to-db."""
        cmd = [
            sys.executable, "-u",
            str(SCRIPTS_DIR / "analyze_single_cif_nearby_mutations.py"),
            cif_path, "--append-to-db",
        ]
        if uniprot:
            cmd.extend(["--uniprot", uniprot])
        self._q("log", "Appending results to proximity database...")
        proc = subprocess.Popen(
            cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
            text=True, encoding="utf-8", errors="replace",
        )
        for line in proc.stdout:
            self._q("log", line.rstrip())
        proc.wait()
        if proc.returncode == 0:
            self._q("log", "Results appended successfully.")
        else:
            self._q("log", "Failed to append results.")

    # ── Cache management ─────────────────────────────────────────────────────

    def _manage_cache(self):
        import shutil

        dialog = ctk.CTkToplevel(self)
        dialog.title("Manage Cache")
        dialog.resizable(False, False)
        dialog.grab_set()

        ctk.CTkLabel(
            dialog, text="Manage Cache",
            font=ctk.CTkFont(size=15, weight="bold"),
        ).grid(row=0, column=0, columnspan=4, padx=20, pady=(16, 2), sticky="w")
        ctk.CTkLabel(
            dialog,
            text="Clear cached API results to force re-fetching on the next run.",
            text_color=_GRAY,
        ).grid(row=1, column=0, columnspan=4, padx=20, pady=(0, 10), sticky="w")

        for col, txt in enumerate(["Step", "Cache", "Entries", ""]):
            ctk.CTkLabel(
                dialog, text=txt, font=ctk.CTkFont(weight="bold"), text_color=_GRAY,
            ).grid(row=2, column=col, padx=(20 if col == 0 else 8, 8), pady=(0, 4), sticky="w")

        count_labels: list = []
        last_step = None

        def _clear(path: Path, is_dir: bool, lbl):
            if is_dir:
                shutil.rmtree(path, ignore_errors=True)
            else:
                path.unlink(missing_ok=True)
            lbl.configure(text="empty")

        def _clear_all():
            for _, _, path, is_dir in _CACHE_ITEMS:
                if is_dir:
                    shutil.rmtree(path, ignore_errors=True)
                else:
                    path.unlink(missing_ok=True)
            for lbl in count_labels:
                lbl.configure(text="empty")

        for i, (step, name, path, is_dir) in enumerate(_CACHE_ITEMS):
            row = 3 + i
            step_text = step if step != last_step else ""
            last_step = step

            ctk.CTkLabel(dialog, text=step_text, text_color=_GRAY).grid(
                row=row, column=0, padx=(20, 8), pady=3, sticky="w")
            ctk.CTkLabel(dialog, text=name).grid(
                row=row, column=1, padx=8, pady=3, sticky="w")

            count_lbl = ctk.CTkLabel(
                dialog, text=_cache_entry_count(path, is_dir),
                text_color=_GRAY, width=90, anchor="e",
            )
            count_lbl.grid(row=row, column=2, padx=8, pady=3, sticky="e")
            count_labels.append(count_lbl)

            ctk.CTkButton(
                dialog, text="Clear", width=64,
                fg_color="gray30", hover_color=_RED,
                command=lambda p=path, d=is_dir, lbl=count_lbl: _clear(p, d, lbl),
            ).grid(row=row, column=3, padx=(0, 20), pady=3)

        bottom_row = 3 + len(_CACHE_ITEMS)
        ctk.CTkButton(
            dialog, text="Clear All Caches",
            fg_color="gray30", hover_color=_RED,
            command=_clear_all,
        ).grid(row=bottom_row, column=0, columnspan=4, padx=20, pady=(10, 16), sticky="ew")

        dialog.after(50, dialog.lift)

    # ── Output folder ────────────────────────────────────────────────────────

    def _open_output_folder(self):
        out = self._output_dir
        out.mkdir(parents=True, exist_ok=True)
        if platform.system() == "Windows":
            subprocess.run(["explorer", str(out)], check=False)
        elif platform.system() == "Darwin":
            subprocess.run(["open", str(out)], check=False)
        else:
            subprocess.run(["xdg-open", str(out)], check=False)
