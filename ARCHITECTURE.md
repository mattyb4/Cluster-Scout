# Cluster-Scout — Architecture & Developer Guide

This document is for **developers** who will maintain or extend Cluster-Scout,
not for end users (see `README.md` for installation/usage and `docs/help.md`
for the in-app user documentation). It explains how the codebase is organized,
why it's structured the way it is, and where to make common kinds of changes.

---

## 1. What the app does, in one paragraph

Cluster-Scout finds recurrent cancer mutations that cluster near — or directly
disrupt — post-translational modification (PTM) sites in 3D protein structure.
It merges mutation data (COSMIC) with PTM data (PTMD), fetches predicted 3D
structures (AlphaFold DB), computes which mutations fall near PTM sites in
folded space, and annotates the results with several functional-prediction
sources. It ships as a desktop GUI so non-developer researchers can run the
whole thing without a command line, but every stage is also runnable headless.

---

## 2. The two-layer design

The codebase has a hard separation between **the pipeline** (the science) and
**the GUI** (the interface). This separation is the most important thing to
understand, because it's what lets the two be developed and tested
independently.

```
┌─────────────────────────────────────────────────┐
│  GUI layer  (app.py + ui/*.py)                    │
│  - CustomTkinter desktop app                      │
│  - Runs pipeline scripts as SUBPROCESSES          │
│  - Never imports pipeline logic to run it inline  │
└───────────────────┬─────────────────────────────┘
                    │ subprocess + stdout parsing
┌───────────────────▼─────────────────────────────┐
│  Pipeline layer  (main.py + scripts/*.py)         │
│  - Pure Python, no GUI dependency                 │
│  - Runnable standalone from the command line      │
│  - Reads/writes TSV files under data/ and Output/ │
└─────────────────────────────────────────────────┘
```

**Why subprocesses instead of importing the pipeline functions directly?**
It keeps the pipeline completely GUI-agnostic (it can run on a headless
server), it isolates long-running scientific work from the UI event loop so a
crash in step 3 can't take down the window, and it makes the GUI and CLI two
thin front-ends over the exact same scripts rather than two code paths that can
drift apart. The cost is that communication happens over stdout text rather
than function returns — see §6.

---

## 3. The pipeline layer

### 3.1 The numbered scripts

The pipeline is a sequence of numbered scripts under `scripts/`, run in order.
The numbers ARE the data flow — each step reads what the previous step wrote:

| Script | Does | Reads | Writes |
|---|---|---|---|
| `1_filter.py` | Merge + filter COSMIC and PTMD; map genes→UniProt | `data/input/` files | `data/steps/*_hotspots_by_protein.tsv` |
| `2_download_structures.py` | Fetch AlphaFold CIF + PAE per protein | step-1 TSV | `cif_models/{uniprot}/` |
| `3_find_nearby_mutations.py` | Compute 3D distances PTM↔mutation | step-1 TSV + CIFs | `Output/*_db.tsv` |
| `4_annotate.py` | Add 14-3-3, PolyPhen, kinase, AIUPred, InterPro | step-3 output | annotated `Output/*_db.tsv` |

Step 4 only runs in `ptm-proximity` mode (see §3.3).

### 3.2 `main.py` — the CLI orchestrator

`main.py` is a thin driver: it decides which steps to run for the chosen mode,
builds each step's command-line invocation, and runs them in order as
subprocesses, timing each. It contains no scientific logic itself. The GUI does
**not** call `main.py` — it orchestrates the same steps itself (see §6) so it
can stream progress. `main.py` is the headless equivalent of that
orchestration. If you add a pipeline step, both `main.py` and the GUI runner
need to know about it.

### 3.3 Pipeline modes

There are two modes, defined as `PTM_PROXIMITY_STEPS` and
`MUTATION_CLUSTERING_STEPS` in `pipeline_utils.py`:

- **ptm-proximity** (default) — the full 4-step pipeline; finds mutations near
  PTM sites.
- **mutation-clustering** — steps 1–3 only, no PTM requirement; finds mutations
  that cluster near each other.

The two modes write to **separate** step-1 output files
(`hotspots_tsv_path()`), because their schemas differ and sharing a path caused
one mode to silently overwrite the other's data. Keep this separation if you
add mode-specific columns.

### 3.4 `pipeline_utils.py` — shared backend utilities

Central module imported by every pipeline script. Holds shared constants
(input-folder names, required-column lists, the step definitions), file
resolution/validation helpers, CIF loading, and small formatters. If you find
yourself copy-pasting a helper between two scripts, it belongs here instead.

### 3.5 Standalone analysis scripts

Not part of the numbered pipeline; each is an independent tool that reuses
`pipeline_utils`:

- `analyze_single_cif_nearby_mutations.py` — run step-3-style analysis on a
  user-supplied CIF (for proteins AlphaFold DB doesn't cover).
- `export_ca_coordinates.py` — export CA coordinates + a ChimeraX heatmap of
  per-residue mutation burden.
- `radius_sweep.py` — sweep the distance cutoff to justify the 10 Å threshold.
- `cif_variance.py` — compare multiple CIF predictions of one protein
  (per-residue positional variance + pLDDT).
- `prototype_hotspot_significance.py` — permutation-test prototype for spatial
  clustering significance.

---

## 4. The GUI layer

### 4.1 `app.py` — window + mixin assembly

`app.py` is deliberately small. The `App` class is composed from one **mixin
per tab**, each in its own `ui/*.py` file:

```python
class App(PipelineTabMixin, PipelineRunnerMixin, ResultsTabMixin,
          VisualizationTabMixin, AnalysisToolsTabMixin, HelpTabMixin, ctk.CTk):
```

**Why mixins?** A single-file GUI for an app this size would be thousands of
lines. Splitting each tab into a mixin keeps each file focused on one tab's
widgets and logic, while still sharing one `self` (so tabs can read each
other's state — e.g. Visualization reading the results dataframe). The tradeoff
is that all mixins share one namespace, so `self._`-attribute names must not
collide across tabs; instance attributes are all initialized in `App.__init__`
so there's one place to see the full shared state.

### 4.2 The CustomTkinter monkeypatches

`app.py` applies four monkeypatches to CustomTkinter at import time. **These are
intentional and each is documented inline** with the bug it fixes and why it's
safe. They address macOS Tcl/Tk lookup, a textbox-scrollbar teardown crash,
DPI/appearance pollers that conflict with the app's tab-switching optimization,
and slow Cocoa tab switches. Don't remove them without reading the docstrings —
they encode real bugs found in use.

### 4.3 `ui/common.py` — shared UI utilities

The GUI-side counterpart to `pipeline_utils.py`: shared constants, the tooltip
widget, colour helpers, and the input-folder registry. Imported by every tab
mixin.

---

## 5. The five tabs

| Tab | Mixin file | Responsibility |
|---|---|---|
| Pipeline | `pipeline_panels.py` (+ `pipeline_runner.py`) | Input selection, settings, run/pause/resume |
| Results | `results_tab.py` | Sortable PTM-site table + per-mutation detail |
| Visualization | `visualization_tab.py` | Lollipop plots per PTM site |
| Analysis Tools | `analysis_tools_tab.py` | Radius sweep, CIF variance, exports |
| Help | `help_tab.py` | Renders `docs/help.md` in-app |

`pipeline_panels.py` builds the Pipeline tab's widgets; `pipeline_runner.py`
holds the execution engine (the largest UI module). They're split because
"what the tab looks like" and "how a run is driven" are separable concerns.

---

## 6. How the GUI runs the pipeline (the critical bridge)

This is the trickiest part of the codebase and worth understanding before
touching `pipeline_runner.py`.

1. When the user clicks Run, the GUI spawns a **background `threading.Thread`**
   (daemon) so the UI event loop stays responsive.
2. That thread runs each pipeline step as a **`subprocess.Popen`**, reading the
   child's stdout line by line.
3. Pipeline scripts print `tqdm` progress bars to stdout; the runner parses
   those lines to update the GUI progress bar.
4. The worker thread can't touch Tk widgets directly (Tk isn't thread-safe), so
   it pushes messages onto a **`queue.Queue`** (`self._queue`).
5. The main thread drains that queue every 100 ms via **`_poll_queue()`**
   (scheduled with `self.after`), and applies updates to widgets there.

```
User clicks Run
   → daemon Thread
       → subprocess.Popen(step N)  ──stdout──▶ parse tqdm
                                                   │
                                        self._queue.put(progress)
                                                   │
   main thread: _poll_queue() every 100ms ◀────────┘
       → update widgets  (only the main thread touches Tk)
```

The same queue is also how background analysis tools (radius sweep, CIF
variance) hand their finished results back to the main thread — look for the
`"viz_data"` branch in `_poll_queue`.

**Rule of thumb:** anything on a worker thread communicates with the UI *only*
through the queue. Never call a widget method from a worker thread.

---

## 7. Data & output layout

```
data/
  input/            user-supplied source files (one per subfolder)
    cosmic/         COSMIC Mutant Census TSV
    ptmd/           PTMD disease-associated PTMs TSV
    1433_interactors/  bundled 14-3-3 interactor file (shipped, not user-supplied)
  steps/            intermediate step-1 output (mode-specific TSVs)
cif_models/         downloaded AlphaFold CIF + PAE, one folder per UniProt ID
Output/
  *_db.tsv          final results
  logs/             skipped-protein / skipped-PTM logs with reasons
  coordinates/      export_ca_coordinates output (incl. ChimeraX files)
  cif_variance/     cif_variance output
```

Output TSVs are written UTF-16 (chosen so Excel opens the `Å` character
correctly on double-click). Note this if adding any tool that re-reads them —
pass `encoding="utf-16"` or they'll appear corrupt.

---

## 8. How to make common changes

**Add a new annotation to step 4:** `4_annotate.py` currently does five
annotations in one file. Add yours as a new function following the existing
ones, wire it into the annotation pass, add its output column(s), and add a
test file under `tests/` mirroring `test_4_annotate_*.py`. (If you have room to
improve things: this file is a good candidate to split into an `annotators/`
subpackage, one module per source.)

**Add a new pipeline step:** create `scripts/N_name.py`, add it to the step
list in `pipeline_utils.py`, and teach BOTH `main.py` and
`pipeline_runner.py` to invoke it. Add a test file.

**Add a new tab:** create `ui/newtab_tab.py` with a `NewTabMixin`, add it to the
`App` base-class list in `app.py`, add a `self._build_newtab_tab()` call in
`_build_ui`, and initialize any new `self._` attributes in `App.__init__`.

**Add a standalone analysis tool:** model it on `radius_sweep.py` — import from
`pipeline_utils`, take argparse args, write to `Output/`. Surface it in the
Analysis Tools tab if it should be GUI-accessible.

---

## 9. Testing

Tests live in `tests/`, one file per pipeline step / tool, run with `pytest`
(configured in `pyproject.toml`). The suite is large and is the best
documentation of intended behavior — when changing a step, run its test file
first to see what contract you're expected to preserve, and add cases for new
behavior. There is no GUI test layer; the tests cover the pipeline and
standalone scripts.

---

## 10. Conventions & gotchas

- **Two shared-utility modules**, mirroring the two layers: `pipeline_utils.py`
  (backend) and `ui/common.py` (frontend). Don't import GUI code into the
  pipeline.
- **UTF-16 output** — see §7.
- **Mode-specific step-1 files** — see §3.3; don't collapse them back to one.
- **Worker threads never touch Tk** — see §6.
- **Mixin attribute names share one namespace** — see §4.1.
- **The monkeypatches are load-bearing** — see §4.2.
- Skipped proteins/PTMs are logged with reasons under `Output/logs/` rather
  than silently dropped — preserve that when adding steps that can skip inputs.
