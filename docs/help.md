# Mutation Cluster Proximity Pipeline — Help

## Getting Started

### Input Files

The pipeline requires two input data files. Use the **Browse** buttons on the Pipeline tab to select each one — the app copies the file into the correct folder automatically, validates that it actually has the expected columns (not just the right file extension), and never overwrites your existing file until the new one is confirmed valid:

| Input | Description | Source |
|-------|-------------|--------|
| **COSMIC** | Mutant Census TSV (~600 MB) | [COSMIC](https://cancer.sanger.ac.uk/cosmic) |
| **PTMD** | Disease-associated PTMs TSV | [PTMD 2.0](https://ptmd.biocuckoo.cn/download.php) |

Each input folder should contain exactly **one** file. Browsing a new file replaces the previous one.

A third input — the **14-3-3 confirmed interactors** spreadsheet — is bundled with the app rather than something you provide. It's small and rarely updated, so it isn't part of the Input Files section.

Before a run starts, the app also checks that your chosen output files aren't locked (e.g. open in Excel) and that your COSMIC/PTMD files pass the same content checks as the Browse dialog — so problems are caught upfront rather than partway through a run.

---

## Pipeline Modes

### PTM Proximity (default)

Finds recurrent cancer mutations that cluster in 3D space near disease-associated post-translational modification (PTM) sites. Runs all 4 steps.

### Mutation Clustering

Finds recurrent cancer mutations that cluster together in 3D space, with no PTM requirement. Runs all 4 steps; step 4 annotates with PolyPhen-2, AIUPred, and InterPro only (14-3-3 and kinase predictions are PTM-site-specific and don't apply here).

### Single Protein

Analyze a single protein by selecting its CIF structure file. The UniProt ID is auto-detected from the file, or you can enter one directly if it isn't. Results can be appended to an existing output database.

Accepts the same **Cutoff**, **Min pLDDT**, and **Max PAE** settings as the main pipeline. **Min samples** is also available, but can only tighten the hotspot threshold already applied when the input TSV was built — mutations below that original threshold aren't in the data to filter in the first place. If you are planning on appending this data to the output database, it is recommended you do not change these variables so they are consistent with the rest of the database.

### CA Coordinates

Export alpha-carbon coordinates for every residue of one or more proteins — each added by **gene symbol** or **UniProt accession** to a list (**+ Add** button or Enter), then exported as a batch with the same options applied to every one. One protein failing (no AlphaFold model, unresolvable gene, etc.) is logged and skipped rather than stopping the rest of the batch.

Each protein's output goes in its own `Output/coordinates/{UniProt}/` folder: `all_ca.tsv` (every residue's coordinates, plus a nearby-patient-count column) and `mutation_ca.tsv` (coordinates only at COSMIC mutation positions). If the AlphaFold structure for a protein is only a fragment (very large proteins are split by AlphaFold into multiple fragments), the app warns you upfront rather than silently analyzing an incomplete structure, and the ChimeraX outputs below are skipped for it.

**Heatmaps** — for single-fragment proteins, also produces ChimeraX scripts you can open directly to reproduce the view with no manual steps:

- **Mutation heatmap** (on by default) — colors the structure by COSMIC patient count near each residue, a red heatmap. Two sub-options apply only to this heatmap: **Log-scale** (compresses heavily skewed patient counts so lower-count regions stay visible instead of being crushed toward one flat color) and **Dim low-confidence residues** (fades each residue in proportion to how low its AlphaFold confidence is, so a hotspot in a poorly-modeled region reads as less certain than an equally hot one in a well-modeled region).
- **pLDDT heatmap** — colors the structure by AlphaFold's own per-residue confidence score, using the same color scheme AlphaFold DB itself uses.

Both heatmaps can be on at once — they're written as separate scripts, since ChimeraX can only show one coloring at a time on a single open structure.

**Markers** — layered on top of whichever heatmap(s) are generated (or a plain, uncolored structure if neither is), without overwriting the heatmap's own color at that residue:

- **Mark PTM sites** — a small green sphere at each known PTM site's position. Requires PTM Proximity mode's Step 1 to have been run at least once, since that's where PTM position data comes from.
- **Show mutation markers** — an orange stick at each COSMIC mutation position.

---

## Pipeline Steps

### Step 1: Filter and Merge Data

Merges PTMD disease-associated PTM sites with COSMIC recurrent mutations. A mutation must appear in at least a minimum number of distinct samples (3 by default — configurable via the **Min samples** field) with confirmed somatic status to be included.

### Step 2: Download Structures

Fetches AlphaFold CIF structure models and PAE (Predicted Aligned Error) files for each protein. Downloaded files are cached in `cif_models/` and reused on subsequent runs; the pipeline checks AlphaFold DB for a newer model version each time and re-downloads only what's changed.

### Step 3: Find Nearby Mutations

Computes 3D distances between PTM sites and mutation hotspots using the AlphaFold structures. The default distance cutoff is **10 Ångströms** (configurable via the **Cutoff** field), with optional **Min pLDDT** and **Max PAE** filters to exclude low-confidence structural regions and low-confidence residue pairs.

### Step 4: Annotate Results

In PTM Proximity mode, adds five types of annotations to each PTM site:

- **14-3-3 binding predictions** — Queries the 14-3-3-Pred API and cross-references experimentally confirmed interactors (Ser/Thr sites only)
- **PolyPhen-2 scores** — Queries myvariant.info for pathogenicity predictions on each mutation
- **Kinase predictions** — Uses the Kinase Library to predict the top 5 upstream kinases for each phosphorylation site (Ser/Thr/Tyr sites only)
- **AIUPred disorder predictions** — Predicts intrinsic disorder and disordered-binding-region propensity, both for the PTM residue and for each nearby mutation's residue
- **InterPro functional domains** — Queries the InterPro REST API for curated domain/family/site entries on each protein, then reports which entry (if any) contains the PTM site's or mutation's specific residue position

Mutation Clustering mode also runs this step, applying the PolyPhen-2, AIUPred, and InterPro annotations above to the anchor mutation and its nearby mutations — 14-3-3 and kinase predictions are skipped, since both require a curated PTM site that this mode has no concept of.

---

## Understanding the Output

### File Encoding

All output TSVs (`ptm_mutation_proximity_db.tsv`, `mutation_cluster_db.tsv`, and their `_long.tsv` companions) are UTF-16 encoded. This is deliberate: it's the encoding Excel reliably auto-detects as tab-delimited when a file is opened directly (double-clicked), unlike plain UTF-8 with or without a BOM. If you're reading these files programmatically instead (pandas, R, etc.), specify the encoding explicitly, or you'll get a decode error or garbled columns: `pd.read_csv(path, sep="\t", encoding="utf-16")` in pandas, or `read.delim(path, fileEncoding="UTF-16")` in R.

### Results Tab

A **PTM Proximity / Mutation Clusters** toggle at the top switches which mode's results are shown. Each shows two linked tables: PTM Proximity mode shows **PTM Sites** (one row per PTM site) and **Mutation Details** (one row per nearby mutation, for whichever PTM site is selected above); Mutation Clustering mode shows the equivalent **Anchor Mutations** and **Nearby Mutations** tables. All four tables have far more columns than are shown by default; click **Columns** on any table to show/hide columns, and hover the **?** badge next to any column name for an explanation of exactly what it means and how it's computed. That in-app reference is the authoritative, up-to-date column list — it isn't duplicated here since column definitions change more often than this document does.

The **Export** button on the Mutation Details table writes exactly what's currently shown (respecting any active search/filter and sort order) to a TSV in the output folder, using the same UTF-16 encoding as the pipeline's own output files.

### Mutation Tags

Mutations shown in the PTM Sites table's raw mutation-list columns include inline tags:

- **(isoform?)** — The reference amino acid in COSMIC doesn't match the AlphaFold structure at this position, possibly due to isoform differences
- **(PP:D,0.999)** — PolyPhen-2 prediction: **D** = Probably Damaging, **P** = Possibly Damaging, **B** = Benign. The number is the confidence score (0-1)
- **(PAE:2.1)** — AlphaFold's Predicted Aligned Error for the residue pair. Lower = higher structural confidence

### Kinase Predictions

Format: `KINASE(log2_score, percentile%)`

- **Log2 score** — Raw motif match strength (higher = better match, can be negative)
- **Percentile** — How the score ranks against a background phosphoproteome (e.g. 95% means better than 95% of all known phosphosites for that kinase)

Only phosphorylation sites (Ser/Thr/Tyr) receive kinase predictions. Other PTM types will have a blank kinase predictions column.

### AIUPred Disorder Scores

A 0-1 score; above 0.5 is treated as "yes" for the corresponding Disordered?/Binding? column. **General** disorder is intrinsic disorder propensity; **binding** disorder specifically flags regions predicted to be disordered in isolation but become ordered upon binding a partner protein.

### Visualization Tab

Selecting a PTM site (or anchor mutation, in Mutation Clustering mode) in the Results tab and clicking **Visualize** draws a lollipop (needle) plot of its nearby mutations here, colored by PolyPhen-2 classification — or pick one directly via the tab's own search/dropdown. A **Data** toggle switches between PTM Proximity and Mutation Clusters results, independent of whichever mode the Results tab is currently showing.

- **View: Single PTM / Whole protein** — Single PTM shows one site's lollipop plot; Whole protein stacks every PTM site/anchor mutation on that protein into one scrollable view, each with its own domain map.
- **Show: All mutations / Unique per position** — All mutations lists every substitution individually; Unique per position collapses same-residue substitutions into one merged lollipop (colored by the most severe PolyPhen class among them).
- Each plotted cluster shows its total and unique mutation counts on the left for quick reference, and a domain map (from InterPro) above the lollipop marks functional domains and the PTM/anchor position.
- **Save PNG** exports the current plot to the output folder.

---

## Analysis Tools

A separate tab for two standalone structural-analysis tools, independent of the main pipeline modes above. A toggle at the top switches between them — each keeps its own settings and plot when you switch away and back. Set parameters, click **Run**, then **Save PNG** to export the resulting plot to the output folder. **Show Details** reveals a log panel with the tool's console output (elbow points, warnings, etc.).

### Radius Sweep

Tests a range of distance cutoffs (not just one fixed value) for one or more genes/proteins, to see how the set of nearby mutations changes as the cutoff changes. At each radius it computes the average number of nearby mutations per PTM site, compares that against a random-placement baseline (mutations shuffled across the same protein), and marks the "elbow" of each curve — the point of diminishing returns.

- **Genes / UniProt IDs** — add one or more proteins, by gene symbol or UniProt accession (**+ Add** button or Enter). Each entry is validated upfront: it must already have hotspot mutation data (from PTM Proximity or Mutation Clustering mode's Step 1) and a downloaded AlphaFold structure (Step 2), or it's rejected with an explanation instead of failing later at Run. A protein spanning multiple AlphaFold fragments is accepted but flagged — only fragment 1 is analyzed, so results may be incomplete for it.
- **Radius range (Å)** — start, stop, and step size for the sweep (default 4–20, step 1). The proximity search re-runs at every radius in this range.
- **Min samples** — minimum distinct COSMIC samples for a mutation to count as a hotspot (default 3). This is computed live from the raw COSMIC file and is independent of whatever **Min samples** value the main pipeline's Step 1 used to build its intermediate file — it can only be tightened or loosened within what that file already contains.
- **Include unfiltered COSMIC comparison** — also sweeps every COSMIC missense mutation for the same genes, not just the recurrent hotspot ones, so hotspot-filtered results can be compared against the full mutation set at each radius. Produces a 4-panel plot instead of 2.

Output is written to `Output/radius_sweep.png` (matching the in-app plot) and a companion `Output/radius_sweep.tsv` with the raw per-radius data; elbow points and the average optimal radius are also printed to the run log.

### CIF Variance

Compares multiple AlphaFold CIF files for the *same* protein — e.g. different model versions, seeds, or predicted fragments — by aligning the structures to an iteratively-refined average reference and computing per-residue positional variance and pLDDT agreement, to see which regions of a prediction are most consistent or most uncertain across models.

- **Input folder** — directory containing two or more `.cif` files for the same protein (default `data/cif_comparison/`), set via **Browse**. A live count below the field confirms how many CIFs were found there.
- **Top N residues** — number of most structurally-variable residues to report (default 10).
- **Report range** — restrict the reported output to a residue range (e.g. start `50`, end `630`). Leave both blank to report every residue.
- **Align range** — use only this residue range for structural alignment, rather than the whole protein; defaults to the Report range if left blank. Useful for excluding disordered regions from alignment while still reporting their variance.
- **UniProt override** — UniProt accession to use for PTM/mutation cross-referencing; auto-detected from the CIF file if left blank.
- **Gene override** — gene symbol used to look up the UniProt ID from the pipeline's intermediate data, only consulted if the UniProt override above is left blank.

Output (written to `Output/cif_variance/`):

- **`variance_plot.png`** — per-residue positional variance and pLDDT (mean ± std) across structures, with PTM and mutation sites marked; the same plot shown in-app
- **`variance_data.tsv`** — per-residue variance, pLDDT stats, and PTM/mutation flags
- **`pairwise_rmsd.tsv`** — RMSD matrix between every pair of input structures

**Generate AlphaFold Seeds JSON** — writes a batch job file for [AlphaFold Server](https://alphafoldserver.com) requesting 10 separate predictions of a protein's sequence, one per seed (seeds 1 through 10), using the same UniProt/gene resolution as the fields above. Upload it at alphafoldserver.com, then drop the resulting CIFs into the input folder above to compare their variance with this tool.

---

## Controls

### Stop / Resume / Cancel

- **Stop** — Freezes the pipeline immediately (mid-step). The process is suspended, not killed.
- **Resume** — Continues exactly where it was frozen. No data is lost or re-processed.
- **Cancel** — Kills the pipeline and restores the previous output file from a backup.

### Output Folder

Use the **Change** button to select a custom output directory. Click **Reset** to return to the default `Output/` folder.

---

## Caching

The pipeline caches data to speed up subsequent runs:

| Cache | Location | Purpose |
|-------|----------|---------|
| Gene → UniProt mappings | `data/cache/uniprot_gene_mapping.tsv`, `data/cache/gene_to_uniprot_mapping.tsv` | Gene symbol / UniProt accession lookups |
| Isoform safe lengths | `data/cache/isoform_safe_lengths.tsv` | Detects when COSMIC's numbering diverges from the canonical AlphaFold sequence |
| CIF structures | `cif_models/` | AlphaFold structure and PAE files |
| 14-3-3 predictions | `data/cache/1433pred/` | Per-protein API responses |
| PolyPhen-2 scores | `data/cache/polyphen.tsv` | Per-mutation pathogenicity |
| Kinase predictions | `data/cache/kinase_predictions.tsv` | Per-sequence-window kinase scores |
| AIUPred disorder | `data/cache/aiupred_disorder.tsv` | Per-residue disorder/binding scores |
| InterPro domains | `data/cache/interpro_domains.tsv` | Per-protein functional domain/family/site entries |

All caches are automatically populated on first run and reused on subsequent runs. Use **Manage Cache** on the Pipeline tab to clear individual caches (or all of them) and force fresh lookups.
