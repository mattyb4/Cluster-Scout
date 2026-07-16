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

Finds recurrent cancer mutations that cluster together in 3D space, with no PTM requirement. Runs steps 1-3 only.

### Single Protein

Analyze a single protein by selecting its CIF structure file. The UniProt ID is auto-detected from the file, or you can enter one directly if it isn't. Results can be appended to an existing output database.

Accepts the same **Cutoff**, **Min pLDDT**, and **Max PAE** settings as the main pipeline. **Min samples** is also available, but can only tighten the hotspot threshold already applied when the input TSV was built — mutations below that original threshold aren't in the data to filter in the first place.

### CA Coordinates

Export alpha-carbon coordinates for every residue of one protein — by UniProt ID or gene symbol — along with its COSMIC missense mutation positions, for use in external visualization tools.

Genes/proteins can be entered as either a **gene symbol** or a **UniProt accession** wherever the app asks for one. If the AlphaFold structure for a protein is only a fragment (very large proteins are split by AlphaFold into multiple fragments), the app warns you upfront rather than silently analyzing an incomplete structure.

---

## Pipeline Steps

### Step 1: Filter and Merge Data

Merges PTMD disease-associated PTM sites with COSMIC recurrent mutations. A mutation must appear in at least a minimum number of distinct samples (3 by default — configurable via the **Min samples** field) with confirmed somatic status to be included.

### Step 2: Download Structures

Fetches AlphaFold CIF structure models and PAE (Predicted Aligned Error) files for each protein. Downloaded files are cached in `cif_models/` and reused on subsequent runs; the pipeline checks AlphaFold DB for a newer model version each time and re-downloads only what's changed.

### Step 3: Find Nearby Mutations

Computes 3D distances between PTM sites and mutation hotspots using the AlphaFold structures. The default distance cutoff is **10 Ångströms** (configurable via the **Cutoff** field), with optional **Min pLDDT** and **Max PAE** filters to exclude low-confidence structural regions and low-confidence residue pairs.

### Step 4: Annotate Results

Adds four types of annotations to each PTM site:

- **14-3-3 binding predictions** — Queries the 14-3-3-Pred API and cross-references experimentally confirmed interactors (Ser/Thr sites only)
- **PolyPhen-2 scores** — Queries myvariant.info for pathogenicity predictions on each mutation
- **Kinase predictions** — Uses the Kinase Library to predict the top 5 upstream kinases for each phosphorylation site (Ser/Thr/Tyr sites only)
- **AIUPred disorder predictions** — Predicts intrinsic disorder and disordered-binding-region propensity, both for the PTM residue and for each nearby mutation's residue

---

## Understanding the Output

### Results Tab

The Results tab shows two linked tables — **PTM Sites** (one row per PTM site) and **Mutation Details** (one row per nearby mutation, for whichever PTM site is selected above). Both tables have far more columns than are shown by default; click **Columns** on either table to show/hide columns, and hover the **?** badge next to any column name for an explanation of exactly what it means and how it's computed. That in-app reference is the authoritative, up-to-date column list — it isn't duplicated here since column definitions change more often than this document does.

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

Selecting a PTM site in the Results tab and clicking **Visualize** draws a lollipop (needle) plot of its nearby mutations on the Visualization tab, colored by PolyPhen-2 classification.

---

## Analysis Tools

A separate tab for two standalone structural-analysis tools, independent of the main pipeline modes above:

### Radius Sweep

Tests a range of distance cutoffs (not just one fixed value) for one or more genes/proteins, to see how the set of nearby mutations changes as the cutoff changes. Optionally also compares against every COSMIC mutation for those genes, not just recurrent hotspots.

### CIF Variance

Compares multiple AlphaFold CIF files for the *same* protein — e.g. different model versions or predicted fragments — by aligning the structures and computing per-residue positional variance, to see which regions of a prediction are most consistent or most uncertain across models.

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

All caches are automatically populated on first run and reused on subsequent runs. Use **Manage Cache** on the Pipeline tab to clear individual caches (or all of them) and force fresh lookups.
