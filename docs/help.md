# Mutation Cluster Proximity Pipeline — Help

## Getting Started

### Input Files

The pipeline requires three input data files. Use the **Browse** buttons on the main screen to select each one:

| Input | Description | Source |
|-------|-------------|--------|
| **COSMIC** | Mutant Census TSV (~600 MB) | [COSMIC](https://cancer.sanger.ac.uk/cosmic) |
| **PTMD** | Disease-associated PTMs TSV | [PTMD 2.0](https://ptmd.biocuckoo.cn/download.php) |
| **14-3-3 Interactors** | Confirmed interactors Excel (optional) | Provided with this tool |

Each input folder should contain exactly **one** file. Uploading a new file replaces the previous one.

---

## Pipeline Modes

### PTM Proximity (default)

Finds recurrent cancer mutations that cluster in 3D space near disease-associated post-translational modification (PTM) sites. Runs all 4 steps.

### Mutation Clustering

Finds recurrent cancer mutations that cluster together in 3D space, with no PTM requirement. Runs steps 1-3 only.

### Single Protein

Analyze a single protein by selecting its CIF structure file. The UniProt ID is auto-detected from the file. Results can be appended to an existing output database.

---

## Pipeline Steps

### Step 1: Filter and Merge Data

Merges PTMD disease-associated PTM sites with COSMIC recurrent mutations. A mutation must appear in at least **3 distinct samples** (with confirmed somatic status) to be included.

### Step 2: Download Structures

Fetches AlphaFold CIF structure models and PAE (Predicted Aligned Error) files for each protein. Downloaded files are cached in `cif_models/` and reused on subsequent runs.

### Step 3: Find Nearby Mutations

Computes 3D distances between PTM sites and mutation hotspots using the AlphaFold structures. The default distance cutoff is **10 Angstroms**.

### Step 4: Annotate Results

Adds three types of annotations to each PTM site:

- **14-3-3 binding predictions** — Queries the 14-3-3-Pred API and cross-references experimentally confirmed interactors
- **PolyPhen-2 scores** — Queries myvariant.info for pathogenicity predictions on each mutation
- **Kinase predictions** — Uses the Kinase Library to predict the top 5 upstream kinases for each phosphorylation site

---

## Understanding the Output

### Output Columns

| Column | Description |
|--------|-------------|
| **UniProt** | UniProt accession ID |
| **gene** | Gene symbol |
| **ptm_site** | PTM position (e.g. S337) |
| **ptm_type** | Type of modification (e.g. Phosphorylation) |
| **mutations_within_5_positions** | Mutations within 5 residues of the PTM site, with 3D distances |
| **mutations_more_than_5_positions** | Mutations beyond 5 residues but within 10 Angstroms in 3D space |
| **mutation_at_ptm_site** | Whether the PTM site itself is a mutation hotspot |
| **ptm_diseases** | Disease associations from PTMD |
| **total_cosmic_missense_patients** | Total patients with any missense mutation in this gene |
| **1433pred_binding_site** | "Yes" if predicted 14-3-3 binding site (Ser/Thr only) |
| **1433pred_consensus** | Raw 14-3-3-Pred consensus score |
| **1433_confirmed_site** | "Yes" if experimentally confirmed 14-3-3 binding site |
| **kinase_predictions** | Top 5 predicted kinases with scores |

### Mutation Tags

Mutations in the output include inline tags:

- **(isoform?)** — The reference amino acid in COSMIC doesn't match the AlphaFold structure at this position, possibly due to isoform differences
- **(PP:D,0.999)** — PolyPhen-2 prediction: **D** = Probably Damaging, **P** = Possibly Damaging, **B** = Benign. The number is the confidence score (0-1)
- **(PAE:2.1)** — AlphaFold's Predicted Aligned Error for the residue pair. Lower = higher structural confidence

### Kinase Predictions

Format: `KINASE(log2_score, percentile%)`

- **Log2 score** — Raw motif match strength (higher = better match, can be negative)
- **Percentile** — How the score ranks against a background phosphoproteome (e.g. 95% means better than 95% of all known phosphosites for that kinase)

Only phosphorylation sites (Ser/Thr/Tyr) receive kinase predictions. Other PTM types will have a blank kinase_predictions column.

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
| UniProt gene mappings | `data/cache/uniprot_gene_mapping.tsv` | Gene symbol lookups |
| CIF structures | `cif_models/` | AlphaFold structure files |
| 14-3-3 predictions | `data/cache/1433pred/` | Per-protein API responses |
| PolyPhen-2 scores | `data/cache/polyphen.tsv` | Per-mutation pathogenicity |
| Kinase predictions | `data/cache/kinase_predictions.tsv` | Per-sequence-window kinase scores |

All caches are automatically populated on first run and reused on subsequent runs. Delete a cache file to force a fresh lookup.
