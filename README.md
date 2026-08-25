# Cluster-Scout

Created by Matt Banks

Based on BYU capstone project created by Matt Banks, Jaden Searle, Tyler Plauche, and Alissa Moulder - 
https://github.com/mattyb4/Bio465Capstone

Owned by Josh Andersen Lab at the Huntsman Cancer Institute

## Introduction

Cluster-Scout began as a 2026 Senior Bioinformatics Capstone project at Brigham Young University, in collaboration with the Huntsman Cancer Institute. That capstone has since concluded, but the project has grown beyond it into a standalone desktop application for finding recurrent cancer mutations that cluster near or directly disrupt post-translational modification (PTM) sites in 3D protein structure.

The app wraps the full analysis pipeline (data filtering, AlphaFold structure lookup, 3D distance calculation, and annotation) along with tools to browse and visualize the results, so no command-line usage is required for day-to-day use. A CLI is still available underneath for scripting or headless runs.

---

## First: Downloading data

The pipeline requires three input data files. Each goes in its own folder under `data/input/`:

| Folder | File | Source |
|---|---|---|
| `data/input/cosmic/` | COSMIC Mutant Census TSV (600+ MB) | [COSMIC](https://cancer.sanger.ac.uk/cosmic) |
| `data/input/ptmd/` | PTMD disease-associated PTMs TSV | [PTMD 2.0](https://ptmd.biocuckoo.cn/download.php) |
| `data/input/1433_interactors/` | 14-3-3 confirmed interactors Excel | Provided in this repository |

**Using the desktop app:** Click the **Browse** button next to each input file to select it. The app copies it into the correct folder automatically.

**Manual setup:** Download each file and place it in the corresponding folder above. Each folder should contain exactly one file — the pipeline will error if it finds multiple files or none.

## Getting Started

### Clone the Repository

First, clone this repository to your local machine and navigate into the project directory:

```bash
git clone https://github.com/mattyb4/Cluster-Scout.git
cd Cluster-Scout
```

### Requirements

Install [uv](https://docs.astral.sh/uv/getting-started/installation/) (handles Python and all dependencies automatically):

**macOS/Linux:**
```bash
curl -LsSf https://astral.sh/uv/install.sh | sh
```

**Windows (PowerShell):**
```powershell
powershell -ExecutionPolicy ByPass -c "irm https://astral.sh/uv/install.ps1 | iex"
```

*If you get the error "uv: command not found", see troubleshooting steps below.*

### Launch the desktop app (recommended)

```bash
uv run app.py
```

This opens Cluster-Scout, a five-tab desktop application:

- **Pipeline** — select input files, choose a mode (PTM Proximity, Mutation Clustering, Single Protein, or CA Coordinates), configure settings (distance cutoff, minimum samples, PolyPhen filters, etc.), and run the analysis with live progress, pause/resume, and cache management.
- **Results** — browse every PTM site found, with sortable columns (mutation counts, unique mutated positions, patient counts, 14-3-3/PolyPhen/disruption flags). Selecting a PTM shows its individual nearby mutations in a detail table below. A **PTM Proximity / Mutation Clusters** toggle switches between that output and the Mutation Clustering mode's own Anchor/Nearby Mutations tables. The Mutation Details table can be exported to its own TSV.
- **Visualization** — generate a lollipop (needle) plot for any PTM site and its nearby mutations directly from the Results tab (via the **📈 Visualize** button or double-clicking a row) or by picking one from the search/dropdown on the tab itself. Mutations are colored by PolyPhen-2 classification and sized by patient count; a broken axis splits mutations within the local sequence window from ones that are spatially close but sequence-distant. Above the lollipop, a domain map draws that protein's InterPro functional domains/families/sites (color-coded by type, in lanes to keep overlapping entries readable) with an arrow marking the PTM site or anchor mutation's own position. A **View: Single PTM / Whole protein** toggle switches between one site's plot and a scrollable stack of every PTM site/anchor mutation on that protein, each with its own domain map. A **Show: All mutations / Unique per position** toggle switches between listing every substitution individually or collapsing same-residue substitutions into one merged lollipop, and each plotted cluster shows its total and unique mutation counts for quick reference. Plots can be exported as PNG.
- **Analysis Tools** — two standalone structural-analysis tools, independent of the main pipeline modes: **Radius Sweep** (test a range of distance cutoffs for one or more genes/proteins to see how the nearby-mutation set changes) and **CIF Variance** (compare multiple AlphaFold CIF predictions of the same protein for structural confidence, and generate an AlphaFold Server batch JSON to request new seeded predictions to compare). Both are also available as standalone CLI scripts — see below.
- **Help / Documentation** — an in-app copy of this project's usage docs (`docs/help.md`).

### Run the pipeline via the command line

The same pipeline can be run headlessly, which is useful for scripting or servers without a display:

The pipeline has two modes:

**PTM-proximity mode** (default) — finds recurrent cancer mutations that cluster in 3D space near disease-associated PTM sites. Runs all four steps and outputs `Output/ptm_mutation_proximity_db.tsv`.

```bash
uv run main.py
```
or equivalently:
```bash
uv run main.py --mode ptm-proximity
```

**Mutation-clustering mode** — finds recurrent cancer mutations that cluster together in 3D space, with no PTM requirement. Runs all four steps (step 4's annotations are PolyPhen-2, AIUPred, and InterPro only — 14-3-3 and kinase predictions are PTM-site-specific and don't apply here) and outputs `Output/mutation_cluster_db.tsv`.

```bash
uv run main.py --mode mutation-clustering
```

**PTM-proximity steps:**

1. **Filter** — merges and filters the PTMD and COSMIC datasets. A mutation must show up in a minimum of 3 distinct samples (and have a confirmed/reported-somatic status) for it to be added to the filtered dataset. This threshold can be changed by editing HOTSPOT_MIN_AFFECTED_CASES near the top of scripts/1_filter.py, or via the "Min samples" setting in the app.
2. **Download structures** — fetches AlphaFold CIF models and PAE files for each protein (will be a little over 2gb) by iterating over all UniProt IDs found in the tsv file generated by step 1. The AlphaFold DB does not seem to have .cif files for proteins multiple-thousand residues long, so some may not be found and will also skip a protein if a canonical sequence is not obvious in the DB. These situations will be logged in Output/logs/download_errors.tsv. You can manually upload these sequences to alphafoldserver.com and generate .cif files there, to then be individually analyzed by analyze_single_cif_nearby_mutations.py. See "Analyzing individual .cif models" below for instructions on how to do that.
3. **Find nearby mutations** — computes 3D distances between PTM sites and nearby cancer mutations. If a PTM residue from the input does not match up with the residue found in the .cif file, it will be skipped and logged in Output/logs/ptm_skipped.tsv
4. **Annotate results** — annotates each PTM site and nearby mutation with 14-3-3 binding predictions (14-3-3-Pred API plus experimentally confirmed interactors), PolyPhen-2 pathogenicity scores (myvariant.info), predicted upstream kinases (Kinase Library), AIUPred intrinsic disorder / binding-region scores, and InterPro functional domains (curated domain/family/site entries covering the PTM site's or mutation's specific residue position, from the InterPro REST API).

**Mutation-clustering steps:**

1. **Filter** — filters the COSMIC dataset for recurrent hotspot mutations and maps gene names to UniProt IDs
2. **Download structures** — same as above (previously downloaded files are reused automatically)
3. **Find mutation clusters** — computes pairwise 3D distances between all recurrent mutations on each protein; outputs mutations that cluster within 10 Å of at least one other mutation
4. **Annotate results** — annotates each anchor mutation and its nearby mutations with PolyPhen-2 pathogenicity scores, AIUPred intrinsic disorder / binding-region scores, and InterPro functional domains

Each mode's Step 1 writes its own intermediate file under `data/steps/PTMD_COSMIC_hotspots_by_protein.tsv` for PTM-proximity mode, `COSMIC_hotspots_by_protein.tsv` for mutation-clustering mode, so running one mode's Step 1 never overwrites the other's data. Tools that are inherently PTM-based (Radius Sweep, `analyze_single_cif_nearby_mutations.py`, the CA Coordinates "Mark PTM sites" option) always read the PTM-proximity file specifically, regardless of which mode you ran most recently.

The main output for PTM-proximity mode is **`Output/ptm_mutation_proximity_db.tsv`** — a table of PTM sites, their nearby COSMIC mutations, 3D distances, 14-3-3 binding predictions, kinase predictions, disorder scores, PolyPhen-2 pathogenicity scores, and InterPro functional domains.

PTM-proximity mode also always produces **`Output/ptm_mutation_proximity_long.tsv`**, with one row per PTM/mutation pair instead of one row per PTM site. This is what powers the per-mutation detail table on the Results tab and the patient-count-aware Visualization plots.

The main output for mutation-clustering mode is **`Output/mutation_cluster_db.tsv`** — a table of recurrent mutations and other mutations clustering within 10 Å of them in 3D space. It also always produces **`Output/mutation_cluster_long.tsv`**, the per-mutation-pair equivalent of `ptm_mutation_proximity_long.tsv` above.

**A note on file encoding:** all of these output TSVs are UTF-16 encoded, chosen deliberately because it's the encoding Excel reliably auto-detects as tab-delimited when a file is opened directly (double-clicked) rather than imported through Data > From Text/CSV — plain UTF-8, with or without a BOM, does not behave the same way in that workflow. If you're reading these files programmatically instead of opening them in Excel, specify the encoding explicitly or you'll get a decode error or garbled columns: `pd.read_csv(path, sep="\t", encoding="utf-16")` in pandas, or `read.delim(path, fileEncoding="UTF-16")` in R.

---
## Interpreting the Data

### Output Database
The main output of this pipeline is ptm_mutation_proximity_db.tsv, found in the Output folder. This tsv file has the following columns:  

**UniProt** - the UniProt ID  
**gene** - the gene the protein is associated with  
**ptm_site** - position within protein sequence where PTM is  
**ptm_type** - the type of PTM  
**mutations_within_5_positions** - list of all mutation hotspots within 5 residues of PTM site. The formatting is initial amino acid, location, AA it mutates to, then optional tags and distance. Tags include `(isoform?)` if the reference residue doesn't match the AlphaFold model, `(PP:D,0.999)` for PolyPhen-2 predictions (D=Damaging, P=Possibly Damaging, B=Benign with score), and the PAE score* in parentheses.  
**mutation_count_within_5_positions** - sum of total mutation hotspots in previous column  
**unique_mutation_position_count_within_5_positions** - count of distinct residue positions represented in mutations_within_5_positions (multiple substitutions at the same residue count once)  
**nearby_muts_total_patient_count** - total distinct patients across all mutations in mutations_within_5_positions  
**mutations_more_than_5_positions** - list of all mutation hotsposts further than 5 residues of PTM site  
**mutation_count_more_than_5_positions** - sum of total mutation hostspots in previous column  
**unique_mutation_position_count_more_than_5_positions** - count of distinct residue positions represented in mutations_more_than_5_positions  
**distant_muts_total_patient_count** - total distinct patients across all mutations in mutations_more_than_5_positions  
**morethan5_linear_distance** - list of distances on linear amino acid sequence for all mutation hotspots in mutations_more_than_5_positions. This allows for easily seeing entries with mutations that are far on the linear sequence but fold close to PTM site in 3D space  
**mutation_at_ptm_site** - indicates if the PTM site itself is a mutation hotspot  
**confirmed_disrupting_mutations** - mutations experimentally shown to disrupt this PTM (from PTMD)  
**ptm_diseases** - lists diseases PTM is associated with according to PTMD 2.0  
**total_cosmic_missense_patients** - total distinct patients with any missense mutation in this gene across COSMIC  
**1433pred_binding_site** - "Yes" if the 14-3-3-Pred consensus score > 0, "No" if ≤ 0, blank for non-Ser/Thr sites  
**1433pred_consensus** - raw 14-3-3-Pred consensus score  
**1433_confirmed_site** - "Yes" if the site appears in the experimentally confirmed 14-3-3 interactors dataset  
**1433_confirmed_pmid** - PubMed ID of the paper that confirmed the 14-3-3 binding site  
**kinase_predictions** - top predicted upstream kinases for the PTM site (phosphorylation sites only), formatted as `KINASE(log2_score, percentile%)`  
**ptm_aiupred_general** - AIUPred general intrinsic disorder score (0-1) at the PTM residue  
**ptm_aiupred_binding** - AIUPred binding-region disorder score (0-1) at the PTM residue  
**ptm_is_disordered** - "yes"/"no", thresholded from ptm_aiupred_general at > 0.5  
**ptm_is_binding** - "yes"/"no", thresholded from ptm_aiupred_binding at > 0.5  
**ptm_domain** - InterPro functional domain/family/site entry (or entries, semicolon-separated) whose residue range contains the PTM site, formatted as `name (type, start-end)`. Curated entries only (InterPro's own cross-database consensus, not every individual member-database hit), so nested/overlapping calls (e.g. a domain within a broader superfamily) can both appear. Blank if the position falls in no annotated entry. The `_long.tsv` companion additionally has a per-mutation **mutation_domain** column, the same lookup for each nearby mutation's own position.  



*Predicted Alignment Error (PAE) score is how confident AlphaFold is that those residues are at that position. Lower score = higher confidence

## Error logging
The pipeline also generates logs found in Output/logs to record any issues where the pipeline was unable to download a file for a certain protein from AlphaFold or unable to run calculations for a PTM and why. For more information, see skipped_ptm_summary.md in Output/logs 

## Analyzing Individual .cif Models

The app's **Single Protein** mode (on the Pipeline tab) runs this same analysis through the GUI: browse to a `.cif` file, and the UniProt ID is auto-detected from it.

To do the same from the command line: if you would like to manually generate the .cif for a skipped protein from AlphaFold and run analysis on it, create a folder in cif_models that is named the exact Uniprot ID for the protein, then put your .cif file in it. Run the following command:

```bash
uv run scripts/analyze_single_cif_nearby_mutations.py <uniprotID goes here>/<.cif file name goes here>
```
Example:
```bash
uv run scripts/analyze_single_cif_nearby_mutations.py P35222/AF-P35222-F1-model_v6.cif
```

By default it prints nearby mutations to the terminal. To also append this new data to the proximity database, add --append-to-db to the end like this:

```bash
uv run scripts/analyze_single_cif_nearby_mutations.py P35222/AF-P35222-F1-model_v6.cif --append-to-db
```
If you would like to output it to a new tsv file instead, run it like this:

```bash
uv run scripts/analyze_single_cif_nearby_mutations.py P35222/AF-P35222-F1-model_v6.cif --append-to-db --output-db Output/outputfilename.tsv
```

Keep in mind that this will still be running analyses based on the input data from PTMD_COSMIC_hotspots_by_protein.tsv generated during the pipeline.

Whichever protein you are running analysis on, in order for it to work, the UniProt ID in "PTMD_COSMIC_hotspots_by_protein.tsv needs to match the name of the folder the .cif file is put in (within the cif_models directory) exactly.

## Exporting Alpha-Carbon Coordinates

The **CA Coordinates** mode (Pipeline tab) exports the 3D coordinates of alpha-carbon atoms for one or more proteins at once, along with ready-to-open ChimeraX visualization scripts. It's also available as a standalone script for CLI/scripting use.

**Batch proteins:** add gene symbols and/or UniProt accessions to a list (each auto-detected, `+ Add` button or Enter), then run — every protein in the list is exported in turn with the same options applied to all of them. One protein failing (no AlphaFold model, unresolvable gene, etc.) is logged and skipped rather than stopping the batch. COSMIC is scanned once and reused for the whole batch rather than re-read per protein.

If a protein's CIF file hasn't been downloaded yet, it's fetched automatically from the AlphaFold DB. Each protein's outputs go in their own `Output/coordinates/{UniProt}/` folder:

- **`all_ca.tsv`** — x/y/z coordinates for every residue, plus a `patients_within_10A` column (total COSMIC patient count summed across all missense mutations within 10 Å of that residue)
- **`mutation_ca.tsv`** — coordinates only at positions with confirmed somatic missense mutations in COSMIC, plus the mutation labels and patient counts

**ChimeraX heatmaps and markers** (single-fragment proteins only — skipped, with a warning, for proteins AlphaFold split into multiple fragments):

| Option | Output file(s) | Effect |
|---|---|---|
| Mutation heatmap (on by default) | `mutations.defattr`, `mutations_view.cxc` | Colors the cartoon by `patients_within_10A`, sequential red palette |
| &nbsp;&nbsp;↳ Log-scale | — | Colors by `log1p(patients_within_10A)` instead of the raw count, for heavily right-skewed data |
| &nbsp;&nbsp;↳ Dim low-confidence residues | — | Fades each residue in proportion to how low its pLDDT is (`100 - pLDDT` percent transparent); switches that script's lighting from ChimeraX's `soft` preset to `simple`, since soft's ambient shadows render incorrectly once part of a model is transparent |
| pLDDT heatmap | `plddt_view.cxc` | Colors the cartoon by AlphaFold's own per-residue confidence, using ChimeraX's built-in `alphafold` palette |
| Mark PTM sites | — | Marks each known PTM site with a small green sphere at its CA coordinate (needs PTM Proximity mode's Step 1 to have been run, for the PTM position data) |
| Show mutation markers | — | Shows each COSMIC mutation position's side chain as an orange stick |

The two heatmaps are independent and can both be on at once — since ChimeraX's `color` command replaces rather than layers, they're written as **separate** `.cxc` scripts rather than combined into one. Markers are layered into whichever script(s) get generated (or their own plain-cartoon script, `markers_view.cxc`, if no heatmap is selected) as independent geometry/atom-display commands, not a recoloring, so they never overwrite a heatmap's own color at that residue. Open any `.cxc` file directly in ChimeraX to reproduce that view with no manual steps.

**CLI usage** (one or more positional tokens, each a gene symbol or UniProt accession):
```bash
uv run scripts/export_ca_coordinates.py P04637
uv run scripts/export_ca_coordinates.py TP53 EGFR P04637
```

**Options** (mirroring the GUI checkboxes above): `--no-mutation-heatmap`, `--plddt-heatmap`, `--log-scale`, `--dim-low-confidence`, `--mark-ptm-sites`, `--mark-mutations`, `--cosmic path/to/COSMIC.tsv`.

## Finding the Optimal Distance Cutoff

Also available from the **Analysis Tools** tab (a gene/UniProt chip list, the same parameters as below, and the resulting plot shown in-app). `scripts/radius_sweep.py` tests a range of distance cutoffs (default 4-20 Å) to help choose the PTM-to-mutation distance threshold used by the pipeline. For a set of genes, it computes the average number of nearby mutations per PTM site at each radius, compares against a random-placement baseline (mutations shuffled across the same protein), and detects the "elbow" of each curve — the point of diminishing returns — using `kneed`.

**Basic usage** (uses a default gene panel: EGFR, TP53, VHL, CANT1, DDR2, PTPN11, LZTR1, CDK12):
```bash
uv run scripts/radius_sweep.py
```

**Choose specific genes:**
```bash
uv run scripts/radius_sweep.py --genes EGFR TP53 VHL
```

**Choose a custom radius range** (start, stop, step in Å):
```bash
uv run scripts/radius_sweep.py --radii 4 25 1
```

**Choose a custom hotspot threshold** (minimum distinct COSMIC samples for a mutation to count as a hotspot, default 3). This is computed live from the raw COSMIC file and is independent of whatever "Min samples" value the main pipeline's step 1 used:
```bash
uv run scripts/radius_sweep.py --min-samples 5
```

**Compare against unfiltered COSMIC mutations** (not just hotspot-filtered ones) with `--unfiltered`. This requires pipeline step 1 to have already been run, since PTM site positions are read from the intermediate `PTMD_COSMIC_hotspots_by_protein.tsv` file, and requires CIF files for the target genes to already be downloaded.

Output is written to `Output/radius_sweep.png` (a 2- or 4-panel plot depending on whether `--unfiltered` was used) and a matching `Output/radius_sweep.tsv` with the raw per-radius data. Elbow points and average optimal radius are also printed to the console.

## Comparing Structural Variance Across CIF Models

Also available from the **Analysis Tools** tab. `scripts/cif_variance.py` compares multiple AlphaFold CIF predictions of the same protein (e.g., different seeds or model versions) to assess structural confidence. It aligns the structures to an iteratively-refined average reference, then reports per-residue positional variance and pLDDT agreement, cross-referenced against PTM and mutation sites from the pipeline's intermediate data.

Place two or more `.cif` files for the same protein in `data/cif_comparison/`, then run:
```bash
uv run scripts/cif_variance.py
```

**Options:**
- `--input-dir` — directory containing the CIF files to compare (default: `data/cif_comparison`)
- `--output-dir` — where results are written (default: `Output/cif_variance`)
- `--top N` — number of top-variance residues to print (default: 10)
- `--range START END` — restrict reported output to a residue range
- `--align-range START END` — use only this residue range for structural alignment (defaults to `--range`); useful for excluding disordered regions from alignment while still reporting their variance
- `--uniprot` / `--gene` — UniProt ID or gene symbol for PTM/mutation cross-referencing (auto-detected from the CIF file when possible)

Output (in `Output/cif_variance/`):
- **`variance_plot.png`** — per-residue positional variance and pLDDT (mean ± std) across structures, with PTM and mutation sites marked
- **`variance_data.tsv`** — per-residue variance, pLDDT stats, and PTM/mutation flags
- **`pairwise_rmsd.tsv`** — RMSD matrix between every pair of input structures

**Generating more seeds to compare:** both the CLI (`--generate-seed-json`) and the Analysis Tools tab (a "Generate AlphaFold Seeds JSON" button) can write a batch JSON file for [AlphaFold Server](https://alphafoldserver.com) requesting 10 separate jobs — one per seed, seeds 1 through 10 — for a protein's canonical sequence, resolved the same way (`--uniprot`/`--gene`/CIF metadata) as the comparison above. Upload it at alphafoldserver.com, then drop the resulting CIFs into `data/cif_comparison/` to compare them here.
```bash
uv run scripts/cif_variance.py --generate-seed-json --uniprot P04637
```

---

### Troubleshooting: `uv: command not found`

**macOS/Linux:** After installing, your shell session needs to reload its PATH. Run:

```bash
source "$HOME/.local/bin/env"
```

Then open a new terminal and `uv` should work. If you use conda, ensure `~/.local/bin` is on your PATH by adding this to your shell profile (e.g. `~/.zshrc` or `~/.bash_profile`) and restarting your terminal:

```bash
export PATH="$HOME/.local/bin:$PATH"
```

**Windows:** After installing, close and reopen PowerShell. If `uv` is still not found, add it to your PATH manually:
1. Search for **"Edit the system environment variables"** in the Start menu.
2. Under **User variables**, select `Path` and click **Edit**.
3. Add `%USERPROFILE%\.local\bin`.
4. Click OK and reopen your terminal.

---

## Notes

**UniProt gene mapping** is fetched live from the UniProt REST API (Step 1). The release version used is printed to the console during Step 1 (`Using UniProt release: ...`).

**AlphaFold structures** are downloaded from AlphaFold DB. The model version is encoded in each downloaded filename (e.g., `AF-P12345-F1-model_v6.cif`). AlphaFold DB v6 covers the full human proteome.

**Input data files** (`PTMD_disease_associated_ptms.tsv`, `Cosmic_MutantCensus_v104_GRCh38.tsv`) are static files downloaded from PTMD and COSMIC.

**Kinase predictions** are generated locally using the Kinase Library package and only computed for phosphorylation sites (Ser/Thr/Tyr).

**AIUPred disorder scores** are computed locally using AIUPred, once per protein, for both general intrinsic disorder and binding-region disorder.

**`ptm_diseases` is pan-cancer:** The `ptm_diseases` column in the output reflects which diseases the PTM site is associated with in PTMD. The nearby COSMIC mutations are pan-cancer and were not filtered by cancer type, so a nearby mutation appearing in the output does not imply it co-occurs in the same cancer type as the PTM disease association.

**404 / Isoforms Only:** Proteins without available AlphaFold structures or lacking canonical models were excluded from structural analysis.
