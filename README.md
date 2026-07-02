# Cluster-Scout

### Identifying Cancer-Driving Mutations that Target Post-Translational Modifications (PTMs)
Created by Matt Banks

Based on BYU capstone project created by Matt Banks, Jaden Searle, Tyler Plauche, and Alissa Moulder - 
https://github.com/mattyb4/Bio465Capstone

Mentored by Dr. Josh Andersen

## Introduction

Cluster-Scout began as a 2026 Senior Bioinformatics Capstone project at Brigham Young University, in collaboration with the Huntsman Cancer Institute. That capstone has since concluded, but the project has grown beyond it into a standalone desktop application for finding recurrent cancer mutations that cluster near — or directly disrupt — post-translational modification (PTM) sites in 3D protein structure.

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

This opens Cluster-Scout, a four-tab desktop application:

- **Pipeline** — select input files, choose a mode, configure settings (distance cutoff, minimum samples, PolyPhen filters, etc.), and run the analysis with live progress, pause/resume, and cache management.
- **Results** — browse every PTM site found, with sortable columns (mutation counts, unique mutated positions, patient counts, 14-3-3/PolyPhen/disruption flags). Selecting a PTM shows its individual nearby mutations in a detail table below.
- **Visualization** — generate a lollipop (needle) plot for any PTM site and its nearby mutations directly from the Results tab (via the **📈 Visualize** button or double-clicking a row) or by picking one from the search/dropdown on the tab itself. Mutations are colored by PolyPhen-2 classification and sized by patient count; a broken axis splits mutations within the local sequence window from ones that are spatially close but sequence-distant. A **Show: All mutations / Unique per position** toggle switches between listing every substitution individually or collapsing same-residue substitutions into one merged lollipop. Plots can be exported as PNG.
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

**Mutation-clustering mode** — finds recurrent cancer mutations that cluster together in 3D space, with no PTM requirement. Runs steps 1-3 only and outputs `Output/mutation_cluster_db.tsv`.

```bash
uv run main.py --mode mutation-clustering
```

**PTM-proximity steps:**

1. **Filter** — merges and filters the PTMD and COSMIC datasets. A mutation must show up in a minimum of 3 distinct samples (and have a confirmed/reported-somatic status) for it to be added to the filtered dataset. This threshold can be changed by editing HOTSPOT_MIN_AFFECTED_CASES near the top of scripts/1_filter.py, or via the "Min samples" setting in the app.
2. **Download structures** — fetches AlphaFold CIF models and PAE files for each protein (will be a little over 2gb) by iterating over all UniProt IDs found in the tsv file generated by step 1. The AlphaFold DB does not seem to have .cif files for proteins multiple-thousand residues long, so some may not be found and will also skip a protein if a canonical sequence is not obvious in the DB. These situations will be logged in Output/logs/download_errors.tsv. You can manually upload these sequences to alphafoldserver.com and generate .cif files there, to then be individually analyzed by analyze_single_cif_nearby_mutations.py. See "Analyzing individual .cif models" below for instructions on how to do that.
3. **Find nearby mutations** — computes 3D distances between PTM sites and nearby cancer mutations. If a PTM residue from the input does not match up with the residue found in the .cif file, it will be skipped and logged in Output/logs/ptm_skipped.tsv
4. **Annotate results** — annotates each PTM site and nearby mutation with 14-3-3 binding predictions (14-3-3-Pred API plus experimentally confirmed interactors), PolyPhen-2 pathogenicity scores (myvariant.info), predicted upstream kinases (Kinase Library), and AIUPred intrinsic disorder / binding-region scores.

**Mutation-clustering steps:**

1. **Filter** — filters the COSMIC dataset for recurrent hotspot mutations and maps gene names to UniProt IDs
2. **Download structures** — same as above (previously downloaded files are reused automatically)
3. **Find mutation clusters** — computes pairwise 3D distances between all recurrent mutations on each protein; outputs mutations that cluster within 10 Å of at least one other mutation

The main output for PTM-proximity mode is **`Output/ptm_mutation_proximity_db.tsv`** — a table of PTM sites, their nearby COSMIC mutations, 3D distances, 14-3-3 binding predictions, kinase predictions, disorder scores, and PolyPhen-2 pathogenicity scores.

Checking **"Long format output"** on the Pipeline tab additionally produces **`Output/ptm_mutation_proximity_long.tsv`**, with one row per PTM/mutation pair instead of one row per PTM site. This is what powers the per-mutation detail table on the Results tab and the patient-count-aware Visualization plots — without it, those views fall back to parsing the wide-format summary columns, which don't carry per-mutation patient counts.

The main output for mutation-clustering mode is **`Output/mutation_cluster_db.tsv`** — a table of recurrent mutations and other mutations clustering within 10 Å of them in 3D space.

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

Keep in mind that this will still be running analyses based on the input data from PTMD_TCGA_hotspots_by_protein.tsv generated during the pipeline.

Whichever protein you are running analysis on, in order for it to work, the UniProt ID in "PTMD_TCGA_hotspots_by_protein.tsv needs to match the name of the folder the .cif file is put in (within the cif_models directory) exactly.

## Exporting Alpha-Carbon Coordinates

A standalone script is available to export the 3D coordinates of alpha-carbon atoms for any protein. It produces two files:

- **`{UniProt}_all_ca.tsv`** — x/y/z coordinates for every residue in the protein
- **`{UniProt}_mutation_ca.tsv`** — coordinates only at positions with confirmed somatic missense mutations in COSMIC, plus the mutation labels and patient counts

If the CIF file has not been downloaded yet, the script will automatically fetch it from the AlphaFold DB. The gene symbol is looked up from the UniProt API (or the local cache if the pipeline has been run before).

**Basic usage:**
```bash
uv run scripts/export_ca_coordinates.py P04637
```

**Provide the gene symbol directly** (skips the UniProt API lookup):
```bash
uv run scripts/export_ca_coordinates.py P04637 --gene TP53
```

**Use a different COSMIC file:**
```bash
uv run scripts/export_ca_coordinates.py P04637 --cosmic path/to/your/COSMIC.tsv
```

Output files are saved to `Output/coordinates/`.

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
