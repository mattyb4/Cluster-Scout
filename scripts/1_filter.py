import argparse
import re
import requests
import pandas as pd
from pathlib import Path
from tqdm import tqdm

SCRIPT_DIR = Path(__file__).parent
PROJECT_ROOT = SCRIPT_DIR.parent
HOTSPOT_MIN_AFFECTED_CASES = 3


def clean_str_list(values):
    cleaned = []
    seen = set()

    for value in values.dropna():
        text = str(value).strip()
        if text and text not in seen:
            seen.add(text)
            cleaned.append(text)

    return "; ".join(cleaned)


def extract_gene_from_protein_change(protein_change):
    if pd.isna(protein_change):
        return None

    text = str(protein_change).strip()
    parts = text.split()

    return parts[0] if parts else None


def extract_aa_change(protein_change):
    if pd.isna(protein_change):
        return None

    text = str(protein_change).strip()
    parts = text.split()

    if len(parts) < 2:
        return None

    return parts[1]


def is_simple_substitution(change):
    if pd.isna(change):
        return False

    text = str(change).strip()
    return bool(re.fullmatch(r"[A-Z]\d+[A-Z]", text))


def build_ptm_site(row):
    residue = str(row["Residue"]).strip() if pd.notna(row["Residue"]) else ""

    position = ""
    if pd.notna(row["Position"]):
        try:
            position = str(int(float(row["Position"])))
        except Exception:
            position = str(row["Position"]).strip()

    site = f"{residue}{position}".strip()
    ptm_type = str(row["Type"]).strip() if pd.notna(row["Type"]) else ""

    return f"{site}:{ptm_type}" if ptm_type else site


def format_mutation_with_count(row):
    return f'{row["mutation"]} ({int(row["affected_cases"])})'


def find_case_count_column(df):
    candidates = [
        "num_cohort_ssm_affected_cases",
        "cohort_ssm_affected_cases",
        "affected_cases",
        "cases",
        "count"
    ]

    for col in candidates:
        if col in df.columns:
            return col

    raise ValueError(
        "Could not find a case-count column in the TCGA file. "
        "Expected one of: num_cohort_ssm_affected_cases, "
        "cohort_ssm_affected_cases, affected_cases, cases, count"
    )


def fetch_uniprot_gene_mapping(uniprot_ids, batch_size=100):
    """Fetch UniProt accession -> primary gene symbol via the UniProt REST API."""
    # Strip variant suffixes (e.g. Q16613_VAR_A129T -> Q16613) — AlphaFold models canonical sequences
    ids = list({uid.split("_")[0] for uid in set(uniprot_ids)})
    if not ids:
        return pd.DataFrame(columns=["UniProt", "gene"])

    rows = []
    uniprot_release = None
    total_batches = (len(ids) + batch_size - 1) // batch_size
    for i in tqdm(range(0, len(ids), batch_size), desc="Fetching UniProt gene names", total=total_batches):
        batch = ids[i : i + batch_size]
        query = " OR ".join(f"accession:{uid}" for uid in batch)
        url = "https://rest.uniprot.org/uniprotkb/search"
        params = {"query": query, "fields": "accession,gene_names", "format": "tsv", "size": batch_size}

        while url:
            resp = requests.get(url, params=params)
            resp.raise_for_status()
            if uniprot_release is None:
                uniprot_release = resp.headers.get("X-UniProt-Release")
            lines = resp.text.strip().split("\n")
            for line in lines[1:]:
                parts = line.split("\t")
                if len(parts) >= 2:
                    accession = parts[0].strip()
                    primary_gene = parts[1].strip().split()[0] if parts[1].strip() else None
                    if primary_gene:
                        rows.append({"UniProt": accession, "gene": primary_gene})

            link_header = resp.headers.get("Link", "")
            match = re.search(r'<([^>]+)>; rel="next"', link_header)
            url = match.group(1) if match else None
            params = None

    if uniprot_release:
        print(f"Using UniProt release: {uniprot_release}")
    return pd.DataFrame(rows).drop_duplicates(subset=["UniProt"])


def fetch_gene_to_uniprot_mapping(gene_names, batch_size=20):
    """Fetch primary gene symbol -> reviewed human UniProt accession via the UniProt REST API."""
    genes = list(set(gene_names))
    gene_set = set(genes)
    rows = []
    total_batches = (len(genes) + batch_size - 1) // batch_size
    for i in tqdm(range(0, len(genes), batch_size), desc="Fetching UniProt IDs for genes", total=total_batches):
        batch = genes[i : i + batch_size]
        gene_query = " OR ".join(f"gene_exact:{g}" for g in batch)
        query = f"({gene_query}) AND organism_id:9606 AND reviewed:true"
        url = "https://rest.uniprot.org/uniprotkb/search"
        params = {
            "query": query,
            "fields": "accession,gene_names",
            "format": "tsv",
            "size": min(batch_size * 3, 500),
        }

        while url:
            resp = requests.get(url, params=params)
            resp.raise_for_status()
            lines = resp.text.strip().split("\n")
            for line in lines[1:]:
                parts = line.split("\t")
                if len(parts) >= 2:
                    accession = parts[0].strip()
                    gene_field = parts[1].strip()
                    primary_gene = gene_field.split()[0] if gene_field else None
                    if primary_gene and primary_gene in gene_set:
                        rows.append({"gene": primary_gene, "UniProt": accession})

            link_header = resp.headers.get("Link", "")
            match = re.search(r'<([^>]+)>; rel="next"', link_header)
            url = match.group(1) if match else None
            params = None

    if not rows:
        return pd.DataFrame(columns=["gene", "UniProt"])
    return pd.DataFrame(rows).drop_duplicates(subset=["gene"])


def _load_and_filter_tcga(tcga_file):
    """Shared TCGA loading and filtering logic used by both pipeline modes."""
    tcga = pd.read_csv(tcga_file, sep="\t")

    tcga["gene"] = tcga["protein_change"].apply(extract_gene_from_protein_change)
    tcga["aa_change"] = tcga["protein_change"].apply(extract_aa_change)

    tcga = tcga[tcga["gene"].notna()].copy()
    tcga = tcga[tcga["aa_change"].apply(is_simple_substitution)].copy()

    case_count_col = find_case_count_column(tcga)
    tcga["affected_cases"] = pd.to_numeric(tcga[case_count_col], errors="coerce")
    tcga = tcga[tcga["affected_cases"].notna()].copy()
    tcga = tcga[tcga["affected_cases"] >= HOTSPOT_MIN_AFFECTED_CASES].copy()

    tcga["mutation"] = tcga["aa_change"]
    tcga["mutation_with_count"] = tcga.apply(format_mutation_with_count, axis=1)

    return tcga, case_count_col


def _run_ptm_proximity_filter(output_file):
    ptmd_file = PROJECT_ROOT / "data" / "PTMD_disease_associated_ptms.tsv"
    tcga_file = PROJECT_ROOT / "data" / "TCGA_frequent_mutations.tsv"

    print("Loading PTMD and TCGA files...")
    ptmd = pd.read_csv(ptmd_file, sep="\t", low_memory=False)
    tcga, case_count_col = _load_and_filter_tcga(tcga_file)

    # -----------------------
    # Filter PTMD disruptions
    # -----------------------
    ptmd = ptmd[ptmd["State"] == "N"].copy()

    # Normalize variant UniProt IDs to canonical accession (e.g. Q16613_VAR_A129T -> Q16613)
    ptmd["UniProt"] = ptmd["UniProt"].str.split("_").str[0]

    # -----------------------
    # Map UniProt -> gene via UniProt REST API
    # -----------------------
    uniprot_ids = ptmd["UniProt"].dropna().unique().tolist()
    print(f"Mapping {len(uniprot_ids)} UniProt IDs to gene names via UniProt API...")
    idmap = fetch_uniprot_gene_mapping(uniprot_ids)

    ptmd = ptmd.merge(idmap, on="UniProt", how="left")

    if "Gene name" in ptmd.columns:
        ptmd["gene"] = ptmd["gene"].fillna(ptmd["Gene name"])

    ptmd = ptmd[ptmd["gene"].notna()].copy()

    # -----------------------
    # Build PTM site and PTM-disease pair
    # -----------------------
    ptmd["ptm_site"] = ptmd.apply(build_ptm_site, axis=1)
    ptmd["ptm_disease_pair"] = ptmd["ptm_site"] + " | " + ptmd["Disease"].astype(str)

    print("Filtering TCGA mutations and aggregating by gene...")
    # -----------------------
    # Aggregate PTMs by gene
    # -----------------------
    ptmd_grouped = (
        ptmd.groupby("gene", as_index=False)
        .agg(
            uniprot_id=("UniProt", "first"),
            ptms_on_protein=("ptm_site", clean_str_list),
            ptm_disease_pairs=("ptm_disease_pair", clean_str_list),
        )
    )

    # -----------------------
    # Aggregate hotspot mutations by gene
    # -----------------------
    tcga_grouped = (
        tcga.groupby("gene", as_index=False)
        .agg(
            mutations_on_protein=("mutation_with_count", clean_str_list),
        )
    )

    # -----------------------
    # Merge datasets
    # -----------------------
    merged = ptmd_grouped.merge(tcga_grouped, on="gene", how="left")

    merged = merged[
        [
            "uniprot_id",
            "gene",
            "ptms_on_protein",
            "mutations_on_protein",
            "ptm_disease_pairs",
        ]
    ]

    merged = merged[merged["mutations_on_protein"].notna()].copy()

    # -----------------------
    # Save result
    # -----------------------
    print("Saving output...")
    merged.to_csv(output_file, sep="\t", index=False)

    print("Done.")
    print(f"Using case count column: {case_count_col}")
    print(f"Hotspot minimum affected cases: {HOTSPOT_MIN_AFFECTED_CASES}")
    print(f"PTMD disruption genes: {ptmd['gene'].nunique()}")
    print(f"TCGA hotspot genes: {tcga['gene'].nunique()}")
    print(f"Final merged proteins: {len(merged)}")
    print(f"Output saved to: {output_file}")


def _run_mutation_clustering_filter(output_file):
    tcga_file = PROJECT_ROOT / "data" / "TCGA_frequent_mutations.tsv"

    print("Loading TCGA file...")
    tcga, case_count_col = _load_and_filter_tcga(tcga_file)

    print("Aggregating hotspot mutations by gene...")
    tcga_grouped = (
        tcga.groupby("gene", as_index=False)
        .agg(mutations_on_protein=("mutation_with_count", clean_str_list))
    )

    gene_names = tcga_grouped["gene"].tolist()
    print(f"Mapping {len(gene_names)} genes to UniProt IDs via UniProt API...")
    gene_map = fetch_gene_to_uniprot_mapping(gene_names)

    result = tcga_grouped.merge(gene_map, on="gene", how="left")
    unmapped = result["UniProt"].isna().sum()
    result = result[result["UniProt"].notna()].copy()
    result = result.rename(columns={"UniProt": "uniprot_id"})
    result = result[["uniprot_id", "gene", "mutations_on_protein"]]

    print("Saving output...")
    result.to_csv(output_file, sep="\t", index=False)

    print("Done.")
    print(f"Using case count column: {case_count_col}")
    print(f"Hotspot minimum affected cases: {HOTSPOT_MIN_AFFECTED_CASES}")
    print(f"TCGA hotspot genes: {tcga['gene'].nunique()}")
    print(f"Genes mapped to UniProt: {len(result)}")
    print(f"Genes not mapped to UniProt (excluded): {unmapped}")
    print(f"Output saved to: {output_file}")


def main():
    parser = argparse.ArgumentParser(description="Filter and prepare input data for the pipeline.")
    parser.add_argument(
        "--mode",
        choices=["ptm-proximity", "mutation-clustering"],
        default="ptm-proximity",
        help=(
            "'ptm-proximity' merges PTMD + TCGA and keeps only genes with both PTMs and mutations. "
            "'mutation-clustering' keeps all recurrent TCGA hotspot mutations regardless of PTMs."
        ),
    )
    args = parser.parse_args()

    output_file = PROJECT_ROOT / "data" / "steps" / "PTMD_TCGA_hotspots_by_protein.tsv"
    output_file.parent.mkdir(parents=True, exist_ok=True)

    if args.mode == "mutation-clustering":
        _run_mutation_clustering_filter(output_file)
    else:
        _run_ptm_proximity_filter(output_file)


if __name__ == "__main__":
    main()
