import argparse
import ast
import re
import requests
import pandas as pd
from pathlib import Path
from tqdm import tqdm

SCRIPT_DIR = Path(__file__).parent
PROJECT_ROOT = SCRIPT_DIR.parent
HOTSPOT_MIN_AFFECTED_CASES = 3

UNMATCHED_GENES_LOG = PROJECT_ROOT / "Output" / "logs" / "ptm_genes_without_cosmic_mutations.tsv"


def clean_str_list(values):
    cleaned = []
    seen = set()

    for value in values.dropna():
        text = str(value).strip()
        if text and text not in seen:
            seen.add(text)
            cleaned.append(text)

    return "; ".join(cleaned)


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


def parse_mutation_site(val):
    """Extract mutations from a PTMD MutationSite cell (e.g. \"['D120N', 'E127D']\")."""
    if pd.isna(val):
        return []
    text = str(val).strip()
    if text in ("", "[]", "nan"):
        return []
    try:
        result = ast.literal_eval(text)
        if isinstance(result, list):
            return [str(m).strip() for m in result if str(m).strip()]
        return [str(result).strip()]
    except (ValueError, SyntaxError):
        return re.findall(r"[A-Z]\d+[A-Z*]", text)


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


COSMIC_SOMATIC_STATUSES = {
    "Confirmed somatic variant",
    "Reported in another cancer sample as somatic",
}


def _load_and_filter_cosmic(cosmic_file):
    """Shared COSMIC Mutant Census loading and filtering logic used by both pipeline modes.

    The Mutant Census is one row per (mutation, sample) occurrence rather than
    pre-aggregated hotspots, so affected-case counts are computed here by counting
    distinct samples per (gene, amino-acid change).
    """
    cols = ["GENE_SYMBOL", "MUTATION_AA", "COSMIC_SAMPLE_ID", "MUTATION_SOMATIC_STATUS"]
    cosmic = pd.read_csv(cosmic_file, sep="\t", usecols=cols, low_memory=False)

    cosmic = cosmic[cosmic["MUTATION_SOMATIC_STATUS"].isin(COSMIC_SOMATIC_STATUSES)].copy()

    cosmic["aa_change"] = cosmic["MUTATION_AA"].str.replace(r"^p\.", "", regex=True)
    cosmic = cosmic[cosmic["aa_change"].apply(is_simple_substitution)].copy()

    cosmic = (
        cosmic.groupby(["GENE_SYMBOL", "aa_change"])["COSMIC_SAMPLE_ID"]
        .nunique()
        .reset_index(name="affected_cases")
        .rename(columns={"GENE_SYMBOL": "gene"})
    )
    cosmic = cosmic[cosmic["affected_cases"] >= HOTSPOT_MIN_AFFECTED_CASES].copy()

    cosmic["mutation"] = cosmic["aa_change"]
    cosmic["mutation_with_count"] = cosmic.apply(format_mutation_with_count, axis=1)

    return cosmic


def _run_ptm_proximity_filter(output_file):
    ptmd_file = PROJECT_ROOT / "data" / "PTMD_disease_associated_ptms.tsv"
    cosmic_file = PROJECT_ROOT / "data" / "Cosmic_MutantCensus_v104_GRCh38.tsv"

    print("Loading PTMD and COSMIC files...")
    ptmd = pd.read_csv(ptmd_file, sep="\t", low_memory=False)
    cosmic = _load_and_filter_cosmic(cosmic_file)

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

    # -----------------------
    # Build known disrupting mutations per PTM site
    # MutationSite records the specific mutations documented in literature as disrupting each PTM.
    # -----------------------
    ptmd["parsed_mutations"] = ptmd["MutationSite"].apply(parse_mutation_site)
    ptmd_with_muts = ptmd[ptmd["parsed_mutations"].map(len) > 0].copy()

    if not ptmd_with_muts.empty:
        ptmd_exploded = ptmd_with_muts.explode("parsed_mutations")
        ptmd_exploded = ptmd_exploded[ptmd_exploded["parsed_mutations"].notna()]
        ptmd_exploded = ptmd_exploded[ptmd_exploded["parsed_mutations"] != ""]
        disruption_map = (
            ptmd_exploded
            .groupby(["UniProt", "ptm_site"])["parsed_mutations"]
            .apply(lambda x: ",".join(sorted(set(x.dropna()))))
            .reset_index()
        )
        disruption_map["disruption_entry"] = (
            disruption_map["ptm_site"] + ">" + disruption_map["parsed_mutations"]
        )
        disruptions_grouped = (
            disruption_map
            .groupby("UniProt")["disruption_entry"]
            .apply(lambda x: "; ".join(x))
            .reset_index()
            .rename(columns={"UniProt": "uniprot_id", "disruption_entry": "ptm_known_disruptions"})
        )
    else:
        disruptions_grouped = pd.DataFrame(columns=["uniprot_id", "ptm_known_disruptions"])

    print("Filtering COSMIC mutations and aggregating by gene...")
    # -----------------------
    # Aggregate PTMs by UniProt so each isoform is kept separate.
    # Grouping by gene previously collapsed all isoforms under one UniProt ID,
    # causing PTM positions from one isoform to be checked against another's structure.
    # -----------------------
    ptmd_grouped = (
        ptmd.groupby("UniProt", as_index=False)
        .agg(
            gene=("gene", "first"),
            ptms_on_protein=("ptm_site", clean_str_list),
            ptm_disease_pairs=("ptm_disease_pair", clean_str_list),
        )
        .rename(columns={"UniProt": "uniprot_id"})
    )
    ptmd_grouped = ptmd_grouped.merge(disruptions_grouped, on="uniprot_id", how="left")

    # -----------------------
    # Aggregate hotspot mutations by gene
    # -----------------------
    cosmic_grouped = (
        cosmic.groupby("gene", as_index=False)
        .agg(
            mutations_on_protein=("mutation_with_count", clean_str_list),
        )
    )

    # -----------------------
    # Merge datasets
    # -----------------------
    merged = ptmd_grouped.merge(cosmic_grouped, on="gene", how="left")

    merged = merged[
        [
            "uniprot_id",
            "gene",
            "ptms_on_protein",
            "mutations_on_protein",
            "ptm_disease_pairs",
            "ptm_known_disruptions",
        ]
    ]

    # -----------------------
    # Log PTMD proteins with no matching COSMIC hotspot mutations for their gene
    # (gene-name mismatch between PTMD/UniProt and COSMIC, or no mutation met
    # HOTSPOT_MIN_AFFECTED_CASES for that gene) before dropping them.
    # -----------------------
    unmatched = merged[merged["mutations_on_protein"].isna()].copy()
    UNMATCHED_GENES_LOG.parent.mkdir(parents=True, exist_ok=True)
    unmatched[["uniprot_id", "gene", "ptms_on_protein"]].to_csv(
        UNMATCHED_GENES_LOG, sep="\t", index=False, encoding="utf-16"
    )
    print(f"Wrote {len(unmatched)} proteins with no matching COSMIC hotspot mutations to: {UNMATCHED_GENES_LOG}")

    merged = merged[merged["mutations_on_protein"].notna()].copy()

    # -----------------------
    # Save result
    # -----------------------
    print("Saving output...")
    merged.to_csv(output_file, sep="\t", index=False)

    print("Done.")
    print(f"Hotspot minimum affected cases: {HOTSPOT_MIN_AFFECTED_CASES}")
    print(f"PTMD disruption genes: {ptmd['gene'].nunique()}")
    print(f"COSMIC hotspot genes: {cosmic['gene'].nunique()}")
    print(f"Final merged proteins: {len(merged)}")
    print(f"Output saved to: {output_file}")


def _run_mutation_clustering_filter(output_file):
    cosmic_file = PROJECT_ROOT / "data" / "Cosmic_MutantCensus_v104_GRCh38.tsv"

    print("Loading COSMIC file...")
    cosmic = _load_and_filter_cosmic(cosmic_file)

    print("Aggregating hotspot mutations by gene...")
    cosmic_grouped = (
        cosmic.groupby("gene", as_index=False)
        .agg(mutations_on_protein=("mutation_with_count", clean_str_list))
    )

    gene_names = cosmic_grouped["gene"].tolist()
    print(f"Mapping {len(gene_names)} genes to UniProt IDs via UniProt API...")
    gene_map = fetch_gene_to_uniprot_mapping(gene_names)

    result = cosmic_grouped.merge(gene_map, on="gene", how="left")
    unmapped = result["UniProt"].isna().sum()
    result = result[result["UniProt"].notna()].copy()
    result = result.rename(columns={"UniProt": "uniprot_id"})
    result = result[["uniprot_id", "gene", "mutations_on_protein"]]

    print("Saving output...")
    result.to_csv(output_file, sep="\t", index=False)

    print("Done.")
    print(f"Hotspot minimum affected cases: {HOTSPOT_MIN_AFFECTED_CASES}")
    print(f"COSMIC hotspot genes: {cosmic['gene'].nunique()}")
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
            "'ptm-proximity' merges PTMD + COSMIC and keeps only genes with both PTMs and mutations. "
            "'mutation-clustering' keeps all recurrent COSMIC hotspot mutations regardless of PTMs."
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
