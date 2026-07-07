import argparse
import ast
import re
import sys
import time
import requests
import pandas as pd
from pathlib import Path
from tqdm import tqdm

sys.path.insert(0, str(Path(__file__).resolve().parent))
from pipeline_utils import (  # noqa: E402
    project_root, COSMIC_SOMATIC_STATUSES,
    input_dir, resolve_input_file, COSMIC_INPUT_DIR, PTMD_INPUT_DIR,
)

PROJECT_ROOT = project_root(__file__)
HOTSPOT_MIN_AFFECTED_CASES = 3

UNMATCHED_GENES_LOG = PROJECT_ROOT / "Output" / "logs" / "ptm_genes_without_cosmic_mutations.tsv"

CACHE_DIR = PROJECT_ROOT / "data" / "cache"

# Gene/UniProt mapping and isoform-mismatch checking are reported to the app as one
# continuous progress bar (rather than two bars that complete and restart) by treating
# them as two equal-weighted phases of a single overall percentage.
_NUM_PHASES = 2


def _emit_progress(phase: int, phase_pct: float, desc: str) -> None:
    """Print overall progress for the app to parse. phase is 0-indexed, phase_pct is 0-100."""
    overall = int((phase * 100 + phase_pct) / _NUM_PHASES)
    print(f"\r##PROGRESS## {overall} {desc}", end="", flush=True)


def _load_cache(filename, columns):
    """Load a TSV cache file into a dict keyed by its first column.

    Values are tuples of the remaining columns as strings ('' if blank).
    Returns {} if the cache file doesn't exist yet.
    """
    path = CACHE_DIR / filename
    if not path.exists():
        return {}
    df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    key_col, *value_cols = columns
    return {row[key_col]: tuple(row[c] for c in value_cols) for _, row in df.iterrows()}


def _save_cache(filename, cache, columns):
    """Write a dict keyed by the first column back to a TSV cache file."""
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    key_col, *value_cols = columns
    rows = [{key_col: key, **dict(zip(value_cols, values))} for key, values in cache.items()]
    pd.DataFrame(rows, columns=columns).to_csv(CACHE_DIR / filename, sep="\t", index=False)


def clean_str_list(values):
    """Deduplicate and join a Series of strings into a semicolon-separated list, preserving order."""
    cleaned = []
    seen = set()

    for value in values.dropna():
        text = str(value).strip()
        if text and text not in seen:
            seen.add(text)
            cleaned.append(text)

    return "; ".join(cleaned)


def is_simple_substitution(change):
    """Return True if the amino-acid change is a single-residue missense substitution (e.g. 'V600E').

    Requires the reference and alternate amino acid to differ so that synonymous
    variants stored by COSMIC as 'p.P100P' (same letter before and after) are
    excluded.  Stop-codon variants ('p.R175*') are also excluded because the
    final character must be a letter.
    """
    if pd.isna(change):
        return False

    text = str(change).strip()
    m = re.fullmatch(r"([A-Z])(\d+)([A-Z])", text)
    return bool(m) and m.group(1) != m.group(3)


def build_ptm_site(row):
    """Build a compact PTM site label like 'S473:Phosphorylation' from a PTMD row's residue, position, and type."""
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
    """Format a mutation and its affected-case count as 'V600E (42)' for display in output columns."""
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


UNIPROT_GENE_CACHE_FILE = "uniprot_gene_mapping.tsv"


def fetch_uniprot_gene_mapping(uniprot_ids, batch_size=100):
    """Fetch UniProt accession -> primary gene symbol via the UniProt REST API.

    Results are cached in data/cache/uniprot_gene_mapping.tsv (including accessions
    with no gene name, recorded as ''), so subsequent runs only query accessions
    that haven't been looked up before.
    """
    # Strip variant suffixes (e.g. Q16613_VAR_A129T -> Q16613) — AlphaFold models canonical sequences
    ids = list({uid.split("_")[0] for uid in set(uniprot_ids)})
    if not ids:
        return pd.DataFrame(columns=["UniProt", "gene"])

    cache = _load_cache(UNIPROT_GENE_CACHE_FILE, ["UniProt", "gene"])
    missing = [uid for uid in ids if uid not in cache]
    print(f"{len(ids) - len(missing)}/{len(ids)} UniProt accessions found in cache; fetching {len(missing)} new...")

    if missing:
        uniprot_release = None
        total_batches = (len(missing) + batch_size - 1) // batch_size
        for batch_num, i in enumerate(
            tqdm(range(0, len(missing), batch_size), desc="Fetching UniProt gene names", total=total_batches), 1,
        ):
            batch = missing[i : i + batch_size]
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
                        primary_gene = parts[1].strip().split()[0] if parts[1].strip() else ""
                        cache[accession] = (primary_gene,)

                link_header = resp.headers.get("Link", "")
                match = re.search(r'<([^>]+)>; rel="next"', link_header)
                url = match.group(1) if match else None
                params = None

            for uid in batch:
                cache.setdefault(uid, ("",))

            _emit_progress(0, batch_num / total_batches * 100,
                           f"Fetching UniProt gene names: batch {batch_num}/{total_batches}")

        if uniprot_release:
            print(f"Using UniProt release: {uniprot_release}")
        _save_cache(UNIPROT_GENE_CACHE_FILE, cache, ["UniProt", "gene"])
    else:
        _emit_progress(0, 100, "UniProt gene names: all cached")

    rows = [{"UniProt": uid, "gene": cache[uid][0]} for uid in ids if cache.get(uid, ("",))[0]]
    return pd.DataFrame(rows, columns=["UniProt", "gene"])


GENE_TO_UNIPROT_CACHE_FILE = "gene_to_uniprot_mapping.tsv"


def fetch_gene_to_uniprot_mapping(gene_names, batch_size=20):
    """Fetch primary gene symbol -> reviewed human UniProt accession via the UniProt REST API.

    Results are cached in data/cache/gene_to_uniprot_mapping.tsv (including genes
    with no reviewed match, recorded as ''), so subsequent runs only query genes
    that haven't been looked up before.
    """
    genes = list(set(gene_names))
    if not genes:
        return pd.DataFrame(columns=["gene", "UniProt"])

    cache = _load_cache(GENE_TO_UNIPROT_CACHE_FILE, ["gene", "UniProt"])
    missing = [g for g in genes if g not in cache]
    print(f"{len(genes) - len(missing)}/{len(genes)} genes found in cache; fetching {len(missing)} new...")

    if missing:
        missing_set = set(missing)
        total_batches = (len(missing) + batch_size - 1) // batch_size
        for batch_num, i in enumerate(
            tqdm(range(0, len(missing), batch_size), desc="Fetching UniProt IDs for genes", total=total_batches), 1,
        ):
            batch = missing[i : i + batch_size]
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
                        if primary_gene and primary_gene in missing_set:
                            cache[primary_gene] = (accession,)

                link_header = resp.headers.get("Link", "")
                match = re.search(r'<([^>]+)>; rel="next"', link_header)
                url = match.group(1) if match else None
                params = None

            for g in batch:
                cache.setdefault(g, ("",))

            _emit_progress(0, batch_num / total_batches * 100,
                           f"Fetching UniProt IDs for genes: batch {batch_num}/{total_batches}")

        _save_cache(GENE_TO_UNIPROT_CACHE_FILE, cache, ["gene", "UniProt"])
    else:
        _emit_progress(0, 100, "UniProt IDs for genes: all cached")

    rows = [{"gene": g, "UniProt": cache[g][0]} for g in genes if cache.get(g, ("",))[0]]
    if not rows:
        return pd.DataFrame(columns=["gene", "UniProt"])
    return pd.DataFrame(rows, columns=["gene", "UniProt"])



def _load_and_filter_cosmic(cosmic_file):
    """Shared COSMIC Mutant Census loading and filtering logic used by both pipeline modes.

    The Mutant Census is one row per (mutation, sample) occurrence rather than
    pre-aggregated hotspots, so affected-case counts are computed here by counting
    distinct samples per (gene, amino-acid change).

    Also returns a gene -> Ensembl transcript accession mapping (one per gene),
    used to detect when COSMIC's mutation numbering follows a different UniProt
    isoform than the canonical AlphaFold-modeled sequence.
    """
    cols = ["GENE_SYMBOL", "MUTATION_AA", "COSMIC_SAMPLE_ID", "MUTATION_SOMATIC_STATUS", "TRANSCRIPT_ACCESSION"]
    cosmic = pd.read_csv(cosmic_file, sep="\t", usecols=cols, low_memory=False)

    cosmic = cosmic[cosmic["MUTATION_SOMATIC_STATUS"].isin(COSMIC_SOMATIC_STATUSES)].copy()

    cosmic["aa_change"] = cosmic["MUTATION_AA"].str.replace(r"^p\.", "", regex=True)
    cosmic = cosmic[cosmic["aa_change"].apply(is_simple_substitution)].copy()

    gene_to_transcript = cosmic.groupby("GENE_SYMBOL")["TRANSCRIPT_ACCESSION"].first().to_dict()

    # Total distinct patients with any missense mutation in each gene, regardless of
    # the hotspot recurrence threshold below. Used as a baseline for comparing
    # nearby/distant mutation patient counts in step 3.
    gene_to_total_missense_patients = (
        cosmic.groupby("GENE_SYMBOL")["COSMIC_SAMPLE_ID"].nunique().to_dict()
    )

    cosmic = (
        cosmic.groupby(["GENE_SYMBOL", "aa_change"])["COSMIC_SAMPLE_ID"]
        .nunique()
        .reset_index(name="affected_cases")
        .rename(columns={"GENE_SYMBOL": "gene"})
    )
    cosmic = cosmic[cosmic["affected_cases"] >= HOTSPOT_MIN_AFFECTED_CASES].copy()

    cosmic["mutation"] = cosmic["aa_change"]
    cosmic["mutation_with_count"] = cosmic.apply(format_mutation_with_count, axis=1)

    return cosmic, gene_to_transcript, gene_to_total_missense_patients


def _fetch_uniprot_sequence(accession):
    """Fetch a UniProt entry's sequence (canonical or a specific isoform) as a plain string."""
    resp = requests.get(f"https://rest.uniprot.org/uniprotkb/{accession}.fasta")
    if resp.status_code != 200:
        return None
    lines = resp.text.strip().split("\n")
    return "".join(lines[1:])


ENSEMBL_XREF_RE = re.compile(r"\.\d+\s*\[([A-Za-z0-9\-]+)\]")


ISOFORM_SAFE_LENGTH_CACHE_FILE = "isoform_safe_lengths.tsv"


def compute_isoform_safe_lengths(gene_to_transcript, gene_to_uniprot, batch_size=10):
    """For genes whose COSMIC transcript is annotated against a UniProt isoform with a
    different sequence than the canonical AlphaFold-modeled accession, compute the
    length of the longest shared prefix between the two sequences. Mutation positions
    beyond this point cannot be reliably mapped onto the canonical structure and should
    be flagged as isoform mismatches in step 3, regardless of whether the residue
    happens to match.

    Results are cached in data/cache/isoform_safe_lengths.tsv, keyed by gene and the
    COSMIC transcript accession used to compute them. A gene is only re-checked if it
    wasn't cached before or its COSMIC transcript accession has changed.

    Returns a DataFrame with columns: gene, isoform_safe_length. Genes not present
    have no restriction (COSMIC numbering matches the canonical sequence).
    """
    genes = [g for g in gene_to_uniprot if g in gene_to_transcript]

    cache = _load_cache(ISOFORM_SAFE_LENGTH_CACHE_FILE, ["gene", "transcript_accession", "isoform_safe_length"])
    missing = [g for g in genes if cache.get(g, (None, None))[0] != gene_to_transcript[g]]
    print(f"{len(genes) - len(missing)}/{len(genes)} genes found in cache; checking {len(missing)} new...")

    if missing:
        total_batches = (len(missing) + batch_size - 1) // batch_size

        for batch_num, i in enumerate(
            tqdm(range(0, len(missing), batch_size), desc="Checking COSMIC transcripts for isoform mismatches", total=total_batches), 1,
        ):
            batch = missing[i : i + batch_size]
            enst_noversion = [gene_to_transcript[g].split(".")[0] for g in batch]
            query = " OR ".join(f"xref:ensembl-{e}" for e in enst_noversion)
            params = {"query": query, "fields": "accession,xref_ensembl,sequence", "format": "tsv", "size": batch_size * 5}
            resp = requests.get("https://rest.uniprot.org/uniprotkb/search", params=params)
            resp.raise_for_status()
            lines = resp.text.strip().split("\n")
            result_rows = [line.split("\t") for line in lines[1:] if line.strip()]

            for gene, enst in zip(batch, enst_noversion):
                canonical_acc = gene_to_uniprot[gene]
                isoform_acc = None
                xref_acc = None
                xref_seq = None
                for parts in result_rows:
                    if len(parts) < 3:
                        continue
                    m = re.search(rf"{re.escape(enst)}{ENSEMBL_XREF_RE.pattern}", parts[1])
                    if m:
                        isoform_acc = m.group(1)
                        xref_acc, xref_seq = parts[0], parts[2]
                        break

                isoform_safe_length = ""
                if isoform_acc is not None and isoform_acc != canonical_acc:
                    canonical_seq = xref_seq if xref_acc == canonical_acc else _fetch_uniprot_sequence(canonical_acc)
                    isoform_seq = _fetch_uniprot_sequence(isoform_acc)
                    if canonical_seq and isoform_seq and canonical_seq != isoform_seq:
                        lcp = 0
                        for a, b in zip(canonical_seq, isoform_seq):
                            if a != b:
                                break
                            lcp += 1
                        isoform_safe_length = str(lcp)

                cache[gene] = (gene_to_transcript[gene], isoform_safe_length)

            _emit_progress(1, batch_num / total_batches * 100,
                           f"Checking isoform mismatches: batch {batch_num}/{total_batches}")

        _save_cache(ISOFORM_SAFE_LENGTH_CACHE_FILE, cache, ["gene", "transcript_accession", "isoform_safe_length"])
    else:
        _emit_progress(1, 100, "Isoform mismatches: all cached")

    rows = [
        {"gene": g, "isoform_safe_length": int(cache[g][1])}
        for g in genes
        if cache.get(g, (None, ""))[1]
    ]
    return pd.DataFrame(rows, columns=["gene", "isoform_safe_length"])


def _run_ptm_proximity_filter(output_file):
    """Run the PTM-proximity pipeline mode: merge PTMD disease PTMs with COSMIC hotspots, keeping only genes with both."""
    ptmd_file = resolve_input_file(input_dir(PROJECT_ROOT, PTMD_INPUT_DIR), (".tsv",))
    cosmic_file = resolve_input_file(input_dir(PROJECT_ROOT, COSMIC_INPUT_DIR), (".tsv",))

    print(f"PTMD file:   {ptmd_file.name}")
    print(f"COSMIC file: {cosmic_file.name}")
    ptmd = pd.read_csv(ptmd_file, sep="\t", low_memory=False)
    cosmic, gene_to_transcript, gene_to_total_missense_patients = _load_and_filter_cosmic(cosmic_file)

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
    t0 = time.time()
    idmap = fetch_uniprot_gene_mapping(uniprot_ids)
    print(f"  UniProt gene mapping completed in {time.time() - t0:.1f}s")

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
    # Flag genes where COSMIC's mutation numbering follows a different UniProt
    # isoform than the canonical AlphaFold-modeled sequence.
    # -----------------------
    print("Checking for COSMIC/canonical isoform mismatches...")
    gene_to_uniprot = dict(zip(merged["gene"], merged["uniprot_id"]))
    t0 = time.time()
    isoform_lengths = compute_isoform_safe_lengths(gene_to_transcript, gene_to_uniprot)
    print(f"  Isoform mismatch check completed in {time.time() - t0:.1f}s")
    merged = merged.merge(isoform_lengths, on="gene", how="left")

    # -----------------------
    # Total COSMIC missense patients per gene, for comparison against the
    # nearby/distant mutation patient counts computed in step 3.
    # -----------------------
    merged["total_cosmic_missense_patients"] = merged["gene"].map(gene_to_total_missense_patients)

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
    """Run the mutation-clustering pipeline mode: keep all recurrent COSMIC hotspots mapped to UniProt, regardless of PTMs."""
    cosmic_file = resolve_input_file(input_dir(PROJECT_ROOT, COSMIC_INPUT_DIR), (".tsv",))
    print(f"COSMIC file: {cosmic_file.name}")

    print("Loading COSMIC file...")
    cosmic, gene_to_transcript, gene_to_total_missense_patients = _load_and_filter_cosmic(cosmic_file)

    print("Aggregating hotspot mutations by gene...")
    cosmic_grouped = (
        cosmic.groupby("gene", as_index=False)
        .agg(mutations_on_protein=("mutation_with_count", clean_str_list))
    )

    gene_names = cosmic_grouped["gene"].tolist()
    print(f"Mapping {len(gene_names)} genes to UniProt IDs via UniProt API...")
    t0 = time.time()
    gene_map = fetch_gene_to_uniprot_mapping(gene_names)
    print(f"  UniProt ID mapping completed in {time.time() - t0:.1f}s")

    result = cosmic_grouped.merge(gene_map, on="gene", how="left")
    unmapped = result["UniProt"].isna().sum()
    result = result[result["UniProt"].notna()].copy()
    result = result.rename(columns={"UniProt": "uniprot_id"})
    result = result[["uniprot_id", "gene", "mutations_on_protein"]]

    # -----------------------
    # Flag genes where COSMIC's mutation numbering follows a different UniProt
    # isoform than the canonical AlphaFold-modeled sequence.
    # -----------------------
    print("Checking for COSMIC/canonical isoform mismatches...")
    gene_to_uniprot = dict(zip(result["gene"], result["uniprot_id"]))
    t0 = time.time()
    isoform_lengths = compute_isoform_safe_lengths(gene_to_transcript, gene_to_uniprot)
    print(f"  Isoform mismatch check completed in {time.time() - t0:.1f}s")
    result = result.merge(isoform_lengths, on="gene", how="left")

    # -----------------------
    # Total COSMIC missense patients per gene, for comparison against the
    # nearby/distant mutation patient counts computed in step 3.
    # -----------------------
    result["total_cosmic_missense_patients"] = result["gene"].map(gene_to_total_missense_patients)

    print("Saving output...")
    result.to_csv(output_file, sep="\t", index=False)

    print("Done.")
    print(f"Hotspot minimum affected cases: {HOTSPOT_MIN_AFFECTED_CASES}")
    print(f"COSMIC hotspot genes: {cosmic['gene'].nunique()}")
    print(f"Genes mapped to UniProt: {len(result)}")
    print(f"Genes not mapped to UniProt (excluded): {unmapped}")
    print(f"Output saved to: {output_file}")


def main():
    """Parse CLI arguments and dispatch to the selected pipeline filter mode (ptm-proximity or mutation-clustering)."""
    global HOTSPOT_MIN_AFFECTED_CASES

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
    parser.add_argument(
        "--min-samples",
        type=int,
        default=HOTSPOT_MIN_AFFECTED_CASES,
        help=f"Minimum distinct COSMIC samples for a mutation to be a hotspot (default: {HOTSPOT_MIN_AFFECTED_CASES})",
    )
    args = parser.parse_args()
    HOTSPOT_MIN_AFFECTED_CASES = args.min_samples

    output_file = PROJECT_ROOT / "data" / "steps" / "PTMD_TCGA_hotspots_by_protein.tsv"
    output_file.parent.mkdir(parents=True, exist_ok=True)

    if args.mode == "mutation-clustering":
        _run_mutation_clustering_filter(output_file)
    else:
        _run_ptm_proximity_filter(output_file)


if __name__ == "__main__":
    main()
