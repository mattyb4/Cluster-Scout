"""Pipeline Step 1: filter and merge the raw input data into the per-protein
hotspot table the rest of the pipeline reads.

--mode ptm-proximity (default) merges PTMD's disease-associated PTM sites
with COSMIC's recurrent missense mutations (>= HOTSPOT_MIN_AFFECTED_CASES
distinct samples, confirmed/reported-somatic only), keeping only genes that
have both a PTM site and at least one qualifying mutation. Writes
data/steps/PTMD_COSMIC_hotspots_by_protein.tsv. With --ptm-source psp, the PTM
sites come from a PhosphoSitePlus download instead (every human site, with its
LTP/HTP evidence counts and any disease associations -- see psp_input.py), and
the output is data/steps/PSP_COSMIC_hotspots_by_protein.tsv.

--mode mutation-clustering keeps every COSMIC hotspot mutation mapped to a
reviewed human UniProt accession, with no PTM requirement. Writes
data/steps/COSMIC_hotspots_by_protein.tsv.

Both modes also map UniProt accessions to gene symbols (or vice versa) via
the UniProt REST API, and put COSMIC's mutation positions onto the canonical
AlphaFold-modeled sequence when COSMIC numbers a protein by a different
isoform (see cosmic_numbering.py) -- both cached under data/cache/ so repeat
runs only look up whatever wasn't already resolved.
"""
import argparse
import ast
import re
import sys
import time
from pathlib import Path

import pandas as pd
import requests
from tqdm import tqdm

sys.path.insert(0, str(Path(__file__).resolve().parent))
from cosmic_numbering import fetch_isoform_sequences, number_protein  # noqa: E402
from pipeline_utils import (  # noqa: E402
    COSMIC_INPUT_DIR,
    COSMIC_SOMATIC_STATUSES,
    DEFAULT_PTM_SOURCE,
    PSP_INPUT_DIR,
    PTM_SOURCES,
    PTMD_INPUT_DIR,
    fmt_time,
    hotspots_tsv_path,
    input_dir,
    project_root,
    ptm_output_paths,
    resolve_input_file,
)
from psp_input import PspInputError, filter_psp_sites, load_psp_diseases, load_psp_sites  # noqa: E402

PROJECT_ROOT = project_root(__file__)
HOTSPOT_MIN_AFFECTED_CASES = 3

UNMATCHED_GENES_LOG = ptm_output_paths(PROJECT_ROOT / "Output", "ptmd")["unmatched"]
PSP_UNMATCHED_GENES_LOG = ptm_output_paths(PROJECT_ROOT / "Output", "psp")["unmatched"]

CACHE_DIR = PROJECT_ROOT / "data" / "cache"

# Two equal-weighted phases (gene mapping, COSMIC numbering) reported as one continuous progress bar
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
    if not site:
        return ""

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
    """Load and filter the COSMIC Mutant Census, shared by both pipeline modes.

    Aggregates the raw (mutation, sample) rows into per-(gene, amino-acid
    change) affected-case counts. Also returns a gene -> Ensembl transcript
    mapping (COSMIC's transcript for each gene).
    """
    cols = ["GENE_SYMBOL", "MUTATION_AA", "COSMIC_SAMPLE_ID", "MUTATION_SOMATIC_STATUS", "TRANSCRIPT_ACCESSION"]
    cosmic = pd.read_csv(cosmic_file, sep="\t", usecols=cols, low_memory=False)

    cosmic = cosmic[cosmic["MUTATION_SOMATIC_STATUS"].isin(COSMIC_SOMATIC_STATUSES)].copy()

    cosmic["aa_change"] = cosmic["MUTATION_AA"].str.replace(r"^p\.", "", regex=True)
    cosmic = cosmic[cosmic["aa_change"].apply(is_simple_substitution)].copy()

    gene_to_transcript = cosmic.groupby("GENE_SYMBOL")["TRANSCRIPT_ACCESSION"].first().to_dict()

    # Total distinct patients per gene, regardless of hotspot threshold below.
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
    # .apply(axis=1) on an empty frame returns a DataFrame, not a Series -- guard
    # against that so a "zero rows survived filtering" input doesn't crash here.
    cosmic["mutation_with_count"] = (
        cosmic.apply(format_mutation_with_count, axis=1) if not cosmic.empty
        else pd.Series(dtype=str)
    )

    return cosmic, gene_to_transcript, gene_to_total_missense_patients


ISOFORM_SEQUENCE_CACHE_FILE = "uniprot_isoform_sequences.tsv"
MUTATION_WITH_COUNT_RE = re.compile(r"^(\S+) \((\d+)\)$")


def apply_canonical_numbering(table):
    """Rewrite each protein's mutations_on_protein ('R248Q (5); ...') in the
    canonical sequence's numbering, and add the columns recording what changed:

      cosmic_numbering_isoform  UniProt isoform COSMIC numbers this protein by
                                (blank = canonical, nothing moved)
      cosmic_mutation_labels    'R249Q=R248Q; ...' -- output label = COSMIC's
                                label, for every mutation that moved
      unmapped_mutations        COSMIC labels with no canonical position (in
                                isoform-only sequence, or no isoform explains
                                the reference residue); step 3 tags these
                                (isoform?)
    """
    isoforms = fetch_isoform_sequences(
        table["uniprot_id"].tolist(), CACHE_DIR / ISOFORM_SEQUENCE_CACHE_FILE,
        progress=lambda pct, desc: _emit_progress(1, pct, desc),
    )
    rewritten, isoform_col, labels_col, unmapped_col = [], [], [], []
    n_moved = n_unmapped = 0
    for uniprot, field in zip(table["uniprot_id"], table["mutations_on_protein"]):
        entries = []
        for token in str(field).split(";"):
            m = MUTATION_WITH_COUNT_RE.match(token.strip())
            if m:
                entries.append((m.group(1), m.group(2)))
        numbering = number_protein(uniprot, isoforms.get(uniprot, {}), [label for label, _ in entries])
        rewritten.append("; ".join(f"{numbering.labels[label]} ({count})" for label, count in entries))
        moved = [(numbering.labels[label], label) for label, _ in entries if numbering.labels[label] != label]
        isoform_col.append(numbering.isoform)
        labels_col.append("; ".join(f"{new}={old}" for new, old in moved))
        unmapped_col.append("; ".join(numbering.unmapped))
        n_moved += len(moved)
        n_unmapped += len(numbering.unmapped)

    table = table.copy()
    table["mutations_on_protein"] = rewritten
    table["cosmic_numbering_isoform"] = isoform_col
    table["cosmic_mutation_labels"] = labels_col
    table["unmapped_mutations"] = unmapped_col
    n_renumbered = sum(1 for iso in isoform_col if iso)
    print(f"  {n_renumbered} proteins use a non-canonical isoform's numbering in COSMIC: "
          f"{n_moved} mutations moved to canonical positions, {n_unmapped} couldn't be placed "
          f"(kept with COSMIC's label, tagged (isoform?) in step 3)")
    return table


def _build_ptmd_ptm_table():
    """Read PTMD's disease-associated PTM sites into one row per protein:
    uniprot_id, gene, ptms_on_protein, ptm_disease_pairs, ptm_known_disruptions.
    Returns (table, number of PTMD genes)."""
    ptmd_file = resolve_input_file(input_dir(PROJECT_ROOT, PTMD_INPUT_DIR), (".tsv",))
    print(f"PTMD file:   {ptmd_file.name}")
    ptmd = pd.read_csv(ptmd_file, sep="\t", low_memory=False)

    # Filter PTMD disruptions
    ptmd = ptmd[ptmd["State"] == "N"].copy()

    # Normalize variant UniProt IDs to canonical accession (e.g. Q16613_VAR_A129T -> Q16613)
    ptmd["UniProt"] = ptmd["UniProt"].str.split("_").str[0]

    # Map UniProt -> gene via UniProt REST API
    uniprot_ids = ptmd["UniProt"].dropna().unique().tolist()
    print(f"Mapping {len(uniprot_ids)} UniProt IDs to gene names via UniProt API...")
    t0 = time.time()
    idmap = fetch_uniprot_gene_mapping(uniprot_ids)
    print(f"  UniProt gene mapping completed in {fmt_time(time.time() - t0)}")

    ptmd = ptmd.merge(idmap, on="UniProt", how="left")

    if "Gene name" in ptmd.columns:
        ptmd["gene"] = ptmd["gene"].fillna(ptmd["Gene name"])

    ptmd = ptmd[ptmd["gene"].notna()].copy()

    # Build PTM site and PTM-disease pair
    ptmd["ptm_site"] = ptmd.apply(build_ptm_site, axis=1)

    # Drop PTMD rows with no usable residue+position (e.g. both fields blank in
    # the source data) -- build_ptm_site returns "" for these; keeping them would
    # add meaningless entries like ":Phosphorylation" to ptms_on_protein.
    missing_site = (ptmd["ptm_site"] == "").sum()
    if missing_site:
        print(f"Dropping {missing_site} PTMD row(s) with no residue/position")
    ptmd = ptmd[ptmd["ptm_site"] != ""].copy()

    ptmd["ptm_disease_pair"] = ptmd["ptm_site"] + " | " + ptmd["Disease"].astype(str)

    # Build known disrupting mutations per PTM site (MutationSite = mutations
    # documented in literature as disrupting that PTM)
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

    # Aggregate PTMs by UniProt, not gene -- grouping by gene collapses all
    # isoforms under one UniProt ID, so PTM positions from one isoform get
    # checked against another isoform's structure.
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
    return ptmd_grouped, ptmd["gene"].nunique()


def format_psp_site_scores(row):
    """Encode a PSP site's evidence counts for the step-1 table as
    'S473:Phosphorylation=LTP/HTP/CST', e.g. 'S473:Phosphorylation=12/40/3'."""
    return f'{row["ptm_site"]}={row["ltp"]}/{row["htp"]}/{row["cst"]}'


def _build_psp_ptm_table(psp_filters=None):
    """Read the human PhosphoSitePlus sites passing *psp_filters* (keyword
    arguments for psp_input.filter_psp_sites; none = every site) into one row
    per protein: uniprot_id, gene, ptms_on_protein, ptm_disease_pairs,
    psp_site_scores. Returns (table, number of PSP genes)."""
    folder = input_dir(PROJECT_ROOT, PSP_INPUT_DIR)
    print(f"PhosphoSitePlus folder: {folder}")
    sites = load_psp_sites(folder)
    diseases = load_psp_diseases(folder)
    sites = filter_psp_sites(sites, diseases, **(psp_filters or {}))
    if sites.empty:
        raise PspInputError("No PhosphoSitePlus sites pass the LTP/HTP/disease filters -- loosen them")
    print(f"{len(sites):,} human PTM sites on {sites['UniProt'].nunique():,} proteins")

    # Same accession -> gene mapping as the PTMD path, so gene symbols match
    # COSMIC's the same way; PSP's own GENE column is the fallback.
    uniprot_ids = sites["UniProt"].unique().tolist()
    print(f"Mapping {len(uniprot_ids)} UniProt IDs to gene names via UniProt API...")
    t0 = time.time()
    idmap = fetch_uniprot_gene_mapping(uniprot_ids)
    print(f"  UniProt gene mapping completed in {fmt_time(time.time() - t0)}")
    sites = sites.merge(idmap, on="UniProt", how="left")
    sites["gene"] = sites["gene"].fillna(sites["psp_gene"].where(sites["psp_gene"] != ""))
    sites = sites[sites["gene"].notna()].copy()

    sites["ptm_site"] = (sites["residue"] + sites["position"].astype(str)
                         + ":" + sites["ptm_type"])
    sites["psp_site_score"] = sites.apply(format_psp_site_scores, axis=1)

    # Disease pairs use PTMD's "site | disease" format, so step 3 reads both
    # sources the same way; the disease text carries PSP's alteration.
    if not diseases.empty:
        diseases["ptm_site"] = (diseases["residue"] + diseases["position"].astype(str)
                                + ":" + diseases["ptm_type"])
        diseases["ptm_disease_pair"] = diseases["ptm_site"] + " | " + diseases["disease"]
        disease_pairs = (
            diseases.groupby("UniProt")["ptm_disease_pair"].apply(clean_str_list)
            .reset_index().rename(columns={"UniProt": "uniprot_id",
                                           "ptm_disease_pair": "ptm_disease_pairs"})
        )
    else:
        disease_pairs = pd.DataFrame(columns=["uniprot_id", "ptm_disease_pairs"])

    sites = sites.sort_values(["UniProt", "position", "ptm_type"])
    psp_grouped = (
        sites.groupby("UniProt", as_index=False)
        .agg(
            gene=("gene", "first"),
            ptms_on_protein=("ptm_site", clean_str_list),
            psp_site_scores=("psp_site_score", clean_str_list),
        )
        .rename(columns={"UniProt": "uniprot_id"})
    )
    psp_grouped = psp_grouped.merge(disease_pairs, on="uniprot_id", how="left")
    psp_grouped = psp_grouped[["uniprot_id", "gene", "ptms_on_protein",
                               "ptm_disease_pairs", "psp_site_scores"]]
    return psp_grouped, sites["gene"].nunique()


def _canonical_numbering_step(table):
    """COSMIC sometimes numbers a protein by a different isoform than the
    canonical AlphaFold-modeled sequence -- move those mutations into
    canonical numbering so they're measured from the right residue."""
    print("Putting COSMIC mutations in canonical-sequence numbering...")
    t0 = time.time()
    table = apply_canonical_numbering(table)
    print(f"  Canonical numbering completed in {fmt_time(time.time() - t0)}")
    return table


def _run_ptm_proximity_filter(output_file, ptm_source=DEFAULT_PTM_SOURCE, psp_filters=None):
    """Run the PTM-proximity pipeline mode: merge PTM sites (from PTMD, or
    PhosphoSitePlus with ptm_source="psp", optionally narrowed by
    *psp_filters*) with COSMIC hotspots, keeping only genes with both."""
    cosmic_file = resolve_input_file(input_dir(PROJECT_ROOT, COSMIC_INPUT_DIR), (".tsv",))
    if ptm_source == "psp":
        source_name, unmatched_log = "PhosphoSitePlus", PSP_UNMATCHED_GENES_LOG
        ptm_grouped, n_ptm_genes = _build_psp_ptm_table(psp_filters)
    else:
        source_name, unmatched_log = "PTMD", UNMATCHED_GENES_LOG
        ptm_grouped, n_ptm_genes = _build_ptmd_ptm_table()

    print(f"COSMIC file: {cosmic_file.name}")
    cosmic, _gene_to_transcript, gene_to_total_missense_patients = _load_and_filter_cosmic(cosmic_file)

    print("Filtering COSMIC mutations and aggregating by gene...")
    # Aggregate hotspot mutations by gene
    cosmic_grouped = (
        cosmic.groupby("gene", as_index=False)
        .agg(
            mutations_on_protein=("mutation_with_count", clean_str_list),
        )
    )

    merged = ptm_grouped.merge(cosmic_grouped, on="gene", how="left")

    # Source-specific columns (disease pairs, disruptions / PSP scores) follow
    # the shared ones, in the order the source's table lists them
    source_cols = [c for c in ptm_grouped.columns
                   if c not in ("uniprot_id", "gene", "ptms_on_protein")]
    merged = merged[["uniprot_id", "gene", "ptms_on_protein", "mutations_on_protein", *source_cols]]

    # Log PTM-source proteins with no matching COSMIC hotspot mutations for their gene
    # (gene-name mismatch, or no mutation met HOTSPOT_MIN_AFFECTED_CASES) before dropping them
    unmatched = merged[merged["mutations_on_protein"].isna()].copy()
    unmatched_log.parent.mkdir(parents=True, exist_ok=True)
    unmatched[["uniprot_id", "gene", "ptms_on_protein"]].to_csv(
        unmatched_log, sep="\t", index=False, encoding="utf-16"
    )
    print(f"Wrote {len(unmatched)} proteins with no matching COSMIC hotspot mutations to: {unmatched_log}")

    merged = merged[merged["mutations_on_protein"].notna()].copy()
    merged = _canonical_numbering_step(merged)

    # For comparison against nearby/distant mutation patient counts computed in step 3
    merged["total_cosmic_missense_patients"] = merged["gene"].map(gene_to_total_missense_patients)

    print("Saving output...")
    merged.to_csv(output_file, sep="\t", index=False)

    print("Done.")
    print(f"Hotspot minimum affected cases: {HOTSPOT_MIN_AFFECTED_CASES}")
    print(f"{source_name} PTM genes: {n_ptm_genes}")
    print(f"COSMIC hotspot genes: {cosmic['gene'].nunique()}")
    print(f"Final merged proteins: {len(merged)}")
    print(f"Output saved to: {output_file}")


def _run_mutation_clustering_filter(output_file):
    """Run the mutation-clustering pipeline mode: keep all recurrent COSMIC hotspots mapped to UniProt,
    regardless of PTMs."""
    cosmic_file = resolve_input_file(input_dir(PROJECT_ROOT, COSMIC_INPUT_DIR), (".tsv",))
    print(f"COSMIC file: {cosmic_file.name}")

    print("Loading COSMIC file...")
    cosmic, _gene_to_transcript, gene_to_total_missense_patients = _load_and_filter_cosmic(cosmic_file)

    print("Aggregating hotspot mutations by gene...")
    cosmic_grouped = (
        cosmic.groupby("gene", as_index=False)
        .agg(mutations_on_protein=("mutation_with_count", clean_str_list))
    )

    gene_names = cosmic_grouped["gene"].tolist()
    print(f"Mapping {len(gene_names)} genes to UniProt IDs via UniProt API...")
    t0 = time.time()
    gene_map = fetch_gene_to_uniprot_mapping(gene_names)
    print(f"  UniProt ID mapping completed in {fmt_time(time.time() - t0)}")

    result = cosmic_grouped.merge(gene_map, on="gene", how="left")
    unmapped = result["UniProt"].isna().sum()
    result = result[result["UniProt"].notna()].copy()
    result = result.rename(columns={"UniProt": "uniprot_id"})
    result = result[["uniprot_id", "gene", "mutations_on_protein"]]

    result = _canonical_numbering_step(result)

    # For comparison against nearby/distant mutation patient counts computed in step 3
    result["total_cosmic_missense_patients"] = result["gene"].map(gene_to_total_missense_patients)

    print("Saving output...")
    result.to_csv(output_file, sep="\t", index=False)

    print("Done.")
    print(f"Hotspot minimum affected cases: {HOTSPOT_MIN_AFFECTED_CASES}")
    print(f"COSMIC hotspot genes: {cosmic['gene'].nunique()}")
    print(f"Genes mapped to UniProt: {len(result)}")
    print(f"Genes not mapped to UniProt (excluded): {unmapped}")
    print(f"Output saved to: {output_file}")


def _non_negative_int(text):
    """argparse type for the PSP evidence-count bounds."""
    value = int(text)
    if value < 0:
        raise argparse.ArgumentTypeError(f"must be 0 or more, got {value}")
    return value


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
        "--ptm-source",
        choices=PTM_SOURCES,
        default=DEFAULT_PTM_SOURCE,
        help="PTM site data for ptm-proximity mode: 'ptmd' (default) or 'psp' (PhosphoSitePlus)",
    )
    parser.add_argument("--psp-min-ltp", type=_non_negative_int, default=0,
                        help="PhosphoSitePlus source only: minimum LTP (low-throughput literature) count")
    parser.add_argument("--psp-max-ltp", type=_non_negative_int, default=None,
                        help="PhosphoSitePlus source only: maximum LTP count (default: no limit)")
    parser.add_argument("--psp-min-htp", type=_non_negative_int, default=0,
                        help="PhosphoSitePlus source only: minimum HTP (high-throughput mass-spec) count")
    parser.add_argument("--psp-max-htp", type=_non_negative_int, default=None,
                        help="PhosphoSitePlus source only: maximum HTP count (default: no limit)")
    parser.add_argument("--psp-disease-only", action="store_true",
                        help="PhosphoSitePlus source only: keep only sites with a PSP disease association")
    parser.add_argument(
        "--min-samples",
        type=int,
        default=HOTSPOT_MIN_AFFECTED_CASES,
        help=f"Minimum distinct COSMIC samples for a mutation to be a hotspot (default: {HOTSPOT_MIN_AFFECTED_CASES})",
    )
    args = parser.parse_args()
    HOTSPOT_MIN_AFFECTED_CASES = args.min_samples
    for kind in ("ltp", "htp"):
        lo, hi = getattr(args, f"psp_min_{kind}"), getattr(args, f"psp_max_{kind}")
        if hi is not None and hi < lo:
            parser.error(f"--psp-max-{kind} ({hi}) is below --psp-min-{kind} ({lo})")
    psp_filters = {
        "min_ltp": args.psp_min_ltp, "max_ltp": args.psp_max_ltp,
        "min_htp": args.psp_min_htp, "max_htp": args.psp_max_htp,
        "disease_only": args.psp_disease_only,
    }

    output_file = hotspots_tsv_path(PROJECT_ROOT, args.mode, args.ptm_source)
    output_file.parent.mkdir(parents=True, exist_ok=True)

    if args.mode == "mutation-clustering":
        _run_mutation_clustering_filter(output_file)
    else:
        try:
            _run_ptm_proximity_filter(output_file, args.ptm_source, psp_filters)
        except PspInputError as exc:
            print(f"\nPhosphoSitePlus input problem:\n{exc}", file=sys.stderr)
            sys.exit(1)


if __name__ == "__main__":
    main()
