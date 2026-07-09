import sys
import numpy as np
from pathlib import Path
import re
import csv
import json
import argparse
from typing import Any
from tqdm import tqdm

sys.path.insert(0, str(Path(__file__).resolve().parent))
from pipeline_utils import (  # noqa: E402
    project_root, AA3TO1, MUT_RE, SITE_RE,
    find_canonical_cif, load_first_chain, load_pae_matrix, get_plddt_map,
)

PROJECT_ROOT = project_root(__file__)
MODELS_ROOT = PROJECT_ROOT / "cif_models"
DEFAULT_OUTPUT_DIR = PROJECT_ROOT / "Output"
PTM_TSV_PATH = PROJECT_ROOT / "data" / "steps" / "PTMD_COSMIC_hotspots_by_protein.tsv"

_PTM_ROWS: list[dict[str, Any]] | None = None

DISTANCE_CUTOFF = 10.0  # Angstroms, adjust as needed

def get_ptm_rows():
    """Lazy-load and cache PTM rows from the intermediate TSV to avoid repeated file I/O."""
    global _PTM_ROWS
    if _PTM_ROWS is None:
        with PTM_TSV_PATH.open("r", encoding="utf-8", newline="") as handle:
            _PTM_ROWS = list(csv.DictReader(handle, delimiter="\t"))
    return _PTM_ROWS

def get_ca_coord(chain, residue_number):
    """Return the alpha-carbon coordinate for a residue number, or None if not found."""
    mask = (chain.res_id == residue_number) & (chain.atom_name == "CA")
    if not np.any(mask):
        return None
    return chain.coord[mask][0]

def compute_distance(coord1, coord2):
    """Compute the Euclidean distance between two 3D coordinate arrays."""
    return np.linalg.norm(coord1 - coord2)

def find_nearby_mutations(chain, ptm_pos, mutation_entries, pae_matrix=None, cutoff=DISTANCE_CUTOFF, max_pae=None):
    """Find mutations within a distance cutoff of a PTM site, with optional PAE filtering."""
    results = []

    ptm_coord = get_ca_coord(chain, ptm_pos)

    if ptm_coord is None:
        return results

    for mutation, mut_pos in mutation_entries:
        mut_coord = get_ca_coord(chain, mut_pos)

        if mut_coord is None:
            continue

        distance = compute_distance(ptm_coord, mut_coord)

        if distance <= cutoff:
            pae = None
            if pae_matrix is not None:
                i, j = ptm_pos - 1, mut_pos - 1
                if 0 <= i < pae_matrix.shape[0] and 0 <= j < pae_matrix.shape[1]:
                    pae = (pae_matrix[i, j] + pae_matrix[j, i]) / 2
            if max_pae is not None and pae is not None and pae > max_pae:
                continue
            results.append({
                "mutation": mutation,
                "mutation_pos": mut_pos,
                "distance": distance,
                "pae": pae,
            })

    return results


def find_mutation_clusters(chain, mutation_entries, pae_matrix=None, cutoff=DISTANCE_CUTOFF, max_pae=None):
    """For each mutation, find other mutations within cutoff Angstroms in 3D space."""
    mut_list = list(mutation_entries)
    results = {}

    for i, (anchor_mut, anchor_pos) in enumerate(mut_list):
        anchor_coord = get_ca_coord(chain, anchor_pos)
        if anchor_coord is None:
            continue

        nearby = []
        for j, (other_mut, other_pos) in enumerate(mut_list):
            if i == j:
                continue
            other_coord = get_ca_coord(chain, other_pos)
            if other_coord is None:
                continue
            distance = compute_distance(anchor_coord, other_coord)
            if distance <= cutoff:
                pae = None
                if pae_matrix is not None:
                    ii, jj = anchor_pos - 1, other_pos - 1
                    if 0 <= ii < pae_matrix.shape[0] and 0 <= jj < pae_matrix.shape[1]:
                        pae = (pae_matrix[ii, jj] + pae_matrix[jj, ii]) / 2
                if max_pae is not None and pae is not None and pae > max_pae:
                    continue
                nearby.append({
                    "mutation": other_mut,
                    "mutation_pos": other_pos,
                    "distance": distance,
                    "pae": pae,
                })

        if nearby:
            results[(anchor_mut, anchor_pos)] = nearby

    return results


PTM_RE = re.compile(r"([A-Z])(\d+)")  # e.g., S557


def parse_ptm_entries(uniprot):
    """Extract deduplicated PTM site positions and modification types for a UniProt ID from the input TSV."""
    entries = {}  # (ptm_site, position) -> ptm_type, deduplicates by site+type

    for row in get_ptm_rows():
        if row.get("uniprot_id") != uniprot:
            continue
        field = row.get("ptms_on_protein", "")
        for token in re.split(r";", field):
            token = token.strip()
            if ":" in token:
                site_part, ptm_type = token.split(":", 1)
                ptm_type = ptm_type.strip()
            else:
                site_part = token
                ptm_type = ""
            match = PTM_RE.search(site_part.strip())
            if match:
                ptm_site = f"{match.group(1)}{match.group(2)}"
                position = int(match.group(2))
                entries[(ptm_site, position)] = ptm_type

    return sorted([(site, pos, mod) for (site, pos), mod in entries.items()], key=lambda x: x[1])


def parse_gene_name(uniprot):
    """Look up the gene symbol for a UniProt ID in the input TSV."""
    for row in get_ptm_rows():
        if row.get("uniprot_id") == uniprot:
            return row.get("gene", "")

    return ""


def parse_isoform_safe_length(uniprot):
    """Return the position past which COSMIC's mutation numbering for this protein
    no longer matches the canonical AlphaFold-modeled sequence, or None if COSMIC's
    numbering matches canonical throughout."""
    for row in get_ptm_rows():
        if row.get("uniprot_id") == uniprot:
            value = row.get("isoform_safe_length", "")
            if value and value.strip():
                return int(float(value))
            return None

    return None

MUT_COUNT_RE = re.compile(r"\((\d+)\)")  # e.g., (5) in "R482H (5)"

def parse_mutation_positions(uniprot=None):
    """Extract all mutation positions and labels for a protein from the input TSV."""
    mutation_entries = set()

    for row in get_ptm_rows():
        if uniprot and row.get("uniprot_id") != uniprot:
            continue

        fields = [row.get("mutations_on_protein", "")]
        for field in fields:
            for token in re.split(r"[;,]", field):
                match = MUT_RE.search(token.strip())
                if match:
                    mutation = f"{match.group(1)}{match.group(2)}{match.group(3)}"
                    mutation_entries.add((mutation, int(match.group(2))))

    return sorted(mutation_entries, key=lambda x: (x[1], x[0]))


def parse_mutation_patient_counts(uniprot):
    """Map each (mutation, position) for this protein to its COSMIC affected-case
    (patient) count, as recorded in mutations_on_protein, e.g. "R482H (5)" -> 5."""
    counts: dict[tuple[str, int], int] = {}

    for row in get_ptm_rows():
        if row.get("uniprot_id") != uniprot:
            continue

        for token in re.split(r"[;,]", row.get("mutations_on_protein", "")):
            token = token.strip()
            mut_match = MUT_RE.search(token)
            count_match = MUT_COUNT_RE.search(token)
            if mut_match and count_match:
                mutation = f"{mut_match.group(1)}{mut_match.group(2)}{mut_match.group(3)}"
                counts[(mutation, int(mut_match.group(2)))] = int(count_match.group(1))

    return counts


def parse_total_cosmic_missense_patients(uniprot):
    """Return the total number of distinct COSMIC patients with any missense
    mutation in this protein's gene, regardless of the hotspot recurrence
    threshold, or None if not available."""
    for row in get_ptm_rows():
        if row.get("uniprot_id") == uniprot:
            value = row.get("total_cosmic_missense_patients", "")
            if value and value.strip():
                return int(float(value))
            return None

    return None


find_model_file = find_canonical_cif


def format_mutations(hits):
    """Format mutation hits as a comma-separated string with distances and PAE scores."""
    if not hits:
        return ""
    parts = []
    for hit in sorted(hits, key=lambda h: (h["mutation_pos"], h["mutation"])):
        entry = f"{hit['mutation']}-{hit['distance']:.2f}Å"
        if hit.get("pae") is not None:
            entry += f"(PAE:{hit['pae']:.1f})"
        parts.append(entry)
    return ", ".join(parts)


def linear_distances(hits, ptm_pos):
    """Compute linear (sequence) distances between unique mutation positions and a PTM site."""
    if not hits:
        return ""
    seen = set()
    distances = []
    for hit in sorted(hits, key=lambda h: (h["mutation_pos"], h["mutation"])):
        pos = hit["mutation_pos"]
        if pos not in seen:
            seen.add(pos)
            distances.append(abs(pos - int(ptm_pos)))
    return ",".join(str(d) for d in distances)


def unique_mutation_position_count(hits):
    """Count distinct mutation positions in a list of hits."""
    return len({hit["mutation_pos"] for hit in hits})

def mutation_at_ptm_site(hits, ptm_pos):
    """Check whether any mutation hit is located at the PTM site itself."""
    return 'yes' if any(hit['mutation_pos'] == int(ptm_pos) for hit in hits) else 'no'

def total_patient_count(hits, patient_counts):
    """Sum COSMIC patient counts across all hits, looking each up by its
    (mutation, position) with any "(isoform?)" tag stripped."""
    total = 0
    for hit in hits:
        mutation = hit["mutation"].replace("(isoform?)", "")
        total += patient_counts.get((mutation, hit["mutation_pos"]), 0)
    return total


def parse_ptm_diseases(uniprot, ptm_site, ptm_type):
    """Extract cancer-related disease associations for a PTM site from PTMD data."""
    CANCER_KEYWORDS = {
        "cancer", "carcinoma", "sarcoma", "lymphoma", "leukemia", "leukaemia",
        "melanoma", "glioma", "glioblastoma", "myeloma", "blastoma", "tumor",
        "tumour", "neoplasm", "mesothelioma", "neuroblastoma", "adenoma",
    }
    diseases = []
    for row in get_ptm_rows():
        if row.get("uniprot_id") != uniprot:
            continue
        for entry in row.get("ptm_disease_pairs", "").split(";"):
            entry = entry.strip()
            if " | " not in entry:
                continue
            site_type, disease = entry.split(" | ", 1)
            if site_type.strip() == f"{ptm_site}:{ptm_type}":
                disease = disease.strip()
                if disease and disease not in diseases:
                    if any(kw in disease.lower() for kw in CANCER_KEYWORDS):
                        diseases.append(disease)
    return "; ".join(diseases)

def parse_ptm_known_disruptions(uniprot):
    """Return dict mapping 'S516:Phosphorylation' -> set of known disrupting mutations for this protein.

    Reads the ptm_known_disruptions column produced by step 1, which encodes entries as
    'S516:Phosphorylation>D120N,E127D; K43:Ubiquitination>R40Q'.
    """
    result = {}
    for row in get_ptm_rows():
        if row.get("uniprot_id") != uniprot:
            continue
        field = row.get("ptm_known_disruptions", "") or ""
        if not field or field.strip() in ("", "nan"):
            continue
        for entry in field.split(";"):
            entry = entry.strip()
            if ">" not in entry:
                continue
            site_full, muts_str = entry.split(">", 1)
            mutations = {m.strip() for m in muts_str.split(",") if m.strip()}
            key = site_full.strip()
            result[key] = result.get(key, set()) | mutations
    return result




def main():
    """Main entry point: parse CLI args and run the PTM-proximity or mutation-clustering pipeline."""
    global DISTANCE_CUTOFF

    parser = argparse.ArgumentParser(description="Scan AFDB models for nearby mutations.")
    parser.add_argument("--uniprot", help="Limit processing to a single UniProt ID.")
    parser.add_argument(
        "--mode",
        choices=["ptm-proximity", "mutation-clustering"],
        default="ptm-proximity",
        help=(
            "'ptm-proximity' (default) finds cancer mutations clustering near PTM sites. "
            "'mutation-clustering' finds recurrent mutations that cluster together in 3D "
            "space, with no PTM anchor."
        ),
    )
    parser.add_argument(
        "--output-dir",
        default=str(DEFAULT_OUTPUT_DIR),
        help="Directory for output files (default: Output/)",
    )
    parser.add_argument(
        "--cutoff",
        type=float,
        default=DISTANCE_CUTOFF,
        help=f"Distance cutoff in Angstroms (default: {DISTANCE_CUTOFF})",
    )
    parser.add_argument(
        "--min-plddt",
        type=float,
        default=None,
        help="Exclude positions with pLDDT below this threshold (default: disabled)",
    )
    parser.add_argument(
        "--max-pae",
        type=float,
        default=None,
        help="Exclude mutation pairs with PAE above this threshold (default: disabled)",
    )
    args = parser.parse_args()
    DISTANCE_CUTOFF = args.cutoff
    MIN_PLDDT = args.min_plddt or 0
    MAX_PAE = args.max_pae
    if MIN_PLDDT > 0:
        print(f"pLDDT filter: excluding positions below {MIN_PLDDT}")
    if MAX_PAE is not None:
        print(f"PAE filter: excluding pairs above {MAX_PAE}")

    output_dir = Path(args.output_dir)
    OUTPUT_PATH = output_dir / "ptm_mutation_proximity_db.tsv"
    LONG_OUTPUT_PATH = output_dir / "ptm_mutation_proximity_long.tsv"
    SKIPPED_PATH = output_dir / "logs" / "ptm_skipped.tsv"
    CLUSTER_OUTPUT_PATH = output_dir / "mutation_cluster_db.tsv"
    CLUSTER_SKIPPED_PATH = output_dir / "logs" / "mutation_cluster_skipped.tsv"
    OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    SKIPPED_PATH.parent.mkdir(parents=True, exist_ok=True)

    SKIP_HEADER = ["UniProt", "gene", "ptm_site", "ptm_type", "skip_reason", "detail"]

    def write_skips(skip_writer, uniprot, gene, ptm_entries, reason, detail):
        """Log a skip reason for every PTM entry of a protein."""
        for ptm_site, _, ptm_type in ptm_entries:
            skip_writer.writerow([uniprot, gene, ptm_site, ptm_type, reason, detail])

    def write_skip(skip_writer, uniprot, gene, ptm_site, ptm_type, reason, detail):
        """Log a skip reason for a single PTM site."""
        skip_writer.writerow([uniprot, gene, ptm_site, ptm_type, reason, detail])


    dirs_present = {d.name for d in MODELS_ROOT.iterdir() if d.is_dir()}

    # ── Mutation-clustering mode ──────────────────────────────────────────────────
    if args.mode == "mutation-clustering":
        CLUSTER_OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)

        all_cluster_uniprots = {row["uniprot_id"] for row in get_ptm_rows() if row.get("uniprot_id")}

        with CLUSTER_OUTPUT_PATH.open("w", encoding="utf-16", newline="") as handle, \
             CLUSTER_SKIPPED_PATH.open("w", encoding="utf-16", newline="") as skip_handle:

            writer = csv.writer(handle, delimiter="\t")
            skip_writer = csv.writer(skip_handle, delimiter="\t")

            writer.writerow([
                "UniProt",
                "gene",
                "anchor_mutation",
                "nearby_mutations",
                "nearby_mutation_count",
                "unique_nearby_position_count",
            ])
            skip_writer.writerow(["UniProt", "gene", "skip_reason", "detail"])

            # Log proteins in input with no downloaded directory
            for uniprot in sorted(all_cluster_uniprots - dirs_present):
                if args.uniprot and uniprot != args.uniprot:
                    continue
                gene = parse_gene_name(uniprot)
                skip_writer.writerow([uniprot, gene, "no_afdb_directory",
                                       "protein not found in AlphaFold DB (no download directory)"])

            uniprot_dirs = [d for d in sorted(MODELS_ROOT.iterdir()) if d.is_dir()]
            for uniprot_dir in tqdm(uniprot_dirs, desc="Scanning structures"):
                uniprot = uniprot_dir.name
                if args.uniprot and uniprot != args.uniprot:
                    continue
                gene = parse_gene_name(uniprot)
                mutation_entries = parse_mutation_positions(uniprot=uniprot)
                if len(mutation_entries) < 2:
                    continue

                model_file = find_model_file(uniprot_dir)
                if model_file is None:
                    skip_writer.writerow([uniprot, gene, "no_canonical_cif",
                                           "AFDB has only isoform models, no canonical sequence model"])
                    tqdm.write(f"  {uniprot}: no canonical CIF file found")
                    continue

                chain = load_first_chain(model_file)
                if chain is None:
                    continue

                pos_to_aa: dict[int, str] = {}
                for atom in chain:
                    if atom.res_id not in pos_to_aa:
                        pos_to_aa[atom.res_id] = AA3TO1.get(atom.res_name, "?")

                plddt_map = get_plddt_map(chain) if MIN_PLDDT > 0 else {}

                safe_length = parse_isoform_safe_length(uniprot)
                mutation_entries = [
                    (mut + "(isoform?)", pos) if (
                        pos not in pos_to_aa
                        or pos_to_aa[pos] != mut[0]
                        or (safe_length is not None and pos > safe_length)
                    )
                    else (mut, pos)
                    for mut, pos in mutation_entries
                ]
                if MIN_PLDDT > 0:
                    before = len(mutation_entries)
                    mutation_entries = [(m, p) for m, p in mutation_entries
                                        if plddt_map.get(p, 0) >= MIN_PLDDT]
                    if before != len(mutation_entries):
                        tqdm.write(f"  {uniprot}: filtered {before - len(mutation_entries)} "
                                   f"mutations below pLDDT {MIN_PLDDT}")

                if len(mutation_entries) < 2:
                    continue

                pae_matrix = load_pae_matrix(uniprot_dir)
                clusters = find_mutation_clusters(chain, mutation_entries, pae_matrix=pae_matrix, cutoff=DISTANCE_CUTOFF, max_pae=MAX_PAE)

                for (anchor_mut, anchor_pos), nearby in sorted(clusters.items(), key=lambda x: (x[0][1], x[0][0])):
                    writer.writerow([
                        uniprot,
                        gene,
                        anchor_mut,
                        format_mutations(nearby),
                        len(nearby),
                        unique_mutation_position_count(nearby),
                    ])

        print(f"Wrote mutation cluster data to {CLUSTER_OUTPUT_PATH}")
        print(f"Wrote skipped proteins to {CLUSTER_SKIPPED_PATH}")

    # ── PTM-proximity mode (default) ──────────────────────────────────────────────
    else:
        # Collect all UniProt IDs in the PTM TSV to catch proteins with no directory at all
        all_ptm_uniprots = {row["uniprot_id"] for row in get_ptm_rows() if row.get("uniprot_id")}

        _LONG_HEADER = [
            "gene", "uniprot_id", "ptm_position", "ptm_type", "ptm_plddt",
            "mutation", "sequence_distance", "distance_angstrom",
            "polyphen_score", "polyphen_class", "mutation_plddt", "pair_pae",
            "patient_count", "total_nearby_patient_count", "total_cosmic_missense_patients",
            "nearby_mutation_count", "confirmed_disrupting_mutation", "ptm_diseases",
            "1433_predicted", "1433_predicted_consensus", "1433_confirmed", "kinase_predictions",
            "ptm_aiupred_general", "ptm_aiupred_binding", "ptm_is_disordered", "ptm_is_binding",
            "mut_aiupred_general", "mut_aiupred_binding", "mut_is_disordered", "mut_is_binding",
        ]

        long_handle = LONG_OUTPUT_PATH.open("w", encoding="utf-16", newline="")
        long_writer = csv.writer(long_handle, delimiter="\t")
        long_writer.writerow(_LONG_HEADER)

        try:
            with OUTPUT_PATH.open("w", encoding="utf-16", newline="") as handle, \
                 SKIPPED_PATH.open("w", encoding="utf-16", newline="") as skip_handle:

                writer = csv.writer(handle, delimiter="\t")
                skip_writer = csv.writer(skip_handle, delimiter="\t")

                writer.writerow([
                    "UniProt",
                    "gene",
                    "ptm_site",
                    "ptm_type",
                    "mutations_within_5_positions",
                    "mutation_count_within_5_positions",
                    "unique_mutation_position_count_within_5_positions",
                    "nearby_muts_total_patient_count",
                    "mutations_more_than_5_positions",
                    "mutation_count_more_than_5_positions",
                    "unique_mutation_position_count_more_than_5_positions",
                    "distant_muts_total_patient_count",
                    "morethan5_linear_distance",
                    "mutation_at_ptm_site",
                    "confirmed_disrupting_mutations",
                    "ptm_diseases",
                    "total_cosmic_missense_patients",
                ])
                skip_writer.writerow(SKIP_HEADER)

                # Proteins in PTM TSV with no downloaded directory at all (NO_ENTRY from AFDB)
                for uniprot in sorted(all_ptm_uniprots - dirs_present):
                    if args.uniprot and uniprot != args.uniprot:
                        continue
                    gene = parse_gene_name(uniprot)
                    ptm_entries = parse_ptm_entries(uniprot)
                    write_skips(skip_writer, uniprot, gene, ptm_entries, "no_afdb_directory",
                                "protein not found in AlphaFold DB (no download directory)")

                uniprot_dirs = [d for d in sorted(MODELS_ROOT.iterdir()) if d.is_dir()]
                for uniprot_dir in tqdm(uniprot_dirs, desc="Scanning structures"):
                    uniprot = uniprot_dir.name
                    if args.uniprot and uniprot != args.uniprot:
                        continue
                    gene = parse_gene_name(uniprot)
                    ptm_entries = parse_ptm_entries(uniprot)
                    if not ptm_entries:
                        write_skip(skip_writer, uniprot, gene, "", "",
                                   "no_ptm_entries", "no PTM sites parsed from PTMD data for this protein")
                        continue
                    mutation_entries = parse_mutation_positions(uniprot=uniprot)
                    known_disruptions = parse_ptm_known_disruptions(uniprot)
                    if not mutation_entries:
                        write_skips(skip_writer, uniprot, gene, ptm_entries, "no_mutation_entries",
                                    "no COSMIC hotspot mutations parsed for this protein")
                        continue

                    model_file = find_model_file(uniprot_dir)
                    if model_file is None:
                        write_skips(skip_writer, uniprot, gene, ptm_entries, "no_canonical_cif",
                                    "AFDB has only isoform models, no canonical sequence model")
                        tqdm.write(f"  {uniprot}: no canonical CIF file found")
                        continue

                    chain = load_first_chain(model_file)
                    if chain is None:
                        write_skips(skip_writer, uniprot, gene, ptm_entries, "cif_parse_failed",
                                    "canonical CIF file exists but could not be parsed")
                        continue

                    # Build residue -> 1-letter AA map for mismatch checking
                    pos_to_aa: dict[int, str] = {}
                    for atom in chain:
                        if atom.res_id not in pos_to_aa:
                            pos_to_aa[atom.res_id] = AA3TO1.get(atom.res_name, "?")

                    plddt_map = get_plddt_map(chain)

                    # Tag mutations whose reference AA does not match this structure, or whose
                    # position is past the point where COSMIC's isoform diverges from canonical
                    # (isoform_safe_length) even if the residue happens to match by coincidence.
                    # Tagging rather than dropping preserves the data and makes mismatches visible.
                    safe_length = parse_isoform_safe_length(uniprot)
                    mutation_entries = [
                        (mut + "(isoform?)", pos) if (
                            pos not in pos_to_aa
                            or pos_to_aa[pos] != mut[0]
                            or (safe_length is not None and pos > safe_length)
                        )
                        else (mut, pos)
                        for mut, pos in mutation_entries
                    ]

                    if MIN_PLDDT > 0:
                        before_muts = len(mutation_entries)
                        mutation_entries = [(m, p) for m, p in mutation_entries
                                            if plddt_map.get(p, 0) >= MIN_PLDDT]
                        before_ptms = len(ptm_entries)
                        ptm_entries = [(s, p, t) for s, p, t in ptm_entries
                                       if plddt_map.get(p, 0) >= MIN_PLDDT]
                        filtered = (before_muts - len(mutation_entries)) + (before_ptms - len(ptm_entries))
                        if filtered:
                            tqdm.write(f"  {uniprot}: filtered {before_muts - len(mutation_entries)} mutations, "
                                       f"{before_ptms - len(ptm_entries)} PTMs below pLDDT {MIN_PLDDT}")

                    patient_counts = parse_mutation_patient_counts(uniprot)
                    total_missense_patients = parse_total_cosmic_missense_patients(uniprot)

                    pae_matrix = load_pae_matrix(uniprot_dir)

                    for ptm_site, ptm_position, ptm_type in ptm_entries:
                        if ptm_position not in pos_to_aa:
                            write_skip(skip_writer, uniprot, gene, ptm_site, ptm_type,
                                       "position_not_in_structure",
                                       f"position {ptm_position} beyond canonical sequence length {max(pos_to_aa) if pos_to_aa else '?'}")
                            continue

                        struct_aa = pos_to_aa[ptm_position]
                        ptm_aa = ptm_site[0]
                        if struct_aa != ptm_aa:
                            write_skip(skip_writer, uniprot, gene, ptm_site, ptm_type,
                                       "residue_mismatch",
                                       f"PTMD={ptm_aa}{ptm_position} but canonical structure has {struct_aa}{ptm_position}")
                            continue

                        nearby = find_nearby_mutations(chain, ptm_position, mutation_entries, pae_matrix=pae_matrix, cutoff=DISTANCE_CUTOFF, max_pae=MAX_PAE)
                        if not nearby:
                            detail = f"no mutations within {DISTANCE_CUTOFF:.0f}A of position {ptm_position}"
                            if MAX_PAE is not None:
                                detail += f" (or all excluded by the max PAE {MAX_PAE} filter)"
                            write_skip(skip_writer, uniprot, gene, ptm_site, ptm_type,
                                       "no_nearby_mutations", detail)
                            continue
                        within_5 = [hit for hit in nearby if abs(hit["mutation_pos"] - ptm_position) <= 5]
                        beyond_5 = [hit for hit in nearby if abs(hit["mutation_pos"] - ptm_position) > 5]
                        site_key = f"{ptm_site}:{ptm_type}" if ptm_type else ptm_site
                        disrupting_set = known_disruptions.get(site_key, set())
                        confirmed = [hit for hit in nearby if hit["mutation"].replace("(isoform?)", "") in disrupting_set]
                        writer.writerow([
                            uniprot,
                            gene,
                            ptm_site,
                            ptm_type,
                            format_mutations(within_5),
                            len(within_5),
                            unique_mutation_position_count(within_5),
                            total_patient_count(within_5, patient_counts),
                            format_mutations(beyond_5),
                            len(beyond_5),
                            unique_mutation_position_count(beyond_5),
                            total_patient_count(beyond_5, patient_counts),
                            linear_distances(beyond_5, ptm_position),
                            mutation_at_ptm_site(within_5, ptm_position),
                            format_mutations(confirmed),
                            parse_ptm_diseases(uniprot, ptm_site, ptm_type),
                            total_missense_patients,
                        ])

                        if long_writer is not None:
                            ptm_plddt_val = plddt_map.get(ptm_position)
                            ptm_plddt_str = f"{ptm_plddt_val:.1f}" if ptm_plddt_val is not None else ""
                            nearby_total_pts = total_patient_count(nearby, patient_counts)
                            ptm_diseases_str = parse_ptm_diseases(uniprot, ptm_site, ptm_type)
                            m_site = SITE_RE.match(ptm_site) if ptm_site else None
                            is_st = bool(m_site and m_site.group(1) in ("S", "T"))
                            for hit in sorted(nearby, key=lambda h: (h["mutation_pos"], h["mutation"])):
                                mut_clean = hit["mutation"].replace("(isoform?)", "")
                                mut_plddt_val = plddt_map.get(hit["mutation_pos"])
                                long_writer.writerow([
                                    gene,
                                    uniprot,
                                    ptm_site,
                                    ptm_type,
                                    ptm_plddt_str,
                                    hit["mutation"],
                                    abs(hit["mutation_pos"] - ptm_position),
                                    f"{hit['distance']:.2f}",
                                    "",  # polyphen_score (filled by step 4)
                                    "",  # polyphen_class (filled by step 4)
                                    f"{mut_plddt_val:.1f}" if mut_plddt_val is not None else "",
                                    f"{hit['pae']:.1f}" if hit["pae"] is not None else "",
                                    patient_counts.get((mut_clean, hit["mutation_pos"]), 0),
                                    nearby_total_pts,
                                    total_missense_patients if total_missense_patients is not None else "",
                                    len(nearby),
                                    "yes" if mut_clean in disrupting_set else "no",
                                    ptm_diseases_str,
                                    "",  # 1433_predicted (filled by step 4)
                                    "",  # 1433_predicted_consensus (filled by step 4)
                                    "",  # 1433_confirmed (filled by step 4)
                                    "",  # kinase_predictions (filled by step 4)
                                    "",  # ptm_aiupred_general (filled by step 4)
                                    "",  # ptm_aiupred_binding (filled by step 4)
                                    "",  # ptm_is_disordered (filled by step 4)
                                    "",  # ptm_is_binding (filled by step 4)
                                    "",  # mut_aiupred_general (filled by step 4)
                                    "",  # mut_aiupred_binding (filled by step 4)
                                    "",  # mut_is_disordered (filled by step 4)
                                    "",  # mut_is_binding (filled by step 4)
                                ])

            print(f"Wrote nearby mutation data to {OUTPUT_PATH}")
            print(f"Wrote skipped PTMs to {SKIPPED_PATH}")
            print(f"Wrote long-format PTM-mutation pairs to {LONG_OUTPUT_PATH}")
        finally:
            long_handle.close()


if __name__ == "__main__":
    main()
