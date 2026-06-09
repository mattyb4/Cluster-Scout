import numpy as np
from pathlib import Path
import re
import csv
import json
import argparse
from typing import Any
from tqdm import tqdm
from biotite.structure.io.pdbx import CIFFile, get_structure  # type: ignore[import-untyped]

SCRIPT_DIR = Path(__file__).parent
PROJECT_ROOT = SCRIPT_DIR.parent
MODELS_ROOT = PROJECT_ROOT / "cif_models"
OUTPUT_PATH = PROJECT_ROOT / "Output" / "ptm_mutation_proximity_db.tsv"
SKIPPED_PATH = PROJECT_ROOT / "Output" / "logs" / "ptm_skipped.tsv"
CLUSTER_OUTPUT_PATH = PROJECT_ROOT / "Output" / "mutation_cluster_db.tsv"
CLUSTER_SKIPPED_PATH = PROJECT_ROOT / "Output" / "logs" / "mutation_cluster_skipped.tsv"

OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)
SKIPPED_PATH.parent.mkdir(parents=True, exist_ok=True)
PTM_TSV_PATH = PROJECT_ROOT / "data" / "steps" / "PTMD_TCGA_hotspots_by_protein.tsv"

_PTM_ROWS: list[dict[str, Any]] | None = None


def get_ptm_rows():
    global _PTM_ROWS
    if _PTM_ROWS is None:
        with PTM_TSV_PATH.open("r", encoding="utf-8", newline="") as handle:
            _PTM_ROWS = list(csv.DictReader(handle, delimiter="\t"))
    return _PTM_ROWS

def get_ca_coord(chain, residue_number):
    mask = (chain.res_id == residue_number) & (chain.atom_name == "CA")
    if not np.any(mask):
        return None
    return chain.coord[mask][0]

#compute distance between two 3D points
def compute_distance(coord1, coord2):
    return np.linalg.norm(coord1 - coord2)

#find mutations within cutoff distance of PTM site
def find_nearby_mutations(chain, ptm_pos, mutation_entries, pae_matrix=None, cutoff=10.0): #adjust cutoff as needed
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
            results.append({
                "mutation": mutation,
                "mutation_pos": mut_pos,
                "distance": distance,
                "pae": pae,
            })

    return results


def find_mutation_clusters(chain, mutation_entries, pae_matrix=None, cutoff=10.0):
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
                nearby.append({
                    "mutation": other_mut,
                    "mutation_pos": other_pos,
                    "distance": distance,
                    "pae": pae,
                })

        if nearby:
            results[(anchor_mut, anchor_pos)] = nearby

    return results


#extracts PTM positions/types from input db
PTM_RE = re.compile(r"([A-Z])(\d+)")  # e.g., S557


def parse_ptm_entries(uniprot):
    # Returns list of (ptm_site, position, ptm_type) tuples.
    # ptms_on_protein tokens are formatted as "S516:Phosphorylation".
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
    for row in get_ptm_rows():
        if row.get("uniprot_id") == uniprot:
            return row.get("gene", "")

    return ""

MUT_RE = re.compile(r"([A-Z])(\d+)([A-Z*])")  # e.g., R482H
#extracts mutation positions from input db and puts them into a list
def parse_mutation_positions(uniprot=None):
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


def find_model_file(uniprot_dir):
    uid = uniprot_dir.name
    # Only match canonical files: AF-{acc}-F{N}-model_v{ver}.cif
    # Isoform files look like AF-{acc}-{M}-F1-model_v{ver}.cif and are excluded.
    canonical_re = re.compile(rf"^AF-{re.escape(uid)}-F\d+-model_v\d+\.", re.IGNORECASE)
    candidates = [p for p in sorted(uniprot_dir.glob("*.cif")) if canonical_re.match(p.name)]
    return candidates[0] if candidates else None


def load_pae_matrix(uniprot_dir):
    uid = uniprot_dir.name
    canonical_re = re.compile(rf"^AF-{re.escape(uid)}-F\d+-predicted_aligned_error_v\d+\.", re.IGNORECASE)
    candidates = [p for p in sorted(uniprot_dir.glob("*.json")) if canonical_re.match(p.name)]
    if not candidates:
        return None
    with candidates[0].open() as f:
        data = json.load(f)
    if isinstance(data, list):
        data = data[0]
    matrix = data.get("predicted_aligned_error")
    return np.array(matrix) if matrix else None


def load_first_chain(model_file):
    try:
        cif = CIFFile.read(str(model_file))
        structure = get_structure(cif, model=1)
    except Exception as exc:
        tqdm.write(f"  Skipping {model_file.name}: failed to parse ({exc})")
        return None

    if structure is None or len(structure) == 0:
        tqdm.write(f"  Skipping {model_file.name}: no models found")
        return None

    chain_ids = list(dict.fromkeys(structure.chain_id))
    if not chain_ids:
        tqdm.write(f"  Skipping {model_file.name}: no chains found")
        return None

    chain_id = chain_ids[0]
    return structure[structure.chain_id == chain_id]


def format_mutations(hits):
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
    return len({hit["mutation_pos"] for hit in hits})

def mutation_at_ptm_site(hits, ptm_pos):
    return 'yes' if any(hit['mutation_pos'] == int(ptm_pos) for hit in hits) else 'no'


def parse_ptm_diseases(uniprot, ptm_site, ptm_type):
    # ptm_disease_pairs format: "S516:Phosphorylation | Bladder cancer; K43:Ubiquitination | Lung adenocarcinoma"
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

AA3TO1 = {
    "ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C","GLN":"Q","GLU":"E",
    "GLY":"G","HIS":"H","ILE":"I","LEU":"L","LYS":"K","MET":"M","PHE":"F",
    "PRO":"P","SER":"S","THR":"T","TRP":"W","TYR":"Y","VAL":"V","SEC":"U","PYL":"O",
}

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
args = parser.parse_args()

SKIP_HEADER = ["UniProt", "gene", "ptm_site", "ptm_type", "skip_reason", "detail"]

def write_skips(skip_writer, uniprot, gene, ptm_entries, reason, detail):
    for ptm_site, _, ptm_type in ptm_entries:
        skip_writer.writerow([uniprot, gene, ptm_site, ptm_type, reason, detail])

def write_skip(skip_writer, uniprot, gene, ptm_site, ptm_type, reason, detail):
    skip_writer.writerow([uniprot, gene, ptm_site, ptm_type, reason, detail])


dirs_present = {d.name for d in MODELS_ROOT.iterdir() if d.is_dir()}

# ── Mutation-clustering mode ──────────────────────────────────────────────────
if args.mode == "mutation-clustering":
    CLUSTER_OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    CLUSTER_SKIPPED_PATH.parent.mkdir(parents=True, exist_ok=True)

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

            pae_matrix = load_pae_matrix(uniprot_dir)
            clusters = find_mutation_clusters(chain, mutation_entries, pae_matrix=pae_matrix)

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
            "mutations_more_than_5_positions",
            "mutation_count_more_than_5_positions",
            "unique_mutation_position_count_more_than_5_positions",
            "morethan5_linear_distance",
            "mutation_at_ptm_site",
            "ptm_diseases",
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
                continue
            mutation_entries = parse_mutation_positions(uniprot=uniprot)
            if not mutation_entries:
                continue

            model_file = find_model_file(uniprot_dir)
            if model_file is None:
                write_skips(skip_writer, uniprot, gene, ptm_entries, "no_canonical_cif",
                            "AFDB has only isoform models, no canonical sequence model")
                tqdm.write(f"  {uniprot}: no canonical CIF file found")
                continue

            chain = load_first_chain(model_file)
            if chain is None:
                continue

            # Build residue -> 1-letter AA map for mismatch checking
            pos_to_aa: dict[int, str] = {}
            for atom in chain:
                if atom.res_id not in pos_to_aa:
                    pos_to_aa[atom.res_id] = AA3TO1.get(atom.res_name, "?")

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

                nearby = find_nearby_mutations(chain, ptm_position, mutation_entries, pae_matrix=pae_matrix)
                if not nearby:
                    continue
                within_5 = [hit for hit in nearby if abs(hit["mutation_pos"] - ptm_position) <= 5]
                beyond_5 = [hit for hit in nearby if abs(hit["mutation_pos"] - ptm_position) > 5]
                writer.writerow([
                    uniprot,
                    gene,
                    ptm_site,
                    ptm_type,
                    format_mutations(within_5),
                    len(within_5),
                    unique_mutation_position_count(within_5),
                    format_mutations(beyond_5),
                    len(beyond_5),
                    unique_mutation_position_count(beyond_5),
                    linear_distances(beyond_5, ptm_position),
                    mutation_at_ptm_site(within_5, ptm_position),
                    parse_ptm_diseases(uniprot, ptm_site, ptm_type),
                ])

    print(f"Wrote nearby mutation data to {OUTPUT_PATH}")
    print(f"Wrote skipped PTMs to {SKIPPED_PATH}")
