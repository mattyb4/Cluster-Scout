#!/usr/bin/env python3

import sys
from Bio.PDB import MMCIFParser
import numpy as np
from pathlib import Path
import argparse
import re
import csv
from typing import Any

sys.path.insert(0, str(Path(__file__).resolve().parent))
from pipeline_utils import AA3TO1  # noqa: E402

parser = MMCIFParser(QUIET=True)
SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent
MODELS_ROOT = PROJECT_ROOT / "cif_models"
DEFAULT_TSV = PROJECT_ROOT / "data" / "steps" / "PTMD_COSMIC_hotspots_by_protein.tsv"
DEFAULT_OUTPUT_DB = PROJECT_ROOT / "Output" / "ptm_mutation_proximity_db.tsv"

OUTPUT_COLUMNS_STEP3 = [
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
]

def get_ca_coord(chain, residue_number):
    for residue in chain:
        if residue.get_id()[1] == residue_number:
            if "CA" in residue:
                return residue["CA"].get_coord()
    return None

#compute distance between two 3D points
def compute_distance(coord1, coord2):
    return np.linalg.norm(coord1 - coord2)

#find mutations within cutoff distance of PTM site
def find_nearby_mutations(chain, ptm_pos, mutation_map, cutoff=10.0): #adjust cutoff as needed
    results = []

    ptm_coord = get_ca_coord(chain, ptm_pos)

    if ptm_coord is None:
        print(f"PTM residue {ptm_pos} not found in structure.")
        return results

    for mut_pos, mut_labels in mutation_map.items():
        mut_coord = get_ca_coord(chain, mut_pos)

        if mut_coord is None:
            continue

        distance = compute_distance(ptm_coord, mut_coord)

        if distance <= cutoff:
            labels = sorted(mut_labels)
            results.append({
                "mutation_pos": mut_pos,
                "mutation_label": labels[0] if labels else f"{mut_pos}",
                "distance": distance
            })

    return results

MUT_RE = re.compile(r"([A-Z])(\d+)([A-Z*])")  # e.g., R482H, S2054L
POS_RE = re.compile(r"[A-Z](\d+)")
PTM_SITE_RE = re.compile(r"([A-Z])(\d+)")


def get_uniprot_value(row):
    return row.get("UniProt") or row.get("uniprot_id")


def get_gene_value(row):
    return row.get("gene") or row.get("Gene") or ""


def extract_ptm_labels_from_list_field(field_value: str):
    ptm_map = {}
    for token in str(field_value).split(";"):
        match = PTM_SITE_RE.search(token.strip())
        if match:
            label = f"{match.group(1)}{match.group(2)}"
            pos = int(match.group(2))
            ptm_map.setdefault(pos, set()).add(label)
    return ptm_map


def extract_mutation_labels_from_field(field_value: str):
    mutation_map = {}
    for token in str(field_value).split(";"):
        match = MUT_RE.search(token.strip())
        if match:
            label = f"{match.group(1)}{match.group(2)}{match.group(3)}"
            pos = int(match.group(2))
            mutation_map.setdefault(pos, set()).add(label)
    return mutation_map

def parse_ptm_entries(uniprot, tsv_path):
    # Returns list of (ptm_site, position, ptm_type) tuples.
    entries = {}

    with tsv_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            if get_uniprot_value(row) != uniprot:
                continue

            # Legacy schema: explicit ptm_site/ptm_pos columns.
            if row.get("ptm_pos"):
                pos = int(row["ptm_pos"])
                label = row.get("ptm_site")
                if label and PTM_SITE_RE.search(str(label)):
                    site = str(label).strip()
                else:
                    site = f"PTM{pos}"
                ptm_type = (row.get("ptm_type") or "").strip()
                entries[(site, pos)] = ptm_type
                continue

            # New schema: semicolon-separated PTM tokens in ptms_on_protein.
            if row.get("ptms_on_protein"):
                for token in str(row["ptms_on_protein"]).split(";"):
                    token = token.strip()
                    if not token:
                        continue
                    if ":" in token:
                        site_part, ptm_type = token.split(":", 1)
                        ptm_type = ptm_type.strip()
                    else:
                        site_part, ptm_type = token, ""
                    m = PTM_SITE_RE.search(site_part.strip())
                    if not m:
                        continue
                    site = f"{m.group(1)}{m.group(2)}"
                    pos = int(m.group(2))
                    entries[(site, pos)] = ptm_type

    return sorted([(site, pos, ptm_type) for (site, pos), ptm_type in entries.items()], key=lambda x: x[1])

def parse_mutation_positions(ptm_pos, tsv_path, uniprot=None): 
    mutation_map = {}

    with tsv_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            row_uniprot = get_uniprot_value(row)
            if uniprot and row_uniprot != uniprot:
                continue

            # Legacy schema: PTM-specific nearby/far mutation columns.
            if row.get("ptm_pos"):
                if int(row["ptm_pos"]) != int(ptm_pos):
                    continue
                fields = [row.get("near_mutation", ""), row.get("far_mutations_prevalence_filtered", "")]
                for field in fields:
                    for token in str(field).split(","):
                        match = MUT_RE.search(token.strip())
                        if match:
                            label = f"{match.group(1)}{match.group(2)}{match.group(3)}"
                            pos = int(match.group(2))
                            mutation_map.setdefault(pos, set()).add(label)
                continue

            # New schema: protein-level mutation list in mutations_on_protein.
            # These mutations are not PTM-position-specific in this table.
            if row.get("mutations_on_protein"):
                parsed = extract_mutation_labels_from_field(row["mutations_on_protein"])
                for pos, labels in parsed.items():
                    mutation_map.setdefault(pos, set()).update(labels)

    return dict(sorted(mutation_map.items()))


def parse_isoform_safe_length(uniprot, tsv_path):
    """Return the residue position past which COSMIC numbering diverges from canonical, or None."""
    with tsv_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            if get_uniprot_value(row) == uniprot:
                value = row.get("isoform_safe_length", "")
                if value and value.strip():
                    return int(float(value))
                return None
    return None


def build_pos_to_aa(chain) -> dict[int, str]:
    """Build a {residue_number: one_letter_aa} dict from a Bio.PDB chain."""
    pos_to_aa = {}
    for residue in chain:
        res_id = residue.get_id()
        if res_id[0] != " ":
            continue
        pos = res_id[1]
        one_letter = AA3TO1.get(residue.get_resname(), "?")
        pos_to_aa[pos] = one_letter
    return pos_to_aa


def tag_isoform_mutations(mutation_map: dict, pos_to_aa: dict, safe_length) -> dict:
    """Return a new mutation_map with (isoform?) tags on mismatched mutations."""
    tagged = {}
    for pos, labels in mutation_map.items():
        new_labels = set()
        for label in labels:
            m = MUT_RE.match(label)
            if m:
                ref_aa = m.group(1)
                mismatch = (
                    pos not in pos_to_aa
                    or pos_to_aa[pos] != ref_aa
                    or (safe_length is not None and pos > safe_length)
                )
                if mismatch:
                    new_labels.add(label + "(isoform?)")
                else:
                    new_labels.add(label)
            else:
                new_labels.add(label)
        tagged[pos] = new_labels
    return tagged


def parse_gene_name(uniprot, tsv_path):
    with tsv_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            if get_uniprot_value(row) == uniprot:
                return get_gene_value(row)
    return ""


def parse_ptm_diseases(uniprot, ptm_site, ptm_type, tsv_path):
    diseases = []
    with tsv_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            if get_uniprot_value(row) != uniprot:
                continue

            # New schema supports ptm_disease_pairs directly.
            for entry in str(row.get("ptm_disease_pairs", "")).split(";"):
                entry = entry.strip()
                if " | " not in entry:
                    continue
                site_type, disease = entry.split(" | ", 1)
                if site_type.strip() == f"{ptm_site}:{ptm_type}":
                    disease = disease.strip()
                    if disease and disease not in diseases:
                        diseases.append(disease)
    return "; ".join(diseases)


def format_mutations(hits):
    if not hits:
        return ""
    parts = []
    for hit in sorted(hits, key=lambda h: (h["mutation_pos"], h["mutation_label"])):
        parts.append(f"{hit['mutation_label']}-{hit['distance']:.2f}A")
    return ", ".join(parts)


def linear_distances(hits, ptm_pos):
    if not hits:
        return ""
    seen = set()
    distances = []
    for hit in sorted(hits, key=lambda h: (h["mutation_pos"], h["mutation_label"])):
        pos = hit["mutation_pos"]
        if pos not in seen:
            seen.add(pos)
            distances.append(abs(pos - int(ptm_pos)))
    return ",".join(str(d) for d in distances)


def unique_mutation_position_count(hits):
    return len({hit["mutation_pos"] for hit in hits})


def mutation_at_ptm_site(hits, ptm_pos):
    return "yes" if any(hit["mutation_pos"] == int(ptm_pos) for hit in hits) else "no"


def _read_existing_table(path: Path):
    if not path.exists():
        return None, []
    for enc in ("utf-16", "utf-8"):
        try:
            with path.open("r", encoding=enc, newline="") as handle:
                reader = csv.DictReader(handle, delimiter="\t")
                rows = list(reader)
                return (reader.fieldnames or []), rows
        except UnicodeError:
            continue
    raise RuntimeError(f"Could not decode existing output file: {path}")


def append_rows_no_duplicates(output_path: Path, new_rows: list[dict[str, Any]]) -> int:
    output_path.parent.mkdir(parents=True, exist_ok=True)

    existing_header, existing_rows = _read_existing_table(output_path)
    write_header = existing_header if existing_header else OUTPUT_COLUMNS_STEP3

    existing_keys = set()
    for row in existing_rows:
        key = tuple(str(row.get(col, "")) for col in OUTPUT_COLUMNS_STEP3)
        existing_keys.add(key)

    rows_to_write = []
    for row in new_rows:
        key = tuple(str(row.get(col, "")) for col in OUTPUT_COLUMNS_STEP3)
        if key in existing_keys:
            continue
        existing_keys.add(key)
        out_row = {col: str(row.get(col, "")) for col in write_header}
        rows_to_write.append(out_row)

    file_exists = output_path.exists()
    with output_path.open("a", encoding="utf-16", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=write_header, delimiter="\t")
        if not file_exists:
            writer.writeheader()
        writer.writerows(rows_to_write)

    return len(rows_to_write)


def parse_args():
    ap = argparse.ArgumentParser(
        description="Find nearby mutations around PTM sites for a selected CIF model."
    )
    ap.add_argument(
        "cif_file",
        help="Path to .cif file. If relative, it is resolved under cif_models/."
    )
    ap.add_argument(
        "--uniprot",
        default=None,
        help="UniProt accession used to filter PTM rows. Defaults to the CIF parent directory name."
    )
    ap.add_argument(
        "--tsv",
        default=str(DEFAULT_TSV),
        help="Input TSV containing PTM/mutation data."
    )
    ap.add_argument(
        "--cutoff",
        type=float,
        default=10.0,
        help="Distance cutoff in Angstroms (default: 10.0)."
    )
    ap.add_argument(
        "--append-to-db",
        action="store_true",
        help="Append nearby-mutation rows to ptm_mutation_proximity_db.tsv (deduplicated)."
    )
    ap.add_argument(
        "--output-db",
        default=str(DEFAULT_OUTPUT_DB),
        help="Output DB path used with --append-to-db."
    )
    return ap.parse_args()


def resolve_cif_path(cif_arg: str) -> Path:
    candidate = Path(cif_arg)
    if candidate.is_absolute():
        resolved = candidate
    else:
        resolved = MODELS_ROOT / candidate
    return resolved.resolve()


def main():
    args = parse_args()

    model_file = resolve_cif_path(args.cif_file)
    if not model_file.exists():
        raise FileNotFoundError(f"CIF file not found: {model_file}")

    tsv_path = Path(args.tsv).resolve()
    if not tsv_path.exists():
        raise FileNotFoundError(f"TSV file not found: {tsv_path}")

    target_uniprot = args.uniprot or model_file.parent.name
    gene = parse_gene_name(target_uniprot, tsv_path)

    structure = parser.get_structure("protein", model_file)
    model = structure[0]   # First (and only) model
    chain = list(model.get_chains())[0]  # AlphaFold usually has one chain

    pos_to_aa = build_pos_to_aa(chain)
    safe_length = parse_isoform_safe_length(target_uniprot, tsv_path)

    ptm_entries = parse_ptm_entries(target_uniprot, tsv_path)
    if not ptm_entries:
        raise ValueError(f"No PTM positions found in TSV for UniProt {target_uniprot}.")

    rows_for_db = []

    for ptm_label, ptm_position, ptm_type in ptm_entries:
        mutation_positions = parse_mutation_positions(ptm_position, tsv_path, uniprot=target_uniprot)
        mutation_positions = tag_isoform_mutations(mutation_positions, pos_to_aa, safe_length)
        nearby = find_nearby_mutations(chain, ptm_position, mutation_positions, cutoff=args.cutoff)

        print(f"\nPTM {ptm_label} ({target_uniprot}):")
        if nearby:
            for hit in nearby:
                print(f"  {hit['mutation_label']} is {hit['distance']:.2f} Å away")

            within_5 = [hit for hit in nearby if abs(hit["mutation_pos"] - ptm_position) <= 5]
            beyond_5 = [hit for hit in nearby if abs(hit["mutation_pos"] - ptm_position) > 5]
            rows_for_db.append({
                "UniProt": target_uniprot,
                "gene": gene,
                "ptm_site": ptm_label,
                "ptm_type": ptm_type,
                "mutations_within_5_positions": format_mutations(within_5),
                "mutation_count_within_5_positions": str(len(within_5)),
                "unique_mutation_position_count_within_5_positions": str(unique_mutation_position_count(within_5)),
                "mutations_more_than_5_positions": format_mutations(beyond_5),
                "mutation_count_more_than_5_positions": str(len(beyond_5)),
                "unique_mutation_position_count_more_than_5_positions": str(unique_mutation_position_count(beyond_5)),
                "morethan5_linear_distance": linear_distances(beyond_5, ptm_position),
                "mutation_at_ptm_site": mutation_at_ptm_site(within_5, ptm_position),
                "ptm_diseases": parse_ptm_diseases(target_uniprot, ptm_label, ptm_type, tsv_path),
            })
        else:
            print("  No nearby mutations found within cutoff distance.")

    if args.append_to_db:
        added = append_rows_no_duplicates(Path(args.output_db).resolve(), rows_for_db)
        print(f"\nAdded {added} new row(s) to {Path(args.output_db).resolve()}")


if __name__ == "__main__":
    main()
