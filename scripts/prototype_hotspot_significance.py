"""PROTOTYPE: permutation-test + FDR significance for PTM-mutation 3D proximity.

Not part of the production pipeline (steps 1-4) and not wired into the app —
this is a standalone demo to evaluate an alternative to the fixed-cutoff
PTM-proximity filter used in scripts/3_find_nearby_mutations.py. It reuses
the same intermediate data (data/steps/PTMD_COSMIC_hotspots_by_protein.tsv)
and downloaded CIF models the real pipeline already produced.

Statistical method and full citations: see the "Structural-hotspot
significance testing" section of scripts/pipeline_utils.py. In short: for
each PTM site, this asks how often a random placement of the same number of
mutations observed in that protein would land at least as many mutations
within `--cutoff` Angstroms of the site as were actually observed — a
permutation-test null adapted from the HotMAPS method (Tokheim et al.,
Cancer Research 2016) — then FDR-corrects across all sites tested
(Benjamini & Hochberg, 1995).

Usage:
    uv run scripts/prototype_hotspot_significance.py [--limit N] [--permutations K]
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from tqdm import tqdm

sys.path.insert(0, str(Path(__file__).resolve().parent))
from pipeline_utils import (  # noqa: E402
    project_root,
    find_canonical_cif, load_first_chain,
    sample_permutation_indices, permutation_pvalue, benjamini_hochberg,
)

PROJECT_ROOT = project_root(__file__)
MODELS_ROOT = PROJECT_ROOT / "cif_models"
PTM_TSV_PATH = PROJECT_ROOT / "data" / "steps" / "PTMD_COSMIC_hotspots_by_protein.tsv"
DISTANCE_CUTOFF = 10.0

MUT_RE = re.compile(r"([A-Z])(\d+)([A-Z*])")
PTM_RE = re.compile(r"([A-Z])(\d+)")


def parse_mutation_positions(row) -> set[int]:
    """Distinct hotspot-mutated residue positions for one protein row."""
    positions = set()
    for token in re.split(r"[;,]", row.get("mutations_on_protein", "") or ""):
        match = MUT_RE.search(token.strip())
        if match:
            positions.add(int(match.group(2)))
    return positions


def parse_ptm_entries(row) -> dict[int, str]:
    """{position: site_label} for one protein row, e.g. {557: 'S557'}."""
    entries: dict[int, str] = {}
    for token in re.split(r";", row.get("ptms_on_protein", "") or ""):
        token = token.strip()
        site_part = token.split(":", 1)[0] if ":" in token else token
        match = PTM_RE.search(site_part.strip())
        if match:
            entries[int(match.group(2))] = f"{match.group(1)}{match.group(2)}"
    return entries


def get_ca_coords(chain) -> tuple[np.ndarray, np.ndarray]:
    """(positions, coords) for every structurally-resolved residue in chain —
    the "universe" the permutation null samples random mutation sets from."""
    ca_mask = chain.atom_name == "CA"
    ca = chain[ca_mask]
    return np.asarray(ca.res_id, dtype=int), np.asarray(ca.coord, dtype=float)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--limit", type=int, default=25,
                     help="Number of proteins to test (default: 25, for a quick look)")
    ap.add_argument("--permutations", type=int, default=5000)
    ap.add_argument("--cutoff", type=float, default=DISTANCE_CUTOFF)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out", default=str(PROJECT_ROOT / "Output" / "prototype_hotspot_significance.tsv"))
    args = ap.parse_args()

    rng = np.random.default_rng(args.seed)

    df = pd.read_csv(PTM_TSV_PATH, sep="\t", dtype=str)
    uniprot_ids = df["uniprot_id"].dropna().unique().tolist()
    if args.limit:
        uniprot_ids = uniprot_ids[:args.limit]

    rows = []
    for uniprot in tqdm(uniprot_ids, desc="Testing proteins"):
        row = df[df["uniprot_id"] == uniprot].iloc[0]
        gene = row.get("gene", "")
        mutated_positions = parse_mutation_positions(row)
        ptm_positions = parse_ptm_entries(row)
        if not mutated_positions or not ptm_positions:
            continue

        model_file = find_canonical_cif(MODELS_ROOT / uniprot)
        if model_file is None:
            continue
        chain = load_first_chain(model_file)
        if chain is None:
            continue

        positions, coords = get_ca_coords(chain)
        pos_to_idx = {int(p): i for i, p in enumerate(positions)}

        # Only positions actually resolved in this structure could have been
        # "placed" anywhere by the permutation null, so they define
        # n_mutations for this protein.
        resolved_mut_positions = [p for p in mutated_positions if p in pos_to_idx]
        n_mutations = len(resolved_mut_positions)
        if n_mutations == 0:
            continue

        # Shared across every PTM site in this protein: n_mutations doesn't
        # vary by site, only which residues fall within cutoff of each site.
        sampled_idx = sample_permutation_indices(
            n_residues=len(positions), n_mutations=n_mutations,
            n_permutations=args.permutations, rng=rng,
        )

        for ptm_pos, ptm_site in ptm_positions.items():
            if ptm_pos not in pos_to_idx:
                continue
            site_coord = coords[pos_to_idx[ptm_pos]]

            observed_count = sum(
                1 for p in resolved_mut_positions
                if np.linalg.norm(coords[pos_to_idx[p]] - site_coord) <= args.cutoff
            )

            p_value, _null_counts = permutation_pvalue(
                site_coord=site_coord, residue_coords=coords,
                observed_count=observed_count, cutoff=args.cutoff,
                sampled_idx=sampled_idx,
            )

            rows.append({
                "uniprot_id": uniprot,
                "gene": gene,
                "ptm_site": ptm_site,
                "n_residues_resolved": len(positions),
                "n_mutations_in_protein": n_mutations,
                "observed_count_within_cutoff": observed_count,
                "passes_fixed_cutoff": observed_count > 0,
                "permutation_p_value": p_value,
            })

    if not rows:
        print("No PTM sites tested (no proteins with both structure and hotspot mutations in --limit).")
        return

    result = pd.DataFrame(rows)
    result["fdr_q_value"] = benjamini_hochberg(result["permutation_p_value"].to_numpy())
    result["significant_q01"] = result["fdr_q_value"] <= 0.01
    result = result.sort_values("permutation_p_value")

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(out_path, sep="\t", index=False)

    n_sites = len(result)
    n_proteins = result["uniprot_id"].nunique()
    n_cutoff = int(result["passes_fixed_cutoff"].sum())
    n_sig = int(result["significant_q01"].sum())
    agree = result[result["passes_fixed_cutoff"] & result["significant_q01"]]
    cutoff_only = result[result["passes_fixed_cutoff"] & ~result["significant_q01"]]

    print(f"\nTested {n_sites} PTM sites across {n_proteins} proteins.")
    print(f"Passes fixed {args.cutoff:.0f}A cutoff (>=1 nearby mutation): {n_cutoff}")
    print(f"Significant at FDR q<=0.01 (permutation test): {n_sig}")
    print(f"\nAgree (near AND significant): {len(agree)}")
    print(f"Fixed cutoff says near, but NOT significant after FDR correction: {len(cutoff_only)}")
    print(f"\nSaved full results to {out_path}")
    print("\nTop 10 most significant PTM sites:")
    print(result.head(10).to_string(index=False))


if __name__ == "__main__":
    main()
