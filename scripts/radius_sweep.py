"""Parameter sweep: test different distance cutoffs to find the optimal radius.

For a set of proteins, computes how many mutations fall within each tested
radius (4-20 Å in 1 Å steps) of each PTM site, averages across PTM sites
per protein, and plots the result.

Usage:
    uv run scripts/radius_sweep.py
    uv run scripts/radius_sweep.py --genes EGFR TP53 VHL
    uv run scripts/radius_sweep.py --radii 4 25 1
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from pipeline_utils import (  # noqa: E402
    project_root, find_canonical_cif, load_first_chain,
)

PROJECT_ROOT = project_root(__file__)
MODELS_ROOT = PROJECT_ROOT / "cif_models"
PTM_TSV = PROJECT_ROOT / "data" / "steps" / "PTMD_TCGA_hotspots_by_protein.tsv"

DEFAULT_GENES = ["EGFR", "TP53", "VHL", "CANT1", "DDR2", "PTPN11", "LZTR1", "CDK12"]


def get_ca_coord(chain, residue_number):
    """Return the CA coordinate for a residue, or None."""
    mask = (chain.res_id == residue_number) & (chain.atom_name == "CA")
    if not np.any(mask):
        return None
    return chain.coord[mask][0]


def load_protein_data(gene: str, df: pd.DataFrame, chain):
    """Extract PTM positions and mutation positions with their CA coordinates.

    Returns (ptm_coords, mutation_coords) where each is a list of (position, coord).
    """
    rows = df[df["gene"] == gene]
    if rows.empty:
        return [], []

    row = rows.iloc[0]

    # Parse PTM positions
    ptm_coords = []
    for token in str(row.get("ptms_on_protein", "")).split(";"):
        token = token.strip()
        if not token:
            continue
        import re
        m = re.search(r"([A-Z])(\d+)", token)
        if m:
            pos = int(m.group(2))
            coord = get_ca_coord(chain, pos)
            if coord is not None:
                ptm_coords.append((pos, coord))

    # Parse mutation positions
    mutation_coords = []
    seen_positions = set()
    for token in str(row.get("mutations_on_protein", "")).split(";"):
        token = token.strip()
        import re
        m = re.search(r"([A-Z])(\d+)([A-Z*])", token)
        if m:
            pos = int(m.group(2))
            if pos not in seen_positions:
                coord = get_ca_coord(chain, pos)
                if coord is not None:
                    mutation_coords.append((pos, coord))
                    seen_positions.add(pos)

    return ptm_coords, mutation_coords


def get_all_ca_coords(chain) -> list[tuple[int, np.ndarray]]:
    """Return all CA (position, coord) pairs from the chain."""
    ca_mask = chain.atom_name == "CA"
    ca_atoms = chain[ca_mask]
    return [(int(ca_atoms.res_id[i]), ca_atoms.coord[i]) for i in range(len(ca_atoms))]


def sweep_radii(ptm_coords, mutation_coords, radii):
    """For each PTM, compute distances to all mutations once, then count at each radius.

    Returns a dict {radius: average_mutation_count_per_ptm}.
    """
    if not ptm_coords or not mutation_coords:
        return {r: 0.0 for r in radii}

    mut_coord_array = np.array([c for _, c in mutation_coords])

    counts_per_radius = {r: [] for r in radii}

    for _ptm_pos, ptm_coord in ptm_coords:
        distances = np.linalg.norm(mut_coord_array - ptm_coord, axis=1)
        for r in radii:
            counts_per_radius[r].append(int(np.sum(distances <= r)))

    return {r: np.mean(counts) for r, counts in counts_per_radius.items()}


def random_baseline(ptm_coords, all_ca_coords, n_mutations, radii, n_permutations=100):
    """Compute the expected sweep curve if mutations were randomly distributed across the protein."""
    if not ptm_coords or n_mutations == 0 or len(all_ca_coords) < n_mutations:
        return {r: 0.0 for r in radii}

    all_coords = np.array([c for _, c in all_ca_coords])
    ptm_coord_array = np.array([c for _, c in ptm_coords])
    n_ptms = len(ptm_coords)

    totals = {r: 0.0 for r in radii}
    for _ in range(n_permutations):
        indices = np.random.choice(len(all_coords), size=n_mutations, replace=False)
        random_coords = all_coords[indices]
        for pi in range(n_ptms):
            distances = np.linalg.norm(random_coords - ptm_coord_array[pi], axis=1)
            for r in radii:
                totals[r] += int(np.sum(distances <= r))

    return {r: v / (n_permutations * n_ptms) for r, v in totals.items()}


def main():
    parser = argparse.ArgumentParser(
        description="Sweep distance cutoffs to find the optimal mutation-capture radius."
    )
    parser.add_argument(
        "--genes", nargs="+", default=DEFAULT_GENES,
        help=f"Gene symbols to test (default: {' '.join(DEFAULT_GENES)})",
    )
    parser.add_argument(
        "--radii", nargs=3, type=float, default=[4, 20, 1],
        metavar=("START", "STOP", "STEP"),
        help="Radius range as start stop step (default: 4 20 1)",
    )
    parser.add_argument(
        "--output", default=str(PROJECT_ROOT / "Output" / "radius_sweep.png"),
        help="Output plot path (default: Output/radius_sweep.png)",
    )
    args = parser.parse_args()

    radii = list(np.arange(args.radii[0], args.radii[1] + args.radii[2] / 2, args.radii[2]))
    print(f"Testing radii: {radii[0]:.0f}-{radii[-1]:.0f} A in {args.radii[2]:.0f} A steps")

    if not PTM_TSV.exists():
        sys.exit(f"Error: intermediate TSV not found at {PTM_TSV}\n"
                 "Run the pipeline (step 1) first to generate it.")

    df = pd.read_csv(PTM_TSV, sep="\t", dtype=str, keep_default_na=False)

    # Map gene names to UniProt IDs
    gene_to_uid = {}
    for g in args.genes:
        rows = df[df["gene"] == g]
        if rows.empty:
            print(f"  Warning: {g} not found in dataset, skipping")
            continue
        gene_to_uid[g] = rows.iloc[0]["uniprot_id"]

    results = []

    for gene, uid in gene_to_uid.items():
        print(f"\n{gene} ({uid}):")

        cif_dir = MODELS_ROOT / uid
        cif_file = find_canonical_cif(cif_dir) if cif_dir.is_dir() else None
        if cif_file is None:
            print(f"  No CIF file found, skipping")
            continue

        chain = load_first_chain(cif_file)
        if chain is None:
            print(f"  Could not parse CIF, skipping")
            continue

        all_ca = get_all_ca_coords(chain)
        protein_length = len(all_ca)

        ptm_coords, mutation_coords = load_protein_data(gene, df, chain)
        print(f"  {len(ptm_coords)} PTM sites, {len(mutation_coords)} unique mutation positions, "
              f"{protein_length} residues")

        if not ptm_coords:
            print(f"  No PTM coordinates found, skipping")
            continue

        avg_counts = sweep_radii(ptm_coords, mutation_coords, radii)
        print("  Computing random baseline (100 permutations)...")
        baseline = random_baseline(ptm_coords, all_ca, len(mutation_coords), radii)

        for r in radii:
            results.append({
                "protein": gene,
                "radius": r,
                "avg_mutation_count": avg_counts[r],
                "random_baseline": baseline[r],
                "avg_normalized": avg_counts[r] / protein_length * 1000,
                "random_normalized": baseline[r] / protein_length * 1000,
                "protein_length": protein_length,
            })
        print(f"  Sweep complete: {avg_counts[radii[0]]:.1f} avg at {radii[0]:.0f}A "
              f"-> {avg_counts[radii[-1]]:.1f} avg at {radii[-1]:.0f}A")

    if not results:
        sys.exit("No data collected. Check that CIF files are downloaded for the target genes.")

    result_df = pd.DataFrame(results)

    # Elbow detection
    from kneed import KneeLocator

    print("\n-- Elbow Detection --")
    elbows: dict[str, float | None] = {}
    for gene in result_df["protein"].unique():
        gene_data = result_df[result_df["protein"] == gene].sort_values("radius")
        try:
            kn = KneeLocator(
                gene_data["radius"].values,
                gene_data["avg_mutation_count"].values,
                curve="convex",
                direction="increasing",
                interp_method="interp1d",
            )
            elbows[gene] = kn.knee
            if kn.knee is not None:
                print(f"  {gene}: optimal radius = {kn.knee:.0f} A")
            else:
                print(f"  {gene}: no elbow detected")
        except Exception:
            elbows[gene] = None
            print(f"  {gene}: could not compute elbow")

    detected = [v for v in elbows.values() if v is not None]
    if detected:
        avg_elbow = np.mean(detected)
        print(f"\n  Average optimal radius across proteins: {avg_elbow:.1f} A")

    # Save data
    tsv_path = Path(args.output).with_suffix(".tsv")
    result_df.to_csv(tsv_path, sep="\t", index=False)
    print(f"\nData saved to: {tsv_path}")

    # Plot: 2-panel figure
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(11, 10), sharex=True)

    proteins = result_df["protein"].unique()
    colors = {}

    # Top panel: raw counts + random baselines
    for gene in proteins:
        gene_data = result_df[result_df["protein"] == gene]
        line, = ax1.plot(gene_data["radius"], gene_data["avg_mutation_count"],
                         marker="o", markersize=4, linewidth=2, label=gene)
        colors[gene] = line.get_color()

        # Random baseline as dashed line (same color)
        ax1.plot(gene_data["radius"], gene_data["random_baseline"],
                 linestyle="--", linewidth=1, alpha=0.5, color=colors[gene])

        # Elbow diamond
        elbow_r = elbows.get(gene)
        if elbow_r is not None:
            elbow_row = gene_data[gene_data["radius"] == elbow_r]
            if not elbow_row.empty:
                ax1.scatter([elbow_r], [elbow_row.iloc[0]["avg_mutation_count"]],
                            s=120, zorder=6, marker="D", color=colors[gene],
                            edgecolors="black", linewidths=1.0)

    ax1.axvline(x=10, color="gray", linestyle="--", alpha=0.5, label="Current cutoff (10 A)")
    ax1.set_ylabel("Avg Mutations per PTM Site", fontsize=13)
    ax1.set_title("Mutation Capture vs. Distance Cutoff", fontsize=15)
    ax1.legend(loc="upper left", fontsize=9, ncol=2)
    ax1.grid(True, alpha=0.3)

    # Add a single label for the dashed random baselines
    ax1.plot([], [], linestyle="--", color="gray", alpha=0.5, label="Random baseline")
    ax1.legend(loc="upper left", fontsize=9, ncol=2)

    if detected:
        ax1.axvline(x=avg_elbow, color="red", linestyle=":", alpha=0.6)
        ax1.text(avg_elbow + 0.3, ax1.get_ylim()[1] * 0.95,
                 f"Avg elbow: {avg_elbow:.1f} A",
                 color="red", fontsize=10, va="top")

    # Bottom panel: normalized by protein length (per 1000 residues)
    for gene in proteins:
        gene_data = result_df[result_df["protein"] == gene]
        plen = int(gene_data.iloc[0]["protein_length"])
        ax2.plot(gene_data["radius"], gene_data["avg_normalized"],
                 marker="o", markersize=4, linewidth=2,
                 label=f"{gene} ({plen} aa)", color=colors[gene])
        ax2.plot(gene_data["radius"], gene_data["random_normalized"],
                 linestyle="--", linewidth=1, alpha=0.5, color=colors[gene])

    ax2.axvline(x=10, color="gray", linestyle="--", alpha=0.5)
    ax2.set_xlabel("Radius (A)", fontsize=13)
    ax2.set_ylabel("Avg Mutations per PTM (per 1000 residues)", fontsize=13)
    ax2.set_title("Size-Normalized Mutation Capture", fontsize=15)
    ax2.legend(loc="upper left", fontsize=9, ncol=2)
    ax2.grid(True, alpha=0.3)

    if detected:
        ax2.axvline(x=avg_elbow, color="red", linestyle=":", alpha=0.6)

    plt.tight_layout()
    plt.savefig(args.output, dpi=150)
    print(f"Plot saved to: {args.output}")
    plt.show()


if __name__ == "__main__":
    main()
