"""Analyze structural variance across multiple AlphaFold CIF predictions.

Compares multiple CIF files for the same protein (e.g., different AlphaFold
seeds or model versions) by aligning structures and computing per-residue
positional variance and pLDDT comparison.

Place CIF files in data/cif_comparison/ and run:
    uv run scripts/cif_variance.py
    uv run scripts/cif_variance.py --input-dir data/cif_comparison --top 20

Output (in Output/cif_variance/):
    - variance_plot.png: per-residue variance and pLDDT with PTM/mutation markers
    - variance_data.tsv: per-residue stats with PTM/mutation annotations
    - pairwise_rmsd.tsv: RMSD matrix between all structure pairs
"""
from __future__ import annotations

import argparse
import itertools
import re
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from Bio.PDB import MMCIFParser, Superimposer

sys.path.insert(0, str(Path(__file__).resolve().parent))
from pipeline_utils import project_root, AA3TO1, extract_uniprot_from_cif  # noqa: E402

PROJECT_ROOT = project_root(__file__)
DEFAULT_INPUT_DIR = PROJECT_ROOT / "data" / "cif_comparison"
DEFAULT_OUTPUT_DIR = PROJECT_ROOT / "Output" / "cif_variance"
PTM_TSV = PROJECT_ROOT / "data" / "steps" / "PTMD_TCGA_hotspots_by_protein.tsv"

_parser = MMCIFParser(QUIET=True)


def load_ca_data(cif_path: Path) -> tuple[list[int], np.ndarray, np.ndarray, list[str]]:
    """Extract CA atom positions, coordinates, pLDDT scores, and residue names from a CIF.

    Returns (positions, coords, plddts, residue_names) where coords is (N, 3).
    """
    structure = _parser.get_structure(cif_path.stem, str(cif_path))
    chain = list(structure[0].get_chains())[0]

    positions, coords, plddts, res_names = [], [], [], []
    for residue in chain.get_residues():
        if residue.get_id()[0] != " ":
            continue
        if "CA" not in residue:
            continue
        ca = residue["CA"]
        positions.append(residue.get_id()[1])
        coords.append(ca.get_vector().get_array())
        plddts.append(ca.get_bfactor())
        res_names.append(AA3TO1.get(residue.get_resname(), "X"))

    return positions, np.array(coords), np.array(plddts), res_names


def align_to_reference(ref_coords: np.ndarray, mobile_coords: np.ndarray,
                       ref_positions: list[int], mobile_positions: list[int]
                       ) -> tuple[np.ndarray, float]:
    """Align mobile_coords onto ref_coords using shared positions. Returns (aligned_coords, rmsd)."""
    shared = set(ref_positions) & set(mobile_positions)
    if len(shared) < 3:
        return mobile_coords, float("inf")

    ref_idx = [ref_positions.index(p) for p in sorted(shared)]
    mob_idx = [mobile_positions.index(p) for p in sorted(shared)]

    from Bio.PDB.Atom import Atom
    from Bio.PDB.Residue import Residue

    ref_atoms = []
    mob_atoms = []
    for ri, mi in zip(ref_idx, mob_idx):
        ra = Atom("CA", ref_coords[ri], 0, 0, " ", "CA", 0)
        ma = Atom("CA", mobile_coords[mi], 0, 0, " ", "CA", 0)
        ref_atoms.append(ra)
        mob_atoms.append(ma)

    sup = Superimposer()
    sup.set_atoms(ref_atoms, mob_atoms)

    # Apply rotation/translation to ALL mobile coordinates
    all_mob_atoms = []
    for i in range(len(mobile_coords)):
        a = Atom("CA", mobile_coords[i], 0, 0, " ", "CA", 0)
        all_mob_atoms.append(a)
    sup.apply(all_mob_atoms)

    aligned = np.array([a.get_vector().get_array() for a in all_mob_atoms])
    return aligned, sup.rms


def iterative_average_alignment(
    all_coords: list[np.ndarray],
    all_positions: list[list[int]],
    max_iterations: int = 10,
    convergence: float = 1e-4,
) -> list[np.ndarray]:
    """Align all structures to an iteratively refined average reference.

    1. Align everything to the first structure.
    2. Compute the average of the aligned coordinates.
    3. Re-align everything to that average.
    4. Repeat until the average stops moving (converges).
    """
    shared = sorted(set.intersection(*(set(p) for p in all_positions)))
    n_structs = len(all_coords)

    # Index maps for extracting shared positions from each structure
    idx_maps = [{p: j for j, p in enumerate(positions)} for positions in all_positions]

    def extract_shared(coords, idx_map):
        return np.array([coords[idx_map[p]] for p in shared])

    # Initial alignment: everything to structure 0
    aligned = [c.copy() for c in all_coords]
    for i in range(1, n_structs):
        aligned[i], _ = align_to_reference(aligned[0], aligned[i],
                                           all_positions[0], all_positions[i])

    for iteration in range(max_iterations):
        # Compute average of shared positions
        shared_stack = np.array([extract_shared(aligned[i], idx_maps[i]) for i in range(n_structs)])
        avg_coords = shared_stack.mean(axis=0)

        # Re-align each structure to the average
        new_aligned = []
        for i in range(n_structs):
            shared_mobile = extract_shared(aligned[i], idx_maps[i])

            from Bio.PDB.Atom import Atom
            ref_atoms = [Atom("CA", avg_coords[j], 0, 0, " ", "CA", 0) for j in range(len(shared))]
            mob_atoms = [Atom("CA", shared_mobile[j], 0, 0, " ", "CA", 0) for j in range(len(shared))]

            sup = Superimposer()
            sup.set_atoms(ref_atoms, mob_atoms)

            all_mob = [Atom("CA", aligned[i][j], 0, 0, " ", "CA", 0) for j in range(len(aligned[i]))]
            sup.apply(all_mob)
            new_aligned.append(np.array([a.get_vector().get_array() for a in all_mob]))

        # Check convergence: has the average moved?
        new_shared = np.array([extract_shared(new_aligned[i], idx_maps[i]) for i in range(n_structs)])
        new_avg = new_shared.mean(axis=0)
        shift = np.sqrt(((new_avg - avg_coords) ** 2).sum(axis=1).mean())

        aligned = new_aligned
        if shift < convergence:
            print(f"  Converged after {iteration + 1} iteration(s) (shift={shift:.6f} A)")
            break
    else:
        print(f"  Reached max iterations ({max_iterations}), shift={shift:.6f} A")

    return aligned


def compute_pairwise_rmsd(all_coords: list[np.ndarray], all_positions: list[list[int]],
                          names: list[str]) -> pd.DataFrame:
    """Compute RMSD between every pair of structures after alignment."""
    n = len(all_coords)
    rmsd_matrix = np.zeros((n, n))

    for i, j in itertools.combinations(range(n), 2):
        _, rmsd = align_to_reference(all_coords[i], all_coords[j],
                                     all_positions[i], all_positions[j])
        rmsd_matrix[i, j] = rmsd
        rmsd_matrix[j, i] = rmsd

    return pd.DataFrame(rmsd_matrix, index=names, columns=names).round(3)


def load_ptm_and_mutation_positions(uniprot: str) -> tuple[set[int], set[int]]:
    """Load PTM site positions and mutation positions from the pipeline's intermediate TSV."""
    ptm_positions: set[int] = set()
    mutation_positions: set[int] = set()

    if not PTM_TSV.exists():
        return ptm_positions, mutation_positions

    df = pd.read_csv(PTM_TSV, sep="\t", dtype=str, keep_default_na=False)
    rows = df[df["uniprot_id"] == uniprot]
    if rows.empty:
        return ptm_positions, mutation_positions

    row = rows.iloc[0]

    for token in str(row.get("ptms_on_protein", "")).split(";"):
        m = re.search(r"([A-Z])(\d+)", token.strip())
        if m:
            ptm_positions.add(int(m.group(2)))

    for token in str(row.get("mutations_on_protein", "")).split(";"):
        m = re.search(r"([A-Z])(\d+)([A-Z*])", token.strip())
        if m:
            mutation_positions.add(int(m.group(2)))

    return ptm_positions, mutation_positions


def main():
    parser = argparse.ArgumentParser(
        description="Analyze structural variance across multiple AlphaFold CIF predictions."
    )
    parser.add_argument(
        "--input-dir", default=str(DEFAULT_INPUT_DIR),
        help=f"Directory containing CIF files to compare (default: {DEFAULT_INPUT_DIR.name})",
    )
    parser.add_argument(
        "--output-dir", default=str(DEFAULT_OUTPUT_DIR),
        help=f"Output directory (default: {DEFAULT_OUTPUT_DIR.name})",
    )
    parser.add_argument(
        "--top", type=int, default=10,
        help="Number of top variable residues to report (default: 10)",
    )
    parser.add_argument(
        "--range", nargs=2, type=int, default=None, metavar=("START", "END"),
        help="Restrict analysis to residue positions START-END (e.g. --range 0 630)",
    )
    parser.add_argument(
        "--uniprot", default=None,
        help="UniProt accession for PTM/mutation cross-referencing (auto-detected from CIF if possible)",
    )
    parser.add_argument(
        "--gene", default=None,
        help="Gene symbol — used to look up UniProt ID from the pipeline's intermediate data",
    )
    args = parser.parse_args()

    input_dir = Path(args.input_dir)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    cif_files = sorted(input_dir.glob("*.cif"))
    if len(cif_files) < 2:
        sys.exit(f"Error: need at least 2 CIF files in {input_dir}, found {len(cif_files)}")

    print(f"Found {len(cif_files)} CIF files in {input_dir}")

    # Determine UniProt ID: explicit arg > CIF metadata > gene lookup
    uniprot = args.uniprot or ""
    if not uniprot:
        uniprot = extract_uniprot_from_cif(cif_files[0]) or ""

    if not uniprot and args.gene and PTM_TSV.exists():
        df_lookup = pd.read_csv(PTM_TSV, sep="\t", dtype=str, keep_default_na=False)
        gene_rows = df_lookup[df_lookup["gene"].str.upper() == args.gene.upper()]
        if not gene_rows.empty:
            uniprot = gene_rows.iloc[0]["uniprot_id"]

    if uniprot:
        print(f"UniProt ID: {uniprot}")
    else:
        print("Warning: no UniProt ID — PTM/mutation cross-referencing will be skipped.")
        print("  Use --uniprot or --gene to provide it.")

    # Load all structures
    all_positions, all_coords, all_plddts, all_names = [], [], [], []
    file_names = []
    for cif in cif_files:
        print(f"  Loading {cif.name}...")
        positions, coords, plddts, res_names = load_ca_data(cif)
        all_positions.append(positions)
        all_coords.append(coords)
        all_plddts.append(plddts)
        all_names.append(res_names)
        file_names.append(cif.stem)

    # Apply range filter before alignment so all analysis uses the same subset
    if args.range:
        rng_start, rng_end = args.range
        print(f"\nFiltering to residue range {rng_start}-{rng_end}")
        for i in range(len(all_positions)):
            mask = [(rng_start <= p <= rng_end) for p in all_positions[i]]
            all_positions[i] = [p for p, keep in zip(all_positions[i], mask) if keep]
            all_coords[i] = all_coords[i][mask]
            all_plddts[i] = all_plddts[i][mask]
            all_names[i] = [n for n, keep in zip(all_names[i], mask) if keep]

    # Iterative alignment to average structure
    print("\nAligning all structures to iterative average reference...")
    aligned_coords = iterative_average_alignment(all_coords, all_positions)

    # Compute pairwise RMSD (using aligned coordinates)
    print("\nComputing pairwise RMSD matrix...")
    rmsd_df = compute_pairwise_rmsd(aligned_coords, all_positions, file_names)
    rmsd_path = output_dir / "pairwise_rmsd.tsv"
    rmsd_df.to_csv(rmsd_path, sep="\t")
    print(f"  Saved to {rmsd_path}")
    print(rmsd_df.to_string())

    # Compute per-residue variance and pLDDT stats
    shared_positions = sorted(set.intersection(*(set(p) for p in all_positions)))
    print(f"\n{len(shared_positions)} shared residue positions"
          + (f" in range {rng_start}-{rng_end}" if args.range else ""))

    # Build aligned coordinate arrays for shared positions
    coord_stack = []  # (n_structures, n_shared, 3)
    plddt_stack = []  # (n_structures, n_shared)
    res_name_list = []

    for i in range(len(aligned_coords)):
        idx_map = {p: j for j, p in enumerate(all_positions[i])}
        c = np.array([aligned_coords[i][idx_map[p]] for p in shared_positions])
        l = np.array([all_plddts[i][idx_map[p]] for p in shared_positions])
        coord_stack.append(c)
        plddt_stack.append(l)

    # Residue names from reference
    first_idx_map = {p: j for j, p in enumerate(all_positions[0])}
    res_name_list = [all_names[0][first_idx_map[p]] for p in shared_positions]

    coord_stack = np.array(coord_stack)   # (n_structures, n_residues, 3)
    plddt_stack = np.array(plddt_stack)   # (n_structures, n_residues)

    # Per-residue positional variance = mean squared deviation of each residue's position
    mean_coords = coord_stack.mean(axis=0)  # (n_residues, 3)
    deviations = coord_stack - mean_coords  # (n_structures, n_residues, 3)
    per_residue_variance = (deviations ** 2).sum(axis=2).mean(axis=0)  # (n_residues,)

    # pLDDT stats
    plddt_mean = plddt_stack.mean(axis=0)
    plddt_std = plddt_stack.std(axis=0)

    # Cross-reference PTMs and mutations
    ptm_positions, mutation_positions = load_ptm_and_mutation_positions(uniprot)
    if ptm_positions:
        print(f"Found {len(ptm_positions)} PTM sites and {len(mutation_positions)} mutation sites")
    else:
        print("No PTM/mutation data found (run pipeline step 1 first for cross-referencing)")

    # Build per-residue data table
    rows = []
    for idx, pos in enumerate(shared_positions):
        rows.append({
            "position": pos,
            "residue": res_name_list[idx],
            "positional_variance": round(float(per_residue_variance[idx]), 4),
            "plddt_mean": round(float(plddt_mean[idx]), 2),
            "plddt_std": round(float(plddt_std[idx]), 2),
            "is_ptm_site": "Yes" if pos in ptm_positions else "",
            "is_mutation_site": "Yes" if pos in mutation_positions else "",
        })
    data_df = pd.DataFrame(rows)

    # Save data
    data_path = output_dir / "variance_data.tsv"
    data_df.to_csv(data_path, sep="\t", index=False)
    print(f"\nPer-residue data saved to {data_path}")

    # Summary stats
    avg_variance = per_residue_variance.mean()
    print(f"\nProtein-wide average positional variance: {avg_variance:.4f} A^2")

    top_n = min(args.top, len(data_df))
    top_var = data_df.nlargest(top_n, "positional_variance")
    print(f"\nTop {top_n} most variable residues:")
    print(top_var.to_string(index=False))

    # ── Plot ──
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(14, 8))

    positions_arr = np.array(shared_positions)

    # Top panel: positional variance
    ax1.plot(positions_arr, per_residue_variance, linewidth=0.8, color="#3a86ff", alpha=0.8)
    ax1.fill_between(positions_arr, per_residue_variance, alpha=0.2, color="#3a86ff")

    # Mark PTM sites
    for pos in ptm_positions:
        if pos in shared_positions:
            idx = shared_positions.index(pos)
            ax1.axvline(x=pos, color="#2ecc71", alpha=0.3, linewidth=0.8)

    # Mark mutation sites
    for pos in mutation_positions:
        if pos in shared_positions:
            idx = shared_positions.index(pos)
            ax1.axvline(x=pos, color="#e74c3c", alpha=0.3, linewidth=0.8)

    # Legend entries for markers
    ax1.axvline(x=-999, color="#2ecc71", alpha=0.5, linewidth=1.5, label="PTM site")
    ax1.axvline(x=-999, color="#e74c3c", alpha=0.5, linewidth=1.5, label="Mutation site")

    ax1.set_xlabel("Residue Position", fontsize=12)
    ax1.set_ylabel("Positional Variance (A^2)", fontsize=12)
    ax1.set_xlim(0, max(shared_positions) + 1)
    ax1.set_title(f"Per-Residue Structural Variance ({len(cif_files)} structures"
                  + (f", {uniprot})" if uniprot else ")"), fontsize=14)
    ax1.legend(loc="upper right", fontsize=10)
    ax1.grid(True, alpha=0.3)

    # Bottom panel: pLDDT mean ± std
    ax2.plot(positions_arr, plddt_mean, linewidth=1, color="#f39c12")
    ax2.fill_between(positions_arr, plddt_mean - plddt_std, plddt_mean + plddt_std,
                     alpha=0.25, color="#f39c12")

    for pos in ptm_positions:
        if pos in shared_positions:
            ax2.axvline(x=pos, color="#2ecc71", alpha=0.3, linewidth=0.8)
    for pos in mutation_positions:
        if pos in shared_positions:
            ax2.axvline(x=pos, color="#e74c3c", alpha=0.3, linewidth=0.8)

    ax2.set_xlabel("Residue Position", fontsize=12)
    ax2.set_ylabel("pLDDT (mean +/- std)", fontsize=12)
    ax2.set_xlim(0, max(shared_positions) + 1)
    ax2.set_title("AlphaFold Confidence (pLDDT)", fontsize=14)
    ax2.set_ylim(0, 100)
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()
    plot_path = output_dir / "variance_plot.png"
    plt.savefig(plot_path, dpi=150)
    print(f"Plot saved to {plot_path}")
    plt.show()


if __name__ == "__main__":
    main()
