"""Parameter sweep: test different distance cutoffs to find the optimal radius.

For a set of proteins, computes how many mutations fall within each tested
radius (4-20 Å in 1 Å steps) of each PTM site, averages across PTM sites
per protein, and plots the result.

Usage:
    uv run scripts/radius_sweep.py
    uv run scripts/radius_sweep.py --genes EGFR TP53 VHL
    uv run scripts/radius_sweep.py --genes P04637 Q06124
    uv run scripts/radius_sweep.py --radii 4 25 1

--genes accepts either gene symbols or UniProt accessions, auto-detected
by format (see looks_like_uniprot_id) — useful when a protein's gene
symbol is missing or ambiguous in the dataset.
"""
from __future__ import annotations

import argparse
import re
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from pipeline_utils import (  # noqa: E402
    project_root, find_canonical_cif, find_canonical_cifs, load_first_chain,
    COSMIC_SOMATIC_STATUSES, input_dir, resolve_input_file, COSMIC_INPUT_DIR,
)

PROJECT_ROOT = project_root(__file__)
MODELS_ROOT = PROJECT_ROOT / "cif_models"
PTM_TSV = PROJECT_ROOT / "data" / "steps" / "PTMD_COSMIC_hotspots_by_protein.tsv"

DEFAULT_GENES = ["EGFR", "TP53", "VHL", "CANT1", "DDR2", "PTPN11", "LZTR1", "CDK12"]
DEFAULT_MIN_CASES = 3

# Standard UniProt accession formats: 6-char (e.g. P04637) and 10-char
# (e.g. A0A099Z4Y8), per UniProt's own published pattern.
_UNIPROT_RE = re.compile(
    r"^([A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2}|[OPQ][0-9][A-Z0-9]{3}[0-9])$",
    re.IGNORECASE,
)


def looks_like_uniprot_id(token: str) -> bool:
    """Heuristic: does *token* look like a UniProt accession (e.g. P04637)
    rather than a gene symbol (e.g. TP53)?"""
    return bool(_UNIPROT_RE.match(token.strip()))


_ptm_tsv_index_cache: pd.DataFrame | None = None


def _load_ptm_tsv_index() -> pd.DataFrame:
    """Load just the gene/uniprot_id columns of PTM_TSV, cached after the
    first read — used for fast gene/UniProt-ID validation without needing
    the full sweep machinery. Returns an empty DataFrame if PTM_TSV doesn't
    exist yet (not cached, so it's retried if the file later appears).
    """
    global _ptm_tsv_index_cache
    if _ptm_tsv_index_cache is not None:
        return _ptm_tsv_index_cache
    if not PTM_TSV.exists():
        return pd.DataFrame(columns=["gene", "uniprot_id"])
    df = pd.read_csv(PTM_TSV, sep="\t", dtype=str, keep_default_na=False,
                      usecols=["gene", "uniprot_id"])
    _ptm_tsv_index_cache = df
    return df


def load_known_genes() -> set[str]:
    """Return the set of gene symbols available in PTM_TSV. Returns an empty
    set if PTM_TSV doesn't exist yet — callers should treat that as "can't
    validate" rather than "nothing is valid", since the missing-file case is
    already reported separately when a sweep actually runs.
    """
    return set(_load_ptm_tsv_index()["gene"].unique())


def resolve_gene_token(token: str) -> tuple[str, str] | None:
    """Resolve *token* — a gene symbol OR a UniProt accession — against
    PTM_TSV, returning (gene_symbol, uniprot_id) if found, else None.
    Accepting either identifier lets a gene with an ambiguous or unmapped
    symbol still be specified directly by its UniProt accession.
    """
    df = _load_ptm_tsv_index()
    if df.empty:
        return None
    token = token.strip()
    if looks_like_uniprot_id(token):
        rows = df[df["uniprot_id"].str.upper() == token.upper()]
    else:
        rows = df[df["gene"].str.upper() == token.upper()]
    if rows.empty:
        return None
    row = rows.iloc[0]
    return row["gene"], row["uniprot_id"]


def has_cif(uniprot_id: str) -> bool:
    """True if a canonical AlphaFold CIF is already downloaded for *uniprot_id*."""
    cif_dir = MODELS_ROOT / uniprot_id
    return cif_dir.is_dir() and find_canonical_cif(cif_dir) is not None


def has_multiple_fragments(uniprot_id: str) -> bool:
    """True if AlphaFold split this protein into multiple structural fragments
    (only happens for very large proteins, roughly >2700 residues). Radius
    Sweep only ever loads fragment 1, so PTM sites/mutations in fragment 2+
    are silently excluded from the analysis.
    """
    cif_dir = MODELS_ROOT / uniprot_id
    if not cif_dir.is_dir():
        return False
    return len(find_canonical_cifs(cif_dir)) > 1


@dataclass
class SweepResult:
    """Everything needed to build the sweep figure, without recomputing anything."""
    result_df: pd.DataFrame
    radii: list[float]
    min_cases: int = DEFAULT_MIN_CASES
    elbows: dict[str, float | None] = field(default_factory=dict)
    uf_elbows: dict[str, float | None] = field(default_factory=dict)
    has_unfiltered: bool = False


def get_ca_coord(chain, residue_number):
    """Return the CA coordinate for a residue, or None."""
    mask = (chain.res_id == residue_number) & (chain.atom_name == "CA")
    if not np.any(mask):
        return None
    return chain.coord[mask][0]


def load_ptm_positions(gene: str, df: pd.DataFrame, chain):
    """Extract PTM site positions with their CA coordinates.

    Returns a list of (position, coord) — PTM sites come from PTMD data, which
    is independent of any COSMIC recurrence threshold (see load_hotspot_mutations
    for the mutation side, which the threshold does apply to).
    """
    rows = df[df["gene"] == gene]
    if rows.empty:
        return []

    row = rows.iloc[0]

    ptm_coords = []
    for token in str(row.get("ptms_on_protein", "")).split(";"):
        token = token.strip()
        if not token:
            continue
        m = re.search(r"([A-Z])(\d+)", token)
        if m:
            pos = int(m.group(2))
            coord = get_ca_coord(chain, pos)
            if coord is not None:
                ptm_coords.append((pos, coord))

    return ptm_coords


_cosmic_cache: pd.DataFrame | None = None
_cosmic_counts_cache: pd.DataFrame | None = None


def _load_cosmic_df(log_cb: Callable[[str], None] = print) -> pd.DataFrame:
    """Load and filter the raw COSMIC file once, caching it for reuse across genes."""
    global _cosmic_cache
    if _cosmic_cache is not None:
        return _cosmic_cache
    cosmic_file = resolve_input_file(input_dir(PROJECT_ROOT, COSMIC_INPUT_DIR), (".tsv",))
    log_cb(f"  Loading COSMIC file: {cosmic_file.name}")
    cols = ["GENE_SYMBOL", "MUTATION_AA", "MUTATION_SOMATIC_STATUS"]
    cosmic = pd.read_csv(cosmic_file, sep="\t", usecols=cols, low_memory=False)
    cosmic = cosmic[cosmic["MUTATION_SOMATIC_STATUS"].isin(COSMIC_SOMATIC_STATUSES)].copy()
    cosmic["aa_change"] = cosmic["MUTATION_AA"].str.replace(r"^p\.", "", regex=True)
    cosmic = cosmic[cosmic["aa_change"].str.match(r"^[A-Z]\d+[A-Z]$", na=False)]
    _cosmic_cache = cosmic
    return cosmic


def load_unfiltered_mutations(gene: str, chain, log_cb: Callable[[str], None] = print) -> list[tuple[int, np.ndarray]]:
    """Load ALL somatic missense mutation positions from raw COSMIC for a gene."""
    cosmic = _load_cosmic_df(log_cb)
    gene_cosmic = cosmic[cosmic["GENE_SYMBOL"] == gene]

    positions = set()
    for aa in gene_cosmic["aa_change"].unique():
        m = re.match(r"[A-Z](\d+)[A-Z]", aa)
        if m:
            positions.add(int(m.group(1)))

    mutation_coords = []
    for pos in sorted(positions):
        coord = get_ca_coord(chain, pos)
        if coord is not None:
            mutation_coords.append((pos, coord))
    return mutation_coords


def _load_cosmic_counts_df(log_cb: Callable[[str], None] = print) -> pd.DataFrame:
    """Load raw COSMIC with per-(gene, aa_change) distinct-sample counts
    ("affected_cases"), caching for reuse across genes. Mirrors 1_filter.py's
    own hotspot-counting logic, but lives here so Radius Sweep's threshold is
    entirely independent of whatever threshold the main pipeline last used.
    """
    global _cosmic_counts_cache
    if _cosmic_counts_cache is not None:
        return _cosmic_counts_cache
    cosmic_file = resolve_input_file(input_dir(PROJECT_ROOT, COSMIC_INPUT_DIR), (".tsv",))
    log_cb(f"  Loading COSMIC file: {cosmic_file.name}")
    cols = ["GENE_SYMBOL", "MUTATION_AA", "COSMIC_SAMPLE_ID", "MUTATION_SOMATIC_STATUS"]
    cosmic = pd.read_csv(cosmic_file, sep="\t", usecols=cols, low_memory=False)
    cosmic = cosmic[cosmic["MUTATION_SOMATIC_STATUS"].isin(COSMIC_SOMATIC_STATUSES)].copy()
    cosmic["aa_change"] = cosmic["MUTATION_AA"].str.replace(r"^p\.", "", regex=True)
    cosmic = cosmic[cosmic["aa_change"].str.match(r"^[A-Z]\d+[A-Z]$", na=False)]
    counts = (
        cosmic.groupby(["GENE_SYMBOL", "aa_change"])["COSMIC_SAMPLE_ID"]
        .nunique()
        .reset_index(name="affected_cases")
        .rename(columns={"GENE_SYMBOL": "gene"})
    )
    _cosmic_counts_cache = counts
    return counts


def load_hotspot_mutations(
    gene: str, chain, min_cases: int, log_cb: Callable[[str], None] = print,
) -> list[tuple[int, np.ndarray]]:
    """Load COSMIC mutation positions for *gene* recurring in >= min_cases
    distinct samples — computed directly from raw COSMIC here, independent of
    whatever threshold the main pipeline's step 1 used to build PTM_TSV.
    """
    counts = _load_cosmic_counts_df(log_cb)
    gene_counts = counts[(counts["gene"] == gene) & (counts["affected_cases"] >= min_cases)]

    positions = set()
    for aa in gene_counts["aa_change"].unique():
        m = re.match(r"[A-Z](\d+)[A-Z]", aa)
        if m:
            positions.add(int(m.group(1)))

    mutation_coords = []
    for pos in sorted(positions):
        coord = get_ca_coord(chain, pos)
        if coord is not None:
            mutation_coords.append((pos, coord))
    return mutation_coords


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


def _detect_elbows(df: pd.DataFrame, proteins, log_cb: Callable[[str], None]) -> dict[str, float | None]:
    """Run kneed's KneeLocator per protein on an avg_mutation_count-by-radius curve."""
    from kneed import KneeLocator

    elbows: dict[str, float | None] = {}
    for gene in proteins:
        gene_data = df[df["protein"] == gene].sort_values("radius")
        if gene_data.empty:
            continue
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
                log_cb(f"  {gene}: optimal radius = {kn.knee:.0f} A")
            else:
                log_cb(f"  {gene}: no elbow detected")
        except Exception:
            elbows[gene] = None
            log_cb(f"  {gene}: could not compute elbow")
    return elbows


def run_sweep(
    genes: list[str],
    radii: list[float],
    min_cases: int = DEFAULT_MIN_CASES,
    unfiltered: bool = False,
    output_tsv_path: Path | None = None,
    log_cb: Callable[[str], None] = print,
) -> SweepResult:
    """Run the radius sweep for a set of genes and return everything needed to plot it.

    *min_cases* is the minimum number of distinct COSMIC samples a mutation must
    recur in to count as a "hotspot" — computed live from raw COSMIC here, fully
    independent of whatever threshold the main pipeline's step 1 last used.

    Raises FileNotFoundError if the pipeline's intermediate TSV is missing, and
    ValueError if no data could be collected for any of the requested genes.
    """
    step = radii[1] - radii[0] if len(radii) > 1 else 0.0
    log_cb(f"Testing radii: {radii[0]:.0f}-{radii[-1]:.0f} A in {step:.0f} A steps "
           f"(hotspot threshold: >= {min_cases} samples)")

    if not PTM_TSV.exists():
        raise FileNotFoundError(
            f"Intermediate TSV not found at {PTM_TSV}\nRun the pipeline (step 1) first to generate it."
        )

    df = pd.read_csv(PTM_TSV, sep="\t", dtype=str, keep_default_na=False)

    # Map gene names (or UniProt IDs) to (gene, uniprot_id) pairs
    gene_to_uid = {}
    for g in genes:
        resolved = resolve_gene_token(g)
        if resolved is None:
            log_cb(f"  Warning: {g} not found in dataset, skipping")
            continue
        gene, uid = resolved
        gene_to_uid[gene] = uid

    results = []

    for gene, uid in gene_to_uid.items():
        log_cb(f"\n{gene} ({uid}):")

        cif_dir = MODELS_ROOT / uid
        cif_file = find_canonical_cif(cif_dir) if cif_dir.is_dir() else None
        if cif_file is None:
            log_cb(f"  No CIF file found, skipping")
            continue

        if has_multiple_fragments(uid):
            log_cb(f"  Warning: {gene} spans multiple AlphaFold fragments — only "
                    f"fragment 1 is analyzed, so PTM sites/mutations beyond it are excluded")

        chain = load_first_chain(cif_file)
        if chain is None:
            log_cb(f"  Could not parse CIF, skipping")
            continue

        all_ca = get_all_ca_coords(chain)
        protein_length = len(all_ca)

        ptm_coords = load_ptm_positions(gene, df, chain)
        mutation_coords = load_hotspot_mutations(gene, chain, min_cases, log_cb)
        log_cb(f"  {len(ptm_coords)} PTM sites, {len(mutation_coords)} unique mutation positions "
               f"(>= {min_cases} samples), {protein_length} residues")

        if not ptm_coords:
            log_cb(f"  No PTM coordinates found, skipping")
            continue

        avg_counts = sweep_radii(ptm_coords, mutation_coords, radii)
        log_cb("  Computing random baseline (100 permutations)...")
        baseline = random_baseline(ptm_coords, all_ca, len(mutation_coords), radii)

        for r in radii:
            results.append({
                "protein": gene,
                "radius": r,
                "dataset": "hotspot",
                "avg_mutation_count": avg_counts[r],
                "random_baseline": baseline[r],
                "avg_normalized": avg_counts[r] / protein_length * 1000,
                "random_normalized": baseline[r] / protein_length * 1000,
                "protein_length": protein_length,
                "n_mutations": len(mutation_coords),
            })
        log_cb(f"  Hotspot: {avg_counts[radii[0]]:.1f} avg at {radii[0]:.0f}A "
               f"-> {avg_counts[radii[-1]]:.1f} avg at {radii[-1]:.0f}A")

        if unfiltered:
            log_cb(f"  Loading unfiltered COSMIC mutations for {gene}...")
            unfiltered_coords = load_unfiltered_mutations(gene, chain, log_cb)
            log_cb(f"  {len(unfiltered_coords)} unfiltered mutation positions")

            uf_counts = sweep_radii(ptm_coords, unfiltered_coords, radii)
            log_cb("  Computing unfiltered random baseline...")
            uf_baseline = random_baseline(ptm_coords, all_ca, len(unfiltered_coords), radii)

            for r in radii:
                results.append({
                    "protein": gene,
                    "radius": r,
                    "dataset": "unfiltered",
                    "avg_mutation_count": uf_counts[r],
                    "random_baseline": uf_baseline[r],
                    "avg_normalized": uf_counts[r] / protein_length * 1000,
                    "random_normalized": uf_baseline[r] / protein_length * 1000,
                    "protein_length": protein_length,
                    "n_mutations": len(unfiltered_coords),
                })
            log_cb(f"  Unfiltered: {uf_counts[radii[0]]:.1f} avg at {radii[0]:.0f}A "
                   f"-> {uf_counts[radii[-1]]:.1f} avg at {radii[-1]:.0f}A")

    if not results:
        raise ValueError("No data collected. Check that CIF files are downloaded for the target genes.")

    result_df = pd.DataFrame(results)

    hotspot_df = result_df[result_df["dataset"] == "hotspot"]
    proteins = hotspot_df["protein"].unique()

    log_cb("\n-- Elbow Detection (hotspot-filtered) --")
    elbows = _detect_elbows(hotspot_df, proteins, log_cb)

    detected = [v for v in elbows.values() if v is not None]
    if detected:
        avg_elbow = np.mean(detected)
        log_cb(f"\n  Average optimal radius across proteins: {avg_elbow:.1f} A")

    has_unfiltered = "unfiltered" in result_df["dataset"].values
    uf_elbows: dict[str, float | None] = {}
    if has_unfiltered:
        unfiltered_df = result_df[result_df["dataset"] == "unfiltered"]
        log_cb("\n-- Elbow Detection (unfiltered) --")
        uf_elbows = _detect_elbows(unfiltered_df, proteins, log_cb)

        uf_detected = [v for v in uf_elbows.values() if v is not None]
        if uf_detected:
            uf_avg_elbow = np.mean(uf_detected)
            log_cb(f"\n  Average optimal radius (unfiltered): {uf_avg_elbow:.1f} A")

    if output_tsv_path is not None:
        result_df.to_csv(output_tsv_path, sep="\t", index=False)
        log_cb(f"\nData saved to: {output_tsv_path}")

    return SweepResult(
        result_df=result_df,
        radii=radii,
        min_cases=min_cases,
        elbows=elbows,
        uf_elbows=uf_elbows,
        has_unfiltered=has_unfiltered,
    )


def build_sweep_figure(result: SweepResult, fig=None):
    """Build the sweep plot from a SweepResult. Uses an injected Figure if given
    (for GUI embedding), otherwise creates one via plt.subplots (for standalone use).
    """
    result_df = result.result_df
    hotspot_df = result_df[result_df["dataset"] == "hotspot"]
    proteins = hotspot_df["protein"].unique()
    elbows = result.elbows
    detected = [v for v in elbows.values() if v is not None]
    avg_elbow = np.mean(detected) if detected else None

    has_unfiltered = result.has_unfiltered
    uf_elbows = result.uf_elbows
    uf_detected = [v for v in uf_elbows.values() if v is not None]
    uf_avg_elbow = np.mean(uf_detected) if uf_detected else None

    if fig is None:
        if has_unfiltered:
            fig, axes = plt.subplots(2, 2, figsize=(16, 10), sharex=False)
            (ax1, ax3), (ax2, ax4) = axes
        else:
            fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(11, 10), sharex=False)
    else:
        if has_unfiltered:
            axes = fig.subplots(2, 2, sharex=False)
            (ax1, ax3), (ax2, ax4) = axes
        else:
            ax1, ax2 = fig.subplots(2, 1, sharex=False)

    colors = {}

    # ── Left/top panel: hotspot raw counts + random baselines ──
    for gene in proteins:
        gene_data = hotspot_df[hotspot_df["protein"] == gene]
        n_muts = int(gene_data.iloc[0]["n_mutations"])
        line, = ax1.plot(gene_data["radius"], gene_data["avg_mutation_count"],
                         marker="o", markersize=4, linewidth=2,
                         label=f"{gene} ({n_muts} muts)")
        colors[gene] = line.get_color()

        ax1.plot(gene_data["radius"], gene_data["random_baseline"],
                 linestyle="--", linewidth=1, alpha=0.5, color=colors[gene])

        elbow_r = elbows.get(gene)
        if elbow_r is not None:
            elbow_row = gene_data[gene_data["radius"] == elbow_r]
            if not elbow_row.empty:
                ax1.scatter([elbow_r], [elbow_row.iloc[0]["avg_mutation_count"]],
                            s=120, zorder=6, marker="D", color=colors[gene],
                            edgecolors="black", linewidths=1.0)

    ax1.axvline(x=10, color="gray", linestyle="--", alpha=0.5)
    ax1.plot([], [], linestyle="--", color="gray", alpha=0.5, label="Random baseline")
    ax1.set_xlabel("Radius (A)", fontsize=12)
    ax1.set_ylabel("Avg Mutations per PTM Site", fontsize=12)
    ax1.set_title(f"Hotspot-Filtered Mutations (>={result.min_cases} samples)", fontsize=14)
    ax1.legend(loc="upper left", fontsize=8, ncol=2)
    ax1.grid(True, alpha=0.3)
    if detected:
        ax1.axvline(x=avg_elbow, color="red", linestyle=":", alpha=0.6)
        ax1.text(avg_elbow + 0.3, ax1.get_ylim()[1] * 0.95,
                 f"Avg elbow: {avg_elbow:.1f} A", color="red", fontsize=10, va="top")

    # ── Left/bottom panel: hotspot normalized ──
    for gene in proteins:
        gene_data = hotspot_df[hotspot_df["protein"] == gene]
        plen = int(gene_data.iloc[0]["protein_length"])
        ax2.plot(gene_data["radius"], gene_data["avg_normalized"],
                 marker="o", markersize=4, linewidth=2,
                 label=f"{gene} ({plen} aa)", color=colors[gene])
        ax2.plot(gene_data["radius"], gene_data["random_normalized"],
                 linestyle="--", linewidth=1, alpha=0.5, color=colors[gene])

    ax2.axvline(x=10, color="gray", linestyle="--", alpha=0.5)
    ax2.set_xlabel("Radius (A)", fontsize=12)
    ax2.set_ylabel("Avg Muts per PTM (per 1000 res)", fontsize=12)
    ax2.set_title("Hotspot - Size Normalized", fontsize=14)
    ax2.legend(loc="upper left", fontsize=8, ncol=2)
    ax2.grid(True, alpha=0.3)
    if detected:
        ax2.axvline(x=avg_elbow, color="red", linestyle=":", alpha=0.6)

    # ── Right panels: unfiltered (if present) ──
    if has_unfiltered:
        unfiltered_df = result_df[result_df["dataset"] == "unfiltered"]

        for gene in proteins:
            gene_data = unfiltered_df[unfiltered_df["protein"] == gene]
            if gene_data.empty:
                continue
            n_muts = int(gene_data.iloc[0]["n_mutations"])
            ax3.plot(gene_data["radius"], gene_data["avg_mutation_count"],
                     marker="o", markersize=4, linewidth=2,
                     label=f"{gene} ({n_muts} muts)", color=colors[gene])
            ax3.plot(gene_data["radius"], gene_data["random_baseline"],
                     linestyle="--", linewidth=1, alpha=0.5, color=colors[gene])

            elbow_r = uf_elbows.get(gene)
            if elbow_r is not None:
                elbow_row = gene_data[gene_data["radius"] == elbow_r]
                if not elbow_row.empty:
                    ax3.scatter([elbow_r], [elbow_row.iloc[0]["avg_mutation_count"]],
                                s=120, zorder=6, marker="D", color=colors[gene],
                                edgecolors="black", linewidths=1.0)

            plen = int(gene_data.iloc[0]["protein_length"])
            ax4.plot(gene_data["radius"], gene_data["avg_normalized"],
                     marker="o", markersize=4, linewidth=2,
                     label=f"{gene} ({plen} aa)", color=colors[gene])
            ax4.plot(gene_data["radius"], gene_data["random_normalized"],
                     linestyle="--", linewidth=1, alpha=0.5, color=colors[gene])

        ax3.axvline(x=10, color="gray", linestyle="--", alpha=0.5)
        ax3.plot([], [], linestyle="--", color="gray", alpha=0.5, label="Random baseline")
        ax3.set_xlabel("Radius (A)", fontsize=12)
        ax3.set_ylabel("Avg Mutations per PTM Site", fontsize=12)
        ax3.set_title("All COSMIC Mutations (unfiltered)", fontsize=14)
        ax3.legend(loc="upper left", fontsize=8, ncol=2)
        ax3.grid(True, alpha=0.3)
        if uf_detected:
            ax3.axvline(x=uf_avg_elbow, color="red", linestyle=":", alpha=0.6)
            ax3.text(uf_avg_elbow + 0.3, ax3.get_ylim()[1] * 0.95,
                     f"Avg elbow: {uf_avg_elbow:.1f} A", color="red", fontsize=10, va="top")

        ax4.axvline(x=10, color="gray", linestyle="--", alpha=0.5)
        ax4.set_xlabel("Radius (A)", fontsize=12)
        ax4.set_ylabel("Avg Muts per PTM (per 1000 res)", fontsize=12)
        ax4.set_title("Unfiltered - Size Normalized", fontsize=14)
        ax4.legend(loc="upper left", fontsize=8, ncol=2)
        ax4.grid(True, alpha=0.3)
        if uf_detected:
            ax4.axvline(x=uf_avg_elbow, color="red", linestyle=":", alpha=0.6)

    fig.tight_layout()
    return fig


def main():
    parser = argparse.ArgumentParser(
        description="Sweep distance cutoffs to find the optimal mutation-capture radius."
    )
    parser.add_argument(
        "--genes", nargs="+", default=DEFAULT_GENES,
        help=f"Gene symbols or UniProt accessions to test "
             f"(default: {' '.join(DEFAULT_GENES)})",
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
    parser.add_argument(
        "--unfiltered", action="store_true",
        help="Also run the sweep with ALL COSMIC mutations (not just hotspots) for comparison",
    )
    parser.add_argument(
        "--min-samples", type=int, default=DEFAULT_MIN_CASES,
        help=f"Minimum distinct COSMIC samples for a mutation to count as a hotspot "
             f"(default: {DEFAULT_MIN_CASES}) — independent of the main pipeline's own threshold",
    )
    args = parser.parse_args()

    radii = list(np.arange(args.radii[0], args.radii[1] + args.radii[2] / 2, args.radii[2]))

    try:
        result = run_sweep(
            args.genes, radii, min_cases=args.min_samples, unfiltered=args.unfiltered,
            output_tsv_path=Path(args.output).with_suffix(".tsv"),
        )
    except (FileNotFoundError, ValueError) as exc:
        sys.exit(str(exc))

    fig = build_sweep_figure(result)
    fig.savefig(args.output, dpi=150)
    print(f"Plot saved to: {args.output}")
    plt.show()


if __name__ == "__main__":
    main()
