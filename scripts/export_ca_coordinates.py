"""Export alpha-carbon coordinates for a protein from its AlphaFold CIF.

Every run's output goes into its own Output/coordinates/{UniProt}/ folder, so
files for different proteins never mix together. Produces two TSV files,
each with a patients_within_10A column giving the total COSMIC patient count
summed across all missense mutations whose CA coordinate is within 10
Angstroms:

  all_ca.tsv       — CA coordinates for every residue
  mutation_ca.tsv  — CA coordinates only at COSMIC missense-mutation positions

For single-fragment proteins, also optionally produces ChimeraX-ready heatmap
scripts (skipped, with a warning, for multi-fragment proteins):

  mutations.defattr    — per-residue patients_within_10A value as a ChimeraX
                          attribute-assignment file (mutation heatmap only)
  mutations_view.cxc   — opens the CIF, loads the attribute file, and colors
                          the cartoon by it (a sequential "Reds" palette,
                          auto-scaled to the attribute's true min/max, or
                          log1p-scaled if requested)
  plddt_view.cxc     — opens the CIF and colors the cartoon by AlphaFold's
                       own per-residue confidence (pLDDT), using ChimeraX's
                       built-in "alphafold" palette
  markers_view.cxc   — opens the CIF with a plain (uncolored) cartoon; only
                       produced when mark_ptm_sites/mark_mutations are on
                       and neither heatmap is

Each heatmap is independently opt-in (mutation_heatmap/plddt_heatmap); open
either .cxc file directly in ChimeraX to see that heatmap with no manual steps.

mark_ptm_sites additionally marks each known PTM site with a small green
sphere at its CA coordinate -- an independent marker model, not a
recoloring. mark_mutations similarly shows each COSMIC mutation position's
side chain as an orange stick. Both are layered into whichever .cxc file(s)
above get written without overwriting that heatmap's own color at the
marked residue.

Accepts multiple proteins in one run -- run_batch_export() (used for both the
positional CLI arguments below and the app's batch UI) runs run_export() once
per protein, in the same output folder, applying the same options to each.
COSMIC is read and cached once for the whole batch (see
_load_cosmic_dataframe) rather than re-scanned per protein.

Usage:
    uv run scripts/export_ca_coordinates.py P04637
    uv run scripts/export_ca_coordinates.py TP53 EGFR P04637
    uv run scripts/export_ca_coordinates.py P04637 --cosmic path/to/COSMIC.tsv
    uv run scripts/export_ca_coordinates.py P04637 --plddt-heatmap --no-mutation-heatmap
    uv run scripts/export_ca_coordinates.py P04637 --mark-ptm-sites --mark-mutations
"""
from __future__ import annotations

import argparse
import os
import re
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Callable

import numpy as np
import pandas as pd
import requests

sys.path.insert(0, str(Path(__file__).resolve().parent))
from pipeline_utils import (  # noqa: E402
    project_root, AA3TO1, COSMIC_SOMATIC_STATUSES,
    find_canonical_cifs, load_first_chain, get_plddt_map,
    input_dir, resolve_input_file, COSMIC_INPUT_DIR,
    hotspots_tsv_path, looks_like_uniprot_id,
)

PROJECT_ROOT = project_root(__file__)
MODELS_ROOT = PROJECT_ROOT / "cif_models"
GENE_CACHE = PROJECT_ROOT / "data" / "cache" / "uniprot_gene_mapping.tsv"
OUTPUT_DIR = PROJECT_ROOT / "Output" / "coordinates"
# PTM-site marking always needs the ptm-proximity file specifically (its
# ptms_on_protein column) -- this tool has no --mode of its own.
PTM_TSV = hotspots_tsv_path(PROJECT_ROOT, "ptm-proximity")

_AF_API = "https://alphafold.ebi.ac.uk/api/prediction/{uid}"
NEARBY_PATIENT_RADIUS_A = 10.0


@dataclass
class ExportResult:
    """Everything produced by a CA-coordinate export run."""
    uid: str
    gene: str
    all_ca_df: pd.DataFrame
    mut_ca_df: pd.DataFrame
    all_out: Path
    mut_out: Path
    mutation_defattr_out: Path | None = None
    mutation_chimerax_script_out: Path | None = None
    plddt_chimerax_script_out: Path | None = None
    plain_chimerax_script_out: Path | None = None


def _download_cif(uid: str, log_cb: Callable[[str], None] = print) -> list[Path]:
    """Fetch CIF file(s) for *uid* from the AlphaFold DB and save to cif_models/{uid}/."""
    log_cb(f"Querying AlphaFold DB for {uid} ...")
    try:
        resp = requests.get(_AF_API.format(uid=uid), timeout=30)
    except requests.RequestException as exc:
        raise RuntimeError(f"AlphaFold API request failed: {exc}") from exc

    if resp.status_code == 404:
        raise ValueError(f"{uid} has no AlphaFold DB entry (404). Check the UniProt accession.")
    resp.raise_for_status()

    records = resp.json()
    if isinstance(records, dict):
        records = [records]

    # Keep only canonical records — isoforms have uniprotAccession like "P11362-9"
    canonical = [r for r in records if r.get("uniprotAccession") == uid]
    if not canonical:
        raise ValueError(f"AlphaFold DB returned no canonical model for {uid} (isoform-only).")

    out_dir = MODELS_ROOT / uid
    out_dir.mkdir(parents=True, exist_ok=True)

    downloaded: list[Path] = []
    with requests.Session() as session:
        session.headers.update({"User-Agent": "export-ca-coordinates/1.0"})
        for record in canonical:
            cif_url = record.get("cifUrl") or record.get("cif_url", "")
            if not cif_url:
                # Fall back to scanning all string values for a .cif URL
                cif_url = next(
                    (v for v in record.values() if isinstance(v, str) and ".cif" in v.lower()),
                    "",
                )
            if not cif_url:
                continue

            filename = cif_url.split("/")[-1]
            dest = out_dir / filename
            if dest.exists() and dest.stat().st_size > 0:
                log_cb(f"  Already downloaded: {filename}")
                downloaded.append(dest)
                continue

            log_cb(f"  Downloading {filename} ...")
            backoff = 1.6
            for attempt in range(4):
                with session.get(cif_url, stream=True, timeout=90) as r:
                    if r.status_code == 200:
                        tmp = dest.with_suffix(dest.suffix + ".part")
                        with open(tmp, "wb") as f:
                            for chunk in r.iter_content(chunk_size=256 * 1024):
                                if chunk:
                                    f.write(chunk)
                        os.replace(tmp, dest)
                        downloaded.append(dest)
                        break
                    time.sleep(backoff ** attempt)
            else:
                log_cb(f"  Warning: failed to download {cif_url}")

    return downloaded


def _load_ca_from_cif(cif_file: Path) -> list[dict]:
    """Extract alpha-carbon coordinates from a CIF file as a list of {residue, position, x, y, z} dicts."""
    chain = load_first_chain(cif_file)
    if chain is None:
        return []

    ca_mask = chain.atom_name == "CA"
    ca_atoms = chain[ca_mask]

    rows = []
    for i in range(len(ca_atoms)):
        one_letter = AA3TO1.get(str(ca_atoms.res_name[i]), "X")
        x, y, z = ca_atoms.coord[i]
        rows.append({
            "residue": one_letter,
            "position": int(ca_atoms.res_id[i]),
            "x": round(float(x), 3),
            "y": round(float(y), 3),
            "z": round(float(z), 3),
        })
    return rows


def _lookup_gene(uniprot_id: str, log_cb: Callable[[str], None] = print) -> str | None:
    """Return gene symbol for *uniprot_id*, checking the local cache first."""
    if GENE_CACHE.exists():
        df = pd.read_csv(GENE_CACHE, sep="\t", dtype=str, keep_default_na=False)
        id_col = "UniProt" if "UniProt" in df.columns else "uniprot_id"
        hits = df[df[id_col] == uniprot_id]
        if not hits.empty:
            gene = hits.iloc[0]["gene"]
            if gene:
                return gene

    log_cb(f"Gene not found in cache — querying UniProt API for {uniprot_id}...")
    try:
        resp = requests.get(
            f"https://rest.uniprot.org/uniprotkb/{uniprot_id}",
            params={"format": "tsv", "fields": "gene_names,protein_name"},
            timeout=15,
        )
        resp.raise_for_status()
        lines = [l for l in resp.text.strip().splitlines() if l]
        if len(lines) >= 2:
            fields = lines[1].split("\t")
            # Detect deleted/merged entries
            protein_name = fields[1].strip() if len(fields) > 1 else ""
            if protein_name.lower() == "deleted":
                raise ValueError(
                    f"UniProt entry {uniprot_id} has been deleted from the database. "
                    "Check whether it was merged into another accession at https://www.uniprot.org"
                )
            gene_field = fields[0].strip()
            gene = gene_field.split()[0] if gene_field else None
            if gene:
                return gene
            raise ValueError(
                f"UniProt entry {uniprot_id} has no gene symbol. "
                "Provide the gene name directly."
            )
    except ValueError:
        raise
    except Exception as exc:
        log_cb(f"  UniProt API error: {exc}")
    return None


def _lookup_uniprot_from_gene(gene: str, log_cb: Callable[[str], None] = print) -> str | None:
    """Return the reviewed human UniProt accession for *gene*, or None if not
    found. The reverse of _lookup_gene — lets a gene symbol be given on its
    own, with the UniProt accession resolved automatically.
    """
    log_cb(f"Looking up UniProt accession for gene {gene} ...")
    try:
        resp = requests.get(
            "https://rest.uniprot.org/uniprotkb/search",
            params={
                "query": f"gene_exact:{gene} AND organism_id:9606 AND reviewed:true",
                "fields": "accession",
                "format": "tsv",
                "size": 1,
            },
            timeout=15,
        )
        resp.raise_for_status()
        lines = [l for l in resp.text.strip().splitlines() if l]
        if len(lines) >= 2:
            accession = lines[1].split("\t")[0].strip()
            if accession:
                return accession
    except Exception as exc:
        log_cb(f"  UniProt API error: {exc}")
    return None


_cosmic_df_cache: dict[str, pd.DataFrame] = {}


def _load_cosmic_dataframe(cosmic_file: Path) -> pd.DataFrame:
    """Load COSMIC's relevant columns, cached by resolved path.

    COSMIC's Mutant Census file is huge (hundreds of MB); a batch export
    scanning many genes should only pay that read-and-parse cost once,
    rather than once per protein.
    """
    key = str(Path(cosmic_file).resolve())
    if key not in _cosmic_df_cache:
        cols = ["GENE_SYMBOL", "MUTATION_AA", "COSMIC_SAMPLE_ID", "MUTATION_SOMATIC_STATUS"]
        _cosmic_df_cache[key] = pd.read_csv(cosmic_file, sep="\t", usecols=cols, low_memory=False)
    return _cosmic_df_cache[key]


def _load_cosmic_mutations(
    gene: str, cosmic_file: Path, log_cb: Callable[[str], None] = print,
) -> tuple[dict[int, list[str]], dict[int, int]]:
    """Load somatic missense mutations from COSMIC for a single gene.

    Returns position-level dicts: {pos: [mutations]} and {pos: patient_count}.
    """
    log_cb(f"Scanning COSMIC for gene {gene} ...")
    df = _load_cosmic_dataframe(cosmic_file)
    df = df[df["GENE_SYMBOL"] == gene].copy()
    df = df[df["MUTATION_SOMATIC_STATUS"].isin(COSMIC_SOMATIC_STATUSES)].copy()
    df["aa_change"] = df["MUTATION_AA"].str.replace(r"^p\.", "", regex=True)
    df = df[df["aa_change"].str.match(r"^[A-Z]\d+[A-Z]$", na=False)].copy()

    agg = (
        df.groupby("aa_change")["COSMIC_SAMPLE_ID"]
        .nunique()
        .reset_index(name="patients")
    )

    pos_mutations: dict[int, list[str]] = {}
    pos_patients: dict[int, int] = {}
    for _, row in agg.iterrows():
        mut = str(row["aa_change"])
        m = re.match(r"[A-Z](\d+)[A-Z]", mut)
        if not m:
            continue
        pos = int(m.group(1))
        pos_mutations.setdefault(pos, []).append(mut)
        pos_patients[pos] = pos_patients.get(pos, 0) + int(row["patients"])

    return pos_mutations, pos_patients


def _compute_patients_within_radius(
    ca_df: pd.DataFrame,
    pos_patients: dict[int, int],
    radius: float = NEARBY_PATIENT_RADIUS_A,
) -> dict[int, int]:
    """For every residue in *ca_df*, sum COSMIC patient counts across all mutation
    positions in *pos_patients* whose CA coordinate is within *radius* Angstroms
    (inclusive), matching the <= cutoff convention used elsewhere in the pipeline.
    A mutation at the residue's own position (distance 0) counts toward its own total.
    """
    if not pos_patients:
        return {int(p): 0 for p in ca_df["position"]}

    mut_rows = ca_df[ca_df["position"].isin(pos_patients)]
    if mut_rows.empty:
        return {int(p): 0 for p in ca_df["position"]}

    mut_coords = mut_rows[["x", "y", "z"]].to_numpy()
    mut_counts = np.array([pos_patients[int(p)] for p in mut_rows["position"]])

    all_coords = ca_df[["x", "y", "z"]].to_numpy()
    dists = np.linalg.norm(all_coords[:, None, :] - mut_coords[None, :, :], axis=2)
    totals = (dists <= radius) @ mut_counts

    return {int(pos): int(total) for pos, total in zip(ca_df["position"], totals)}


def load_ptm_positions(uniprot: str) -> list[tuple[int, str]]:
    """Load (position, raw PTM token) pairs from the pipeline's intermediate
    PTM/mutation-hotspot TSV for *uniprot*, e.g. [(15, "S15:Phosphorylation"), ...].

    Returns an empty list if the TSV doesn't exist yet (pipeline step 1 hasn't
    been run) or has no row for this protein.
    """
    if not PTM_TSV.exists():
        return []
    df = pd.read_csv(PTM_TSV, sep="\t", dtype=str, keep_default_na=False)
    rows = df[df["uniprot_id"] == uniprot]
    if rows.empty:
        return []

    row = rows.iloc[0]
    positions: list[tuple[int, str]] = []
    for token in str(row.get("ptms_on_protein", "")).split(";"):
        token = token.strip()
        m = re.search(r"[A-Z](\d+)", token)
        if m:
            positions.append((int(m.group(1)), token))
    return positions


def build_ptm_marker_lines(
    ca_df: pd.DataFrame, ptm_positions: list[tuple[int, str]],
    radius: float = 1.2, color: str = "green",
) -> list[str]:
    """Build one ChimeraX `shape sphere` command per PTM position, marking its
    CA coordinate with an independent sphere model rather than recoloring the
    residue -- so it can be layered on top of any heatmap (or a plain
    cartoon) without overwriting that heatmap's own color at that residue.

    PTM positions with no matching CA coordinate in *ca_df* (e.g. outside the
    exported fragment) are silently skipped.
    """
    coords = {int(row["position"]): (row["x"], row["y"], row["z"]) for _, row in ca_df.iterrows()}
    lines = []
    for pos, token in ptm_positions:
        coord = coords.get(pos)
        if coord is None:
            continue
        x, y, z = coord
        name = re.sub(r"[^\w]+", "_", token).strip("_") or f"ptm_{pos}"
        lines.append(f"shape sphere radius {radius:g} center {x:g},{y:g},{z:g} color {color} name {name}")
    return lines


def build_mutation_marker_lines(
    positions: list[int], chain_id: str = "A", color: str = "orange",
) -> list[str]:
    """Build ChimeraX commands that show each mutation position's side chain
    as a colored stick.

    Unlike PTM markers (independent `shape sphere` models with no connection
    to the real structure), this reveals and restyles the residue's actual
    side-chain atoms -- simpler since ChimeraX resolves the geometry itself
    from the residue spec, no CA-coordinate math needed. `target ab`
    restricts the coloring to atoms/bonds only ("c"/"r" would be cartoons),
    so it never overwrites a heatmap's own cartoon color at that residue.
    "sidechain" (not "sideonly") is used specifically because it includes the
    CA atom, which keeps the stick visually connected to the backbone.
    """
    lines = []
    for pos in sorted(set(int(p) for p in positions)):
        spec = f"/{chain_id}:{pos} & sidechain"
        lines.append(f"show {spec} atoms")
        lines.append(f"style {spec} stick")
        lines.append(f"color {spec} {color} target ab")
    return lines


def build_confidence_dim_lines(plddt_map: dict[int, float], chain_id: str = "A") -> list[str]:
    """Build ChimeraX `transparency` commands that dim each residue's cartoon
    in proportion to how low its pLDDT confidence is: transparency percent =
    100 - pLDDT, so a low-confidence hotspot still shows its mutation-heatmap
    color (still "how much"), but fades toward invisible rather than being
    trusted at face value the way a fully-opaque residue would be.

    There's no ChimeraX "transparency byattribute" command (unlike `color
    byattribute`), so this sets it explicitly per residue. `target c`
    restricts it to the cartoon only, leaving any atom-level markers (PTM
    spheres, mutation sticks) fully opaque and unaffected.
    """
    lines = []
    for pos, plddt in sorted(plddt_map.items()):
        pct = max(0, min(100, round(100 - plddt)))
        lines.append(f"transparency /{chain_id}:{pos} {pct} target c")
    return lines


def write_defattr_file(
    ca_df: pd.DataFrame, out_path: Path, chain_id: str = "A",
    attr_name: str = "patients_within_10A",
) -> Path:
    """Write a ChimeraX attribute-assignment file (.defattr) with one
    per-residue value: https://www.cgl.ucsf.edu/chimerax/docs/user/formats/defattr.html

    Each data line needs a leading tab before the residue spec (ChimeraX's
    parser rejects a line missing it), and the residue spec needs a leading
    "/" before the chain letter (bare "A:1" is rejected as "Bad atom specifier").
    """
    lines = [f"attribute: {attr_name}", "recipient: residues", "#"]
    for _, row in ca_df.iterrows():
        lines.append(f"\t/{chain_id}:{int(row['position'])}\t{row[attr_name]}")
    # newline="\n" forces LF-only line endings (matching ChimeraX's own
    # shipped .defattr files) instead of write_text()'s default platform
    # translation, which would write CRLF on Windows.
    out_path.write_text("\n".join(lines) + "\n", encoding="utf-8", newline="\n")
    return out_path


def write_chimerax_script(
    cif_path: Path, defattr_path: Path, out_path: Path,
    attr_name: str = "patients_within_10A",
    value_range: tuple[float, float] | None = None,
    palette: str = "Reds",
    extra_lines: list[str] = (),
    lighting: str = "soft",
) -> Path:
    """Write a ChimeraX command script (.cxc) that opens *cif_path*, loads the
    attribute data from *defattr_path*, and colors the cartoon by it as a
    heatmap. Open directly in ChimeraX to reproduce the view with no manual steps.

    *palette* defaults to "Reds" (sequential, light-to-dark) rather than a
    diverging scale, since mutation density has no meaningful zero-midpoint --
    every value is "how much", not "which direction".

    *value_range*, if given, clamps ChimeraX's `range` instead of auto-scaling
    to the true min/max: COSMIC patient counts are heavily right-skewed, so an
    unclamped range lets one hotspot outlier crush every other residue into
    the lightest color.

    *extra_lines*, if given, are inserted after the coloring command and
    before the final lighting command -- used to layer independent
    (non-recoloring) markers like PTM-site spheres on top of the heatmap.

    *lighting* defaults to "soft" (ambient-only, depends entirely on 64-way
    ambient shadowing for depth/edge definition), but callers that also apply
    per-residue transparency (see build_confidence_dim_lines) should pass
    "simple" instead: ChimeraX's ambient shadow computation doesn't handle
    transparent geometry correctly, and once any part of the model is
    transparent the *whole* model's shadows can break, leaving even opaque
    residues flatly lit at full ambient intensity with no edge definition
    ("blinding" and hard to read). "simple" uses real directional key/fill
    lights instead, which don't have this failure mode.
    """
    range_clause = f" range {value_range[0]:g},{value_range[1]:g}" if value_range else ""
    lines = [
        f'open "{cif_path}"',
        f'open "{defattr_path}"',
        "hide atoms",
        "cartoon",
        f"color byattribute r:{attr_name} #1 palette {palette} target c noValueColor gray{range_clause}",
        *extra_lines,
        f"lighting {lighting}",
    ]
    out_path.write_text("\n".join(lines) + "\n", encoding="utf-8", newline="\n")
    return out_path


def write_plddt_chimerax_script(cif_path: Path, out_path: Path, extra_lines: list[str] = ()) -> Path:
    """Write a ChimeraX command script (.cxc) that opens *cif_path* and colors
    the cartoon by pLDDT confidence.

    No defattr file is needed: AlphaFold CIFs already carry the per-residue
    pLDDT score in the standard B-factor field, and ChimeraX's own built-in
    "alphafold" palette (the same scheme the AlphaFold DB itself uses) reads
    it directly via `color bfactor`.

    *extra_lines*, if given, are inserted after the coloring command and
    before the final lighting command -- see write_chimerax_script.
    """
    lines = [
        f'open "{cif_path}"',
        "hide atoms",
        "cartoon",
        "color bfactor #1 palette alphafold",
        *extra_lines,
        "lighting soft",
    ]
    out_path.write_text("\n".join(lines) + "\n", encoding="utf-8", newline="\n")
    return out_path


def write_plain_chimerax_script(cif_path: Path, out_path: Path, extra_lines: list[str] = ()) -> Path:
    """Write a ChimeraX command script (.cxc) that opens *cif_path* and shows
    a plain cartoon with no heatmap coloring -- used when PTM-site markers
    are requested but neither heatmap is.
    """
    lines = [
        f'open "{cif_path}"',
        "hide atoms",
        "cartoon",
        *extra_lines,
        "lighting soft",
    ]
    out_path.write_text("\n".join(lines) + "\n", encoding="utf-8", newline="\n")
    return out_path


def run_export(
    uniprot: str | None = None,
    gene: str | None = None,
    cosmic_file: Path | None = None,
    output_dir: Path = OUTPUT_DIR,
    mutation_heatmap: bool = True,
    plddt_heatmap: bool = False,
    mark_ptm_sites: bool = False,
    mark_mutations: bool = False,
    log_scale: bool = False,
    dim_low_confidence: bool = False,
    log_cb: Callable[[str], None] = print,
) -> ExportResult:
    """Export CA coordinates (all residues + COSMIC mutation positions) for a protein.

    Either *uniprot* or *gene* must be given; if *uniprot* is omitted, the
    UniProt accession is resolved from *gene* via a live UniProt API lookup.

    *mutation_heatmap*/*plddt_heatmap* independently control which ChimeraX
    heatmap script(s) get written (single-fragment proteins only -- see step
    8 below). *log_scale* only affects the mutation heatmap: if True, it's
    colored by log1p(patients_within_10A) instead of the raw count, under a
    separate "patients_within_10A_log" attribute name -- see
    write_chimerax_script's docstring for why this can help with heavily
    right-skewed patient counts.

    *mark_ptm_sites*, if True, marks each known PTM site (from the pipeline's
    intermediate PTM/mutation-hotspot TSV) with a small green sphere at its CA
    coordinate in whichever heatmap script(s) get written -- or its own plain
    (uncolored) script if neither heatmap is requested. See
    build_ptm_marker_lines's docstring for why this is a separate marker
    rather than a recoloring.

    *mark_mutations*, if True, similarly shows each COSMIC mutation position's
    side chain as an orange stick -- see build_mutation_marker_lines's
    docstring. Independent of *mark_ptm_sites*; both can be on at once.

    *dim_low_confidence*, if True, dims each residue's mutation-heatmap color
    in proportion to how low its pLDDT confidence is -- see
    build_confidence_dim_lines's docstring. Only affects the mutation
    heatmap; has no effect if *mutation_heatmap* is False.

    Raises ValueError if neither uniprot nor gene is given, if a gene-only
    lookup can't be resolved to a UniProt accession, or if no AlphaFold
    structure, no CA atoms, or no gene symbol could be resolved. Raises
    FileNotFoundError if the COSMIC file is missing.
    """
    uniprot = (uniprot or "").strip()
    gene = (gene or "").strip() or None
    if not uniprot and not gene:
        raise ValueError("Provide a UniProt accession, a gene symbol, or both.")

    if uniprot:
        uid = uniprot.upper()
    else:
        resolved_uid = _lookup_uniprot_from_gene(gene, log_cb)
        if resolved_uid is None:
            raise ValueError(
                f"Could not find a reviewed human UniProt accession for gene '{gene}'. "
                "Provide the UniProt accession directly."
            )
        uid = resolved_uid.upper()
        log_cb(f"Resolved {gene} -> {uid}")

    if cosmic_file is None:
        cosmic_file = resolve_input_file(input_dir(PROJECT_ROOT, COSMIC_INPUT_DIR), (".tsv",))
    cosmic_file = Path(cosmic_file)

    # ── 1. Locate CIF files (download from AlphaFold if not present) ─────────
    uniprot_dir = MODELS_ROOT / uid
    cif_files = find_canonical_cifs(uniprot_dir) if uniprot_dir.is_dir() else []

    if not cif_files:
        _download_cif(uid, log_cb)
        cif_files = find_canonical_cifs(uniprot_dir)

    if not cif_files:
        raise ValueError(f"No canonical AlphaFold CIF files found in {uniprot_dir}")

    log_cb(f"CIF fragment(s): {[f.name for f in cif_files]}")

    # ── 2. Extract CA coordinates from all fragments ──────────────────────────
    all_records: list[dict] = []
    for cf in cif_files:
        records = _load_ca_from_cif(cf)
        log_cb(f"  {cf.name}: {len(records)} CA atoms")
        all_records.extend(records)

    if not all_records:
        raise ValueError("No CA atoms could be extracted")

    all_ca_df = (
        pd.DataFrame(all_records, columns=["residue", "position", "x", "y", "z"])
        .drop_duplicates(subset=["position"], keep="first")  # deduplicate overlapping fragments
        .sort_values("position")
        .reset_index(drop=True)
    )
    log_cb(f"Total unique CA atoms: {len(all_ca_df)}")

    # ── 3. Gene symbol ────────────────────────────────────────────────────────
    resolved_gene = gene or _lookup_gene(uid, log_cb)
    if resolved_gene is None:
        raise ValueError(
            f"Could not determine gene symbol for {uid}. Provide the gene symbol directly."
        )
    log_cb(f"Gene: {resolved_gene}")

    # ── 4. COSMIC missense mutations ──────────────────────────────────────────
    if not cosmic_file.exists():
        raise FileNotFoundError(f"COSMIC file not found: {cosmic_file}")

    pos_mutations, pos_patients = _load_cosmic_mutations(resolved_gene, cosmic_file, log_cb)
    log_cb(f"Missense mutation positions in COSMIC: {len(pos_mutations)}")

    # ── 5. Patient counts within radius of each coordinate ────────────────────
    patients_within = _compute_patients_within_radius(all_ca_df, pos_patients)
    all_ca_df["patients_within_10A"] = all_ca_df["position"].map(patients_within).astype(int)

    # ── 6. Filter to mutation positions ───────────────────────────────────────
    mut_rows = []
    for _, row in all_ca_df.iterrows():
        pos = int(row["position"])
        if pos not in pos_mutations:
            continue
        mut_rows.append({
            "residue": row["residue"],
            "position": pos,
            "x": row["x"],
            "y": row["y"],
            "z": row["z"],
            "mutations": "; ".join(sorted(pos_mutations[pos])),
            "total_patients": pos_patients[pos],
            "patients_within_10A": patients_within.get(pos, 0),
        })

    mut_ca_df = pd.DataFrame(
        mut_rows,
        columns=["residue", "position", "x", "y", "z", "mutations", "total_patients", "patients_within_10A"],
    )
    log_cb(f"CA atoms at mutation positions: {len(mut_ca_df)}")

    # ── 7. Write outputs ──────────────────────────────────────────────────────
    # Everything for this protein goes in its own {uid} subfolder, so a
    # second export (or a different protein) never mixes files together.
    output_dir = Path(output_dir) / uid
    output_dir.mkdir(parents=True, exist_ok=True)
    all_out = output_dir / "all_ca.tsv"
    mut_out = output_dir / "mutation_ca.tsv"

    all_ca_df.to_csv(all_out, sep="\t", index=False)
    mut_ca_df.to_csv(mut_out, sep="\t", index=False)

    log_cb("")
    log_cb("Done.")
    log_cb(f"  All CA coordinates : {all_out}  ({len(all_ca_df)} rows)")
    log_cb(f"  Mutation CA coords : {mut_out}  ({len(mut_ca_df)} rows)")

    # ── 8. ChimeraX heatmap/marker files (single-fragment proteins only) ──────
    mutation_defattr_out = mutation_chimerax_script_out = None
    plddt_chimerax_script_out = plain_chimerax_script_out = None
    if not (mutation_heatmap or plddt_heatmap or mark_ptm_sites or mark_mutations):
        log_cb("  No heatmaps or markers selected, skipping ChimeraX files.")
    elif len(cif_files) > 1:
        log_cb(
            f"  Skipping ChimeraX files: {uid} spans {len(cif_files)} AlphaFold "
            f"fragments, and only fragment 1's residues were exported above."
        )
    else:
        cif_path = cif_files[0].resolve()

        marker_lines: list[str] = []
        if mark_ptm_sites:
            ptm_positions = load_ptm_positions(uid)
            if ptm_positions:
                ptm_marker_lines = build_ptm_marker_lines(all_ca_df, ptm_positions)
                marker_lines += ptm_marker_lines
                log_cb(f"  Marking {len(ptm_marker_lines)} PTM site(s) as green spheres")
            else:
                log_cb(f"  No PTM site data found for {uid} in the pipeline's "
                       f"intermediate data — nothing to mark.")

        if mark_mutations:
            if not mut_ca_df.empty:
                mutation_marker_lines = build_mutation_marker_lines(mut_ca_df["position"].tolist())
                marker_lines += mutation_marker_lines
                log_cb(f"  Marking {len(mut_ca_df)} mutation position(s) as orange sticks")
            else:
                log_cb(f"  No COSMIC mutation positions found for {uid} — nothing to mark.")

        if mutation_heatmap:
            mutation_defattr_out = output_dir / "mutations.defattr"
            mutation_chimerax_script_out = output_dir / "mutations_view.cxc"

            # A separate copy (not all_ca_df itself) so the log-scaled column
            # never leaks into the returned ExportResult or the TSVs already
            # written above -- it exists only for this heatmap.
            heatmap_attr = "patients_within_10A"
            heatmap_df = all_ca_df
            if log_scale:
                heatmap_attr = "patients_within_10A_log"
                heatmap_df = all_ca_df.copy()
                heatmap_df[heatmap_attr] = np.log1p(heatmap_df["patients_within_10A"])
                log_cb("  Log-scaling mutation heatmap (log1p of patients_within_10A)")

            dim_lines: list[str] = []
            if dim_low_confidence:
                chain = load_first_chain(cif_path)
                if chain is not None:
                    dim_lines = build_confidence_dim_lines(get_plddt_map(chain))
                    log_cb(f"  Dimming {len(dim_lines)} residue(s) by confidence on the mutation heatmap")

            write_defattr_file(heatmap_df, mutation_defattr_out, attr_name=heatmap_attr)
            write_chimerax_script(
                cif_path, mutation_defattr_out.resolve(), mutation_chimerax_script_out,
                attr_name=heatmap_attr, extra_lines=dim_lines + marker_lines,
                # "soft" lighting's depth cues come entirely from ambient
                # shadowing, which breaks once any part of the model is
                # transparent (see write_chimerax_script's docstring) --
                # only switch away from it when dim_lines actually added
                # per-residue transparency.
                lighting="simple" if dim_lines else "soft",
            )
            log_cb(f"  Mutation heatmap attribute file : {mutation_defattr_out}")
            log_cb(f"  Mutation heatmap script (open this in ChimeraX) : {mutation_chimerax_script_out}")

        if plddt_heatmap:
            plddt_chimerax_script_out = output_dir / "plddt_view.cxc"
            write_plddt_chimerax_script(cif_path, plddt_chimerax_script_out, extra_lines=marker_lines)
            log_cb(f"  pLDDT heatmap script (open this in ChimeraX) : {plddt_chimerax_script_out}")

        if (mark_ptm_sites or mark_mutations) and not mutation_heatmap and not plddt_heatmap:
            plain_chimerax_script_out = output_dir / "markers_view.cxc"
            write_plain_chimerax_script(cif_path, plain_chimerax_script_out, extra_lines=marker_lines)
            log_cb(f"  Marker script (open this in ChimeraX) : {plain_chimerax_script_out}")

    return ExportResult(
        uid=uid, gene=resolved_gene,
        all_ca_df=all_ca_df, mut_ca_df=mut_ca_df,
        all_out=all_out, mut_out=mut_out,
        mutation_defattr_out=mutation_defattr_out,
        mutation_chimerax_script_out=mutation_chimerax_script_out,
        plddt_chimerax_script_out=plddt_chimerax_script_out,
        plain_chimerax_script_out=plain_chimerax_script_out,
    )


@dataclass
class BatchExportItem:
    """One protein's outcome within a run_batch_export() batch."""
    token: str
    result: ExportResult | None = None
    error: str | None = None


def run_batch_export(
    tokens: list[str],
    cosmic_file: Path | None = None,
    output_dir: Path = OUTPUT_DIR,
    mutation_heatmap: bool = True,
    plddt_heatmap: bool = False,
    mark_ptm_sites: bool = False,
    mark_mutations: bool = False,
    log_scale: bool = False,
    dim_low_confidence: bool = False,
    progress_cb: Callable[[int, int, str], None] | None = None,
    log_cb: Callable[[str], None] = print,
) -> list[BatchExportItem]:
    """Run run_export() once per entry in *tokens* (each a gene symbol or a
    UniProt accession, auto-detected via looks_like_uniprot_id), applying the
    same heatmap/marker options to every protein.

    Each protein's own outputs land in their own Output/coordinates/{uid}/
    subfolder (run_export's own design), so a batch never mixes proteins'
    files together. A failure on one protein (bad token, no AlphaFold model,
    no gene symbol resolvable, etc.) is logged and skipped rather than
    aborting the rest of the batch -- check each BatchExportItem.error to see
    which, if any, failed.

    *progress_cb*, if given, is called as progress_cb(index, total, token)
    (1-based index) right before each protein starts, for a caller that wants
    to show "N/total" progress independent of the log stream.

    COSMIC is read once and cached across the whole batch (see
    _load_cosmic_dataframe) -- without that, scanning it fresh per protein
    would make a large batch prohibitively slow.
    """
    if cosmic_file is None:
        cosmic_file = resolve_input_file(input_dir(PROJECT_ROOT, COSMIC_INPUT_DIR), (".tsv",))
    cosmic_file = Path(cosmic_file)

    total = len(tokens)
    items: list[BatchExportItem] = []
    for i, raw_token in enumerate(tokens, 1):
        token = raw_token.strip()
        if progress_cb is not None:
            progress_cb(i, total, token)
        log_cb("")
        log_cb(f"── [{i}/{total}] {token} " + "─" * max(0, 40 - len(token)))

        kwargs = dict(
            cosmic_file=cosmic_file, output_dir=output_dir,
            mutation_heatmap=mutation_heatmap, plddt_heatmap=plddt_heatmap,
            mark_ptm_sites=mark_ptm_sites, mark_mutations=mark_mutations,
            log_scale=log_scale, dim_low_confidence=dim_low_confidence,
            log_cb=log_cb,
        )
        try:
            if looks_like_uniprot_id(token):
                result = run_export(uniprot=token, **kwargs)
            else:
                result = run_export(gene=token, **kwargs)
        except (FileNotFoundError, ValueError, RuntimeError) as exc:
            log_cb(f"  Error: {exc}")
            items.append(BatchExportItem(token=token, error=str(exc)))
        else:
            items.append(BatchExportItem(token=token, result=result))

    n_ok = sum(1 for item in items if item.result is not None)
    log_cb("")
    log_cb(f"Batch complete: {n_ok}/{total} succeeded.")
    return items


def main() -> None:
    """Export all alpha-carbon coordinates and mutation-site coordinates for a given UniProt protein."""
    parser = argparse.ArgumentParser(
        description=(
            "Export alpha-carbon coordinates for all residues and COSMIC missense-mutation sites."
        )
    )
    parser.add_argument(
        "proteins", nargs="*", default=[],
        help="One or more gene symbols and/or UniProt accessions (e.g. "
             "TP53 P04637 EGFR), auto-detected and each exported in turn. "
             "A UniProt accession's gene symbol, or a gene symbol's UniProt "
             "accession, is resolved automatically as needed.",
    )
    parser.add_argument(
        "--gene",
        help="Deprecated alias for a single positional token -- kept for "
             "backward compatibility. Added to the proteins list above if given.",
    )
    parser.add_argument(
        "--cosmic",
        default=None,
        help="Path to COSMIC Mutant Census TSV (default: auto-detected from data/input/cosmic/)",
    )
    parser.add_argument(
        "--mutation-heatmap", action=argparse.BooleanOptionalAction, default=True,
        help="Generate the mutation (patients-within-10A) ChimeraX heatmap. "
             "Enabled by default; use --no-mutation-heatmap to skip it.",
    )
    parser.add_argument(
        "--plddt-heatmap", action="store_true",
        help="Also generate a ChimeraX heatmap script colored by AlphaFold's "
             "per-residue pLDDT confidence, using ChimeraX's built-in "
             "AlphaFold palette.",
    )
    parser.add_argument(
        "--log-scale", action="store_true",
        help="Color the mutation heatmap by log1p(patients_within_10A) instead "
             "of the raw count -- helps when patient counts are heavily "
             "right-skewed and a linear scale would crush most residues into "
             "one flat color. Has no effect on the pLDDT heatmap.",
    )
    parser.add_argument(
        "--mark-ptm-sites", action="store_true",
        help="Mark each known PTM site with a small green sphere at its CA "
             "coordinate, layered on top of whichever heatmap(s) are "
             "generated (or a plain cartoon if neither is).",
    )
    parser.add_argument(
        "--mark-mutations", action="store_true",
        help="Show each COSMIC mutation position's side chain as an orange "
             "stick, layered on top of whichever heatmap(s) are generated "
             "(or a plain cartoon if neither is). Independent of "
             "--mark-ptm-sites; both can be used together.",
    )
    parser.add_argument(
        "--dim-low-confidence", action="store_true",
        help="Dim each residue's mutation-heatmap color in proportion to how "
             "low its pLDDT confidence is (100 - pLDDT percent transparent). "
             "Only affects the mutation heatmap, and switches its lighting "
             "from 'soft' to 'simple' (ChimeraX's ambient shadows don't "
             "render correctly with transparent geometry present).",
    )
    args = parser.parse_args()

    tokens = list(args.proteins)
    if args.gene:
        tokens.append(args.gene)
    if not tokens:
        sys.exit("Error: Provide at least one gene symbol or UniProt accession.")

    try:
        items = run_batch_export(
            tokens,
            cosmic_file=Path(args.cosmic) if args.cosmic else None,
            mutation_heatmap=args.mutation_heatmap,
            plddt_heatmap=args.plddt_heatmap,
            mark_ptm_sites=args.mark_ptm_sites,
            mark_mutations=args.mark_mutations,
            log_scale=args.log_scale,
            dim_low_confidence=args.dim_low_confidence,
        )
    except (FileNotFoundError, ValueError, RuntimeError) as exc:
        sys.exit(f"Error: {exc}")

    if not any(item.result is not None for item in items):
        sys.exit(1)


if __name__ == "__main__":
    main()
