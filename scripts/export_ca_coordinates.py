"""Export alpha-carbon coordinates for a protein from its AlphaFold CIF.

Every run's output goes into its own Output/coordinates/{UniProt}/ folder, so
files for different proteins never mix together. Produces two TSV files,
each with a patients_within_10A column giving the total COSMIC patient count
summed across all missense mutations whose CA coordinate is within 10
Angstroms:

  all_ca.tsv       — CA coordinates for every residue
  mutation_ca.tsv  — CA coordinates only at COSMIC missense-mutation positions

For single-fragment proteins, also produces a ChimeraX-ready pair (skipped,
with a warning, for multi-fragment proteins — see write_chimerax_files):

  mutations.defattr — per-residue heatmap value as a ChimeraX
                       attribute-assignment file — log1p(patients_within_10A)
                       by default (log_scale=True), or the raw count if
                       log-scaling is turned off
  view.cxc           — a ChimeraX command script that opens the CIF, loads
                        the attribute file, and colors the cartoon by it as
                        a heatmap (a sequential "Reds" palette, capped at the
                        99th/90th percentile for log/linear mode
                        respectively, to stay legible despite how skewed
                        mutation density usually is). Open this file
                        directly in ChimeraX to see the result.

Usage:
    uv run scripts/export_ca_coordinates.py P04637
    uv run scripts/export_ca_coordinates.py P04637 --gene TP53
    uv run scripts/export_ca_coordinates.py P04637 --cosmic path/to/COSMIC.tsv
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
    find_canonical_cifs, load_first_chain,
    input_dir, resolve_input_file, COSMIC_INPUT_DIR,
)

PROJECT_ROOT = project_root(__file__)
MODELS_ROOT = PROJECT_ROOT / "cif_models"
GENE_CACHE = PROJECT_ROOT / "data" / "cache" / "uniprot_gene_mapping.tsv"
OUTPUT_DIR = PROJECT_ROOT / "Output" / "coordinates"

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
    defattr_out: Path | None = None
    chimerax_script_out: Path | None = None


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


def _load_cosmic_mutations(
    gene: str, cosmic_file: Path, log_cb: Callable[[str], None] = print,
) -> tuple[dict[int, list[str]], dict[int, int]]:
    """Load somatic missense mutations from COSMIC for a single gene.

    Returns position-level dicts: {pos: [mutations]} and {pos: patient_count}.
    """
    cols = ["GENE_SYMBOL", "MUTATION_AA", "COSMIC_SAMPLE_ID", "MUTATION_SOMATIC_STATUS"]
    log_cb(f"Scanning COSMIC for gene {gene} ...")
    df = pd.read_csv(cosmic_file, sep="\t", usecols=cols, low_memory=False)
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


def write_defattr_file(
    ca_df: pd.DataFrame, out_path: Path, chain_id: str = "A",
    attr_name: str = "patients_within_10A",
) -> Path:
    """Write a ChimeraX attribute-assignment file (.defattr) with one
    per-residue value, using the format documented at
    https://www.cgl.ucsf.edu/chimerax/docs/user/formats/defattr.html

    Each data line needs a LEADING tab before the residue spec (confirmed
    against ChimeraX's own shipped example files via raw byte inspection —
    the rendered docs example is easy to misread as spec+tab+value only,
    but ChimeraX's parser rejects a data line missing that initial tab,
    reporting it as "not of the form 'name: value'"), and the residue spec
    itself needs a leading "/" before the chain letter (ChimeraX's atom-spec
    grammar requires it for a chain specifier — bare "A:1" is rejected as
    "Bad atom specifier"; it must be "/A:1").
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
) -> Path:
    """Write a ChimeraX command script (.cxc) that opens *cif_path*, loads
    the attribute data from *defattr_path*, and colors the cartoon by it as
    a heatmap. Open this file directly in ChimeraX (File > Open, or drag
    onto the window) to reproduce the view with no manual steps.

    *palette* defaults to "Reds", a built-in ColorBrewer sequential palette
    (light-to-dark, one hue) — a diverging blue-white-red scale implies a
    meaningful midpoint/"neutral" value, which mutation density doesn't
    have; every value is "how much", not "which direction from zero".

    *value_range*, if given, is passed to ChimeraX's `range` option instead
    of letting it auto-scale to the attribute's true min/max. COSMIC patient
    counts are heavily right-skewed — one true hotspot residue can be 10-100x
    any other position — so an unclamped auto-range stretches the whole
    gradient to fit that single outlier, crushing every other residue down
    into the lightest color. Values above the given upper bound are simply
    clamped to the top color rather than further stretching the scale.
    """
    range_clause = f" range {value_range[0]:g},{value_range[1]:g}" if value_range else ""
    lines = [
        f'open "{cif_path}"',
        f'open "{defattr_path}"',
        "hide atoms",
        "cartoon",
        f"color byattribute r:{attr_name} #1 palette {palette} target c noValueColor gray{range_clause}",
        "lighting soft",
    ]
    out_path.write_text("\n".join(lines) + "\n", encoding="utf-8", newline="\n")
    return out_path


def run_export(
    uniprot: str,
    gene: str | None = None,
    cosmic_file: Path | None = None,
    output_dir: Path = OUTPUT_DIR,
    log_scale: bool = True,
    log_cb: Callable[[str], None] = print,
) -> ExportResult:
    """Export CA coordinates (all residues + COSMIC mutation positions) for a protein.

    Raises FileNotFoundError if the COSMIC file is missing, and ValueError if no
    AlphaFold structure, no CA atoms, or no gene symbol could be resolved.

    *log_scale* controls the ChimeraX heatmap's coloring, not the TSV outputs
    (those always report the raw patients_within_10A count). Mutation density
    is heavily right-skewed, so log-scaling (the default) spreads out the
    color gradient across far more of the protein than a linear scale would;
    the color-scale cap is also mode-dependent — the 99th percentile of
    log1p(patients_within_10A) when log-scaled, or the 90th percentile of the
    raw count on a linear scale, since linear's wider color steps need a
    tighter cap to keep the bulk of the structure differentiated.
    """
    uid = uniprot.strip().upper()
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

    # ── 8. ChimeraX heatmap files (single-fragment proteins only) ─────────────
    defattr_out = chimerax_script_out = None
    if len(cif_files) > 1:
        log_cb(
            f"  Skipping ChimeraX files: {uid} spans {len(cif_files)} AlphaFold "
            f"fragments, and only fragment 1's residues were exported above."
        )
    else:
        defattr_out = output_dir / "mutations.defattr"
        chimerax_script_out = output_dir / "view.cxc"

        if log_scale:
            heatmap_attr = "log_patients_within_10A"
            all_ca_df[heatmap_attr] = np.log1p(all_ca_df["patients_within_10A"])
            cap_percentile = 0.99
        else:
            heatmap_attr = "patients_within_10A"
            cap_percentile = 0.90

        write_defattr_file(all_ca_df, defattr_out, attr_name=heatmap_attr)

        # Cap the color scale at a percentile rather than the true max:
        # patients_within_10A is heavily right-skewed (one hotspot residue can
        # dwarf every other position), so an uncapped range crushes almost the
        # whole structure into the lightest color. Values above the cap just
        # render as the most intense color instead of stretching the scale
        # further. Linear mode uses a tighter cap (90th vs. 99th percentile)
        # since its wider color steps need a smaller span to stay legible.
        cap = float(all_ca_df[heatmap_attr].quantile(cap_percentile))
        value_range = (0.0, cap) if cap > 0 else None

        write_chimerax_script(
            cif_files[0].resolve(), defattr_out.resolve(), chimerax_script_out,
            attr_name=heatmap_attr, value_range=value_range,
        )
        log_cb(f"  ChimeraX attribute file ({'log-scaled' if log_scale else 'linear'}) : {defattr_out}")
        log_cb(f"  ChimeraX script (open this in ChimeraX) : {chimerax_script_out}")

    return ExportResult(
        uid=uid, gene=resolved_gene,
        all_ca_df=all_ca_df, mut_ca_df=mut_ca_df,
        all_out=all_out, mut_out=mut_out,
        defattr_out=defattr_out, chimerax_script_out=chimerax_script_out,
    )


def main() -> None:
    """Export all alpha-carbon coordinates and mutation-site coordinates for a given UniProt protein."""
    parser = argparse.ArgumentParser(
        description=(
            "Export alpha-carbon coordinates for all residues and COSMIC missense-mutation sites."
        )
    )
    parser.add_argument("uniprot", help="UniProt accession (e.g. P04637 for TP53)")
    parser.add_argument("--gene", help="Gene symbol — skips the UniProt API gene lookup")
    parser.add_argument(
        "--cosmic",
        default=None,
        help="Path to COSMIC Mutant Census TSV (default: auto-detected from data/input/cosmic/)",
    )
    parser.add_argument(
        "--linear-scale", action="store_true",
        help="Color the ChimeraX heatmap by raw patient count instead of log1p(count) "
             "(default: log-scaled, since mutation density is usually heavily skewed)",
    )
    args = parser.parse_args()

    try:
        run_export(
            args.uniprot,
            gene=args.gene,
            cosmic_file=Path(args.cosmic) if args.cosmic else None,
            log_scale=not args.linear_scale,
        )
    except (FileNotFoundError, ValueError, RuntimeError) as exc:
        sys.exit(f"Error: {exc}")


if __name__ == "__main__":
    main()
