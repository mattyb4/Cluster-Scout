"""Export alpha-carbon coordinates for a protein from its AlphaFold CIF.

Produces two TSV files in Output/coordinates/:

  {UniProt}_all_ca.tsv       — CA coordinates for every residue
  {UniProt}_mutation_ca.tsv  — CA coordinates only at COSMIC missense-mutation positions

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


@dataclass
class ExportResult:
    """Everything produced by a CA-coordinate export run."""
    uid: str
    gene: str
    all_ca_df: pd.DataFrame
    mut_ca_df: pd.DataFrame
    all_out: Path
    mut_out: Path


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


def run_export(
    uniprot: str,
    gene: str | None = None,
    cosmic_file: Path | None = None,
    output_dir: Path = OUTPUT_DIR,
    log_cb: Callable[[str], None] = print,
) -> ExportResult:
    """Export CA coordinates (all residues + COSMIC mutation positions) for a protein.

    Raises FileNotFoundError if the COSMIC file is missing, and ValueError if no
    AlphaFold structure, no CA atoms, or no gene symbol could be resolved.
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

    # ── 5. Filter to mutation positions ───────────────────────────────────────
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
        })

    mut_ca_df = pd.DataFrame(
        mut_rows,
        columns=["residue", "position", "x", "y", "z", "mutations", "total_patients"],
    )
    log_cb(f"CA atoms at mutation positions: {len(mut_ca_df)}")

    # ── 6. Write outputs ──────────────────────────────────────────────────────
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    all_out = output_dir / f"{uid}_all_ca.tsv"
    mut_out = output_dir / f"{uid}_mutation_ca.tsv"

    all_ca_df.to_csv(all_out, sep="\t", index=False)
    mut_ca_df.to_csv(mut_out, sep="\t", index=False)

    log_cb("")
    log_cb("Done.")
    log_cb(f"  All CA coordinates : {all_out}  ({len(all_ca_df)} rows)")
    log_cb(f"  Mutation CA coords : {mut_out}  ({len(mut_ca_df)} rows)")

    return ExportResult(
        uid=uid, gene=resolved_gene,
        all_ca_df=all_ca_df, mut_ca_df=mut_ca_df,
        all_out=all_out, mut_out=mut_out,
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
    args = parser.parse_args()

    try:
        run_export(
            args.uniprot,
            gene=args.gene,
            cosmic_file=Path(args.cosmic) if args.cosmic else None,
        )
    except (FileNotFoundError, ValueError, RuntimeError) as exc:
        sys.exit(f"Error: {exc}")


if __name__ == "__main__":
    main()
