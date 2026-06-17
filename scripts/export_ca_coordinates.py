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
from pathlib import Path

import pandas as pd
import requests
from biotite.structure.io.pdbx import CIFFile, get_structure  # type: ignore[import-untyped]

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent
MODELS_ROOT = PROJECT_ROOT / "cif_models"
DEFAULT_COSMIC = PROJECT_ROOT / "data" / "Cosmic_MutantCensus_v104_GRCh38.tsv"
GENE_CACHE = PROJECT_ROOT / "data" / "cache" / "uniprot_gene_mapping.tsv"
OUTPUT_DIR = PROJECT_ROOT / "Output" / "coordinates"

_AF_API = "https://alphafold.ebi.ac.uk/api/prediction/{uid}"

_SOMATIC_STATUSES = {
    "Confirmed somatic variant",
    "Reported in another cancer sample as somatic",
}

# Standard three-letter → one-letter amino acid codes
_AA3TO1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    "SEC": "U", "PYL": "O",
}


def _download_cif(uid: str) -> list[Path]:
    """Fetch CIF file(s) for *uid* from the AlphaFold DB and save to cif_models/{uid}/."""
    print(f"Querying AlphaFold DB for {uid} ...")
    try:
        resp = requests.get(_AF_API.format(uid=uid), timeout=30)
    except requests.RequestException as exc:
        sys.exit(f"Error: AlphaFold API request failed: {exc}")

    if resp.status_code == 404:
        sys.exit(f"Error: {uid} has no AlphaFold DB entry (404). Check the UniProt accession.")
    resp.raise_for_status()

    records = resp.json()
    if isinstance(records, dict):
        records = [records]

    # Keep only canonical records — isoforms have uniprotAccession like "P11362-9"
    canonical = [r for r in records if r.get("uniprotAccession") == uid]
    if not canonical:
        sys.exit(f"Error: AlphaFold DB returned no canonical model for {uid} (isoform-only).")

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
                print(f"  Already downloaded: {filename}")
                downloaded.append(dest)
                continue

            print(f"  Downloading {filename} ...")
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
                print(f"  Warning: failed to download {cif_url}", file=sys.stderr)

    return downloaded


def _find_cif_files(uniprot_dir: Path) -> list[Path]:
    uid = uniprot_dir.name
    pat = re.compile(rf"^AF-{re.escape(uid)}-F(\d+)-model_v\d+\.", re.IGNORECASE)
    hits = [(int(pat.match(p.name).group(1)), p) for p in uniprot_dir.glob("*.cif") if pat.match(p.name)]
    return [p for _, p in sorted(hits)]


def _load_ca_from_cif(cif_file: Path) -> list[dict]:
    try:
        cif = CIFFile.read(str(cif_file))
        structure = get_structure(cif, model=1)
    except Exception as exc:
        print(f"  Warning: could not parse {cif_file.name}: {exc}", file=sys.stderr)
        return []

    if structure is None or len(structure) == 0:
        return []

    chain_ids = list(dict.fromkeys(structure.chain_id))
    if not chain_ids:
        return []

    chain = structure[structure.chain_id == chain_ids[0]]
    ca_mask = chain.atom_name == "CA"
    ca_atoms = chain[ca_mask]

    rows = []
    for i in range(len(ca_atoms)):
        one_letter = _AA3TO1.get(str(ca_atoms.res_name[i]), "X")
        x, y, z = ca_atoms.coord[i]
        rows.append({
            "residue": one_letter,
            "position": int(ca_atoms.res_id[i]),
            "x": round(float(x), 3),
            "y": round(float(y), 3),
            "z": round(float(z), 3),
        })
    return rows


def _lookup_gene(uniprot_id: str) -> str | None:
    """Return gene symbol for *uniprot_id*, checking the local cache first."""
    if GENE_CACHE.exists():
        df = pd.read_csv(GENE_CACHE, sep="\t", dtype=str, keep_default_na=False)
        id_col = "UniProt" if "UniProt" in df.columns else "uniprot_id"
        hits = df[df[id_col] == uniprot_id]
        if not hits.empty:
            gene = hits.iloc[0]["gene"]
            if gene:
                return gene

    print(f"Gene not found in cache — querying UniProt API for {uniprot_id}...")
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
                sys.exit(
                    f"Error: UniProt entry {uniprot_id} has been deleted from the database.\n"
                    "Check whether it was merged into another accession at https://www.uniprot.org"
                )
            gene_field = fields[0].strip()
            gene = gene_field.split()[0] if gene_field else None
            if gene:
                return gene
            sys.exit(
                f"Error: UniProt entry {uniprot_id} has no gene symbol.\n"
                "Use --gene SYMBOL to provide the gene name directly."
            )
    except SystemExit:
        raise
    except Exception as exc:
        print(f"  UniProt API error: {exc}", file=sys.stderr)
    return None


def _load_cosmic_mutations(
    gene: str, cosmic_file: Path
) -> tuple[dict[int, list[str]], dict[int, int]]:
    """Return position-level dicts: {pos: [mutations]} and {pos: patient_count}."""
    cols = ["GENE_SYMBOL", "MUTATION_AA", "COSMIC_SAMPLE_ID", "MUTATION_SOMATIC_STATUS"]
    print(f"Scanning COSMIC for gene {gene} ...")
    df = pd.read_csv(cosmic_file, sep="\t", usecols=cols, low_memory=False)
    df = df[df["GENE_SYMBOL"] == gene].copy()
    df = df[df["MUTATION_SOMATIC_STATUS"].isin(_SOMATIC_STATUSES)].copy()
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


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Export alpha-carbon coordinates for all residues and COSMIC missense-mutation sites."
        )
    )
    parser.add_argument("uniprot", help="UniProt accession (e.g. P04637 for TP53)")
    parser.add_argument("--gene", help="Gene symbol — skips the UniProt API gene lookup")
    parser.add_argument(
        "--cosmic",
        default=str(DEFAULT_COSMIC),
        help=f"Path to COSMIC Mutant Census TSV (default: {DEFAULT_COSMIC.name})",
    )
    args = parser.parse_args()

    uid = args.uniprot.strip().upper()
    cosmic_file = Path(args.cosmic)

    # ── 1. Locate CIF files (download from AlphaFold if not present) ─────────
    uniprot_dir = MODELS_ROOT / uid
    cif_files = _find_cif_files(uniprot_dir) if uniprot_dir.is_dir() else []

    if not cif_files:
        _download_cif(uid)
        cif_files = _find_cif_files(uniprot_dir)

    if not cif_files:
        sys.exit(f"Error: no canonical AlphaFold CIF files found in {uniprot_dir}")

    print(f"CIF fragment(s): {[f.name for f in cif_files]}")

    # ── 2. Extract CA coordinates from all fragments ──────────────────────────
    all_records: list[dict] = []
    for cf in cif_files:
        records = _load_ca_from_cif(cf)
        print(f"  {cf.name}: {len(records)} CA atoms")
        all_records.extend(records)

    if not all_records:
        sys.exit("Error: no CA atoms could be extracted")

    all_ca_df = (
        pd.DataFrame(all_records, columns=["residue", "position", "x", "y", "z"])
        .drop_duplicates(subset=["position"], keep="first")  # deduplicate overlapping fragments
        .sort_values("position")
        .reset_index(drop=True)
    )
    print(f"Total unique CA atoms: {len(all_ca_df)}")

    # ── 3. Gene symbol ────────────────────────────────────────────────────────
    gene = args.gene or _lookup_gene(uid)
    print(f"Gene: {gene}")

    # ── 4. COSMIC missense mutations ──────────────────────────────────────────
    if not cosmic_file.exists():
        sys.exit(f"Error: COSMIC file not found: {cosmic_file}")

    pos_mutations, pos_patients = _load_cosmic_mutations(gene, cosmic_file)
    print(f"Missense mutation positions in COSMIC: {len(pos_mutations)}")

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
    print(f"CA atoms at mutation positions: {len(mut_ca_df)}")

    # ── 6. Write outputs ──────────────────────────────────────────────────────
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    all_out = OUTPUT_DIR / f"{uid}_all_ca.tsv"
    mut_out = OUTPUT_DIR / f"{uid}_mutation_ca.tsv"

    all_ca_df.to_csv(all_out, sep="\t", index=False)
    mut_ca_df.to_csv(mut_out, sep="\t", index=False)

    print(f"\nDone.")
    print(f"  All CA coordinates : {all_out}  ({len(all_ca_df)} rows)")
    print(f"  Mutation CA coords : {mut_out}  ({len(mut_ca_df)} rows)")


if __name__ == "__main__":
    main()
