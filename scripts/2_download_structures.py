from __future__ import annotations

import argparse
import os
import re
import time
from pathlib import Path
from typing import Optional, Tuple, Dict, Any, List

import pandas as pd
import requests
from tqdm import tqdm

AF_PRED_ENDPOINT = "https://alphafold.ebi.ac.uk/api/prediction/{}"
ACC_RE = re.compile(r"^(?:[A-NR-Z0-9][A-Z0-9]{5}|[OPQ][0-9][A-Z0-9]{3}[0-9])(?:-\d+)?$")


def clean_accession(cell: Any) -> Optional[str]:
    """Extract the first valid UniProt accession from a cell that may contain delimiters or junk."""
    if cell is None:
        return None
    s = str(cell).strip()
    if not s:
        return None
    parts = re.split(r"[;,\s|]+", s)
    for p in parts:
        p = p.strip()
        if ACC_RE.match(p):
            return p
    return None


def pick_urls(record: dict, prefer: str = "cif") -> dict[str, str]:
    """Select the best structure and PAE download URLs from an AlphaFold prediction record."""
    urls = [v for v in record.values() if isinstance(v, str) and v.startswith("http")]

    def pick_by_priority(priorities):
        best_url = ""
        best_score = -1
        for u in urls:
            lu = u.lower()
            for endings, score in priorities:
                if any(lu.endswith(e) for e in endings):
                    s = score + (1 if lu.endswith(".gz") else 0)
                    if s > best_score:
                        best_score = s
                        best_url = u
        return best_url

    if prefer == "cif":
        structure_url = pick_by_priority([
            ((".cif.gz", ".cif"), 300),
            ((".bcif.gz", ".bcif"), 200),
            ((".pdb.gz", ".pdb"), 100),
        ])
    else:
        structure_url = pick_by_priority([
            ((".pdb.gz", ".pdb"), 300),
            ((".cif.gz", ".cif"), 200),
            ((".bcif.gz", ".bcif"), 100),
        ])
    pae_url = record.get("paeDocUrl", "") or pick_by_priority([
        ((".pae.json.gz", ".pae.json"), 300),
    ])
    return {"structure_url": structure_url, "pae_url": pae_url}

def fetch_prediction(acc: str, session: requests.Session, retries: int = 4) -> Optional[Any]:
    """Query the AlphaFold DB API for prediction metadata, with retry and backoff."""
    url = AF_PRED_ENDPOINT.format(acc)
    backoff = 1.6
    for attempt in range(retries):
        r = session.get(url, timeout=30)
        if r.status_code == 200:
            return r.json()
        if r.status_code == 404:
            return None
        time.sleep(backoff ** attempt)
    r.raise_for_status()
    return None

def download(url: str, outpath: Path, session: requests.Session, retries: int = 4) -> None:
    """Download a file via streaming GET with retry, using an atomic temp-file rename."""
    outpath.parent.mkdir(parents=True, exist_ok=True)
    if outpath.exists() and outpath.stat().st_size > 0:
        return
    backoff = 1.6
    for attempt in range(retries):
        with session.get(url, stream=True, timeout=90) as r:
            if r.status_code == 200:
                tmp = outpath.with_suffix(outpath.suffix + ".part")
                with open(tmp, "wb") as f:
                    for chunk in r.iter_content(chunk_size=1024 * 256):
                        if chunk:
                            f.write(chunk)
                os.replace(tmp, outpath)
                return
        time.sleep(backoff ** attempt)
    raise RuntimeError(f"Failed download after retries: {url}")


def _is_cached(acc_dir: Path, acc: str, also_pae: bool) -> bool:
    """Return True when all needed files for this accession are already present locally.

    AlphaFold encodes the model version in the filename (e.g. model_v6.cif), so this
    check is implicitly version-aware: a version bump produces a new filename that won't
    exist yet, causing the cache check to fail and triggering a fresh download.
    """
    if not acc_dir.is_dir():
        return False
    struct_re = re.compile(rf"^AF-{re.escape(acc)}-F\d+-model_v\d+\.", re.IGNORECASE)
    pae_re    = re.compile(rf"^AF-{re.escape(acc)}-F\d+-predicted_aligned_error_v\d+\.", re.IGNORECASE)
    files = list(acc_dir.iterdir())
    if not any(struct_re.match(f.name) for f in files):
        return False
    if also_pae and not any(pae_re.match(f.name) for f in files):
        return False
    return True


def read_table(path: str) -> pd.DataFrame:
    """Read a tabular input file, auto-detecting format from the file extension."""
    p = Path(path)
    if p.suffix.lower() in (".xlsx", ".xls"):
        return pd.read_excel(p)
    # default TSV/CSV guess
    if p.suffix.lower() == ".csv":
        return pd.read_csv(p)
    return pd.read_csv(p, sep="\t")


SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent
DEFAULT_LOGS_DIR = PROJECT_ROOT / "Output" / "logs"


def main(in_path: str, id_column: str, out_dir: str, prefer: str, also_pae: bool, delay: float, logs_dir: Path) -> None:
    """Download AlphaFold CIF structures (and optionally PAE files) for all UniProt IDs in the input table."""
    df = read_table(in_path)
    if id_column not in df.columns:
        raise ValueError(f"Column '{id_column}' not found. Columns: {list(df.columns)}")

    raw = df[id_column].tolist()
    accs = sorted({a for a in (clean_accession(x) for x in raw) if a})

    out_base = Path(out_dir)
    out_base.mkdir(parents=True, exist_ok=True)
    logs_dir.mkdir(parents=True, exist_ok=True)

    report_rows = []
    cached_count = 0

    with requests.Session() as s:
        s.headers.update({"User-Agent": "bulk-afdb-downloader/1.0"})
        for acc in tqdm(accs, desc="Downloading AlphaFold structures"):
            row = {"UniProt": acc, "status": "", "structure_file": "", "pae_file": "", "note": ""}

            if _is_cached(out_base / acc, acc, also_pae):
                row["status"] = "ALREADY_CACHED"
                row["note"] = "Files already present; skipped API call"
                report_rows.append(row)
                cached_count += 1
                continue

            try:
                meta = fetch_prediction(acc, s)
                if not meta:
                    row["status"] = "NO_ENTRY"
                    row["note"] = "No AlphaFold DB record (404)"
                    report_rows.append(row)
                    tqdm.write(f"  {acc}: no AlphaFold DB entry (404)")
                    continue

                records = meta if isinstance(meta, list) else [meta]
                # Keep only canonical records (uniprotAccession == acc, no isoform dash-suffix).
                # Isoform records have uniprotAccession like "P11362-9" and model different sequences.
                canonical_records = [r for r in records if r.get("uniprotAccession") == acc]
                if not canonical_records:
                    row["status"] = "NO_CANONICAL_MODEL"
                    row["note"] = f"No canonical AFDB model; {len(records)} isoform-only record(s)"
                    report_rows.append(row)
                    tqdm.write(f"  {acc}: no canonical model (isoforms only)")
                    continue
                records = canonical_records
                downloaded = 0
                for record in records:
                    urls = pick_urls(record, prefer=prefer)
                    if not urls.get("structure_url"):
                        continue
                    struct_url = urls["structure_url"]
                    struct_name = struct_url.split("/")[-1]
                    struct_path = out_base / acc / struct_name
                    download(struct_url, struct_path, s)
                    row["structure_file"] = str(struct_path)
                    downloaded += 1
                    if also_pae and urls.get("pae_url"):
                        pae_url = urls["pae_url"]
                        pae_name = pae_url.split("/")[-1]
                        pae_path = out_base / acc / pae_name
                        download(pae_url, pae_path, s)
                        row["pae_file"] = str(pae_path)
                    time.sleep(delay)

                if not downloaded:
                    row["status"] = "NO_STRUCTURE_URL"
                    row["note"] = f"No structure URL in any of {len(records)} record(s)"
                    report_rows.append(row)
                    tqdm.write(f"  {acc}: no structure URL found")
                    continue

                row["status"] = "DOWNLOADED"
                row["note"] = struct_name
                report_rows.append(row)

            except Exception as e:
                row["status"] = "ERROR"
                row["note"] = str(e)
                report_rows.append(row)
                tqdm.write(f"  {acc}: ERROR — {e}")

    rep = pd.DataFrame(report_rows)
    rep_path = logs_dir / "download_report.tsv"
    rep.to_csv(rep_path, sep="\t", index=False)
    print(f"\nWrote report: {rep_path}")

    downloaded_count = (rep["status"] == "DOWNLOADED").sum()
    print(f"Already cached (skipped): {cached_count}  |  Downloaded this run: {downloaded_count}")

    error_rows = rep[~rep["status"].isin({"DOWNLOADED", "ALREADY_CACHED"})]
    if not error_rows.empty:
        err_path = logs_dir / "download_errors.tsv"
        error_rows.to_csv(err_path, sep="\t", index=False)
        print(f"Wrote error log ({len(error_rows)} failed): {err_path}")
    else:
        print("No errors — all accessions downloaded successfully.")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("in_path", help="Input .tsv/.csv/.xlsx file")
    ap.add_argument("--id_column", default="UniProt", help="Column containing UniProt accessions")
    ap.add_argument("--out_dir", default="afdb_models", help="Output directory")
    ap.add_argument("--prefer", choices=["cif", "pdb"], default="cif", help="Prefer mmCIF or PDB when both available")
    ap.add_argument("--also_pae", action="store_true", help="Also download PAE json when url contains 'pae'")
    ap.add_argument("--delay", type=float, default=0.1, help="Polite delay between requests (seconds)")
    ap.add_argument("--logs_dir", default=str(DEFAULT_LOGS_DIR), help="Directory for download_report.tsv and download_errors.tsv")
    args = ap.parse_args()

    main(args.in_path, args.id_column, args.out_dir, args.prefer, args.also_pae, args.delay, Path(args.logs_dir))
