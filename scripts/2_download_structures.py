from __future__ import annotations

import argparse
import os
import re
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Optional, Tuple, Dict, Any, List

import pandas as pd
import requests
from tqdm import tqdm

AF_PRED_ENDPOINT = "https://alphafold.ebi.ac.uk/api/prediction/{}"
ACC_RE = re.compile(r"^(?:[A-NR-Z0-9][A-Z0-9]{5}|[OPQ][0-9][A-Z0-9]{3}[0-9])(?:-\d+)?$")

# Be polite to the shared EBI-hosted API — same modest-concurrency convention
# used elsewhere in this project for other academic servers (kinase_library,
# AIUPred).
_DOWNLOAD_MAX_WORKERS = 6


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


_VERSION_RE = re.compile(r"model_v(\d+)\.", re.IGNORECASE)
_PAE_VERSION_RE = re.compile(r"predicted_aligned_error_v(\d+)\.", re.IGNORECASE)


def _cached_version(acc_dir: Path, acc: str, pattern: str) -> Optional[int]:
    """Return the highest version number already downloaded for *acc* matching
    *pattern* (structure or PAE), or None if nothing is cached. Multi-fragment
    proteins are expected to share one version across all fragments, so the
    max is representative."""
    if not acc_dir.is_dir():
        return None
    pat = re.compile(rf"^AF-{re.escape(acc)}-F\d+-{pattern}_v(\d+)\.", re.IGNORECASE)
    versions = [int(m.group(1)) for f in acc_dir.iterdir() if (m := pat.match(f.name))]
    return max(versions) if versions else None


def _remove_stale_version(acc_dir: Path, acc: str, keep_version: int) -> None:
    """Delete structure/PAE files for *acc* whose version isn't *keep_version*.

    Without this, an update would leave both the old and new version files
    present — and find_canonical_cif's alphabetical glob sort would pick
    whichever sorts first as a string (e.g. "v4" < "v6"), silently using the
    STALE file even after a newer one was downloaded.
    """
    for pattern in ("model", "predicted_aligned_error"):
        pat = re.compile(rf"^AF-{re.escape(acc)}-F\d+-{pattern}_v(\d+)\.", re.IGNORECASE)
        for f in acc_dir.iterdir():
            m = pat.match(f.name)
            if m and int(m.group(1)) != keep_version:
                f.unlink()


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


def _check_and_download(
    acc: str, out_base: Path, prefer: str, also_pae: bool, delay: float,
    session: requests.Session,
) -> dict:
    """Check AlphaFold DB's current version for *acc* against what's cached,
    and (re)download only if it's missing or outdated. Runs in a worker
    thread — returns a report row rather than printing directly, so all
    console output happens back on the main thread.
    """
    row = {"UniProt": acc, "status": "", "structure_file": "", "pae_file": "", "note": ""}
    acc_dir = out_base / acc
    cached_struct_v = _cached_version(acc_dir, acc, "model")
    cached_pae_v = _cached_version(acc_dir, acc, "predicted_aligned_error") if also_pae else None

    try:
        meta = fetch_prediction(acc, session)
        if not meta:
            if cached_struct_v is not None:
                row["status"] = "ALREADY_CACHED"
                row["note"] = f"AlphaFold DB has no record now (404) — kept cached v{cached_struct_v}"
            else:
                row["status"] = "NO_ENTRY"
                row["note"] = "No AlphaFold DB record (404)"
            return row

        records = meta if isinstance(meta, list) else [meta]
        # Keep only canonical records (uniprotAccession == acc, no isoform dash-suffix).
        # Isoform records have uniprotAccession like "P11362-9" and model different sequences.
        canonical_records = [r for r in records if r.get("uniprotAccession") == acc]
        if not canonical_records:
            row["status"] = "NO_CANONICAL_MODEL"
            row["note"] = f"No canonical AFDB model; {len(records)} isoform-only record(s)"
            return row

        record_urls = [pick_urls(r, prefer=prefer) for r in canonical_records]
        latest_struct_v = None
        latest_pae_v = None
        for urls in record_urls:
            m = _VERSION_RE.search(urls.get("structure_url", ""))
            if m:
                latest_struct_v = max(latest_struct_v or 0, int(m.group(1)))
            if also_pae:
                m = _PAE_VERSION_RE.search(urls.get("pae_url", ""))
                if m:
                    latest_pae_v = max(latest_pae_v or 0, int(m.group(1)))

        up_to_date = (
            cached_struct_v is not None and latest_struct_v is not None
            and cached_struct_v >= latest_struct_v
            and (not also_pae or (cached_pae_v is not None and latest_pae_v is not None
                                   and cached_pae_v >= latest_pae_v))
        )
        if up_to_date:
            row["status"] = "ALREADY_CACHED"
            row["note"] = f"Up to date (v{cached_struct_v})"
            return row

        was_cached = cached_struct_v is not None
        downloaded = 0
        struct_name = ""
        for urls in record_urls:
            if not urls.get("structure_url"):
                continue
            struct_url = urls["structure_url"]
            struct_name = struct_url.split("/")[-1]
            struct_path = acc_dir / struct_name
            download(struct_url, struct_path, session)
            row["structure_file"] = str(struct_path)
            downloaded += 1
            if also_pae and urls.get("pae_url"):
                pae_url = urls["pae_url"]
                pae_name = pae_url.split("/")[-1]
                pae_path = acc_dir / pae_name
                download(pae_url, pae_path, session)
                row["pae_file"] = str(pae_path)
            time.sleep(delay)

        if not downloaded:
            row["status"] = "NO_STRUCTURE_URL"
            row["note"] = f"No structure URL in any of {len(canonical_records)} record(s)"
            return row

        if was_cached and latest_struct_v is not None:
            _remove_stale_version(acc_dir, acc, latest_struct_v)
            row["status"] = "UPDATED"
            row["note"] = f"v{cached_struct_v} -> v{latest_struct_v}"
        else:
            row["status"] = "DOWNLOADED"
            row["note"] = struct_name
        return row

    except Exception as e:
        row["status"] = "ERROR"
        row["note"] = str(e)
        return row


def main(in_path: str, id_column: str, out_dir: str, prefer: str, also_pae: bool, delay: float, logs_dir: Path) -> None:
    """Download AlphaFold CIF structures (and optionally PAE files) for all UniProt IDs in the input table.

    Every accession is checked against AlphaFold DB's current version, in
    parallel (_DOWNLOAD_MAX_WORKERS workers) — not just skipped outright when
    something with a matching name already exists locally — so a newer
    release on AlphaFold DB gets picked up automatically.
    """
    df = read_table(in_path)
    if id_column not in df.columns:
        raise ValueError(f"Column '{id_column}' not found. Columns: {list(df.columns)}")

    raw = df[id_column].tolist()
    accs = sorted({a for a in (clean_accession(x) for x in raw) if a})

    out_base = Path(out_dir)
    out_base.mkdir(parents=True, exist_ok=True)
    logs_dir.mkdir(parents=True, exist_ok=True)

    report_rows = []

    with requests.Session() as s:
        s.headers.update({"User-Agent": "bulk-afdb-downloader/1.0"})
        with ThreadPoolExecutor(max_workers=_DOWNLOAD_MAX_WORKERS) as pool:
            futures = {
                pool.submit(_check_and_download, acc, out_base, prefer, also_pae, delay, s): acc
                for acc in accs
            }
            for future in tqdm(as_completed(futures), total=len(futures),
                                desc="Checking/downloading AlphaFold structures"):
                acc = futures[future]
                try:
                    row = future.result()
                except Exception as e:
                    row = {"UniProt": acc, "status": "ERROR", "structure_file": "",
                           "pae_file": "", "note": str(e)}
                report_rows.append(row)
                if row["status"] == "NO_ENTRY":
                    tqdm.write(f"  {acc}: no AlphaFold DB entry (404)")
                elif row["status"] == "NO_CANONICAL_MODEL":
                    tqdm.write(f"  {acc}: no canonical model (isoforms only)")
                elif row["status"] == "NO_STRUCTURE_URL":
                    tqdm.write(f"  {acc}: no structure URL found")
                elif row["status"] == "UPDATED":
                    tqdm.write(f"  {acc}: updated ({row['note']})")
                elif row["status"] == "ERROR":
                    tqdm.write(f"  {acc}: ERROR — {row['note']}")

    rep = pd.DataFrame(report_rows)
    rep_path = logs_dir / "download_report.tsv"
    rep.to_csv(rep_path, sep="\t", index=False)
    print(f"\nWrote report: {rep_path}")

    cached_count = (rep["status"] == "ALREADY_CACHED").sum()
    downloaded_count = (rep["status"] == "DOWNLOADED").sum()
    updated_count = (rep["status"] == "UPDATED").sum()
    print(f"Already up to date: {cached_count}  |  Newly downloaded: {downloaded_count}  |  "
          f"Updated to a newer AlphaFold version: {updated_count}")

    error_rows = rep[~rep["status"].isin({"DOWNLOADED", "ALREADY_CACHED", "UPDATED"})]
    if not error_rows.empty:
        err_path = logs_dir / "download_errors.tsv"
        error_rows.to_csv(err_path, sep="\t", index=False)
        print(f"Wrote error log ({len(error_rows)} failed): {err_path}")
    else:
        print("No errors — all accessions checked successfully.")


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
