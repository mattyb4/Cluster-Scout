"""Pipeline Step 2: download AlphaFold structures for every protein Step 1 found.

Reads UniProt accessions from a column of the given input table (Step 1's
hotspot TSV in the main pipeline, but usable standalone against any table),
and for each one queries the AlphaFold DB API, downloading the canonical
model's CIF (or PDB) file -- and, with --also_pae, its Predicted Aligned
Error JSON -- into out_dir/{accession}/. Every accession is checked against
AlphaFold DB's current version in parallel, even if a file with a matching
name already exists locally, so a newer release is picked up (and the stale
version removed) automatically rather than only downloading what's missing.
Writes a per-accession status report (download_report.tsv) and a filtered
error log (download_errors.tsv) to --logs_dir.
"""
from __future__ import annotations

import argparse
import os
import random
import re
import threading
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Any, Optional

import pandas as pd
import requests
from tqdm import tqdm

AF_PRED_ENDPOINT = "https://alphafold.ebi.ac.uk/api/prediction/{}"
AF_FILES_BASE = "https://alphafold.ebi.ac.uk/files"
ACC_RE = re.compile(r"^(?:[A-NR-Z0-9][A-Z0-9]{5}|[OPQ][0-9][A-Z0-9]{3}[0-9])(?:-\d+)?$")

# Be polite to the shared EBI-hosted API — same modest-concurrency convention
# used elsewhere in this project for other academic servers (kinase_library,
# AIUPred).
_DOWNLOAD_MAX_WORKERS = 2

# AFDB's /files host throttles bursts (HTTP 429, no Retry-After header): ~30
# requests in a few seconds trips it, while a steady ~2 req/s never does. The
# limiter below caps the *combined* request rate of all worker threads.
_DEFAULT_MAX_RPS = 2.0

# Statuses worth retrying. 403 is included because AFDB's CDN answers throttled
# or filtered clients with 403 as well as 429.
_RETRYABLE_STATUSES = frozenset({403, 408, 425, 429, 500, 502, 503, 504})
_THROTTLE_STATUSES = frozenset({403, 429})
# "No such record" statuses: 404 for an unknown accession, 400 for a malformed one.
_NO_RECORD_STATUSES = frozenset({400, 404})


class RateLimiter:
    """Thread-safe global request pacer shared by every worker.

    Spaces requests at least ``interval`` seconds apart. When the server
    throttles us (429/403) every thread pauses until a cool-down passes -- the
    server's Retry-After if sent, else an escalating default -- and the pacing
    interval is widened, then relaxes back toward the baseline as requests
    succeed. ``max_rps <= 0`` disables pacing entirely.
    """

    def __init__(self, max_rps: float = 0.0):
        self.base_interval = 1.0 / max_rps if max_rps > 0 else 0.0
        self.interval = self.base_interval
        self._lock = threading.Lock()
        self._next_slot = 0.0
        self._pause_until = 0.0
        self._strikes = 0

    def wait(self) -> None:
        """Block until this thread may send its next request."""
        with self._lock:
            now = time.monotonic()
            slot = max(now, self._next_slot, self._pause_until)
            self._next_slot = slot + self.interval
        delay = slot - time.monotonic()
        if delay > 0:
            time.sleep(delay)

    def success(self) -> None:
        with self._lock:
            self._strikes = 0
            if self.interval > self.base_interval:
                self.interval = max(self.base_interval, self.interval * 0.9)

    def throttled(self, retry_after: Optional[float] = None) -> float:
        """Record a throttle response; returns the cool-down (s) all threads now observe."""
        with self._lock:
            now = time.monotonic()
            if now < self._pause_until:
                # Several in-flight requests all got 429 from the same burst: they are
                # one throttle episode, so don't escalate (or re-pause) once per thread.
                return self._pause_until - now
            self._strikes += 1
            cooldown = retry_after if retry_after else min(5.0 * 2 ** (self._strikes - 1), 30.0)
            cooldown += random.uniform(0, 1.0)  # de-synchronise threads waking together
            self._pause_until = now + cooldown
            if self.base_interval:
                self.interval = min(max(self.interval, self.base_interval) * 1.25, 2.0)
            return cooldown


# Configured by main(); disabled by default so importing/unit-testing the module
# never sleeps on pacing.
_LIMITER = RateLimiter(0.0)


def _retry_after_seconds(r: Any) -> Optional[float]:
    """Parse a Retry-After header (delta-seconds form) if the response carries one."""
    try:
        return float((getattr(r, "headers", None) or {}).get("Retry-After"))
    except (TypeError, ValueError):
        return None


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

    # Safety net: if AFDB renames/removes the URL fields but still identifies the
    # model, rebuild the URL from its (stable, documented) file naming scheme.
    model_id = _record_model_id(record)
    version = record.get("latestVersion")
    if model_id and isinstance(version, int):
        if not structure_url:
            ext = "cif" if prefer == "cif" else "pdb"
            structure_url = f"{AF_FILES_BASE}/{model_id}-model_v{version}.{ext}"
        if not pae_url:
            pae_url = f"{AF_FILES_BASE}/{model_id}-predicted_aligned_error_v{version}.json"
    return {"structure_url": structure_url, "pae_url": pae_url}


def _record_model_id(record: dict) -> str:
    """AFDB model id (e.g. 'AF-P04637-F1') from whichever field carries it.

    ``entryId`` was deprecated in favour of ``modelEntityId`` (sunset
    2026-06-25); accept both, newest first, so either schema works.
    """
    for key in ("modelEntityId", "entryId"):
        v = record.get(key)
        if isinstance(v, str) and v.startswith("AF-"):
            return v
    return ""


def is_canonical_record(record: dict, acc: str) -> bool:
    """True if *record* models *acc* itself rather than one of its isoforms.

    Isoform records carry accessions like 'P11362-9'. Checks the UniProt
    accession field first, then falls back to the model id ('AF-P11362-F1' vs
    'AF-P11362-9-F1') so a renamed accession field doesn't silently turn every
    protein into 'isoform-only'.
    """
    if record.get("uniprotAccession") == acc:
        return True
    return re.fullmatch(rf"AF-{re.escape(acc)}-F\d+", _record_model_id(record)) is not None


def _get_with_retry(url: str, session: requests.Session, *, stream: bool, timeout: int,
                    retries: int) -> Optional[Any]:
    """GET *url* with global pacing, throttle-aware retry and jittered backoff.

    Returns the open 200 response (caller must close it / use it as a context
    manager), or None for a "no such record" status. Raises RuntimeError once
    retries are exhausted or immediately on a non-retryable error status.
    """
    last = "no attempt made"
    for attempt in range(retries):
        _LIMITER.wait()
        try:
            r = session.get(url, timeout=timeout, stream=stream)
        except requests.RequestException as e:
            last = f"{type(e).__name__}: {e}"
            time.sleep(min(1.6 ** attempt, 30) + random.uniform(0, 0.5))
            continue

        status = r.status_code
        if status == 200:
            _LIMITER.success()
            return r
        if status in _NO_RECORD_STATUSES:
            return None
        retry_after = _retry_after_seconds(r)
        if hasattr(r, "close"):
            r.close()
        last = f"HTTP {status}"
        if status not in _RETRYABLE_STATUSES:
            raise RuntimeError(f"{last} for {url}")
        if status in _THROTTLE_STATUSES:
            _LIMITER.throttled(retry_after)  # next wait() sleeps out the cool-down, for every thread
        else:
            time.sleep(min(1.6 ** attempt, 30) + random.uniform(0, 0.5))
    raise RuntimeError(f"{last} after {retries} attempts: {url}")


def fetch_prediction(acc: str, session: requests.Session, retries: int = 6) -> Optional[Any]:
    """Query the AlphaFold DB API for prediction metadata, with retry and backoff.

    Returns the parsed JSON, or None when AFDB has no record for *acc*.
    """
    r = _get_with_retry(AF_PRED_ENDPOINT.format(acc), session, stream=False, timeout=30,
                        retries=retries)
    if r is None:
        return None
    try:
        return r.json()
    except ValueError as e:
        raise RuntimeError(f"AFDB returned non-JSON for {acc} (API format change?): {e}") from e


def _looks_like_error_page(first_chunk: bytes) -> bool:
    """True if a "successful" body is actually an HTML error/interstitial page."""
    head = first_chunk[:256].lstrip().lower()
    return head.startswith((b"<!doctype", b"<html", b"<?xml"))


def download(url: str, outpath: Path, session: requests.Session, retries: int = 6) -> None:
    """Download a file via streaming GET with retry, using an atomic temp-file rename."""
    outpath.parent.mkdir(parents=True, exist_ok=True)
    if outpath.exists() and outpath.stat().st_size > 0:
        return
    tmp = outpath.with_suffix(outpath.suffix + ".part")
    last = "no attempt made"
    for attempt in range(retries):
        r = _get_with_retry(url, session, stream=True, timeout=90, retries=retries)
        if r is None:
            raise RuntimeError(f"File not found (404): {url}")
        try:
            with r:
                chunks = r.iter_content(chunk_size=1024 * 256)
                first = next((c for c in chunks if c), b"")
                if not first or _looks_like_error_page(first):
                    last = "empty or HTML body instead of a structure file"
                    continue
                with open(tmp, "wb") as f:
                    f.write(first)
                    for chunk in chunks:
                        if chunk:
                            f.write(chunk)
            os.replace(tmp, outpath)
            return
        except (requests.RequestException, OSError) as e:  # connection dropped mid-stream
            last = f"{type(e).__name__}: {e}"
        finally:
            if tmp.exists():
                tmp.unlink(missing_ok=True)
        time.sleep(min(1.6 ** attempt, 30) + random.uniform(0, 0.5))
    raise RuntimeError(f"Failed download after retries ({last}): {url}")


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
        canonical_records = [r for r in records if isinstance(r, dict) and is_canonical_record(r, acc)]
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


def main(in_path: str, id_column: str, out_dir: str, prefer: str, also_pae: bool, delay: float, logs_dir: Path,
         max_rps: float = _DEFAULT_MAX_RPS) -> None:
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

    global _LIMITER
    _LIMITER = RateLimiter(max_rps)

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

    # Canary for upstream API drift: a handful of NO_STRUCTURE_URL / non-JSON
    # results is normal, but a large share means AFDB changed its response shape.
    n_checked = len(rep)
    n_no_url = int((rep["status"] == "NO_STRUCTURE_URL").sum())
    n_bad_json = int(rep["note"].astype(str).str.contains("non-JSON", regex=False).sum())
    if n_checked and (n_no_url + n_bad_json) / n_checked > 0.05:
        print(f"WARNING: {n_no_url + n_bad_json}/{n_checked} accessions returned no usable structure URL or "
              f"non-JSON metadata. The AlphaFold DB API response format may have changed -- "
              f"see https://alphafold.ebi.ac.uk/api-docs and the EBI/PDBe AFDB announcements.")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("in_path", help="Input .tsv/.csv/.xlsx file")
    ap.add_argument("--id_column", default="UniProt", help="Column containing UniProt accessions")
    ap.add_argument("--out_dir", default="afdb_models", help="Output directory")
    ap.add_argument("--prefer", choices=["cif", "pdb"], default="cif", help="Prefer mmCIF or PDB when both available")
    ap.add_argument("--also_pae", action="store_true", help="Also download PAE json when url contains 'pae'")
    ap.add_argument("--delay", type=float, default=0.1, help="Polite delay between requests (seconds)")
    ap.add_argument("--logs_dir", default=str(DEFAULT_LOGS_DIR), help="Directory for download_report.tsv and download_errors.tsv")
    ap.add_argument("--max_rps", type=float, default=_DEFAULT_MAX_RPS,
                    help="Max combined requests/second across all workers (0 = unlimited); AFDB returns 429 on bursts")
    args = ap.parse_args()

    main(args.in_path, args.id_column, args.out_dir, args.prefer, args.also_pae, args.delay, Path(args.logs_dir),
         args.max_rps)
