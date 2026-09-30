"""Read a PhosphoSitePlus (PSP) bulk download as the PTM-site input for
ptm-proximity mode (the alternative to PTMD).

The user drops a PSP download into data/input/phosphositeplus/ as-is. Files
are identified by their COLUMN HEADERS, not their names: browsers rename
downloads, unzipping adds folder levels, and PSP may rename files between
releases, but the columns are what the pipeline actually depends on. Plain
text and .gz files are both read; anything unrecognized is listed and ignored.

Two file kinds are used:
  site dataset   (required, one per PTM type, e.g. Phosphorylation_site_dataset)
                 -- every known site with its LTP/HTP evidence counts
  disease sites  (optional, Disease-associated_sites) -- disease annotations

All site datasets share one format, and a row's PTM type comes from its
MOD_RSD suffix (S473-p -> Phosphorylation), so a new or missing PTM file needs
no code change. Unknown suffixes are counted and reported, never silently
dropped.
"""
from __future__ import annotations

import gzip
import io
import re
from dataclasses import dataclass, field
from pathlib import Path

import pandas as pd
from pipeline_utils import looks_like_uniprot_id

# MOD_RSD suffix -> PTM type name. Names match PTMD's where the two overlap
# (step 4 recognizes phosphosites by "phosphorylation" in the name).
PSP_MOD_TYPES = {
    "p": "Phosphorylation",
    "ub": "Ubiquitination",
    "ac": "Acetylation",
    "sm": "SUMOylation",
    "m1": "Monomethylation",
    "m2": "Dimethylation",
    "m3": "Trimethylation",
    "me": "Methylation",
    "gl": "O-GlcNAc glycosylation",
    "ga": "O-GalNAc glycosylation",
    # Seen only in PSP's annotation files so far, not its site datasets
    "pa": "Palmitoylation",
    "sc": "Succinylation",
    "ng": "N-linked glycosylation",
    "ca": "Caspase cleavage",
}

# Residues each PTM type is expected on. A site outside these is kept but
# counted as a warning -- it points to a data problem, not a pipeline bug.
_EXPECTED_RESIDUES = {
    "p": set("STYH"), "ub": set("K"), "ac": set("K"), "sm": set("K"),
    "m1": set("KR"), "m2": set("KR"), "m3": set("K"), "me": set("KR"),
    "gl": set("ST"), "ga": set("ST"),
}

SITE_COLUMNS = ("ACC_ID", "GENE", "ORGANISM", "MOD_RSD", "LT_LIT", "MS_LIT", "MS_CST")
DISEASE_COLUMNS = ("DISEASE", "ALTERATION", "ACC_ID", "ORGANISM", "MOD_RSD")

# Columns that identify each file kind. A site dataset must NOT have the
# annotation files' distinguishing columns, since those carry the same
# site/score columns too.
_KIND_SIGNATURES = {
    "disease": ({"DISEASE", "ALTERATION", "ACC_ID", "MOD_RSD"}, set()),
    "site": ({"ACC_ID", "MOD_RSD", "ORGANISM", "LT_LIT", "MS_LIT", "MS_CST"},
             {"DISEASE", "ON_FUNCTION", "KINASE", "SUB_MOD_RSD"}),
}

_TEXT_SUFFIXES = {"", ".txt", ".tsv", ".gz"}
_HEADER_SEARCH_LINES = 50
MOD_RSD_RE = re.compile(r"^([A-Z])(\d+)-([a-z0-9]+)$")

# If more than this fraction of human rows fails to parse, the format has
# probably changed -- stop rather than carry on with a chunk of data missing.
MAX_UNPARSEABLE_FRACTION = 0.05


class PspInputError(Exception):
    """The PSP input folder can't be used; the message says why."""


@dataclass
class PspFiles:
    site: list[Path] = field(default_factory=list)
    disease: list[Path] = field(default_factory=list)
    ignored: list[Path] = field(default_factory=list)


def _open_text(path: Path) -> io.StringIO:
    """Open a PSP file as text: gzip-aware, UTF-8 with a Latin-1 fallback
    (PSP files have shipped in both)."""
    raw = gzip.open(path, "rb").read() if path.suffix.lower() == ".gz" else path.read_bytes()
    try:
        text = raw.decode("utf-8")
    except UnicodeDecodeError:
        text = raw.decode("latin-1")
    return io.StringIO(text)


def _read_head(path: Path) -> list[str]:
    """First lines of a file, without decoding the whole thing -- used to
    identify files cheaply."""
    opener = gzip.open if path.suffix.lower() == ".gz" else open
    lines = []
    with opener(path, "rb") as fh:
        for _ in range(_HEADER_SEARCH_LINES):
            line = fh.readline()
            if not line:
                break
            lines.append(line.decode("latin-1").rstrip("\r\n"))
    return lines


def _find_header(lines: list[str]) -> tuple[int, list[str]] | None:
    """Index and column names of the header line: the first tab-separated
    line containing both ACC_ID and MOD_RSD. PSP files start with a date and
    a license paragraph whose length isn't guaranteed."""
    for i, line in enumerate(lines):
        cols = [c.strip() for c in line.split("\t")]
        if "ACC_ID" in cols and "MOD_RSD" in cols:
            return i, cols
    return None


def classify_file(path: Path) -> str | None:
    """Return "site", "disease", or None for a file that isn't one of those."""
    if path.suffix.lower() not in _TEXT_SUFFIXES:
        return None
    try:
        found = _find_header(_read_head(path))
    except (OSError, EOFError, gzip.BadGzipFile):
        return None
    if found is None:
        return None
    cols = set(found[1])
    for kind, (required, excluded) in _KIND_SIGNATURES.items():
        if required <= cols and not (excluded & cols):
            return kind
    return None


def discover_psp_files(folder: Path) -> PspFiles:
    """Classify every file under *folder* (recursively)."""
    files = PspFiles()
    if not folder.is_dir():
        return files
    for path in sorted(p for p in folder.rglob("*") if p.is_file()):
        if "__MACOSX" in path.parts or path.name.startswith("."):
            continue
        kind = classify_file(path)
        if kind == "site":
            files.site.append(path)
        elif kind == "disease":
            files.disease.append(path)
        else:
            files.ignored.append(path)
    return files


def read_psp_table(path: Path, required: tuple[str, ...]) -> pd.DataFrame:
    """Read a PSP file into a DataFrame of strings, starting at its header
    line. Raises PspInputError naming any missing required column."""
    handle = _open_text(path)
    lines = handle.getvalue().splitlines()
    found = _find_header(lines[:_HEADER_SEARCH_LINES])
    if found is None:
        raise PspInputError(
            f"{path.name}: no header line (with ACC_ID and MOD_RSD) in the first "
            f"{_HEADER_SEARCH_LINES} lines -- the file format may have changed"
        )
    header_idx, cols = found
    missing = [c for c in required if c not in cols]
    if missing:
        raise PspInputError(f"{path.name}: missing expected column(s) {', '.join(missing)}")
    handle.seek(0)
    df = pd.read_csv(handle, sep="\t", skiprows=header_idx, dtype=str,
                     keep_default_na=False, quoting=3)
    df.columns = [c.strip() for c in df.columns]
    return df


def _peek_mod_suffixes(path: Path) -> set[str]:
    """MOD_RSD suffixes in a site file's first data rows -- enough to tell
    which PTM type(s) it holds without reading the whole file."""
    lines = _read_head(path)
    found = _find_header(lines)
    if found is None:
        return set()
    header_idx, cols = found
    col = cols.index("MOD_RSD")
    suffixes = set()
    for line in lines[header_idx + 1:]:
        parts = line.split("\t")
        if len(parts) > col:
            m = MOD_RSD_RE.match(parts[col].strip())
            if m:
                suffixes.add(m.group(3))
    return suffixes


def validate_psp_folder(folder: Path) -> list[str]:
    """Folder- and file-level checks, cheap enough for the Pipeline tab's
    status display. Returns problem descriptions; empty means usable."""
    files = discover_psp_files(folder)
    if not files.site:
        found = ", ".join(p.name for p in files.ignored[:8]) or "nothing"
        return [f"No PhosphoSitePlus site dataset found in {folder} (found: {found}). "
                f"Expected files like Phosphorylation_site_dataset, with columns "
                f"{', '.join(SITE_COLUMNS)}."]

    problems = []
    owner: dict[str, Path] = {}
    for path in files.site:
        try:
            header = _find_header(_read_head(path))
        except OSError as exc:
            problems.append(f"{path.name}: could not be read ({exc})")
            continue
        cols = header[1] if header else []
        missing = [c for c in SITE_COLUMNS if c not in cols]
        if missing:
            problems.append(f"{path.name}: missing expected column(s) {', '.join(missing)}")
        for suffix in _peek_mod_suffixes(path):
            if suffix in owner and owner[suffix] != path:
                problems.append(
                    f"{owner[suffix].name} and {path.name} both contain "
                    f"{PSP_MOD_TYPES.get(suffix, suffix)} sites -- remove one "
                    f"(e.g. an older download) so sites aren't counted twice"
                )
            owner.setdefault(suffix, path)
    for path in files.disease:
        header = _find_header(_read_head(path))
        cols = header[1] if header else []
        missing = [c for c in DISEASE_COLUMNS if c not in cols]
        if missing:
            problems.append(f"{path.name}: missing expected column(s) {', '.join(missing)}")
    if len(files.disease) > 1:
        problems.append("More than one disease-associated sites file found: "
                        + ", ".join(p.name for p in files.disease) + " -- remove extras")
    return problems


def psp_folder_summary(folder: Path) -> str:
    """One-line description of a valid folder's contents for status displays."""
    files = discover_psp_files(folder)
    parts = [f"{len(files.site)} site file{'s' if len(files.site) != 1 else ''}"]
    parts.append("disease file" if files.disease else "no disease file")
    return ", ".join(parts)


def _to_count(series: pd.Series) -> tuple[pd.Series, int]:
    """PSP evidence counts: blank means 0. Returns (ints, n_non_numeric)."""
    stripped = series.str.strip()
    numeric = pd.to_numeric(stripped, errors="coerce")
    bad = int((numeric.isna() & (stripped != "")).sum())
    return numeric.fillna(0).astype(int), bad


def load_psp_sites(folder: Path, log=print) -> pd.DataFrame:
    """Load every human site from the folder's site datasets.

    Returns one row per (accession, MOD_RSD) with columns: UniProt, psp_gene,
    residue, position, ptm_type, ltp, htp, cst. Logs a per-file summary and
    skipped-row counts by reason; raises PspInputError on anything that
    makes the input unusable.
    """
    problems = validate_psp_folder(folder)
    if problems:
        raise PspInputError("\n".join(problems))
    files = discover_psp_files(folder)

    frames = []
    for path in files.site:
        df = read_psp_table(path, SITE_COLUMNS)
        human = df[df["ORGANISM"].str.strip().str.lower() == "human"].copy()
        n_human = len(human)
        if n_human == 0:
            log(f"  {path.name}: no human rows -- skipped")
            continue

        parsed = human["MOD_RSD"].str.strip().str.extract(MOD_RSD_RE)
        parsed.columns = ["residue", "position", "suffix"]
        human = human.join(parsed)

        skipped: dict[str, int] = {}
        unparseable = human["suffix"].isna()
        skipped["unparseable MOD_RSD"] = int(unparseable.sum())
        if n_human and unparseable.mean() > MAX_UNPARSEABLE_FRACTION:
            raise PspInputError(
                f"{path.name}: {int(unparseable.sum())} of {n_human} human rows have an "
                f"unrecognized MOD_RSD format (expected e.g. S473-p) -- the file format "
                f"may have changed"
            )
        human = human[~unparseable]

        unknown = ~human["suffix"].isin(PSP_MOD_TYPES)
        for suffix, n in human.loc[unknown, "suffix"].value_counts().items():
            skipped[f"unknown PTM type suffix '-{suffix}'"] = int(n)
        human = human[~unknown]

        # Isoform accessions (P31946-2) number positions on that isoform, which
        # can't be mapped onto the canonical sequence AlphaFold models.
        isoform = human["ACC_ID"].str.strip().str.contains("-", regex=False)
        skipped["isoform accession (not canonical sequence)"] = int(isoform.sum())
        human = human[~isoform]

        # A few rows carry non-UniProt IDs (e.g. GenBank's AAA58698), which
        # have no AlphaFold model and which UniProt's API rejects outright.
        not_uniprot = ~human["ACC_ID"].str.strip().map(looks_like_uniprot_id)
        skipped["not a UniProt accession"] = int(not_uniprot.sum())
        human = human[~not_uniprot]

        n_mismatched = sum(
            int((~group["residue"].isin(_EXPECTED_RESIDUES[suffix])).sum())
            for suffix, group in human.groupby("suffix") if suffix in _EXPECTED_RESIDUES
        )

        counts = {}
        for col, name in (("LT_LIT", "ltp"), ("MS_LIT", "htp"), ("MS_CST", "cst")):
            counts[name], n_bad = _to_count(human[col])
            if n_bad:
                skipped[f"non-numeric {col} (treated as 0)"] = n_bad

        out = pd.DataFrame({
            "UniProt": human["ACC_ID"].str.strip(),
            "psp_gene": human["GENE"].str.strip(),
            "residue": human["residue"],
            "position": human["position"].astype(int),
            "ptm_type": human["suffix"].map(PSP_MOD_TYPES),
            **counts,
        })
        frames.append(out)

        types = ", ".join(f"{t} {n:,}" for t, n in out["ptm_type"].value_counts().items())
        log(f"  {path.name}: {len(out):,} human sites used ({types})")
        for reason, n in skipped.items():
            if n:
                log(f"    skipped {n:,}: {reason}")
        if n_mismatched:
            log(f"    warning: {n_mismatched:,} sites on an unexpected residue for their "
                f"PTM type (kept)")

    if not frames:
        raise PspInputError("No human sites found in any PhosphoSitePlus site dataset")
    sites = pd.concat(frames, ignore_index=True)
    # A site listed twice (e.g. once per file) keeps its highest evidence counts
    sites = (sites.groupby(["UniProt", "residue", "position", "ptm_type"], as_index=False)
             .agg(psp_gene=("psp_gene", "first"), ltp=("ltp", "max"),
                  htp=("htp", "max"), cst=("cst", "max")))
    return sites


def filter_psp_sites(sites: pd.DataFrame, diseases: pd.DataFrame, min_ltp: int = 0,
                     max_ltp: int | None = None, min_htp: int = 0, max_htp: int | None = None,
                     disease_only: bool = False, log=print) -> pd.DataFrame:
    """Keep sites whose LTP and HTP counts fall in [min, max] (None = no upper
    bound) and, with *disease_only*, that PSP lists at least one disease for.
    *sites*/*diseases* are load_psp_sites()/load_psp_diseases() output."""
    keep = sites["ltp"].ge(min_ltp) & sites["htp"].ge(min_htp)
    if max_ltp is not None:
        keep &= sites["ltp"].le(max_ltp)
    if max_htp is not None:
        keep &= sites["htp"].le(max_htp)
    if disease_only:
        if diseases.empty:
            raise PspInputError(
                "'Disease-associated sites only' needs PhosphoSitePlus's Disease-associated_sites "
                "file, which isn't in the input folder"
            )
        site_key = ["UniProt", "residue", "position", "ptm_type"]
        with_disease = pd.MultiIndex.from_frame(diseases[site_key].drop_duplicates())
        keep &= pd.MultiIndex.from_frame(sites[site_key]).isin(with_disease)

    def _range(lo, hi):
        return f"{lo}-{hi}" if hi is not None else f"{lo}+"

    filtered = sites[keep].reset_index(drop=True)
    log(f"PhosphoSitePlus filters: LTP {_range(min_ltp, max_ltp)}, HTP {_range(min_htp, max_htp)}, "
        f"disease-associated only: {'yes' if disease_only else 'no'} -- "
        f"kept {len(filtered):,} of {len(sites):,} sites")
    return filtered


def load_psp_diseases(folder: Path, log=print) -> pd.DataFrame:
    """Load human disease associations from the optional disease-associated
    sites file. Returns columns UniProt, residue, position, ptm_type, disease
    (with PSP's alteration appended, e.g. "Breast cancer (increased)"), or an
    empty frame if the file is absent."""
    empty = pd.DataFrame(columns=["UniProt", "residue", "position", "ptm_type", "disease"])
    files = discover_psp_files(folder)
    if not files.disease:
        log("  No disease-associated sites file -- disease columns will be blank")
        return empty
    path = files.disease[0]
    df = read_psp_table(path, DISEASE_COLUMNS)
    df = df[df["ORGANISM"].str.strip().str.lower() == "human"].copy()
    parsed = df["MOD_RSD"].str.strip().str.extract(MOD_RSD_RE)
    parsed.columns = ["residue", "position", "suffix"]
    df = df.join(parsed)
    df = df[df["suffix"].isin(PSP_MOD_TYPES) & df["DISEASE"].str.strip().ne("")]
    df = df[df["ACC_ID"].str.strip().map(looks_like_uniprot_id)]
    if df.empty:
        return empty
    alteration = df["ALTERATION"].str.strip()
    disease = df["DISEASE"].str.strip().where(
        alteration == "", df["DISEASE"].str.strip() + " (" + alteration + ")"
    )
    out = pd.DataFrame({
        "UniProt": df["ACC_ID"].str.strip(),
        "residue": df["residue"],
        "position": df["position"].astype(int),
        "ptm_type": df["suffix"].map(PSP_MOD_TYPES),
        "disease": disease,
    }).drop_duplicates()
    log(f"  {path.name}: {len(out):,} human site-disease associations "
        f"({out['disease'].nunique():,} distinct)")
    return out
