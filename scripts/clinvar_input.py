"""Read ClinVar's variant_summary file as a mutation source: missense variants
with a clinical classification, counted per (gene, protein change).

variant_summary.txt is several GB (every variant, on both GRCh37 and GRCh38),
so it is read in chunks and reduced once to a small cache of its GRCh38
missense variants (data/cache/clinvar_missense.tsv.gz, rebuilt whenever the
file changes). Choosing different classification classes afterwards reads only
that cache.

A ClinVar variant has one row, unlike COSMIC's one row per sample, so a
mutation's count here is its number of ClinVar submitters, summed over every
variant that gives the same protein change (e.g. two codon changes to the same
amino acid).
"""
from __future__ import annotations

import json
import re
from pathlib import Path

import pandas as pd
from pipeline_utils import AA3TO1

# Selectable classification classes -> the ClinVar classification terms each covers
CLINVAR_CLASSES: dict[str, tuple[str, ...]] = {
    "P": ("Pathogenic",),
    "LP": ("Likely pathogenic",),
    "VUS": ("Uncertain significance",),
    "LB": ("Likely benign",),
    "B": ("Benign",),
    "CONFLICTING": ("Conflicting classifications of pathogenicity",),
}
CLINVAR_CLASS_LABELS = {
    "P": "Pathogenic",
    "LP": "Likely pathogenic",
    "VUS": "Uncertain significance",
    "LB": "Likely benign",
    "B": "Benign",
    "CONFLICTING": "Conflicting classifications",
}
DEFAULT_CLINVAR_CLASSES = ("P", "LP")

CLINVAR_CACHE_FILE = "clinvar_missense.tsv.gz"
CLINVAR_CACHE_META = "clinvar_missense.json"
CLINVAR_CACHE_VERSION = 1
_CHUNK_ROWS = 500_000
_READ_COLUMNS = ["Name", "GeneSymbol", "ClinicalSignificance", "Assembly", "NumberSubmitters"]

# '... (p.Arg175His)' at the end of a variant's Name
_PROTEIN_CHANGE_RE = re.compile(r"\(p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})\)\s*$")
# 'NM_000546.6(TP53):c.524G>A ...' -> TP53 (the gene a multi-gene row's protein change is in)
_NAME_GENE_RE = re.compile(r"^[^(]*\(([^)]+)\):")
# selenocysteine/pyrrolysine aren't in PolyPhen/AlphaFold's alphabet
_STANDARD_AA3 = {three: one for three, one in AA3TO1.items() if one not in ("U", "O")}


class ClinVarInputError(Exception):
    """The ClinVar file is unusable (wrong format, or no variants survive filtering)."""


def parse_protein_changes(names: pd.Series) -> pd.Series:
    """One-letter missense changes ('R175H') from variant Names such as
    'NM_000546.6(TP53):c.524G>A (p.Arg175His)'; '' where the Name has no
    single-residue substitution (synonymous, stop, frameshift, indel...)."""
    parts = names.str.extract(_PROTEIN_CHANGE_RE)
    ref = parts[0].str.upper().map(_STANDARD_AA3)
    alt = parts[2].str.upper().map(_STANDARD_AA3)
    valid = ref.notna() & alt.notna() & (ref != alt)
    return (ref + parts[1] + alt).where(valid, "")


def _gene_symbols(chunk: pd.DataFrame) -> pd.Series:
    """GeneSymbol, except for variants spanning several genes ('A;B'), where
    the gene in the variant's Name is the one its protein change is in."""
    genes = chunk["GeneSymbol"]
    multi = genes.str.contains(";", regex=False)
    if multi.any():
        from_name = chunk.loc[multi, "Name"].str.extract(_NAME_GENE_RE)[0]
        genes = genes.copy()
        genes.loc[multi] = from_name.fillna("")
    return genes


def _read_variant_summary(path: Path) -> pd.DataFrame:
    """GRCh38 missense variants of a variant_summary file, one row per
    (gene, protein change, classification) with its submitter total."""
    pieces = []
    reader = pd.read_csv(path, sep="\t", usecols=_READ_COLUMNS, dtype=str,
                         keep_default_na=False, chunksize=_CHUNK_ROWS)
    for chunk in reader:
        chunk = chunk[chunk["Assembly"] == "GRCh38"]
        chunk = chunk.assign(aa_change=parse_protein_changes(chunk["Name"]))
        chunk = chunk[chunk["aa_change"] != ""]
        if chunk.empty:
            continue
        chunk = chunk.assign(gene=_gene_symbols(chunk))
        chunk = chunk[chunk["gene"] != ""]
        # a variant with no submitters recorded still has been submitted once
        submitters = pd.to_numeric(chunk["NumberSubmitters"], errors="coerce").fillna(1).clip(lower=1)
        pieces.append(chunk.assign(submitters=submitters.astype(int))
                      .rename(columns={"ClinicalSignificance": "significance"})
                      [["gene", "aa_change", "significance", "submitters"]])
    if not pieces:
        raise ClinVarInputError(
            f"No GRCh38 missense variants found in {path.name}. Is this ClinVar's variant_summary.txt?")
    return (pd.concat(pieces)
            .groupby(["gene", "aa_change", "significance"], as_index=False, sort=False)["submitters"].sum())


def load_clinvar_missense(path: Path, cache_dir: Path) -> pd.DataFrame:
    """The reduced missense table for *path*, from the cache when it was built
    from this exact file (same size and modification time)."""
    cache_file, meta_file = cache_dir / CLINVAR_CACHE_FILE, cache_dir / CLINVAR_CACHE_META
    stat = path.stat()
    fingerprint = {"version": CLINVAR_CACHE_VERSION, "name": path.name,
                   "size": stat.st_size, "mtime": int(stat.st_mtime)}
    try:
        if cache_file.exists() and json.loads(meta_file.read_text()) == fingerprint:
            return pd.read_csv(cache_file, sep="\t", dtype={"gene": str, "aa_change": str, "significance": str})
    except (OSError, ValueError):
        pass  # unreadable cache -> rebuild

    print(f"Reading {path.name} ({stat.st_size / 1e9:.1f} GB) -- one-time, cached for next runs...", flush=True)
    table = _read_variant_summary(path)
    cache_dir.mkdir(parents=True, exist_ok=True)
    table.to_csv(cache_file, sep="\t", index=False)
    meta_file.write_text(json.dumps(fingerprint))
    return table


def classes_selected(significance: pd.Series, classes) -> pd.Series:
    """True where a ClinVar classification belongs to one of *classes*
    (keys of CLINVAR_CLASSES). Combined classifications such as
    'Pathogenic/Likely pathogenic' count for each part."""
    wanted = {term for key in classes for term in CLINVAR_CLASSES[key]}
    # ~20 distinct classifications across millions of rows: decide each once
    verdict = {text: any(part in wanted for part in str(text).split("/")) for text in significance.unique()}
    return significance.map(verdict)


def load_clinvar_counts(path: Path, cache_dir: Path, classes=DEFAULT_CLINVAR_CLASSES):
    """Per-(gene, protein change) ClinVar submitter counts for *classes*.

    Returns (counts, gene_totals): counts has columns gene, aa_change,
    affected_cases; gene_totals maps gene -> submitters over all of its
    selected-class missense variants (the analogue of COSMIC's per-gene
    total-patient figure, kept before any hotspot threshold)."""
    unknown = [c for c in classes if c not in CLINVAR_CLASSES]
    if unknown:
        raise ClinVarInputError(f"Unknown ClinVar class(es): {', '.join(unknown)}")
    table = load_clinvar_missense(path, cache_dir)
    table = table[classes_selected(table["significance"], classes)]
    counts = (table.groupby(["gene", "aa_change"], as_index=False)["submitters"].sum()
              .rename(columns={"submitters": "affected_cases"}))
    gene_totals = counts.groupby("gene")["affected_cases"].sum().to_dict()
    return counts, gene_totals
