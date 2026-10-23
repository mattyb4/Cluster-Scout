"""Put COSMIC mutation positions onto the canonical UniProt sequence that
AlphaFold models.

COSMIC numbers each mutation along its own transcript, which is sometimes a
different UniProt isoform than the canonical one -- so "R248Q" in COSMIC can be
residue 249 (or an exon the canonical form doesn't have at all) in the
structure. Transcript cross-references proved unreliable for finding that
isoform (COSMIC's transcript IDs are often retired or re-versioned), so the
mutations themselves decide: the isoform whose residues agree with the most of
COSMIC's reference residues is taken as COSMIC's numbering. That isoform is
then aligned to the canonical sequence, and each mutation is moved to its
aligned canonical position -- but only if the canonical residue there is the
mutation's reference residue. Mutations that can't be placed that way keep
their COSMIC label and are reported as unmapped -- as are mutations whose
position equally good alignments place differently (gaps sliding, repeats) --
and every mutation of a protein that no isoform's numbering explains (>1% of
reference residues unmatched).

The full methodology, with rationale, validation figures and limitations, is
in docs/methods.md -- keep it in sync with any change here.
"""
from __future__ import annotations

import re
import time
from dataclasses import dataclass, field
from itertools import islice
from pathlib import Path

import pandas as pd
import requests
from Bio.Align import PairwiseAligner

MUTATION_RE = re.compile(r"^([A-Z])(\d+)([A-Z*])$")

# Global alignment tuned for isoforms: identical sequence broken up by whole
# inserted/skipped exons, so gaps are cheap to extend but costly to open.
_ALIGNER = PairwiseAligner(mode="global", match_score=2, mismatch_score=-3,
                           open_gap_score=-5, extend_gap_score=-0.1)

# A protein whose best numbering sequence still leaves more than 1% of COSMIC's
# reference residues unexplained is numbered on a transcript UniProt doesn't have
# (docs/methods.md); its mutations can't be placed reliably
MIN_NUMBERING_FIT = 0.99

# Equally-scored alignments examined when looking for positions whose canonical
# partner isn't pinned down (see ambiguous_positions)
_MAX_ALIGNMENTS = 256

_UNIPROT_STREAM_URL = "https://rest.uniprot.org/uniprotkb/stream"
_BATCH_SIZE = 25
CACHE_COLUMNS = ["query_accession", "accession", "sequence"]


@dataclass
class ProteinNumbering:
    """COSMIC -> canonical numbering for one protein's mutations."""
    isoform: str = ""                                   # isoform COSMIC numbers by; "" = canonical
    labels: dict[str, str] = field(default_factory=dict)   # COSMIC label -> label used in outputs
    unmapped: list[str] = field(default_factory=list)       # COSMIC labels with no canonical position
    ambiguous: list[str] = field(default_factory=list)      # subset of unmapped: equal-scoring alignments disagree
    poor_fit: bool = False                                  # no isoform explains the reference residues: all unmapped


def parse_fasta(text: str) -> dict[str, str]:
    """UniProt FASTA -> {accession: sequence} (isoforms keep their -N suffix)."""
    seqs: dict[str, list[str]] = {}
    current = None
    for line in text.splitlines():
        if line.startswith(">"):
            parts = line[1:].split("|")
            current = parts[1] if len(parts) > 1 else None
            if current:
                seqs[current] = []
        elif current:
            seqs[current].append(line.strip())
    return {acc: "".join(lines) for acc, lines in seqs.items()}


def load_isoform_cache(path: Path) -> dict[str, dict[str, str]]:
    """{canonical accession: {accession: sequence}}; an accession UniProt
    returned nothing for is cached with no sequences."""
    if not path.exists():
        return {}
    df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    cache: dict[str, dict[str, str]] = {}
    for query, acc, seq in zip(df["query_accession"], df["accession"], df["sequence"]):
        entry = cache.setdefault(query, {})
        if acc:
            entry[acc] = seq
    return cache


def save_isoform_cache(path: Path, cache: dict[str, dict[str, str]]) -> None:
    rows = []
    for query, seqs in cache.items():
        rows += [(query, acc, seq) for acc, seq in seqs.items()] or [(query, "", "")]
    path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows, columns=CACHE_COLUMNS).to_csv(path, sep="\t", index=False)


def _fetch_fasta(session, accessions: list[str], retries: int, backoff: float) -> str | None:
    """UniProt FASTA (isoforms included) for *accessions*, or None if the
    request still fails after *retries* attempts."""
    query = " OR ".join(f"accession:{a}" for a in accessions)
    for attempt in range(retries):
        try:
            resp = session.get(_UNIPROT_STREAM_URL, params={
                "query": f"({query})", "includeIsoform": "true", "format": "fasta"}, timeout=60)
            if resp.status_code == 200:
                return resp.text
        except requests.RequestException:
            pass
        if attempt < retries - 1:
            time.sleep(backoff * 2 ** attempt)
    return None


def fetch_isoform_sequences(accessions, cache_path: Path, progress=None, session=None,
                            retries: int = 3, backoff: float = 2.0) -> dict[str, dict[str, str]]:
    """Every UniProt isoform sequence (canonical included) for each accession,
    cached so repeat runs only fetch accessions not seen before.

    UniProt's API occasionally fails a request (HTTP 500s have been seen), so
    a failing batch is retried, then fetched one accession at a time. An
    accession that still fails is left out -- not cached, so the next run
    retries it -- and its mutations simply keep COSMIC's numbering.
    """
    cache = load_isoform_cache(cache_path)
    missing = sorted({a for a in accessions if a not in cache})
    print(f"{len(set(accessions)) - len(missing)}/{len(set(accessions))} proteins' isoform sequences "
          f"found in cache; fetching {len(missing)} new...")
    session = session or requests.Session()
    batches = [missing[i:i + _BATCH_SIZE] for i in range(0, len(missing), _BATCH_SIZE)]
    failed: list[str] = []

    def _store(batch, text):
        for acc, seq in parse_fasta(text).items():
            base = acc.split("-")[0]
            if base in batch:
                cache.setdefault(base, {})[acc] = seq
        for acc in batch:
            cache.setdefault(acc, {})

    try:
        for n, batch in enumerate(batches, 1):
            text = _fetch_fasta(session, batch, retries, backoff)
            if text is not None:
                _store(batch, text)
            else:
                for acc in batch:
                    single = _fetch_fasta(session, [acc], retries, backoff)
                    if single is not None:
                        _store([acc], single)
                    else:
                        failed.append(acc)
            if progress:
                progress(n / len(batches) * 100, f"Fetching isoform sequences: batch {n}/{len(batches)}")
    finally:
        if missing:
            save_isoform_cache(cache_path, cache)
    if failed:
        print(f"  Warning: UniProt didn't return isoform sequences for {len(failed)} protein(s) "
              f"({', '.join(failed[:10])}{'...' if len(failed) > 10 else ''}) -- their mutations keep "
              f"COSMIC's numbering this run; they'll be retried next run")
    return cache


def _matches(seq: str, ref: str, pos: int) -> bool:
    return 0 < pos <= len(seq) and seq[pos - 1] == ref


def choose_numbering_isoform(canonical: str, isoforms: dict[str, str], mutations: list[str]) -> str:
    """The isoform whose residues agree with the most mutations' reference
    residues. The canonical sequence wins ties -- another isoform must explain
    strictly more of COSMIC's mutations to be chosen."""
    parsed = [(m.group(1), int(m.group(2))) for m in map(MUTATION_RE.match, mutations) if m]

    def score(seq):
        return sum(_matches(seq, ref, pos) for ref, pos in parsed)

    best, best_score = canonical, score(isoforms.get(canonical, ""))
    for acc, seq in sorted(isoforms.items()):
        if acc != canonical and score(seq) > best_score:
            best, best_score = acc, score(seq)
    return best


def _alignment_map(alignment) -> dict[int, int]:
    mapping = {}
    for (i0, i1), (c0, _c1) in zip(*alignment.aligned):
        for k in range(i1 - i0):
            mapping[i0 + k + 1] = c0 + k + 1
    return mapping


def numbering_fit(seq: str, mutations: list[str]) -> float:
    """Fraction of the mutations' reference residues found at their stated
    positions in *seq* (NaN if no label parses)."""
    parsed = [(m.group(1), int(m.group(2))) for m in map(MUTATION_RE.match, mutations) if m]
    return sum(_matches(seq, ref, pos) for ref, pos in parsed) / len(parsed) if parsed else float("nan")


def numbering_fit(seq: str, mutations: list[str]) -> float:
    """Fraction of the mutations' reference residues found at their stated
    positions in *seq* (NaN if no label parses)."""
    parsed = [(m.group(1), int(m.group(2))) for m in map(MUTATION_RE.match, mutations) if m]
    return sum(_matches(seq, ref, pos) for ref, pos in parsed) / len(parsed) if parsed else float("nan")


def position_map(isoform_seq: str, canonical_seq: str) -> dict[int, int]:
    """1-based isoform position -> 1-based canonical position, for every
    aligned (not gapped) column of the global alignment."""
    return _alignment_map(_ALIGNER.align(isoform_seq, canonical_seq)[0])


def position_map_with_ambiguity(isoform_seq: str, canonical_seq: str) -> tuple[dict[int, int], set[int]]:
    """position_map plus the isoform positions that equally good alignments
    place differently (or place in a gap in some and not others).

    A gap can slide when the residues at its edges match, and in repeats it
    can land in a different unit at the same score. A mutation there passes
    the reference-residue check at either canonical position, so the check
    can't catch a wrong move -- only comparing optimal alignments can. At most
    _MAX_ALIGNMENTS optimal alignments are compared."""
    alignments = _ALIGNER.align(isoform_seq, canonical_seq)
    maps = [_alignment_map(a) for a in islice(alignments, _MAX_ALIGNMENTS)]
    first = maps[0]
    ambiguous = {p for m in maps[1:] for p in first.keys() | m.keys() if first.get(p) != m.get(p)}
    return first, ambiguous


def number_protein(canonical: str, isoforms: dict[str, str], mutations: list[str]) -> ProteinNumbering:
    """Work out where each COSMIC mutation label sits on the canonical sequence."""
    result = ProteinNumbering()
    canonical_seq = isoforms.get(canonical)
    if not canonical_seq:
        # No sequence (e.g. an obsolete accession) -- nothing to check against
        result.labels = {m: m for m in mutations}
        return result

    chosen = choose_numbering_isoform(canonical, isoforms, mutations)
    if numbering_fit(isoforms[chosen], mutations) < MIN_NUMBERING_FIT:
        # Even the best isoform leaves COSMIC's residues unexplained, so COSMIC numbered on
        # a transcript UniProt doesn't hold. A residue can still match by chance at a wrong
        # position, and the residue check can't tell -- keep every label, flag every one
        result.poor_fit = True
        result.labels = {m: m for m in mutations}
        result.unmapped = list(mutations)
        return result
    pmap, ambiguous = (position_map_with_ambiguity(isoforms[chosen], canonical_seq)
                       if chosen != canonical else (None, set()))
    if pmap is not None:
        result.isoform = chosen

    for label in mutations:
        m = MUTATION_RE.match(label)
        if not m:
            result.labels[label] = label
            continue
        ref, pos, alt = m.group(1), int(m.group(2)), m.group(3)
        mapped = pmap.get(pos) if pmap is not None else pos
        if pos in ambiguous:
            # Equal-scoring alignments disagree on where this residue goes, and the
            # residue check can't tell them apart: don't guess, keep COSMIC's label
            result.labels[label] = label
            result.unmapped.append(label)
            result.ambiguous.append(label)
        elif mapped is not None and _matches(canonical_seq, ref, mapped):
            result.labels[label] = f"{ref}{mapped}{alt}"
        elif _matches(canonical_seq, ref, pos):
            # Already right in canonical numbering (the alignment disagrees,
            # but the residue check is the stronger evidence)
            result.labels[label] = label
        else:
            result.labels[label] = label
            result.unmapped.append(label)
    return result
