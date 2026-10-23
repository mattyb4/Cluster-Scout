"""Check which isoform cosmic_numbering picks against COSMIC's own transcript.

The Mutant Census file names the Ensembl protein (ENSP) each mutation was
called on, and COSMIC numbers the mutation along that protein. The pipeline
doesn't use this to choose an isoform (cross-references to UniProt were too
patchy), but it is an independent check on the choice that *is* made: fetch
the ENSP sequence and renumber the same mutations with it as the numbering
sequence, then compare labels with what the pipeline produced.

Ensembl serves only the current version of a protein, so the sequence is first
checked against COSMIC's own reference residues; a protein whose current
sequence no longer fits is reported as unusable rather than counted either way.

Usage:  uv run python scripts/validate_isoform_choice.py [--steps FILE] [--out DIR]
"""
from __future__ import annotations

import argparse
import json
import re
import sys
import time
from collections import Counter
from pathlib import Path

import pandas as pd
import requests

sys.path.insert(0, str(Path(__file__).resolve().parent))
import cosmic_numbering as cn  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
COSMIC_FILE = ROOT / "data" / "input" / "cosmic" / "Cosmic_MutantCensus_v104_GRCh38.tsv"
ISOFORM_CACHE = ROOT / "data" / "cache" / "uniprot_isoform_sequences.tsv"
ENSEMBL_CACHE = ROOT / "data" / "cache" / "ensembl_protein_sequences.json"
_ENSEMBL_URL = "https://rest.ensembl.org/sequence/id"
_BATCH = 50
_ENTRY_RE = re.compile(r"^(\S+) \((\d+)\)$")
MIN_ENSEMBL_FIT = 0.9


def cosmic_proteins(path: Path = COSMIC_FILE) -> pd.DataFrame:
    """Rows of (gene, transcript, ensp, n): the protein COSMIC numbers each gene's mutations on."""
    parts = []
    for chunk in pd.read_csv(path, sep="\t", usecols=["GENE_SYMBOL", "TRANSCRIPT_ACCESSION", "HGVSP"],
                             dtype=str, chunksize=500_000):
        chunk = chunk.dropna(subset=["HGVSP"])
        chunk["ENSP"] = chunk["HGVSP"].str.split(":").str[0]
        parts.append(chunk.groupby(["GENE_SYMBOL", "TRANSCRIPT_ACCESSION", "ENSP"]).size())
    out = pd.concat(parts).groupby(level=[0, 1, 2]).sum().reset_index(name="n")
    return out.rename(columns={"GENE_SYMBOL": "gene", "TRANSCRIPT_ACCESSION": "transcript", "ENSP": "ensp"})


def fetch_ensembl_proteins(ensps, cache_path: Path = ENSEMBL_CACHE, retries: int = 3) -> dict[str, str]:
    """{ENSP (unversioned): current sequence}; ids Ensembl no longer knows are cached as ''."""
    cache = json.loads(cache_path.read_text()) if cache_path.exists() else {}
    ids = sorted({e.split(".")[0] for e in ensps} - cache.keys())
    session = requests.Session()
    for i in range(0, len(ids), _BATCH):
        batch = ids[i:i + _BATCH]
        for attempt in range(retries):
            try:
                r = session.post(_ENSEMBL_URL, json={"ids": batch, "type": "protein"},
                                 headers={"Accept": "application/json"}, timeout=60)
                if r.status_code == 200:
                    got = {d["query"]: d["seq"] for d in r.json()}
                    cache.update({e: got.get(e, "") for e in batch})
                    break
                if r.status_code == 429:
                    time.sleep(float(r.headers.get("Retry-After", 2)))
                    continue
            except requests.RequestException:
                pass
            time.sleep(2 ** attempt)
        time.sleep(0.2)
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    cache_path.write_text(json.dumps(cache))
    return cache


def cosmic_labels(row) -> list[str]:
    """COSMIC's own mutation labels for a step-1 row (undoing any renumbering)."""
    new_to_old = dict(e.strip().split("=") for e in str(row["cosmic_mutation_labels"]).split(";") if "=" in e)
    labels = []
    for token in str(row["mutations_on_protein"]).split(";"):
        m = _ENTRY_RE.match(token.strip())
        if m:
            labels.append(new_to_old.get(m.group(1), m.group(1)))
    return labels


def reference_fit(seq: str, labels: list[str]) -> float:
    """Fraction of *labels* whose reference residue sits at their position in *seq*."""
    return cn.numbering_fit(seq, labels)


def compare_protein(uniprot: str, isoforms: dict[str, str], labels: list[str], ensembl_seq: str) -> dict:
    """Renumber *labels* with the pipeline's isoform choice and with Ensembl's
    sequence as the only candidate; return the label-level comparison."""
    pipeline = cn.number_protein(uniprot, isoforms, labels)
    reference = cn.number_protein(uniprot, {uniprot: isoforms[uniprot], "ENSEMBL": ensembl_seq}, labels)
    counts = Counter()
    disagreements = []
    for label in labels:
        p_unmapped, r_unmapped = label in pipeline.unmapped, label in reference.unmapped
        if p_unmapped and r_unmapped:
            kind = "both_unmapped"
        elif p_unmapped:
            kind = "pipeline_unmapped_only"       # recall loss: Ensembl places it, pipeline doesn't
        elif r_unmapped:
            kind = "ensembl_unmapped_only"
        elif pipeline.labels[label] == reference.labels[label]:
            kind = "same_label"
        else:
            kind = "different_label"              # the dangerous case
            disagreements.append((label, pipeline.labels[label], reference.labels[label]))
        counts[kind] += 1
    chosen_seq = isoforms.get(pipeline.isoform or uniprot, "")
    return {"pipeline_isoform": pipeline.isoform or "canonical", "reference_isoform": reference.isoform or "canonical",
            "chosen_equals_ensembl": chosen_seq == ensembl_seq, "counts": counts, "disagreements": disagreements}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--steps", type=Path, default=ROOT / "data" / "steps" / "COSMIC_hotspots_by_protein.tsv")
    ap.add_argument("--out", type=Path, default=ROOT / "data" / "validation")
    args = ap.parse_args()

    rows = pd.read_csv(args.steps, sep="\t", dtype=str, keep_default_na=False)
    isoforms = cn.load_isoform_cache(ISOFORM_CACHE)
    proteins = cosmic_proteins()
    per_gene = proteins.sort_values("n", ascending=False).drop_duplicates("gene").set_index("gene")["ensp"]
    print(f"{len(rows)} proteins; fetching {per_gene.reindex(rows['gene']).nunique()} Ensembl sequences...")
    ens = fetch_ensembl_proteins(per_gene.reindex(rows["gene"]).dropna().tolist())

    totals, records, disagreements, excluded = Counter(), [], [], Counter()
    for _, row in rows.iterrows():
        ensp = per_gene.get(row["gene"])
        iso = isoforms.get(row["uniprot_id"], {})
        seq = ens.get(ensp.split(".")[0], "") if isinstance(ensp, str) else ""
        labels = cosmic_labels(row)
        if not seq or row["uniprot_id"] not in iso:
            excluded["no Ensembl or UniProt sequence"] += 1
            continue
        fit = reference_fit(seq, labels)
        if fit < MIN_ENSEMBL_FIT:
            excluded["Ensembl sequence doesn't fit COSMIC's residues (version drift)"] += 1
            continue
        res = compare_protein(row["uniprot_id"], iso, labels, seq)
        totals["proteins"] += 1
        totals["chosen_equals_ensembl"] += res["chosen_equals_ensembl"]
        totals["isoform_choice_agrees"] += res["pipeline_isoform"] == res["reference_isoform"]
        totals.update(res["counts"])
        records.append({"gene": row["gene"], "uniprot": row["uniprot_id"], "ensp": ensp,
                        "ensembl_fit": round(fit, 4), "pipeline_isoform": res["pipeline_isoform"],
                        "ensembl_numbering_isoform": res["reference_isoform"],
                        "chosen_equals_ensembl": res["chosen_equals_ensembl"], **res["counts"]})
        disagreements += [(row["gene"], row["uniprot_id"], *d) for d in res["disagreements"]]

    args.out.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(records).to_csv(args.out / "isoform_choice_by_protein.tsv", sep="\t", index=False)
    pd.DataFrame(disagreements, columns=["gene", "uniprot", "cosmic_label", "pipeline_label", "ensembl_label"]
                 ).to_csv(args.out / "isoform_choice_label_disagreements.tsv", sep="\t", index=False)

    n = sum(totals[k] for k in ("same_label", "different_label", "both_unmapped",
                                "pipeline_unmapped_only", "ensembl_unmapped_only"))
    print(f"\nProteins compared: {totals['proteins']}  (excluded: {dict(excluded) or 'none'})")
    print(f"  pipeline's isoform choice == Ensembl-based choice: {totals['isoform_choice_agrees']}")
    print(f"  chosen sequence identical to Ensembl protein:      {totals['chosen_equals_ensembl']}")
    print(f"Mutations compared: {n}")
    for k in ("same_label", "both_unmapped", "pipeline_unmapped_only", "ensembl_unmapped_only", "different_label"):
        print(f"  {k:24s} {totals[k]:7d}  ({100 * totals[k] / n:.3f}%)")
    print(f"Written to {args.out}")


if __name__ == "__main__":
    main()
