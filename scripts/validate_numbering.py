"""Check cosmic_numbering.position_map against UniProt's own isoform annotation.

UniProt records every isoform as a list of "Alternative sequence" (VAR_SEQ)
features -- exact canonical start/end coordinates plus the replacement text.
Applying those edits to the canonical sequence gives an isoform->canonical
position map that doesn't depend on any alignment, so it is an independent
ground truth for the Needleman-Wunsch map used to renumber COSMIC mutations.

Positions inside a replaced block (a substitution rather than a pure
insertion/deletion) have no unambiguous truth and are reported separately.

Usage:  uv run python scripts/validate_numbering.py [--out DIR]
Reads data/cache/uniprot_isoform_sequences.tsv (written by step 1), fetches
VAR_SEQ annotation from UniProt (cached in data/cache/uniprot_var_seq.json)
and writes a summary plus a per-position disagreement table.
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from collections import Counter
from pathlib import Path

import requests

sys.path.insert(0, str(Path(__file__).resolve().parent))
import cosmic_numbering as cn  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
ISOFORM_CACHE = ROOT / "data" / "cache" / "uniprot_isoform_sequences.tsv"
VAR_SEQ_CACHE = ROOT / "data" / "cache" / "uniprot_var_seq.json"
_ENTRY_URL = "https://rest.uniprot.org/uniprotkb/{}.json"


def parse_var_seq(entry: dict) -> dict[str, list[tuple[int, int, str]]]:
    """UniProt JSON entry -> {isoform id: [(start, end, replacement), ...]}.
    Only isoforms built from VAR_SEQ edits are returned."""
    feats = {}
    for f in entry.get("features", []):
        if f.get("type") != "Alternative sequence":
            continue
        alts = f.get("alternativeSequence", {}).get("alternativeSequences", [""])
        feats[f["featureId"]] = (f["location"]["start"]["value"], f["location"]["end"]["value"],
                                 alts[0] if alts else "")
    out = {}
    for comment in entry.get("comments", []):
        for iso in comment.get("isoforms", []) if comment.get("commentType") == "ALTERNATIVE PRODUCTS" else []:
            ids = [i for i in iso.get("sequenceIds", []) if i in feats]
            if ids and iso.get("isoformIds"):
                out[iso["isoformIds"][0]] = sorted(feats[i] for i in ids)
    return out


def build_truth(canonical: str, edits: list[tuple[int, int, str]]) -> tuple[str, dict[int, int], set[int]]:
    """Apply VAR_SEQ *edits* to *canonical*.

    Returns (isoform sequence, {isoform pos: canonical pos} for residues in
    unedited stretches, isoform positions inside replaced blocks)."""
    seq, truth, replaced = [], {}, set()
    cpos = 1
    for start, end, alt in edits:
        for p in range(cpos, start):
            seq.append(canonical[p - 1])
            truth[len(seq)] = p
        for residue in alt:
            seq.append(residue)
            replaced.add(len(seq))
        cpos = end + 1
    for p in range(cpos, len(canonical) + 1):
        seq.append(canonical[p - 1])
        truth[len(seq)] = p
    return "".join(seq), truth, replaced


def compare(canonical: str, isoform: str, truth: dict[int, int], replaced: set[int]) -> dict:
    """Score the aligner's map against the ground truth for one isoform."""
    pmap, ambiguous = cn.position_map_with_ambiguity(isoform, canonical)
    wrong = {p: (c, pmap.get(p)) for p, c in truth.items() if pmap.get(p) != c}
    return {
        "truth_positions": len(truth),
        "agree": len(truth) - len(wrong),
        # The dangerous kind: placed on the wrong canonical residue (vs. left unmapped, which is safe)
        "wrong_mapped": sum(got is not None for _, got in wrong.values()),
        "wrong_mapped_flagged": sum(got is not None and p in ambiguous for p, (_, got) in wrong.items()),
        "ambiguous_total": len(ambiguous & truth.keys()),
        "wrong": wrong,                       # isoform pos -> (true canonical pos, aligner's)
        "replaced_positions": len(replaced),
        "replaced_mapped": sum(p in pmap for p in replaced),
    }


def fetch_var_seq(accessions, cache_path: Path = VAR_SEQ_CACHE, retries: int = 3) -> dict:
    cache = json.loads(cache_path.read_text()) if cache_path.exists() else {}
    session = requests.Session()
    for n, acc in enumerate(a for a in accessions if a not in cache):
        for attempt in range(retries):
            try:
                r = session.get(_ENTRY_URL.format(acc), params={"fields": "ft_var_seq,cc_alternative_products"}, timeout=60)
                if r.status_code == 200:
                    cache[acc] = {k: v for k, v in parse_var_seq(r.json()).items()}
                    break
                if r.status_code == 404:
                    cache[acc] = {}
                    break
            except requests.RequestException:
                pass
            time.sleep(2 ** attempt)
        if n % 50 == 0:
            cache_path.write_text(json.dumps(cache))
    cache_path.write_text(json.dumps(cache))
    return cache


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--out", type=Path, default=ROOT / "data" / "validation")
    args = ap.parse_args()

    seqs = cn.load_isoform_cache(ISOFORM_CACHE)
    proteins = sorted(a for a, s in seqs.items() if a in s and len(s) > 1)
    print(f"{len(proteins)} proteins with isoforms in cache; fetching VAR_SEQ annotation...")
    var_seq = fetch_var_seq(proteins)

    totals, rows, skipped = Counter(), [], Counter()
    for acc in proteins:
        canonical = seqs[acc][acc]
        for iso_id, edits in var_seq.get(acc, {}).items():
            iso_seq = seqs[acc].get(iso_id)
            if iso_seq is None:
                skipped["isoform not in cache"] += 1
                continue
            built, truth, replaced = build_truth(canonical, [tuple(e) for e in edits])
            if built != iso_seq:
                # Annotation and cached sequence disagree (version drift) -- no usable truth
                skipped["VAR_SEQ edits don't rebuild cached sequence"] += 1
                continue
            res = compare(canonical, iso_seq, truth, replaced)
            totals["isoforms"] += 1
            for k in ("truth_positions", "agree", "replaced_positions", "replaced_mapped",
                      "wrong_mapped", "wrong_mapped_flagged", "ambiguous_total"):
                totals[k] += res[k]
            if res["wrong"]:
                totals["isoforms_with_disagreement"] += 1
            for p, (true_c, got) in res["wrong"].items():
                rows.append((acc, iso_id, p, true_c, "" if got is None else got))

    args.out.mkdir(parents=True, exist_ok=True)
    with open(args.out / "position_map_disagreements.tsv", "w") as fh:
        fh.write("canonical_accession\tisoform\tisoform_pos\ttrue_canonical_pos\talignment_canonical_pos\n")
        for r in rows:
            fh.write("\t".join(map(str, r)) + "\n")

    t = totals
    pct = 100 * t["agree"] / t["truth_positions"] if t["truth_positions"] else float("nan")
    print(f"\nIsoforms checked:                 {t['isoforms']}  (skipped: {dict(skipped) or 'none'})")
    print(f"Positions with unambiguous truth: {t['truth_positions']}")
    print(f"  alignment agrees:               {t['agree']}  ({pct:.4f}%)")
    print(f"  alignment disagrees:            {len(rows)}  in {t['isoforms_with_disagreement']} isoforms")
    print(f"  placed on a wrong position:     {t['wrong_mapped']}  (rest were left unmapped, which is safe)")
    print(f"    flagged by the ambiguity check: {t['wrong_mapped_flagged']}  (unflagged = would be moved silently)")
    print(f"Correct positions also flagged:   {t['ambiguous_total'] - t['wrong_mapped_flagged']}  (cost of the check)")
    print(f"Positions in replaced blocks:     {t['replaced_positions']}  "
          f"(alignment placed {t['replaced_mapped']} of them on the canonical sequence)")
    print(f"Disagreements written to {args.out / 'position_map_disagreements.tsv'}")


if __name__ == "__main__":
    main()
