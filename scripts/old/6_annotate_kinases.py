"""Annotate phosphorylation sites with predicted upstream kinases.

Uses the Kinase Library to predict the top 5 most likely kinases for each
Ser/Thr/Tyr phosphorylation site, based on a ±7 residue sequence window
extracted from the AlphaFold CIF structure.

Adds one column to the proximity database:

  kinase_predictions — Top 5 kinases formatted as
                       "KINASE(log2_score,percentile%); ..."
"""
from __future__ import annotations

import sys
import warnings
from pathlib import Path

import pandas as pd
from tqdm import tqdm

sys.path.insert(0, str(Path(__file__).resolve().parent))
from pipeline_utils import (  # noqa: E402
    project_root, AA3TO1, SITE_RE,
    find_canonical_cif, load_first_chain,
)

PROJECT_ROOT = project_root(__file__)
MODELS_ROOT = PROJECT_ROOT / "cif_models"
PROXIMITY_DB = PROJECT_ROOT / "Output" / "ptm_mutation_proximity_db.tsv"

_TOP_K = 5
_WINDOW = 7  # residues on each side of the phosphosite


def extract_sequence(chain) -> dict[int, str]:
    """Build a {position: one_letter_aa} dict from a biotite chain's CA atoms."""
    ca_mask = chain.atom_name == "CA"
    ca_atoms = chain[ca_mask]
    return {
        int(ca_atoms.res_id[i]): AA3TO1.get(str(ca_atoms.res_name[i]), "X")
        for i in range(len(ca_atoms))
    }


def build_window(pos_to_aa: dict[int, str], site_pos: int) -> str | None:
    """Build a 15-mer sequence window centered on site_pos with lowercase phosphosite."""
    residue = pos_to_aa.get(site_pos)
    if not residue or residue not in ("S", "T", "Y"):
        return None

    chars = []
    for offset in range(-_WINDOW, _WINDOW + 1):
        p = site_pos + offset
        if offset == 0:
            chars.append(residue.lower())
        else:
            chars.append(pos_to_aa.get(p, "_"))

    return "".join(chars)


def predict_kinases(window: str) -> str:
    """Run the Kinase Library on a 15-mer window and return a formatted top-5 string."""
    import kinase_library as kl

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        sub = kl.Substrate(window)
        result = sub.predict()

    top = result.head(_TOP_K)
    parts = []
    for kinase, row in top.iterrows():
        parts.append(f"{kinase}({row['Score']:.2f},{row['Percentile']:.1f}%)")
    return "; ".join(parts)


def main() -> None:
    """Annotate phosphorylation sites in the proximity DB with kinase predictions."""
    print(f"Reading proximity DB: {PROXIMITY_DB}")
    df = pd.read_csv(PROXIMITY_DB, sep="\t", encoding="utf-16", dtype=str,
                     keep_default_na=False)

    # Group rows by UniProt to load each CIF only once
    unique_uniprots = df["UniProt"].unique().tolist()
    print(f"{len(unique_uniprots)} unique proteins to process")

    # Build sequence maps from CIF files
    seq_maps: dict[str, dict[int, str]] = {}
    skipped = 0
    for uid in unique_uniprots:
        uniprot_dir = MODELS_ROOT / uid
        cif_file = find_canonical_cif(uniprot_dir) if uniprot_dir.is_dir() else None
        if cif_file is None:
            skipped += 1
            continue
        chain = load_first_chain(cif_file)
        if chain is None:
            skipped += 1
            continue
        seq_maps[uid] = extract_sequence(chain)

    if skipped:
        print(f"  Skipped {skipped} proteins (no CIF or unparseable)")
    print(f"  Loaded sequences for {len(seq_maps)} proteins")

    # Predict kinases for each phosphorylation site
    predictions: list[str] = []
    annotated = 0
    for _, row in tqdm(df.iterrows(), total=len(df), desc="Predicting kinases"):
        uid = row.get("UniProt", "")
        ptm_site = row.get("ptm_site", "")
        ptm_type = row.get("ptm_type", "")

        # Only predict for phosphorylation sites
        if "phosphorylation" not in ptm_type.lower():
            predictions.append("")
            continue

        m = SITE_RE.match(ptm_site.strip()) if ptm_site else None
        if not m:
            predictions.append("")
            continue

        pos = int(m.group(2))
        pos_to_aa = seq_maps.get(uid)
        if pos_to_aa is None:
            predictions.append("")
            continue

        window = build_window(pos_to_aa, pos)
        if window is None:
            predictions.append("")
            continue

        try:
            pred = predict_kinases(window)
            predictions.append(pred)
            annotated += 1
        except Exception:
            predictions.append("")

    df["kinase_predictions"] = predictions

    print(f"Annotated {annotated}/{len(df)} rows with kinase predictions.")
    df.to_csv(PROXIMITY_DB, sep="\t", index=False, encoding="utf-16")
    print(f"Updated proximity DB written to: {PROXIMITY_DB}")


if __name__ == "__main__":
    main()
