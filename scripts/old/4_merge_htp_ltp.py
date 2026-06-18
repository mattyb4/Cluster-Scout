"""
Merge HTP/LTP scores from data/htp_ltp_scores.tsv into
ptm_mutation_proximity_db.tsv, joining on UniProt + ptm_site + ptm_type.

Adds two columns to the proximity db:
  - LTP_score
  - HTP_score

Writes the result back to ptm_mutation_proximity_db.tsv (in-place).
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent
OUTPUT_DIR = PROJECT_ROOT / "Output"

PROXIMITY_DB   = OUTPUT_DIR / "ptm_mutation_proximity_db.tsv"
HTP_LTP_FILE   = PROJECT_ROOT / "data" / "htp_ltp_scores.tsv"


def main() -> None:
    print(f"Reading proximity DB : {PROXIMITY_DB}")
    proximity = pd.read_csv(PROXIMITY_DB, sep="\t", encoding="utf-16")

    print(f"Reading HTP/LTP file : {HTP_LTP_FILE}")
    htp_ltp = pd.read_csv(HTP_LTP_FILE, sep="\t")

    htp_ltp.columns = htp_ltp.columns.str.strip()

    score_map = (
        htp_ltp[["UniProt", "ptm_site", "ptm_type", "LTP Score", "HTP Score"]]
        .dropna(subset=["UniProt", "ptm_site", "ptm_type"])
        .drop_duplicates(subset=["UniProt", "ptm_site", "ptm_type"])
        .set_index(["UniProt", "ptm_site", "ptm_type"])
    )

    def lookup(row, col):
        key = (row["UniProt"], row["ptm_site"], row["ptm_type"])
        if key in score_map.index:
            return score_map.at[key, col]
        return None

    proximity["LTP_score"] = proximity.apply(lambda r: lookup(r, "LTP Score"), axis=1)
    proximity["HTP_score"] = proximity.apply(lambda r: lookup(r, "HTP Score"), axis=1)

    matched = proximity["LTP_score"].notna().sum()
    total = len(proximity)
    print(f"Matched {matched}/{total} rows with HTP/LTP scores.")

    proximity.to_csv(PROXIMITY_DB, sep="\t", index=False, encoding="utf-16")
    print(f"Updated proximity DB written to: {PROXIMITY_DB}")


if __name__ == "__main__":
    main()
