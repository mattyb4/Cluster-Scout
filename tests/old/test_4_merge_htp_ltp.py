"""Unit tests for scripts/4_merge_htp_ltp.py."""
import pandas as pd
import pytest


@pytest.fixture
def mod(merge_module):
    return merge_module


def test_main_merges_scores_on_uniprot_site_and_type(mod, tmp_path, monkeypatch):
    proximity_path = tmp_path / "proximity.tsv"
    htp_ltp_path = tmp_path / "htp_ltp.tsv"

    proximity_df = pd.DataFrame({
        "UniProt": ["P04637", "Q06124"],
        "gene": ["TP53", "PTPN11"],
        "ptm_site": ["S15", "T59"],
        "ptm_type": ["Phosphorylation", "Phosphorylation"],
    })
    proximity_df.to_csv(proximity_path, sep="\t", index=False, encoding="utf-16")

    # HTP/LTP file has stray whitespace in column names, as the real file does.
    htp_ltp_df = pd.DataFrame({
        " UniProt": ["P04637"],
        "ptm_site ": ["S15"],
        " ptm_type": ["Phosphorylation"],
        "LTP Score": [2],
        "HTP Score": [5],
    })
    htp_ltp_df.to_csv(htp_ltp_path, sep="\t", index=False)

    monkeypatch.setattr(mod, "PROXIMITY_DB", proximity_path)
    monkeypatch.setattr(mod, "HTP_LTP_FILE", htp_ltp_path)

    mod.main()

    result = pd.read_csv(proximity_path, sep="\t", encoding="utf-16")
    tp53_row = result.loc[result["UniProt"] == "P04637"].iloc[0]
    ptpn11_row = result.loc[result["UniProt"] == "Q06124"].iloc[0]

    assert tp53_row["LTP_score"] == 2
    assert tp53_row["HTP_score"] == 5
    assert pd.isna(ptpn11_row["LTP_score"])
    assert pd.isna(ptpn11_row["HTP_score"])
