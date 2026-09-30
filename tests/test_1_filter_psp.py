"""Tests for step 1's PhosphoSitePlus source (_build_psp_ptm_table in
scripts/1_filter.py): PSP sites become the same per-protein columns the PTMD
source produces, plus psp_site_scores. The UniProt gene lookup is stubbed."""
import pandas as pd
import pytest
from test_psp_input import DISEASE_HEADER, SITE_HEADER, disease_row, site_row, write


@pytest.fixture
def mod(filter_module):
    return filter_module


@pytest.fixture
def psp_table(mod, tmp_path, monkeypatch):
    write(tmp_path / "phospho", SITE_HEADER, [
        site_row("AKT1", "P31749", "T308-p", ltp="5", htp="20", cst=""),
        site_row("AKT1", "P31749", "S473-p", ltp="12", htp="40", cst="3"),
        site_row("OLDNAME", "Q00001", "S5-p"),
    ])
    write(tmp_path / "ubiq", SITE_HEADER, [site_row("AKT1", "P31749", "K8-ub")])
    write(tmp_path / "diseases", DISEASE_HEADER, [
        disease_row("breast cancer", "increased", "P31749", "S473-p"),
    ])
    monkeypatch.setattr(mod, "input_dir", lambda root, sub: tmp_path)
    # Q00001 has no UniProt gene -- PSP's own GENE column should fill in
    monkeypatch.setattr(mod, "fetch_uniprot_gene_mapping",
                        lambda ids: pd.DataFrame({"UniProt": ["P31749"], "gene": ["AKT1"]}))
    table, n_genes = mod._build_psp_ptm_table()
    return table.set_index("uniprot_id"), n_genes


def test_columns_match_the_ptmd_table_plus_scores(psp_table):
    table, _ = psp_table
    assert list(table.reset_index().columns) == [
        "uniprot_id", "gene", "ptms_on_protein", "ptm_disease_pairs", "psp_site_scores",
    ]


def test_sites_are_in_position_order_with_types(psp_table):
    table, _ = psp_table
    assert table.loc["P31749", "ptms_on_protein"] == (
        "K8:Ubiquitination; T308:Phosphorylation; S473:Phosphorylation"
    ), "sites from every PTM file should be combined per protein, in sequence order"


def test_scores_and_diseases_are_encoded_per_site(psp_table):
    table, _ = psp_table
    assert "S473:Phosphorylation=12/40/3" in table.loc["P31749", "psp_site_scores"]
    assert "T308:Phosphorylation=5/20/0" in table.loc["P31749", "psp_site_scores"], (
        "a blank MS_CST must be written as 0"
    )
    assert table.loc["P31749", "ptm_disease_pairs"] == "S473:Phosphorylation | breast cancer (increased)", (
        "disease pairs must use PTMD's 'site | disease' format so step 3 reads both sources alike"
    )


def test_psp_gene_is_the_fallback(psp_table):
    table, n_genes = psp_table
    assert table.loc["Q00001", "gene"] == "OLDNAME"
    assert n_genes == 2


def test_filters_narrow_the_step1_table(mod, tmp_path, monkeypatch):
    write(tmp_path / "phospho", SITE_HEADER, [
        site_row("AKT1", "P31749", "S473-p", ltp="12", htp="40"),
        site_row("AKT1", "P31749", "T308-p", ltp="0", htp="3"),
    ])
    monkeypatch.setattr(mod, "input_dir", lambda root, sub: tmp_path)
    monkeypatch.setattr(mod, "fetch_uniprot_gene_mapping",
                        lambda ids: pd.DataFrame({"UniProt": ["P31749"], "gene": ["AKT1"]}))
    table, _ = mod._build_psp_ptm_table({"min_ltp": 1})
    assert list(table["ptms_on_protein"]) == ["S473:Phosphorylation"]

    from psp_input import PspInputError
    with pytest.raises(PspInputError, match="No PhosphoSitePlus sites pass"):
        mod._build_psp_ptm_table({"min_ltp": 100})
