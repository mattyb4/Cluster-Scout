"""Unit tests for _load_and_filter_cosmic in scripts/1_filter.py.

Verifies the COSMIC Mutant Census filtering rules:
 - only "confirmed/reported somatic" rows count
 - only simple amino-acid substitutions (e.g. R175H) are kept
 - hotspots require >= HOTSPOT_MIN_AFFECTED_CASES distinct samples
 - gene_to_total_missense_patients counts ALL missense patients per gene,
   regardless of the hotspot threshold
"""
import pandas as pd
import pytest


@pytest.fixture
def mod(filter_module):
    return filter_module


def _make_cosmic_file(tmp_path, rows, name="cosmic.tsv"):
    df = pd.DataFrame(
        rows,
        columns=["GENE_SYMBOL", "MUTATION_AA", "COSMIC_SAMPLE_ID", "MUTATION_SOMATIC_STATUS", "TRANSCRIPT_ACCESSION"],
    )
    path = tmp_path / name
    df.to_csv(path, sep="\t", index=False)
    return path


@pytest.fixture
def cosmic_file(tmp_path):
    rows = [
        # TP53/R175H: 3 confirmed-somatic samples -> meets hotspot threshold
        ("TP53", "p.R175H", "S1", "Confirmed somatic variant", "ENST00000269305.9"),
        ("TP53", "p.R175H", "S2", "Confirmed somatic variant", "ENST00000269305.9"),
        ("TP53", "p.R175H", "S3", "Confirmed somatic variant", "ENST00000269305.9"),
        # TP53/R273H: only 2 confirmed-somatic samples -> below hotspot threshold
        ("TP53", "p.R273H", "S4", "Confirmed somatic variant", "ENST00000269305.9"),
        ("TP53", "p.R273H", "S5", "Reported in another cancer sample as somatic", "ENST00000269305.9"),
        # TP53/E11*: stop-codon, not a simple substitution -> dropped entirely
        ("TP53", "p.E11*", "S6", "Confirmed somatic variant", "ENST00000269305.9"),
        # TP53/R175H reported again but with a non-somatic status -> dropped
        ("TP53", "p.R175H", "S7", "Variant of unknown origin", "ENST00000269305.9"),
        # PTPN11/E76A: 2 confirmed-somatic samples -> below hotspot threshold
        ("PTPN11", "p.E76A", "S1", "Confirmed somatic variant", "ENST00000351677.7"),
        ("PTPN11", "p.E76A", "S2", "Confirmed somatic variant", "ENST00000351677.7"),
    ]
    return _make_cosmic_file(tmp_path, rows)


def test_hotspot_threshold_filtering(mod, cosmic_file):
    cosmic, _, _ = mod._load_and_filter_cosmic(cosmic_file)

    # Only TP53/R175H has >= HOTSPOT_MIN_AFFECTED_CASES (3) confirmed-somatic samples
    assert list(cosmic["gene"]) == ["TP53"], (
        f"only TP53/R175H meets the hotspot threshold -- PTPN11/E76A (2 samples) and "
        f"TP53/R273H (2 samples) should both be filtered out, got genes {list(cosmic['gene'])}"
    )
    assert list(cosmic["mutation"]) == ["R175H"], (
        f"the surviving row's mutation should be the 'p.' prefix stripped, got {list(cosmic['mutation'])}"
    )
    assert list(cosmic["affected_cases"]) == [3], (
        f"affected_cases should count distinct confirmed-somatic samples (S1,S2,S3), "
        f"not the raw row count (which would also include the non-somatic S7 row), "
        f"got {list(cosmic['affected_cases'])}"
    )
    assert list(cosmic["mutation_with_count"]) == ["R175H (3)"], (
        f"mutation_with_count should combine the mutation label and its case count, "
        f"got {list(cosmic['mutation_with_count'])}"
    )


def test_gene_to_transcript_mapping(mod, cosmic_file):
    _, gene_to_transcript, _ = mod._load_and_filter_cosmic(cosmic_file)

    assert gene_to_transcript["TP53"] == "ENST00000269305.9", (
        "gene_to_transcript should map each gene to its COSMIC TRANSCRIPT_ACCESSION"
    )
    assert gene_to_transcript["PTPN11"] == "ENST00000351677.7", (
        "the transcript mapping should be built for every gene present, not just the "
        "one that met the hotspot threshold"
    )


def test_total_missense_patients_ignores_hotspot_threshold(mod, cosmic_file):
    _, _, gene_to_total = mod._load_and_filter_cosmic(cosmic_file)

    # TP53: R175H (S1,S2,S3) + R273H (S4,S5) = 5 distinct patients with a
    # qualifying missense mutation, even though R273H itself didn't meet the
    # hotspot threshold and the non-somatic/stop-codon rows were excluded.
    assert gene_to_total["TP53"] == 5, (
        f"total missense patient count must include R273H's 2 samples even though "
        f"R273H didn't individually meet the hotspot threshold (3+2=5), got {gene_to_total['TP53']}"
    )
    # PTPN11: E76A (S1,S2) = 2 distinct patients, despite not being a hotspot.
    assert gene_to_total["PTPN11"] == 2, (
        f"PTPN11 should still be counted (2 patients) even though its only mutation "
        f"never reached hotspot status, got {gene_to_total['PTPN11']}"
    )


def test_no_rows_survive_filtering_returns_empty_not_error(mod, tmp_path):
    # Every row here fails at least one filter: non-somatic status, or a
    # stop-codon change that isn't a simple substitution.
    rows = [
        ("TP53", "p.R175H", "S1", "Variant of unknown origin", "ENST00000269305.9"),
        ("TP53", "p.E11*", "S2", "Confirmed somatic variant", "ENST00000269305.9"),
    ]
    cosmic_file = _make_cosmic_file(tmp_path, rows)

    cosmic, gene_to_transcript, gene_to_total = mod._load_and_filter_cosmic(cosmic_file)

    assert cosmic.empty, (
        f"if every row is filtered out (non-somatic status, stop-codon changes), the "
        f"result must be an empty DataFrame, not raise -- got {len(cosmic)} row(s)"
    )
    assert gene_to_total == {}, (
        "with no qualifying somatic missense rows for any gene, gene_to_total should "
        f"be an empty dict, got {gene_to_total}"
    )
