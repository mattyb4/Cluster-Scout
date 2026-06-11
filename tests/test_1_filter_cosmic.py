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
    df = pd.DataFrame(
        rows,
        columns=["GENE_SYMBOL", "MUTATION_AA", "COSMIC_SAMPLE_ID", "MUTATION_SOMATIC_STATUS", "TRANSCRIPT_ACCESSION"],
    )
    path = tmp_path / "cosmic.tsv"
    df.to_csv(path, sep="\t", index=False)
    return path


def test_hotspot_threshold_filtering(mod, cosmic_file):
    cosmic, _, _ = mod._load_and_filter_cosmic(cosmic_file)

    # Only TP53/R175H has >= HOTSPOT_MIN_AFFECTED_CASES (3) confirmed-somatic samples
    assert list(cosmic["gene"]) == ["TP53"]
    assert list(cosmic["mutation"]) == ["R175H"]
    assert list(cosmic["affected_cases"]) == [3]
    assert list(cosmic["mutation_with_count"]) == ["R175H (3)"]


def test_gene_to_transcript_mapping(mod, cosmic_file):
    _, gene_to_transcript, _ = mod._load_and_filter_cosmic(cosmic_file)

    assert gene_to_transcript["TP53"] == "ENST00000269305.9"
    assert gene_to_transcript["PTPN11"] == "ENST00000351677.7"


def test_total_missense_patients_ignores_hotspot_threshold(mod, cosmic_file):
    _, _, gene_to_total = mod._load_and_filter_cosmic(cosmic_file)

    # TP53: R175H (S1,S2,S3) + R273H (S4,S5) = 5 distinct patients with a
    # qualifying missense mutation, even though R273H itself didn't meet the
    # hotspot threshold and the non-somatic/stop-codon rows were excluded.
    assert gene_to_total["TP53"] == 5
    # PTPN11: E76A (S1,S2) = 2 distinct patients, despite not being a hotspot.
    assert gene_to_total["PTPN11"] == 2
