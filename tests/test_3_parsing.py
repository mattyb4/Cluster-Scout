"""Unit tests for the PTM-row parsing helpers in scripts/3_find_nearby_mutations.py.

These functions read from the module-level _PTM_ROWS cache (normally populated
from data/steps/PTMD_COSMIC_hotspots_by_protein.tsv via get_ptm_rows()). Tests
inject synthetic rows directly into _PTM_ROWS to avoid touching that file.
"""
import pytest


@pytest.fixture
def mod(nearby_module):
    return nearby_module


PTPN11_ROW = {
    "uniprot_id": "Q06124",
    "gene": "PTPN11",
    "ptms_on_protein": "T59:Phosphorylation; Y62:Phosphorylation",
    "mutations_on_protein": "T59A (281); E76A (50); D61Y (3)",
    "ptm_disease_pairs": (
        "T59:Phosphorylation | Noonan syndrome; "
        "T59:Phosphorylation | Acute myeloid leukemia"
    ),
    "ptm_known_disruptions": "T59:Phosphorylation>T59A,D61Y",
    "isoform_safe_length": "408",
    "total_cosmic_missense_patients": "1690",
}

NO_RESTRICTION_ROW = {
    "uniprot_id": "P04637",
    "gene": "TP53",
    "ptms_on_protein": "S15:Phosphorylation",
    "mutations_on_protein": "R175H (10)",
    "ptm_disease_pairs": "S15:Phosphorylation | Breast cancer",
    "ptm_known_disruptions": "",
    "isoform_safe_length": "",
    "total_cosmic_missense_patients": "",
}


@pytest.fixture
def ptm_rows(mod):
    mod._PTM_ROWS = [PTPN11_ROW, NO_RESTRICTION_ROW]
    yield
    mod._PTM_ROWS = None


class TestParsePtmEntries(object):
    def test_parses_and_sorts_by_position(self, mod, ptm_rows):
        assert mod.parse_ptm_entries("Q06124") == [
            ("T59", 59, "Phosphorylation"),
            ("Y62", 62, "Phosphorylation"),
        ]

    def test_unknown_uniprot_returns_empty(self, mod, ptm_rows):
        assert mod.parse_ptm_entries("Z99999") == []


class TestParseGeneName:
    def test_known_uniprot(self, mod, ptm_rows):
        assert mod.parse_gene_name("Q06124") == "PTPN11"

    def test_unknown_uniprot(self, mod, ptm_rows):
        assert mod.parse_gene_name("Z99999") == ""


class TestParseIsoformSafeLength:
    def test_returns_int_when_present(self, mod, ptm_rows):
        assert mod.parse_isoform_safe_length("Q06124") == 408

    def test_returns_none_when_blank(self, mod, ptm_rows):
        assert mod.parse_isoform_safe_length("P04637") is None

    def test_returns_none_when_unknown_uniprot(self, mod, ptm_rows):
        assert mod.parse_isoform_safe_length("Z99999") is None


class TestParseMutationPositions:
    def test_parses_and_sorts_by_position_then_mutation(self, mod, ptm_rows):
        assert mod.parse_mutation_positions("Q06124") == [
            ("T59A", 59),
            ("D61Y", 61),
            ("E76A", 76),
        ]

    def test_no_uniprot_filter_combines_all_proteins(self, mod, ptm_rows):
        positions = mod.parse_mutation_positions()
        assert ("R175H", 175) in positions
        assert ("T59A", 59) in positions


class TestParseMutationPatientCounts:
    def test_maps_mutation_position_to_patient_count(self, mod, ptm_rows):
        assert mod.parse_mutation_patient_counts("Q06124") == {
            ("T59A", 59): 281,
            ("E76A", 76): 50,
            ("D61Y", 61): 3,
        }

    def test_unknown_uniprot_returns_empty_dict(self, mod, ptm_rows):
        assert mod.parse_mutation_patient_counts("Z99999") == {}


class TestParseTotalCosmicMissensePatients:
    def test_returns_int_when_present(self, mod, ptm_rows):
        assert mod.parse_total_cosmic_missense_patients("Q06124") == 1690

    def test_returns_none_when_blank(self, mod, ptm_rows):
        assert mod.parse_total_cosmic_missense_patients("P04637") is None


class TestParsePtmDiseases:
    def test_filters_to_cancer_associated_diseases(self, mod, ptm_rows):
        # "Noonan syndrome" isn't cancer-associated and should be excluded;
        # "Acute myeloid leukemia" should be kept.
        assert mod.parse_ptm_diseases("Q06124", "T59", "Phosphorylation") == "Acute myeloid leukemia"

    def test_no_match_returns_empty_string(self, mod, ptm_rows):
        assert mod.parse_ptm_diseases("Q06124", "Y62", "Phosphorylation") == ""


class TestParsePtmKnownDisruptions:
    def test_parses_site_to_mutation_set(self, mod, ptm_rows):
        assert mod.parse_ptm_known_disruptions("Q06124") == {
            "T59:Phosphorylation": {"T59A", "D61Y"},
        }

    def test_blank_field_returns_empty_dict(self, mod, ptm_rows):
        assert mod.parse_ptm_known_disruptions("P04637") == {}


class TestFormatMutations:
    def test_empty_hits(self, mod):
        assert mod.format_mutations([]) == ""

    def test_formats_and_sorts_by_position_then_mutation(self, mod):
        hits = [
            {"mutation": "C20D", "mutation_pos": 20, "distance": 3.456, "pae": 1.23},
            {"mutation": "A15B", "mutation_pos": 15, "distance": 5.0, "pae": None},
        ]
        assert mod.format_mutations(hits) == "A15B-5.00Å, C20D-3.46Å(PAE:1.2)"


class TestLinearDistances:
    def test_empty_hits(self, mod):
        assert mod.linear_distances([], 10) == ""

    def test_distances_sorted_and_deduped_by_position(self, mod):
        hits = [
            {"mutation": "A20B", "mutation_pos": 20, "distance": 1, "pae": None},
            {"mutation": "C15D", "mutation_pos": 15, "distance": 1, "pae": None},
            {"mutation": "E15F", "mutation_pos": 15, "distance": 1, "pae": None},  # same position, dedup
        ]
        assert mod.linear_distances(hits, ptm_pos=10) == "5,10"


class TestUniqueMutationPositionCount:
    def test_counts_distinct_positions(self, mod):
        hits = [
            {"mutation_pos": 15},
            {"mutation_pos": 15},
            {"mutation_pos": 20},
        ]
        assert mod.unique_mutation_position_count(hits) == 2


class TestMutationAtPtmSite:
    def test_yes_when_position_matches(self, mod):
        hits = [{"mutation_pos": 59}]
        assert mod.mutation_at_ptm_site(hits, 59) == "yes"

    def test_no_when_position_does_not_match(self, mod):
        hits = [{"mutation_pos": 59}]
        assert mod.mutation_at_ptm_site(hits, 100) == "no"


class TestTotalPatientCount:
    def test_sums_counts_across_hits(self, mod):
        hits = [
            {"mutation": "T59A", "mutation_pos": 59},
            {"mutation": "D61Y", "mutation_pos": 61},
        ]
        patient_counts = {("T59A", 59): 281, ("D61Y", 61): 3}
        assert mod.total_patient_count(hits, patient_counts) == 284

    def test_strips_isoform_tag_before_lookup(self, mod):
        hits = [{"mutation": "D61Y(isoform?)", "mutation_pos": 61}]
        patient_counts = {("D61Y", 61): 3}
        assert mod.total_patient_count(hits, patient_counts) == 3

    def test_unmatched_hit_contributes_zero(self, mod):
        hits = [{"mutation": "Z99Z", "mutation_pos": 1}]
        assert mod.total_patient_count(hits, {}) == 0
