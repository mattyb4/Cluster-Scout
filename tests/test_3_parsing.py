"""Unit tests for the PTM-row parsing helpers in scripts/3_find_nearby_mutations.py.

Most of these functions read from the module-level _PTM_ROWS cache (normally
populated from data/steps/PTMD_COSMIC_hotspots_by_protein.tsv via
get_ptm_rows()). Tests inject synthetic rows directly into _PTM_ROWS to avoid
touching that file -- get_ptm_rows() itself (the lazy-load/file-read/cache
behavior) is tested separately in TestGetPtmRows below, against a real
temp-file path.
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


class TestGetPtmRows:
    """Exercises the actual lazy-load-from-file/cache behavior directly,
    unlike every other test in this file which bypasses it via the ptm_rows
    fixture's direct _PTM_ROWS injection.
    """

    @pytest.fixture(autouse=True)
    def _reset_cache(self, mod):
        # _PTM_ROWS is a module-level global shared across the whole (session-
        # scoped) nearby_module fixture -- must be reset before AND after each
        # test in this class so it can neither leak in stale state from a
        # previous test nor leak out and corrupt tests in other files.
        mod._PTM_ROWS = None
        yield
        mod._PTM_ROWS = None

    def test_loads_rows_from_the_tsv_file(self, mod, tmp_path, monkeypatch):
        path = tmp_path / "hotspots.tsv"
        path.write_text("uniprot_id\tgene\nP04637\tTP53\n", encoding="utf-8")
        monkeypatch.setattr(mod, "PTM_TSV_PATH", path)

        rows = mod.get_ptm_rows()
        assert rows == [{"uniprot_id": "P04637", "gene": "TP53"}], (
            f"get_ptm_rows should parse the TSV into a list of dicts keyed by header, got {rows}"
        )

    def test_second_call_is_served_from_the_in_memory_cache(self, mod, tmp_path, monkeypatch):
        path = tmp_path / "hotspots.tsv"
        path.write_text("uniprot_id\tgene\nP04637\tTP53\n", encoding="utf-8")
        monkeypatch.setattr(mod, "PTM_TSV_PATH", path)

        first = mod.get_ptm_rows()
        # Change the file on disk -- if get_ptm_rows re-read it, this test would
        # see the new content instead of the cached one.
        path.write_text("uniprot_id\tgene\nQ06124\tPTPN11\n", encoding="utf-8")
        second = mod.get_ptm_rows()

        assert second is first, (
            "a second call within the same process must return the cached list "
            "object (no re-read of the file), matching the module's own "
            "'Lazy-load and cache' docstring -- got a different object"
        )
        assert second == [{"uniprot_id": "P04637", "gene": "TP53"}], (
            "the cached (stale) content should still reflect what was on disk at "
            "the time of the FIRST call, proving the file wasn't re-read"
        )


class TestParsePtmEntries(object):
    def test_parses_and_sorts_by_position(self, mod, ptm_rows):
        result = mod.parse_ptm_entries("Q06124")
        assert result == [
            ("T59", 59, "Phosphorylation"),
            ("Y62", 62, "Phosphorylation"),
        ], f"PTM entries should be parsed from 'ptms_on_protein' and sorted by position, got {result}"

    def test_unknown_uniprot_returns_empty(self, mod, ptm_rows):
        result = mod.parse_ptm_entries("Z99999")
        assert result == [], f"a UniProt ID with no matching row should yield no PTM entries, got {result}"


class TestParseGeneName:
    def test_known_uniprot(self, mod, ptm_rows):
        result = mod.parse_gene_name("Q06124")
        assert result == "PTPN11", f"should return the gene symbol for a known accession, got {result!r}"

    def test_unknown_uniprot(self, mod, ptm_rows):
        result = mod.parse_gene_name("Z99999")
        assert result == "", f"an unknown accession should default to an empty string, not raise, got {result!r}"


class TestParseIsoformSafeLength:
    def test_returns_int_when_present(self, mod, ptm_rows):
        result = mod.parse_isoform_safe_length("Q06124")
        assert result == 408, f"a present numeric isoform_safe_length should parse to an int, got {result!r}"

    def test_returns_none_when_blank(self, mod, ptm_rows):
        result = mod.parse_isoform_safe_length("P04637")
        assert result is None, (
            f"a blank isoform_safe_length means 'no restriction' and should return None, got {result!r}"
        )

    def test_returns_none_when_unknown_uniprot(self, mod, ptm_rows):
        result = mod.parse_isoform_safe_length("Z99999")
        assert result is None, f"an unknown accession should return None, not raise, got {result!r}"


class TestParseMutationPositions:
    def test_parses_and_sorts_by_position_then_mutation(self, mod, ptm_rows):
        result = mod.parse_mutation_positions("Q06124")
        assert result == [
            ("T59A", 59),
            ("D61Y", 61),
            ("E76A", 76),
        ], f"mutations should be parsed and sorted by (position, mutation label), got {result}"

    def test_no_uniprot_filter_combines_all_proteins(self, mod, ptm_rows):
        positions = mod.parse_mutation_positions()
        assert ("R175H", 175) in positions, (
            f"with uniprot=None, mutations from EVERY protein row should be included, but "
            f"TP53's R175H is missing from {positions}"
        )
        assert ("T59A", 59) in positions, (
            f"with uniprot=None, PTPN11's T59A should also be included, but is missing from {positions}"
        )

    def test_malformed_mutation_tokens_are_skipped_not_raised(self, mod):
        mod._PTM_ROWS = [{
            "uniprot_id": "X00000", "gene": "JUNK",
            "mutations_on_protein": "not a real mutation; also-not-one; ;",
        }]
        try:
            result = mod.parse_mutation_positions("X00000")
            assert result == [], (
                f"tokens that don't match the mutation regex (MUT_RE) must be silently "
                f"skipped, not raise or produce garbage entries, got {result}"
            )
        finally:
            mod._PTM_ROWS = None


class TestParseMutationPatientCounts:
    def test_maps_mutation_position_to_patient_count(self, mod, ptm_rows):
        result = mod.parse_mutation_patient_counts("Q06124")
        assert result == {
            ("T59A", 59): 281,
            ("E76A", 76): 50,
            ("D61Y", 61): 3,
        }, f"each mutation's '(N)' suffix should map to an integer patient count, got {result}"

    def test_unknown_uniprot_returns_empty_dict(self, mod, ptm_rows):
        result = mod.parse_mutation_patient_counts("Z99999")
        assert result == {}, f"an unknown accession has no counts to report, got {result}"


class TestParseTotalCosmicMissensePatients:
    def test_returns_int_when_present(self, mod, ptm_rows):
        result = mod.parse_total_cosmic_missense_patients("Q06124")
        assert result == 1690, f"a present value should parse to an int, got {result!r}"

    def test_returns_none_when_blank(self, mod, ptm_rows):
        result = mod.parse_total_cosmic_missense_patients("P04637")
        assert result is None, f"a blank value should return None rather than 0 or '', got {result!r}"

    def test_returns_none_when_unknown_uniprot(self, mod, ptm_rows):
        result = mod.parse_total_cosmic_missense_patients("Z99999")
        assert result is None, f"an accession with no row at all should also return None, got {result!r}"


class TestParsePtmDiseases:
    def test_filters_to_cancer_associated_diseases(self, mod, ptm_rows):
        # "Noonan syndrome" isn't cancer-associated and should be excluded;
        # "Acute myeloid leukemia" should be kept.
        result = mod.parse_ptm_diseases("Q06124", "T59", "Phosphorylation")
        assert result == "Acute myeloid leukemia", (
            f"only cancer-keyword-matching diseases should survive the filter ('Noonan "
            f"syndrome' must be excluded, 'Acute myeloid leukemia' kept), got {result!r}"
        )

    def test_no_match_returns_empty_string(self, mod, ptm_rows):
        result = mod.parse_ptm_diseases("Q06124", "Y62", "Phosphorylation")
        assert result == "", (
            f"a PTM site with no disease-pair entries at all should return an empty "
            f"string, not raise or return None, got {result!r}"
        )

    def test_entry_missing_delimiter_is_skipped(self, mod):
        # An entry without ' | ' can't be split into site/disease and must be
        # skipped via the `continue` branch, not raise a ValueError from unpacking.
        mod._PTM_ROWS = [{
            "uniprot_id": "X00000", "gene": "JUNK",
            "ptm_disease_pairs": "this entry has no pipe delimiter at all",
        }]
        try:
            result = mod.parse_ptm_diseases("X00000", "T1", "Phosphorylation")
            assert result == "", (
                f"a malformed disease-pair entry missing the ' | ' delimiter must be "
                f"skipped, not raise, got {result!r}"
            )
        finally:
            mod._PTM_ROWS = None

    def test_duplicate_diseases_are_not_repeated(self, mod):
        mod._PTM_ROWS = [{
            "uniprot_id": "X00000", "gene": "JUNK",
            "ptm_disease_pairs": (
                "T1:Phosphorylation | Breast cancer; T1:Phosphorylation | Breast cancer"
            ),
        }]
        try:
            result = mod.parse_ptm_diseases("X00000", "T1", "Phosphorylation")
            assert result == "Breast cancer", (
                f"the same disease listed twice for one site should be deduplicated, "
                f"not joined as 'Breast cancer; Breast cancer', got {result!r}"
            )
        finally:
            mod._PTM_ROWS = None


class TestParsePtmKnownDisruptions:
    def test_parses_site_to_mutation_set(self, mod, ptm_rows):
        result = mod.parse_ptm_known_disruptions("Q06124")
        assert result == {
            "T59:Phosphorylation": {"T59A", "D61Y"},
        }, f"should parse 'site:type>mut1,mut2' into a dict of sets, got {result}"

    def test_blank_field_returns_empty_dict(self, mod, ptm_rows):
        result = mod.parse_ptm_known_disruptions("P04637")
        assert result == {}, f"a blank ptm_known_disruptions field should give an empty dict, got {result}"

    def test_literal_nan_string_is_treated_as_blank(self, mod):
        # The field is documented to sometimes literally contain the string "nan"
        # (e.g. from a pandas NaN written out via str()), not just be empty.
        mod._PTM_ROWS = [{"uniprot_id": "X00000", "gene": "JUNK", "ptm_known_disruptions": "nan"}]
        try:
            result = mod.parse_ptm_known_disruptions("X00000")
            assert result == {}, (
                f"the literal string 'nan' must be treated the same as a truly blank "
                f"field, not parsed as a real entry, got {result}"
            )
        finally:
            mod._PTM_ROWS = None

    def test_entry_missing_greater_than_is_skipped(self, mod):
        mod._PTM_ROWS = [{
            "uniprot_id": "X00000", "gene": "JUNK",
            "ptm_known_disruptions": "T1:Phosphorylation no separator here",
        }]
        try:
            result = mod.parse_ptm_known_disruptions("X00000")
            assert result == {}, (
                f"an entry missing the '>' separator can't be split into site/mutations "
                f"and must be skipped, not raise, got {result}"
            )
        finally:
            mod._PTM_ROWS = None


class TestFormatMutations:
    def test_empty_hits(self, mod):
        result = mod.format_mutations([])
        assert result == "", f"no hits means nothing to format -- expected '', got {result!r}"

    def test_formats_and_sorts_by_position_then_mutation(self, mod):
        hits = [
            {"mutation": "C20D", "mutation_pos": 20, "distance": 3.456, "pae": 1.23},
            {"mutation": "A15B", "mutation_pos": 15, "distance": 5.0, "pae": None},
        ]
        result = mod.format_mutations(hits)
        assert result == "A15B-5.00Å, C20D-3.46Å(PAE:1.2)", (
            f"hits should sort by position, format distance to 2 decimals, and only "
            f"append a (PAE:...) suffix when pae is not None, got {result!r}"
        )


class TestLinearDistances:
    def test_empty_hits(self, mod):
        result = mod.linear_distances([], 10)
        assert result == "", f"no hits means no distances to report -- expected '', got {result!r}"

    def test_distances_sorted_and_deduped_by_position(self, mod):
        hits = [
            {"mutation": "A20B", "mutation_pos": 20, "distance": 1, "pae": None},
            {"mutation": "C15D", "mutation_pos": 15, "distance": 1, "pae": None},
            {"mutation": "E15F", "mutation_pos": 15, "distance": 1, "pae": None},  # same position, dedup
        ]
        result = mod.linear_distances(hits, ptm_pos=10)
        assert result == "5,10", (
            f"two hits sharing position 15 should be deduplicated to one distance (5), "
            f"sorted ascending with position 20's distance (10), got {result!r}"
        )


class TestUniqueMutationPositionCount:
    def test_counts_distinct_positions(self, mod):
        hits = [
            {"mutation_pos": 15},
            {"mutation_pos": 15},
            {"mutation_pos": 20},
        ]
        result = mod.unique_mutation_position_count(hits)
        assert result == 2, f"position 15 appears twice but should count once -- expected 2 distinct positions, got {result}"

    def test_empty_hits_returns_zero(self, mod):
        result = mod.unique_mutation_position_count([])
        assert result == 0, f"no hits means zero distinct positions, got {result}"


class TestMutationAtPtmSite:
    def test_yes_when_position_matches(self, mod):
        hits = [{"mutation_pos": 59}]
        result = mod.mutation_at_ptm_site(hits, 59)
        assert result == "yes", f"a hit exactly at the PTM position should report 'yes', got {result!r}"

    def test_no_when_position_does_not_match(self, mod):
        hits = [{"mutation_pos": 59}]
        result = mod.mutation_at_ptm_site(hits, 100)
        assert result == "no", f"no hit at the PTM position should report 'no', got {result!r}"


class TestTotalPatientCount:
    def test_sums_counts_across_hits(self, mod):
        hits = [
            {"mutation": "T59A", "mutation_pos": 59},
            {"mutation": "D61Y", "mutation_pos": 61},
        ]
        patient_counts = {("T59A", 59): 281, ("D61Y", 61): 3}
        result = mod.total_patient_count(hits, patient_counts)
        assert result == 284, f"total should sum every hit's individual count (281+3), got {result}"

    def test_strips_isoform_tag_before_lookup(self, mod):
        hits = [{"mutation": "D61Y(isoform?)", "mutation_pos": 61}]
        patient_counts = {("D61Y", 61): 3}
        result = mod.total_patient_count(hits, patient_counts)
        assert result == 3, (
            f"the '(isoform?)' suffix must be stripped before looking up the count "
            f"dict (which is keyed by the plain mutation label), got {result}"
        )

    def test_unmatched_hit_contributes_zero(self, mod):
        hits = [{"mutation": "Z99Z", "mutation_pos": 1}]
        result = mod.total_patient_count(hits, {})
        assert result == 0, (
            f"a hit with no entry in patient_counts should default to contributing 0, "
            f"not raise a KeyError, got {result}"
        )
