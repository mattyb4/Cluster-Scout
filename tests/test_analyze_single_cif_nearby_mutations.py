"""Unit tests for scripts/analyze_single_cif_nearby_mutations.py (the Single
Protein mode's analysis engine) -- gap-filling everything except
append_rows_no_duplicates, which is already covered by
test_single_protein_append.py.
"""
import csv

import numpy as np
import pytest

from conftest import import_script

mod = import_script("analyze_single_cif_nearby_mutations.py")


def _make_residue(resname, pos, coord=(0.0, 0.0, 0.0), bfactor=80.0, hetero=False, has_ca=True):
    """Build a real Bio.PDB Residue (with an optional CA atom) -- mirrors the
    direct-Atom-construction technique already used in test_cif_variance.py,
    since this script's chain-walking code (`for residue in chain`,
    `"CA" in residue`) works the same whether the chain is a real Bio.PDB
    Chain or just a plain list of Residue objects.
    """
    from Bio.PDB.Residue import Residue
    from Bio.PDB.Atom import Atom

    hetfield = "H" if hetero else " "
    res = Residue((hetfield, pos, " "), resname, " ")
    if has_ca:
        atom = Atom("CA", np.array(coord, dtype=float), bfactor, 1.0, " ", "CA", pos, element="C")
        res.add(atom)
    return res


def _write_tsv(path, fieldnames, rows):
    with path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


class TestGetCaCoord:
    def test_returns_coord(self):
        chain = [_make_residue("SER", 10, coord=(1.0, 2.0, 3.0))]
        result = mod.get_ca_coord(chain, 10)
        assert np.array_equal(result, [1.0, 2.0, 3.0]), f"should return residue 10's CA coordinate, got {result}"

    def test_returns_none_for_missing_residue(self):
        chain = [_make_residue("SER", 10)]
        assert mod.get_ca_coord(chain, 999) is None

    def test_returns_none_when_residue_has_no_ca(self):
        chain = [_make_residue("SER", 10, has_ca=False)]
        assert mod.get_ca_coord(chain, 10) is None, (
            "a residue present in the chain but missing its CA atom should return None, not raise"
        )


class TestComputeDistance:
    def test_basic_distance(self):
        result = mod.compute_distance(np.array([0.0, 0.0, 0.0]), np.array([3.0, 4.0, 0.0]))
        assert result == 5.0, f"a 3-4-5 right triangle should give distance 5.0, got {result}"


class TestFindNearbyMutations:
    def _chain(self):
        return [
            _make_residue("SER", 10, coord=(0.0, 0.0, 0.0)),
            _make_residue("ALA", 15, coord=(5.0, 0.0, 0.0)),  # 5A away
            _make_residue("GLY", 30, coord=(50.0, 0.0, 0.0)),  # 50A away
        ]

    def test_only_returns_hits_within_cutoff(self):
        mutation_map = {15: {"A15B"}, 30: {"G30H"}}
        results = mod.find_nearby_mutations(self._chain(), 10, mutation_map, cutoff=10.0)
        assert len(results) == 1, f"only the 5A-away mutation should be within a 10A cutoff, got {len(results)}"
        assert results[0]["mutation_label"] == "A15B"

    def test_returns_empty_if_ptm_position_missing(self):
        results = mod.find_nearby_mutations(self._chain(), 999, {15: {"A15B"}}, cutoff=10.0)
        assert results == [], "a PTM position with no CA atom in the structure should give no results, not raise"

    def test_picks_first_label_alphabetically_when_multiple(self):
        # A position with 2 substitution labels (e.g. two different COSMIC
        # entries at the same residue) -- the sorted-first one is used as the
        # representative label.
        mutation_map = {15: {"A15Z", "A15B"}}
        results = mod.find_nearby_mutations(self._chain(), 10, mutation_map, cutoff=10.0)
        assert results[0]["mutation_label"] == "A15B", (
            f"with multiple labels at one position, the alphabetically-first should be "
            f"chosen as the representative, got {results[0]['mutation_label']!r}"
        )

    def test_max_pae_excludes_hits_above_threshold(self):
        pae_matrix = np.zeros((30, 30))
        pae_matrix[9, 14] = 8.0
        pae_matrix[14, 9] = 8.0
        results = mod.find_nearby_mutations(
            self._chain(), 10, {15: {"A15B"}}, cutoff=10.0, pae_matrix=pae_matrix, max_pae=5.0,
        )
        assert results == [], "a hit whose averaged PAE (8.0) exceeds max_pae (5.0) must be excluded"


class TestGetUniprotValue:
    def test_prefers_uniprot_column(self):
        assert mod.get_uniprot_value({"UniProt": "P04637", "uniprot_id": "Q99999"}) == "P04637"

    def test_falls_back_to_uniprot_id_column(self):
        assert mod.get_uniprot_value({"uniprot_id": "P04637"}) == "P04637"

    def test_missing_both_returns_none(self):
        assert mod.get_uniprot_value({}) is None


class TestGetGeneValue:
    def test_prefers_gene_column(self):
        assert mod.get_gene_value({"gene": "TP53", "Gene": "OTHER"}) == "TP53"

    def test_falls_back_to_capitalized_gene_column(self):
        assert mod.get_gene_value({"Gene": "TP53"}) == "TP53"

    def test_missing_both_returns_empty_string(self):
        assert mod.get_gene_value({}) == "", "with no gene column at all, should default to '' rather than None"


class TestExtractPtmLabelsFromListField:
    def test_parses_semicolon_separated_tokens(self):
        result = mod.extract_ptm_labels_from_list_field("S15:Phosphorylation; T18:Phosphorylation")
        assert result == {15: {"S15"}, 18: {"T18"}}, f"expected positions mapped to label sets, got {result}"

    def test_empty_field_returns_empty_dict(self):
        assert mod.extract_ptm_labels_from_list_field("") == {}


class TestExtractMutationLabelsFromField:
    def test_parses_semicolon_separated_tokens(self):
        result = mod.extract_mutation_labels_from_field("R175H (10); R273H (5)")
        assert result == {175: {"R175H"}, 273: {"R273H"}}, f"expected positions mapped to mutation label sets, got {result}"

    def test_multiple_labels_at_same_position_are_grouped(self):
        result = mod.extract_mutation_labels_from_field("R175H (10); R175C (2)")
        assert result == {175: {"R175H", "R175C"}}, (
            f"two different substitutions at the SAME position should be grouped into "
            f"one set, not overwrite each other, got {result}"
        )


class TestParsePtmEntries:
    def test_legacy_schema(self, tmp_path):
        tsv = tmp_path / "data.tsv"
        _write_tsv(tsv, ["UniProt", "ptm_pos", "ptm_site", "ptm_type"], [
            {"UniProt": "P04637", "ptm_pos": "15", "ptm_site": "S15", "ptm_type": "Phosphorylation"},
        ])
        result = mod.parse_ptm_entries("P04637", tsv)
        assert result == [("S15", 15, "Phosphorylation")], f"legacy-schema row should parse directly, got {result}"

    def test_new_schema(self, tmp_path):
        tsv = tmp_path / "data.tsv"
        _write_tsv(tsv, ["uniprot_id", "ptms_on_protein"], [
            {"uniprot_id": "P04637", "ptms_on_protein": "S15:Phosphorylation; T18:Phosphorylation"},
        ])
        result = mod.parse_ptm_entries("P04637", tsv)
        assert result == [("S15", 15, "Phosphorylation"), ("T18", 18, "Phosphorylation")], (
            f"new-schema semicolon-separated tokens should be parsed and sorted by position, got {result}"
        )

    def test_filters_by_uniprot(self, tmp_path):
        tsv = tmp_path / "data.tsv"
        _write_tsv(tsv, ["uniprot_id", "ptms_on_protein"], [
            {"uniprot_id": "P04637", "ptms_on_protein": "S15:Phosphorylation"},
            {"uniprot_id": "Q99999", "ptms_on_protein": "S99:Phosphorylation"},
        ])
        result = mod.parse_ptm_entries("P04637", tsv)
        assert result == [("S15", 15, "Phosphorylation")], f"rows for other proteins must be excluded, got {result}"


class TestParseMutationPositions:
    def test_legacy_schema(self, tmp_path):
        tsv = tmp_path / "data.tsv"
        _write_tsv(tsv, ["UniProt", "ptm_pos", "near_mutation", "far_mutations_prevalence_filtered"], [
            {"UniProt": "P04637", "ptm_pos": "15", "near_mutation": "A16B", "far_mutations_prevalence_filtered": "C30D"},
        ])
        result = mod.parse_mutation_positions(15, tsv, uniprot="P04637")
        assert result == {16: {"A16B"}, 30: {"C30D"}}, (
            f"legacy schema should combine near_mutation and far_mutations columns, got {result}"
        )

    def test_legacy_schema_ignores_other_ptm_positions(self, tmp_path):
        tsv = tmp_path / "data.tsv"
        _write_tsv(tsv, ["UniProt", "ptm_pos", "near_mutation", "far_mutations_prevalence_filtered"], [
            {"UniProt": "P04637", "ptm_pos": "99", "near_mutation": "A16B", "far_mutations_prevalence_filtered": ""},
        ])
        result = mod.parse_mutation_positions(15, tsv, uniprot="P04637")
        assert result == {}, "a row for a DIFFERENT ptm_pos (99, not the requested 15) must be excluded"

    def test_new_schema_protein_level_not_ptm_specific(self, tmp_path):
        tsv = tmp_path / "data.tsv"
        _write_tsv(tsv, ["uniprot_id", "mutations_on_protein"], [
            {"uniprot_id": "P04637", "mutations_on_protein": "R175H (10)"},
        ])
        # New schema's mutation list isn't tied to a specific ptm_pos -- the same
        # result should come back regardless of which ptm_pos is requested.
        result_a = mod.parse_mutation_positions(15, tsv, uniprot="P04637")
        result_b = mod.parse_mutation_positions(999, tsv, uniprot="P04637")
        assert result_a == result_b == {175: {"R175H"}}, (
            f"new-schema mutations are protein-level, not PTM-position-specific -- "
            f"both lookups should return the same set, got {result_a} vs {result_b}"
        )


class TestParseIsoformSafeLength:
    def test_returns_int_when_present(self, tmp_path):
        tsv = tmp_path / "data.tsv"
        _write_tsv(tsv, ["uniprot_id", "isoform_safe_length"], [{"uniprot_id": "P04637", "isoform_safe_length": "393"}])
        assert mod.parse_isoform_safe_length("P04637", tsv) == 393

    def test_returns_none_when_blank(self, tmp_path):
        tsv = tmp_path / "data.tsv"
        _write_tsv(tsv, ["uniprot_id", "isoform_safe_length"], [{"uniprot_id": "P04637", "isoform_safe_length": ""}])
        assert mod.parse_isoform_safe_length("P04637", tsv) is None

    def test_returns_none_when_unknown_uniprot(self, tmp_path):
        tsv = tmp_path / "data.tsv"
        _write_tsv(tsv, ["uniprot_id", "isoform_safe_length"], [{"uniprot_id": "P04637", "isoform_safe_length": "393"}])
        assert mod.parse_isoform_safe_length("Z99999", tsv) is None


class TestBuildPosToAa:
    def test_maps_position_to_one_letter_code(self):
        chain = [_make_residue("SER", 1), _make_residue("ALA", 2)]
        result = mod.build_pos_to_aa(chain)
        assert result == {1: "S", 2: "A"}, f"expected one-letter codes per position, got {result}"

    def test_hetero_residues_are_skipped(self):
        chain = [_make_residue("SER", 1), _make_residue("HOH", 2, hetero=True)]
        result = mod.build_pos_to_aa(chain)
        assert result == {1: "S"}, f"a hetero residue (e.g. water) must be excluded, got {result}"

    def test_unknown_residue_name_maps_to_question_mark(self):
        chain = [_make_residue("XYZ", 1)]
        result = mod.build_pos_to_aa(chain)
        assert result == {1: "?"}, (
            f"an unrecognized residue name should map to '?' (not raise or 'X' -- this "
            f"file's own fallback differs from pipeline_utils' 'X'), got {result}"
        )


class TestBuildPosToPlddt:
    def test_maps_position_to_bfactor(self):
        chain = [_make_residue("SER", 1, bfactor=55.5)]
        result = mod.build_pos_to_plddt(chain)
        assert result == {1: 55.5}, f"pLDDT should come from the CA atom's B-factor, got {result}"

    def test_residue_without_ca_is_skipped(self):
        chain = [_make_residue("SER", 1, has_ca=False)]
        result = mod.build_pos_to_plddt(chain)
        assert result == {}, "a residue with no CA atom has no pLDDT to report and must be excluded"


class TestParseMutationPatientCounts:
    def test_maps_label_position_to_count(self, tmp_path):
        tsv = tmp_path / "data.tsv"
        _write_tsv(tsv, ["uniprot_id", "mutations_on_protein"], [
            {"uniprot_id": "P04637", "mutations_on_protein": "R175H (10); R273H (5)"},
        ])
        result = mod.parse_mutation_patient_counts("P04637", tsv)
        assert result == {("R175H", 175): 10, ("R273H", 273): 5}, f"expected (label,pos)->count mapping, got {result}"

    def test_legacy_schema_rows_produce_no_counts(self, tmp_path):
        tsv = tmp_path / "data.tsv"
        _write_tsv(tsv, ["uniprot_id", "ptm_pos", "near_mutation"], [
            {"uniprot_id": "P04637", "ptm_pos": "15", "near_mutation": "A16B"},
        ])
        result = mod.parse_mutation_patient_counts("P04637", tsv)
        assert result == {}, (
            "legacy-schema rows have no mutations_on_protein column at all, so no "
            "counts can be extracted -- must be empty, not raise a KeyError"
        )


class TestFilterMutationsByMinSamples:
    def test_drops_labels_below_threshold(self):
        mutation_map = {175: {"R175H"}, 273: {"R273H"}}
        patient_counts = {("R175H", 175): 10, ("R273H", 273): 2}
        result = mod.filter_mutations_by_min_samples(mutation_map, patient_counts, min_samples=3)
        assert result == {175: {"R175H"}}, (
            f"R273H has only 2 known samples (below the threshold of 3) and should be "
            f"dropped; position 273 should disappear entirely since nothing survives, got {result}"
        )

    def test_labels_with_unknown_count_are_kept(self):
        # No entry in patient_counts at all (e.g. legacy-schema data) -- can't tell
        # whether it would pass, so it must be kept rather than assumed to fail.
        mutation_map = {175: {"R175H"}}
        result = mod.filter_mutations_by_min_samples(mutation_map, {}, min_samples=3)
        assert result == {175: {"R175H"}}, (
            f"a label with no known patient count must be kept (can't recover data "
            f"already excluded upstream), got {result}"
        )


class TestTagIsoformMutations:
    def test_tags_mismatched_reference_aa(self):
        mutation_map = {175: {"R175H"}}
        pos_to_aa = {175: "K"}  # structure has Lysine, not the expected Arginine
        result = mod.tag_isoform_mutations(mutation_map, pos_to_aa, safe_length=None)
        assert result == {175: {"R175H(isoform?)"}}, (
            f"a mutation whose reference AA (R) doesn't match the structure's actual "
            f"residue (K) should be tagged, got {result}"
        )

    def test_matching_reference_aa_is_untagged(self):
        mutation_map = {175: {"R175H"}}
        pos_to_aa = {175: "R"}
        result = mod.tag_isoform_mutations(mutation_map, pos_to_aa, safe_length=None)
        assert result == {175: {"R175H"}}, f"a matching reference AA should be left untagged, got {result}"

    def test_position_beyond_safe_length_is_tagged(self):
        mutation_map = {500: {"A500B"}}
        pos_to_aa = {500: "A"}  # reference AA matches...
        result = mod.tag_isoform_mutations(mutation_map, pos_to_aa, safe_length=400)  # ...but beyond safe_length
        assert result == {500: {"A500B(isoform?)"}}, (
            f"even with a matching reference AA, a position past the isoform-safe "
            f"length boundary must still be tagged, got {result}"
        )


class TestParseGeneName:
    def test_known_uniprot(self, tmp_path):
        tsv = tmp_path / "data.tsv"
        _write_tsv(tsv, ["uniprot_id", "gene"], [{"uniprot_id": "P04637", "gene": "TP53"}])
        assert mod.parse_gene_name("P04637", tsv) == "TP53"

    def test_unknown_uniprot_returns_empty_string(self, tmp_path):
        tsv = tmp_path / "data.tsv"
        _write_tsv(tsv, ["uniprot_id", "gene"], [{"uniprot_id": "P04637", "gene": "TP53"}])
        assert mod.parse_gene_name("Z99999", tsv) == ""


class TestParsePtmDiseases:
    def test_matches_site_and_type(self, tmp_path):
        tsv = tmp_path / "data.tsv"
        _write_tsv(tsv, ["uniprot_id", "ptm_disease_pairs"], [
            {"uniprot_id": "P04637", "ptm_disease_pairs": "S15:Phosphorylation | Breast cancer"},
        ])
        result = mod.parse_ptm_diseases("P04637", "S15", "Phosphorylation", tsv)
        assert result == "Breast cancer", f"a matching site:type entry should return its disease, got {result!r}"

    def test_no_match_returns_empty_string(self, tmp_path):
        tsv = tmp_path / "data.tsv"
        _write_tsv(tsv, ["uniprot_id", "ptm_disease_pairs"], [
            {"uniprot_id": "P04637", "ptm_disease_pairs": "S15:Phosphorylation | Breast cancer"},
        ])
        result = mod.parse_ptm_diseases("P04637", "T99", "Phosphorylation", tsv)
        assert result == "", f"a PTM site with no matching disease entries should return '', got {result!r}"


class TestFormatMutations:
    def test_empty_hits(self):
        assert mod.format_mutations([]) == ""

    def test_formats_with_plain_a_suffix_not_angstrom_symbol(self):
        # Unlike 3_find_nearby_mutations.py's format_mutations (which uses the
        # literal "Å" symbol), this script's version uses a plain ASCII "A".
        hits = [{"mutation_pos": 15, "mutation_label": "A15B", "distance": 5.0, "pae": None}]
        result = mod.format_mutations(hits)
        assert result == "A15B-5.00A", f"expected a plain ASCII 'A' suffix (not 'Å'), got {result!r}"

    def test_includes_pae_when_present(self):
        hits = [{"mutation_pos": 15, "mutation_label": "A15B", "distance": 5.0, "pae": 2.34}]
        result = mod.format_mutations(hits)
        assert result == "A15B-5.00A(PAE:2.3)", f"a present PAE should append a (PAE:...) suffix, got {result!r}"


class TestLinearDistances:
    def test_empty_hits(self):
        assert mod.linear_distances([], 10) == ""

    def test_deduped_and_sorted(self):
        hits = [
            {"mutation_pos": 20, "mutation_label": "A20B", "distance": 1, "pae": None},
            {"mutation_pos": 15, "mutation_label": "C15D", "distance": 1, "pae": None},
        ]
        result = mod.linear_distances(hits, ptm_pos=10)
        assert result == "5,10", f"expected sorted-by-position linear distances, got {result!r}"


class TestUniqueMutationPositionCount:
    def test_counts_distinct_positions(self):
        hits = [{"mutation_pos": 15}, {"mutation_pos": 15}, {"mutation_pos": 20}]
        assert mod.unique_mutation_position_count(hits) == 2


class TestMutationAtPtmSite:
    def test_yes_when_position_matches(self):
        assert mod.mutation_at_ptm_site([{"mutation_pos": 59}], 59) == "yes"

    def test_no_when_position_does_not_match(self):
        assert mod.mutation_at_ptm_site([{"mutation_pos": 59}], 100) == "no"

    def test_ptm_pos_as_string_is_coerced(self):
        # ptm_pos sometimes arrives as a string (straight from argparse/CSV) --
        # must still compare correctly against the int mutation_pos.
        assert mod.mutation_at_ptm_site([{"mutation_pos": 59}], "59") == "yes", (
            "mutation_at_ptm_site must coerce ptm_pos to int before comparing, "
            "since callers may pass it through as a string"
        )


class TestReadExistingTable:
    def test_missing_file_returns_none_and_empty_list(self, tmp_path):
        header, rows = mod._read_existing_table(tmp_path / "does_not_exist.tsv")
        assert header is None and rows == [], f"a missing file should give (None, []), got ({header}, {rows})"

    def test_reads_utf16_file(self, tmp_path):
        path = tmp_path / "data.tsv"
        with path.open("w", encoding="utf-16", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=["UniProt", "gene"], delimiter="\t")
            writer.writeheader()
            writer.writerow({"UniProt": "P04637", "gene": "TP53"})

        header, rows = mod._read_existing_table(path)
        assert header == ["UniProt", "gene"] and rows == [{"UniProt": "P04637", "gene": "TP53"}], (
            f"a UTF-16 file (the pipeline's own output encoding) should be read "
            f"correctly, got header={header} rows={rows}"
        )

    def test_reads_utf8_file(self, tmp_path):
        path = tmp_path / "data.tsv"
        with path.open("w", encoding="utf-8", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=["UniProt", "gene"], delimiter="\t")
            writer.writeheader()
            writer.writerow({"UniProt": "P04637", "gene": "TP53"})

        header, rows = mod._read_existing_table(path)
        assert header == ["UniProt", "gene"], (
            f"a UTF-8 file should also be read correctly via the encoding fallback, got {header}"
        )


class TestResolveCifPath:
    def test_relative_path_resolves_under_models_root(self, monkeypatch, tmp_path):
        monkeypatch.setattr(mod, "MODELS_ROOT", tmp_path)
        result = mod.resolve_cif_path("P04637/model.cif")
        assert result == (tmp_path / "P04637" / "model.cif").resolve(), (
            f"a relative path should be resolved relative to MODELS_ROOT, got {result}"
        )

    def test_absolute_path_used_as_is(self, tmp_path):
        absolute = (tmp_path / "somewhere" / "model.cif").resolve()
        result = mod.resolve_cif_path(str(absolute))
        assert result == absolute, f"an absolute path should be used directly, not re-rooted under MODELS_ROOT, got {result}"
