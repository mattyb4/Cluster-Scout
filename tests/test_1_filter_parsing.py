"""Unit tests for the small pure parsing/formatting helpers in scripts/1_filter.py."""
import pandas as pd
import pytest


@pytest.fixture
def mod(filter_module):
    return filter_module


class TestCleanStrList:
    def test_dedupes_preserving_order(self, mod):
        values = pd.Series(["b", "a", "b", "c", "a"])
        assert mod.clean_str_list(values) == "b; a; c"

    def test_drops_nan_and_blank_entries(self, mod):
        values = pd.Series(["x", None, "  ", float("nan")])
        assert mod.clean_str_list(values) == "x"

    def test_empty_series_returns_empty_string(self, mod):
        assert mod.clean_str_list(pd.Series([], dtype=object)) == ""


class TestIsSimpleSubstitution:
    @pytest.mark.parametrize(
        "change,expected",
        [
            ("R482H", True),
            ("E291Q", True),
            ("p.R482H", False),  # 'p.' prefix must already be stripped by the caller
            ("R482*", False),  # stop-codon mutations are not simple substitutions
            ("R213_R214del", False),  # in-frame deletions are excluded
            (None, False),
            (float("nan"), False),
        ],
    )
    def test_matches_expected(self, mod, change, expected):
        assert mod.is_simple_substitution(change) is expected


class TestBuildPtmSite:
    def test_with_type(self, mod):
        row = pd.Series({"Residue": "S", "Position": 516.0, "Type": "Phosphorylation"})
        assert mod.build_ptm_site(row) == "S516:Phosphorylation"

    def test_without_type(self, mod):
        row = pd.Series({"Residue": "S", "Position": 516.0, "Type": float("nan")})
        assert mod.build_ptm_site(row) == "S516"

    def test_position_already_a_string(self, mod):
        row = pd.Series({"Residue": "K", "Position": "43", "Type": "Ubiquitination"})
        assert mod.build_ptm_site(row) == "K43:Ubiquitination"


class TestFormatMutationWithCount:
    def test_basic(self, mod):
        row = pd.Series({"mutation": "R482H", "affected_cases": 5})
        assert mod.format_mutation_with_count(row) == "R482H (5)"


class TestParseMutationSite:
    def test_list_literal(self, mod):
        assert mod.parse_mutation_site("['D120N', 'E127D']") == ["D120N", "E127D"]

    def test_empty_list_literal(self, mod):
        assert mod.parse_mutation_site("[]") == []

    def test_nan(self, mod):
        assert mod.parse_mutation_site(float("nan")) == []

    def test_single_value_string(self, mod):
        assert mod.parse_mutation_site("D120N") == ["D120N"]

    def test_falls_back_to_regex_for_malformed_literal(self, mod):
        assert mod.parse_mutation_site("Mutations D120N and E127D reported") == ["D120N", "E127D"]
