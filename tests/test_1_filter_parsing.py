"""Unit tests for the small pure parsing/formatting helpers in scripts/1_filter.py."""
import pandas as pd
import pytest


@pytest.fixture
def mod(filter_module):
    return filter_module


class TestCleanStrList:
    def test_dedupes_preserving_order(self, mod):
        values = pd.Series(["b", "a", "b", "c", "a"])
        result = mod.clean_str_list(values)
        assert result == "b; a; c", (
            f"duplicates should be dropped while keeping first-seen order, got {result!r}"
        )

    def test_drops_nan_and_blank_entries(self, mod):
        values = pd.Series(["x", None, "  ", float("nan")])
        result = mod.clean_str_list(values)
        assert result == "x", (
            f"None, whitespace-only, and NaN entries must all be filtered out, got {result!r}"
        )

    def test_empty_series_returns_empty_string(self, mod):
        result = mod.clean_str_list(pd.Series([], dtype=object))
        assert result == "", f"an empty Series has nothing to join, so the result should be '', got {result!r}"


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
        result = mod.is_simple_substitution(change)
        assert result is expected, (
            f"is_simple_substitution({change!r}) should be {expected} -- got {result}. "
            "Only single-residue missense changes (letter, digits, letter) count as simple "
            "substitutions; stop-codons, deletions, and non-string inputs must not."
        )


class TestBuildPtmSite:
    def test_with_type(self, mod):
        row = pd.Series({"Residue": "S", "Position": 516.0, "Type": "Phosphorylation"})
        result = mod.build_ptm_site(row)
        assert result == "S516:Phosphorylation", (
            f"residue+position+type should combine as 'S516:Phosphorylation', got {result!r}"
        )

    def test_without_type(self, mod):
        row = pd.Series({"Residue": "S", "Position": 516.0, "Type": float("nan")})
        result = mod.build_ptm_site(row)
        assert result == "S516", (
            f"a NaN Type should be omitted entirely (no trailing ':'), got {result!r}"
        )

    def test_position_already_a_string(self, mod):
        row = pd.Series({"Residue": "K", "Position": "43", "Type": "Ubiquitination"})
        result = mod.build_ptm_site(row)
        assert result == "K43:Ubiquitination", (
            f"a Position already given as a string should work the same as a float, got {result!r}"
        )

    def test_non_numeric_position_falls_back_to_raw_string(self, mod):
        # int(float("abc")) raises -- the except branch must fall back to the raw string
        # rather than propagating the exception or silently dropping the position.
        row = pd.Series({"Residue": "S", "Position": "abc", "Type": "Phosphorylation"})
        result = mod.build_ptm_site(row)
        assert result == "Sabc:Phosphorylation", (
            "a non-numeric Position value can't go through int(float(...)), so build_ptm_site "
            f"must fall back to using it as-is (str().strip()) rather than raising, got {result!r}"
        )

    def test_nan_residue_is_omitted(self, mod):
        row = pd.Series({"Residue": float("nan"), "Position": 100.0, "Type": "Phosphorylation"})
        result = mod.build_ptm_site(row)
        assert result == "100:Phosphorylation", (
            f"a NaN Residue should contribute an empty string, not the literal 'nan', got {result!r}"
        )


class TestFormatMutationWithCount:
    def test_basic(self, mod):
        row = pd.Series({"mutation": "R482H", "affected_cases": 5})
        result = mod.format_mutation_with_count(row)
        assert result == "R482H (5)", f"expected the mutation label with a parenthesized count, got {result!r}"


class TestParseMutationSite:
    def test_list_literal(self, mod):
        result = mod.parse_mutation_site("['D120N', 'E127D']")
        assert result == ["D120N", "E127D"], (
            f"a Python-list-literal string should parse via ast.literal_eval, got {result}"
        )

    def test_empty_list_literal(self, mod):
        result = mod.parse_mutation_site("[]")
        assert result == [], f"an empty list literal should parse to an empty list, got {result}"

    def test_nan(self, mod):
        result = mod.parse_mutation_site(float("nan"))
        assert result == [], f"a NaN cell has no mutations to report -- should return [], got {result}"

    def test_single_value_string(self, mod):
        result = mod.parse_mutation_site("D120N")
        assert result == ["D120N"], (
            f"a bare (non-list-literal) mutation string should still be wrapped in a list, got {result}"
        )

    def test_falls_back_to_regex_for_malformed_literal(self, mod):
        result = mod.parse_mutation_site("Mutations D120N and E127D reported")
        assert result == ["D120N", "E127D"], (
            "text that fails ast.literal_eval (free-form prose, not a Python literal) must "
            f"still extract mutation tokens via the regex fallback, got {result}"
        )
