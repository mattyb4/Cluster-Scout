"""Unit tests for the post-hoc PolyPhen-class filter in scripts/4_annotate.py
(applied to the wide-format proximity DB after Phase 2 has tagged mutations)."""
import pandas as pd
import pytest

from conftest import import_script

mod = import_script("4_annotate.py")


class TestFilterMutStr:
    def test_removes_entries_with_excluded_code(self):
        s = "A100B(PP:D,0.99)-1.00Å, C200D(PP:B,0.10)-2.00Å"
        result = mod._filter_mut_str(s, {"D"})
        assert result == "C200D(PP:B,0.10)-2.00Å", (
            f"the D-tagged (probably damaging) entry should be removed while the "
            f"B-tagged (benign) entry is kept, got {result!r}"
        )

    def test_empty_string_returns_unchanged(self):
        result = mod._filter_mut_str("", {"D"})
        assert result == "", f"an empty mutation string has nothing to filter -- expected '', got {result!r}"

    def test_empty_exclude_codes_returns_unchanged(self):
        s = "A100B(PP:D,0.99)-1.00Å"
        result = mod._filter_mut_str(s, set())
        assert result == s, (
            f"with no classes excluded, the string should pass through completely "
            f"unmodified, got {result!r}"
        )

    def test_unscored_entries_are_always_kept(self):
        # An entry with no (PP:...) tag at all -- e.g. PolyPhen never scored it.
        s = "A100B-1.00Å"
        result = mod._filter_mut_str(s, {"D", "P", "B"})
        assert result == s, (
            f"an entry with no PolyPhen tag has no code to match against exclude_codes "
            f"and must never be filtered out, even excluding all 3 classes, got {result!r}"
        )


class TestPositionsFromStr:
    def test_extracts_positions_in_order(self):
        s = "A100B-1.00Å, C200D-2.00Å"
        result = mod._positions_from_str(s)
        assert result == [100, 200], f"positions should be extracted in the order they appear, got {result}"

    def test_dedupes_repeated_positions(self):
        s = "A100B-1.00Å, C100D-1.50Å"
        result = mod._positions_from_str(s)
        assert result == [100], (
            f"two different substitutions at the SAME position should count once, "
            f"not appear twice, got {result}"
        )

    def test_empty_string_returns_empty_list(self):
        assert mod._positions_from_str("") == [], "an empty string has no positions to extract"

    def test_none_returns_empty_list(self):
        assert mod._positions_from_str(None) == [], (
            "a None value (e.g. an unset pandas cell) should be handled via the "
            "`mutation_str or ''` guard, not raise an AttributeError on .split()"
        )


class TestApplyPolyphenFilter:
    def _row(self, ptm_site, within, beyond, confirmed=""):
        return {
            "ptm_site": ptm_site,
            "mutations_within_5_positions": within,
            "mutations_more_than_5_positions": beyond,
            "confirmed_disrupting_mutations": confirmed,
            "mutation_count_within_5_positions": 0,
            "mutation_count_more_than_5_positions": 0,
            "unique_mutation_position_count_within_5_positions": 0,
            "unique_mutation_position_count_more_than_5_positions": 0,
            "morethan5_linear_distance": "",
            "mutation_at_ptm_site": "no",
        }

    def test_no_op_when_exclude_classes_is_empty(self):
        df = pd.DataFrame([self._row("S100", "A101B(PP:D,0.99)-1.00Å", "")])
        result = mod.apply_polyphen_filter(df, [])
        assert result is df, (
            "with nothing excluded, apply_polyphen_filter should return the SAME "
            "DataFrame object untouched (no-op fast path), not a filtered copy"
        )

    def test_no_op_when_exclude_classes_not_in_code_map(self):
        df = pd.DataFrame([self._row("S100", "A101B(PP:D,0.99)-1.00Å", "")])
        result = mod.apply_polyphen_filter(df, ["not_a_real_class"])
        assert result is df, (
            "an exclude_classes list with no entries recognized in _PP_CODE_MAP should "
            "also take the no-op fast path"
        )

    def test_removes_excluded_mutation_but_keeps_row_with_remaining_mutations(self):
        df = pd.DataFrame([self._row(
            "S100",
            "A101B(PP:D,0.99)-1.00Å, C102D(PP:B,0.10)-2.00Å",
            "",
        )])
        result = mod.apply_polyphen_filter(df, ["probably_damaging"])

        assert len(result) == 1, (
            f"the row still has one qualifying (benign) mutation left after filtering "
            f"-- it must survive, got {len(result)} row(s)"
        )
        assert result.iloc[0]["mutations_within_5_positions"] == "C102D(PP:B,0.10)-2.00Å", (
            f"the D-tagged mutation should be removed from the string, got "
            f"{result.iloc[0]['mutations_within_5_positions']!r}"
        )
        assert result.iloc[0]["mutation_count_within_5_positions"] == 1, (
            f"the count column must be recomputed to reflect the filtered string "
            f"(1 remaining), got {result.iloc[0]['mutation_count_within_5_positions']}"
        )

    def test_drops_row_when_all_mutations_filtered_out(self):
        df = pd.DataFrame([
            self._row("S100", "A101B(PP:D,0.99)-1.00Å, C102D(PP:B,0.10)-2.00Å", ""),  # survives
            self._row("S200", "E201F(PP:D,0.99)-1.00Å", ""),  # only a D mutation -> dropped
        ])
        result = mod.apply_polyphen_filter(df, ["probably_damaging"])

        assert list(result["ptm_site"]) == ["S100"], (
            f"S200's only mutation was entirely filtered out (both within_5 and "
            f"more_than_5 end up empty), so the whole row must be dropped -- expected "
            f"only S100 to remain, got {list(result['ptm_site'])}"
        )

    def test_recomputes_mutation_at_ptm_site_after_filtering(self):
        # The mutation AT the PTM site itself gets filtered out -> mutation_at_ptm_site
        # must flip from "yes" to "no", not keep the stale pre-filter value.
        row = self._row("S100", "S100A(PP:D,0.99)-0.00Å, C102D(PP:B,0.10)-2.00Å", "")
        row["mutation_at_ptm_site"] = "yes"
        df = pd.DataFrame([row])

        result = mod.apply_polyphen_filter(df, ["probably_damaging"])
        assert result.iloc[0]["mutation_at_ptm_site"] == "no", (
            f"after removing the D-tagged mutation at position 100 (the PTM site "
            f"itself), mutation_at_ptm_site must be recomputed to 'no', not left "
            f"stale as 'yes', got {result.iloc[0]['mutation_at_ptm_site']!r}"
        )

    def test_recomputes_linear_distance_after_filtering_beyond_5(self):
        row = self._row(
            "S100",
            "",
            "A120B(PP:D,0.99)-9.00Å, C130D(PP:B,0.10)-9.50Å",
        )
        row["morethan5_linear_distance"] = "20,30"
        df = pd.DataFrame([row])

        result = mod.apply_polyphen_filter(df, ["probably_damaging"])
        assert result.iloc[0]["morethan5_linear_distance"] == "30", (
            f"position 120 (linear distance 20 from PTM site 100) was removed by the "
            f"filter -- morethan5_linear_distance must be recomputed to only include "
            f"the surviving position 130 (distance 30), got "
            f"{result.iloc[0]['morethan5_linear_distance']!r}"
        )

    def test_confirmed_disrupting_mutations_column_is_also_filtered(self):
        df = pd.DataFrame([self._row(
            "S100",
            "A101B(PP:D,0.99)-1.00Å, C102D(PP:B,0.10)-2.00Å",
            "",
            confirmed="A101B(PP:D,0.99)-1.00Å",
        )])
        result = mod.apply_polyphen_filter(df, ["probably_damaging"])
        assert result.iloc[0]["confirmed_disrupting_mutations"] == "", (
            f"the confirmed_disrupting_mutations column must go through the same "
            f"class-based filtering as the main mutation columns, got "
            f"{result.iloc[0]['confirmed_disrupting_mutations']!r}"
        )
