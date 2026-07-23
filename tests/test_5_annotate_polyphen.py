"""Unit tests for PolyPhen-2 annotation functions in scripts/4_annotate.py."""
import pandas as pd
import pytest

from conftest import import_script

mod = import_script("4_annotate.py")


SAMPLE_HITS = [
    {
        "dbnsfp": {
            "polyphen2": {
                "hdiv": {
                    "pred": ["D", "P", "D"],
                    "score": [0.999, 0.632, 0.998],
                }
            }
        }
    }
]


class TestBestPrediction:
    def test_picks_most_severe_across_transcripts(self):
        pred, score = mod._pp_best_prediction(SAMPLE_HITS)
        assert pred == "D", (
            f"among D/P/D across transcripts, the most severe (Probably Damaging = D) "
            f"should win, got {pred!r}"
        )
        assert float(score) == pytest.approx(0.999), (
            f"the score reported should be the one paired with the winning prediction, got {score}"
        )

    def test_empty_hits_returns_blanks(self):
        result = mod._pp_best_prediction([])
        assert result == ("", ""), f"no hits means no prediction to report -- expected ('', ''), got {result}"

    def test_single_string_pred_and_float_score(self):
        hits = [{"dbnsfp": {"polyphen2": {"hdiv": {"pred": "B", "score": 0.012}}}}]
        pred, score = mod._pp_best_prediction(hits)
        assert pred == "B", (
            f"myvariant.info sometimes returns a scalar pred/score instead of a list "
            f"(single-transcript case) -- must handle both shapes, got pred={pred!r}"
        )
        assert float(score) == pytest.approx(0.012), f"expected the scalar score preserved, got {score}"

    def test_polyphen2_as_list(self):
        hits = [{
            "dbnsfp": {
                "polyphen2": [
                    {"hdiv": {"pred": ["B"], "score": [0.01]}},
                    {"hdiv": {"pred": ["D"], "score": [0.995]}},
                ]
            }
        }]
        pred, score = mod._pp_best_prediction(hits)
        assert pred == "D", (
            f"the 'polyphen2' field itself can be a list of per-transcript dicts (not "
            f"just a single dict) -- the more severe D across both must win, got {pred!r}"
        )
        assert float(score) == pytest.approx(0.995), f"expected the D transcript's score, got {score}"

    def test_missing_hdiv_returns_blanks(self):
        hits = [{"dbnsfp": {"polyphen2": {"hvar": {"pred": "D", "score": 1.0}}}}]
        result = mod._pp_best_prediction(hits)
        assert result == ("", ""), (
            f"only the HDIV PolyPhen-2 model is used -- an HVAR-only response has no "
            f"HDIV prediction and should return blanks, not fall back to HVAR, got {result}"
        )


class TestAnnotateMutationString:
    def setup_method(self):
        self.cache = {
            ("TP53", "R175H"): ("D", "0.999"),
            ("TP53", "V143A"): ("B", "0.012"),
        }

    def test_inserts_pp_tag_before_distance(self):
        result = mod.annotate_mutation_string(
            "R175H-3.52Å(PAE:2.1)", "TP53", self.cache
        )
        assert result == "R175H(PP:D,0.999)-3.52Å(PAE:2.1)", (
            f"the (PP:class,score) tag must be inserted right after the mutation label "
            f"and before the distance/PAE suffix, got {result!r}"
        )

    def test_preserves_isoform_tag(self):
        result = mod.annotate_mutation_string(
            "V143A(isoform?)-8.10Å", "TP53", self.cache
        )
        assert result == "V143A(isoform?)(PP:B,0.012)-8.10Å", (
            f"an existing '(isoform?)' tag must be preserved, with the PP tag inserted "
            f"after it (not overwriting or displacing it), got {result!r}"
        )

    def test_skips_already_tagged(self):
        s = "R175H(PP:D,0.999)-3.52Å"
        result = mod.annotate_mutation_string(s, "TP53", self.cache)
        assert result == s, (
            f"re-annotating an already-tagged string must be idempotent (no double "
            f"tagging), got {result!r}"
        )

    def test_unknown_mutation_left_unchanged(self):
        s = "X999Y-1.00Å"
        result = mod.annotate_mutation_string(s, "TP53", self.cache)
        assert result == s, (
            f"a mutation with no entry in the PolyPhen cache should be left untouched, "
            f"not tagged with blank/garbage values, got {result!r}"
        )

    def test_empty_string_passthrough(self):
        result = mod.annotate_mutation_string("", "TP53", self.cache)
        assert result == "", f"an empty input string should return empty, not raise, got {result!r}"

    def test_multiple_entries(self):
        s = "R175H-3.52Å, V143A-8.10Å"
        result = mod.annotate_mutation_string(s, "TP53", self.cache)
        assert "(PP:D,0.999)" in result, f"the first entry's tag should be present, got {result!r}"
        assert "(PP:B,0.012)" in result, f"the second entry's tag should also be present, got {result!r}"


class TestFetchPolyphen:
    def test_caches_and_serves(self, tmp_path, monkeypatch):
        monkeypatch.setattr(mod, "_PP_CACHE_FILE", tmp_path / "pp.tsv")
        calls = []

        def fake_get(url, **kwargs):
            calls.append(url)
            class R:
                status_code = 200
                def raise_for_status(self): pass
                def json(self):
                    return {"hits": SAMPLE_HITS}
            return R()

        monkeypatch.setattr(mod._pp_session, "get", fake_get)

        p1, s1 = mod.fetch_polyphen("TP53", "R175H")
        assert p1 == "D", f"the live-fetched prediction should reflect the mocked API response, got {p1!r}"

        mod._pp_save_cache({("TP53", "R175H"): (p1, s1)})
        cache = mod._pp_load_cache()
        assert ("TP53", "R175H") in cache, (
            f"a saved cache entry should round-trip through _pp_save_cache/_pp_load_cache, "
            f"but the key is missing from {cache}"
        )

    def test_returns_blanks_on_network_error(self, tmp_path, monkeypatch):
        monkeypatch.setattr(mod, "_PP_CACHE_FILE", tmp_path / "pp.tsv")
        monkeypatch.setattr(
            mod._pp_session, "get",
            lambda *a, **k: (_ for _ in ()).throw(mod.requests.RequestException("timeout"))
        )
        result = mod.fetch_polyphen("TP53", "R175H")
        assert result == ("", ""), (
            f"a network exception must be caught and return blanks, not propagate and "
            f"crash the whole annotation phase, got {result}"
        )

    def test_stop_codon_returns_blanks(self, tmp_path, monkeypatch):
        monkeypatch.setattr(mod, "_PP_CACHE_FILE", tmp_path / "pp.tsv")
        result = mod.fetch_polyphen("TP53", "R175*")
        assert result == ("", ""), (
            f"PolyPhen-2 doesn't score stop-codon changes -- must short-circuit to blanks "
            f"without even making an API call, got {result}"
        )

    def test_mutation_not_matching_pattern_returns_blanks_without_request(self, tmp_path, monkeypatch):
        monkeypatch.setattr(mod, "_PP_CACHE_FILE", tmp_path / "pp.tsv")
        calls = []
        monkeypatch.setattr(mod._pp_session, "get", lambda *a, **k: calls.append(1))

        result = mod.fetch_polyphen("TP53", "not-a-valid-mutation-string")
        assert result == ("", ""), (
            f"a mutation string that doesn't match MUT_RE at all (garbage input) must "
            f"return blanks via the early-return guard, got {result}"
        )
        assert calls == [], (
            "the early-return for an unparseable mutation string should happen BEFORE "
            f"any API call is attempted, but got {len(calls)} call(s)"
        )


class TestPpLookupSingle:
    def setup_method(self):
        self.cache = {("TP53", "R175H"): ("D", "0.99")}

    def test_finds_bare_mutation(self):
        result = mod._pp_lookup_single("R175H", "TP53", self.cache)
        assert result == ("D", "0.99"), (
            f"a bare mutation label with no distance suffix (e.g. an anchor_mutation "
            f"cell) should look up directly in the cache, got {result}"
        )

    def test_strips_isoform_tag_before_lookup(self):
        result = mod._pp_lookup_single("R175H(isoform?)", "TP53", self.cache)
        assert result == ("D", "0.99"), (
            f"the '(isoform?)' suffix must be stripped before the cache lookup, got {result}"
        )

    def test_unknown_mutation_returns_blanks(self):
        result = mod._pp_lookup_single("X999Y", "TP53", self.cache)
        assert result == ("", ""), (
            f"a mutation with no cache entry should return blanks, not raise, got {result}"
        )

    def test_none_mutation_returns_blanks(self):
        result = mod._pp_lookup_single(None, "TP53", self.cache)
        assert result == ("", ""), (
            f"a None mutation (e.g. a missing anchor_mutation cell) must not raise on "
            f".replace(), should return blanks, got {result}"
        )


class TestRunPolyphenPhaseMutationCols:
    def test_only_tags_specified_columns(self, tmp_path, monkeypatch):
        monkeypatch.setattr(mod, "_PP_CACHE_FILE", tmp_path / "pp.tsv")
        df = pd.DataFrame([{
            "gene": "TP53",
            "nearby_mutations": "R175H-3.52Å(PAE:2.1)",
            "mutations_within_5_positions": "V143A-1.00Å",
        }])
        # Pre-seed the cache so both mutations are already known -- no network call needed.
        pd.DataFrame([
            {"gene": "TP53", "mutation": "R175H", "pred": "D", "score": "0.99"},
            {"gene": "TP53", "mutation": "V143A", "pred": "B", "score": "0.01"},
        ]).to_csv(mod._PP_CACHE_FILE, sep="\t", index=False)

        mod.run_polyphen_phase(df, mutation_cols=["nearby_mutations"])

        assert "(PP:D,0.99)" in df.iloc[0]["nearby_mutations"], (
            f"nearby_mutations was passed via mutation_cols -- it should be tagged, "
            f"got {df.iloc[0]['nearby_mutations']!r}"
        )
        assert df.iloc[0]["mutations_within_5_positions"] == "V143A-1.00Å", (
            f"mutations_within_5_positions was NOT passed in mutation_cols -- it must "
            f"be left completely untouched, got {df.iloc[0]['mutations_within_5_positions']!r}"
        )

    def test_default_mutation_cols_covers_ptm_proximity_columns(self, tmp_path, monkeypatch):
        monkeypatch.setattr(mod, "_PP_CACHE_FILE", tmp_path / "pp.tsv")
        df = pd.DataFrame([{
            "gene": "TP53",
            "mutations_within_5_positions": "R175H-3.52Å",
            "mutations_more_than_5_positions": "V143A-8.10Å",
        }])
        pd.DataFrame([
            {"gene": "TP53", "mutation": "R175H", "pred": "D", "score": "0.99"},
            {"gene": "TP53", "mutation": "V143A", "pred": "B", "score": "0.01"},
        ]).to_csv(mod._PP_CACHE_FILE, sep="\t", index=False)

        mod.run_polyphen_phase(df)

        assert "(PP:D,0.99)" in df.iloc[0]["mutations_within_5_positions"], (
            "with no mutation_cols override, the default _MUTATION_COLS (ptm-proximity's "
            "within/beyond columns) should still be tagged"
        )
        assert "(PP:B,0.01)" in df.iloc[0]["mutations_more_than_5_positions"], (
            "the second default column should also be tagged"
        )

    def test_bare_mutation_cols_fetch_but_are_not_tagged_inline(self, tmp_path, monkeypatch):
        monkeypatch.setattr(mod, "_PP_CACHE_FILE", tmp_path / "pp.tsv")
        df = pd.DataFrame([{
            "gene": "TP53",
            "anchor_mutation": "R175H",
            "nearby_mutations": "",
        }])
        calls = []

        def fake_fetch(gene, mut):
            calls.append((gene, mut))
            return "D", "0.99"

        monkeypatch.setattr(mod, "fetch_polyphen", fake_fetch)

        cache = mod.run_polyphen_phase(
            df, mutation_cols=["nearby_mutations"], bare_mutation_cols=("anchor_mutation",),
        )

        assert ("TP53", "R175H") in calls, (
            f"anchor_mutation is an uncached bare mutation label -- it must trigger a "
            f"live fetch (via bare_mutation_cols), not silently stay blank forever, "
            f"got fetched pairs: {calls}"
        )
        assert cache.get(("TP53", "R175H")) == ("D", "0.99"), (
            f"the fetched result must land in the returned cache so callers can look "
            f"it up afterward, got {cache.get(('TP53', 'R175H'))}"
        )
        assert df.iloc[0]["anchor_mutation"] == "R175H", (
            "bare_mutation_cols have no '-distance' entry format to tag inline -- "
            f"anchor_mutation must be left as the plain label, got {df.iloc[0]['anchor_mutation']!r}"
        )

    def test_bare_mutation_col_already_cached_triggers_no_fetch(self, tmp_path, monkeypatch):
        monkeypatch.setattr(mod, "_PP_CACHE_FILE", tmp_path / "pp.tsv")
        pd.DataFrame([
            {"gene": "TP53", "mutation": "R175H", "pred": "D", "score": "0.99"},
        ]).to_csv(mod._PP_CACHE_FILE, sep="\t", index=False)
        df = pd.DataFrame([{"gene": "TP53", "anchor_mutation": "R175H", "nearby_mutations": ""}])
        calls = []
        monkeypatch.setattr(mod, "fetch_polyphen", lambda g, m: calls.append((g, m)))

        mod.run_polyphen_phase(
            df, mutation_cols=["nearby_mutations"], bare_mutation_cols=("anchor_mutation",),
        )
        assert calls == [], (
            f"anchor_mutation is already in the on-disk cache -- it must not trigger "
            f"another fetch, got {calls}"
        )
