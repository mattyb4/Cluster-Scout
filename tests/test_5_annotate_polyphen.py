"""Unit tests for PolyPhen-2 annotation functions in scripts/4_annotate.py."""
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
