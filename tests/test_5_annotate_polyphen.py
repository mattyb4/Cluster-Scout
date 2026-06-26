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
        assert pred == "D"
        assert float(score) == pytest.approx(0.999)

    def test_empty_hits_returns_blanks(self):
        assert mod._pp_best_prediction([]) == ("", "")

    def test_single_string_pred_and_float_score(self):
        hits = [{"dbnsfp": {"polyphen2": {"hdiv": {"pred": "B", "score": 0.012}}}}]
        pred, score = mod._pp_best_prediction(hits)
        assert pred == "B"
        assert float(score) == pytest.approx(0.012)

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
        assert pred == "D"
        assert float(score) == pytest.approx(0.995)

    def test_missing_hdiv_returns_blanks(self):
        hits = [{"dbnsfp": {"polyphen2": {"hvar": {"pred": "D", "score": 1.0}}}}]
        assert mod._pp_best_prediction(hits) == ("", "")


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
        assert result == "R175H(PP:D,0.999)-3.52Å(PAE:2.1)"

    def test_preserves_isoform_tag(self):
        result = mod.annotate_mutation_string(
            "V143A(isoform?)-8.10Å", "TP53", self.cache
        )
        assert result == "V143A(isoform?)(PP:B,0.012)-8.10Å"

    def test_skips_already_tagged(self):
        s = "R175H(PP:D,0.999)-3.52Å"
        assert mod.annotate_mutation_string(s, "TP53", self.cache) == s

    def test_unknown_mutation_left_unchanged(self):
        s = "X999Y-1.00Å"
        assert mod.annotate_mutation_string(s, "TP53", self.cache) == s

    def test_empty_string_passthrough(self):
        assert mod.annotate_mutation_string("", "TP53", self.cache) == ""

    def test_multiple_entries(self):
        s = "R175H-3.52Å, V143A-8.10Å"
        result = mod.annotate_mutation_string(s, "TP53", self.cache)
        assert "(PP:D,0.999)" in result
        assert "(PP:B,0.012)" in result


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
        assert p1 == "D"

        mod._pp_save_cache({("TP53", "R175H"): (p1, s1)})
        cache = mod._pp_load_cache()
        assert ("TP53", "R175H") in cache

    def test_returns_blanks_on_network_error(self, tmp_path, monkeypatch):
        monkeypatch.setattr(mod, "_PP_CACHE_FILE", tmp_path / "pp.tsv")
        monkeypatch.setattr(
            mod._pp_session, "get",
            lambda *a, **k: (_ for _ in ()).throw(mod.requests.RequestException("timeout"))
        )
        assert mod.fetch_polyphen("TP53", "R175H") == ("", "")

    def test_stop_codon_returns_blanks(self, tmp_path, monkeypatch):
        monkeypatch.setattr(mod, "_PP_CACHE_FILE", tmp_path / "pp.tsv")
        assert mod.fetch_polyphen("TP53", "R175*") == ("", "")
