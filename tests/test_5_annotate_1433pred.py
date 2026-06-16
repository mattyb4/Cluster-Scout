"""Unit tests for scripts/5_annotate_1433pred.py."""
import json

import pytest

from conftest import import_script

mod = import_script("5_annotate_1433pred.py")


SAMPLE_API_RESPONSE = [
    {"Site": 59, "Peptide": "HIKIQNtGDYY", "SVM": "-0.577", "PSSM": "0.181",
     "ANN": "0.318", "Consensus": "-0.026", "pSer/Thr": "Yes"},
    {"Site": 76, "Peptide": "GEVTEAsIGGS", "SVM": "1.234", "PSSM": "0.900",
     "ANN": "0.750", "Consensus": "1.046", "pSer/Thr": "-"},
    {"Site": 100, "Peptide": "XXXXXXsXXXX", "SVM": "0.0", "PSSM": "0.0",
     "ANN": "0.0", "Consensus": "0.0", "pSer/Thr": "-"},
]


class TestBuildSiteScoreMap:
    def test_maps_position_to_consensus(self):
        scores = mod.build_site_score_map(SAMPLE_API_RESPONSE)
        assert scores[59] == pytest.approx(-0.026)
        assert scores[76] == pytest.approx(1.046)

    def test_skips_malformed_entries(self):
        bad = [{"Site": "x", "Consensus": "1.0"}, {"Site": 5, "Consensus": "bad"}]
        assert mod.build_site_score_map(bad) == {}

    def test_empty_response(self):
        assert mod.build_site_score_map([]) == {}


class TestAnnotateRow:
    def setup_method(self):
        self.scores = mod.build_site_score_map(SAMPLE_API_RESPONSE)

    def test_positive_consensus_is_yes(self):
        binding, consensus = mod.annotate_row("S76", self.scores)
        assert binding == "Yes"
        assert float(consensus) == pytest.approx(1.046)

    def test_negative_consensus_is_no(self):
        binding, consensus = mod.annotate_row("T59", self.scores)
        assert binding == "No"
        assert float(consensus) == pytest.approx(-0.026)

    def test_zero_consensus_is_no(self):
        binding, _ = mod.annotate_row("S100", self.scores)
        assert binding == "No"

    def test_non_ser_thr_is_blank(self):
        binding, consensus = mod.annotate_row("Y62", self.scores)
        assert binding == ""
        assert consensus == ""

    def test_lysine_is_blank(self):
        binding, consensus = mod.annotate_row("K43", self.scores)
        assert binding == ""
        assert consensus == ""

    def test_position_not_in_response_is_blank(self):
        binding, consensus = mod.annotate_row("T999", self.scores)
        assert binding == ""
        assert consensus == ""

    def test_malformed_ptm_site_is_blank(self):
        binding, consensus = mod.annotate_row("", self.scores)
        assert binding == ""
        assert consensus == ""


class TestFetch1433pred:
    def test_caches_response_and_serves_from_cache_on_second_call(self, tmp_path, monkeypatch):
        monkeypatch.setattr(mod, "CACHE_DIR", tmp_path)
        calls = []

        def fake_get(url, timeout=30):
            calls.append(url)
            class R:
                status_code = 200
                def json(self):
                    return SAMPLE_API_RESPONSE
            return R()

        monkeypatch.setattr(mod.requests, "get", fake_get)

        data1 = mod.fetch_1433pred("Q06124")
        data2 = mod.fetch_1433pred("Q06124")

        assert data1 == SAMPLE_API_RESPONSE
        assert data2 == SAMPLE_API_RESPONSE
        assert len(calls) == 1  # second call served from cache

    def test_returns_none_on_non_200(self, tmp_path, monkeypatch):
        monkeypatch.setattr(mod, "CACHE_DIR", tmp_path)

        def fake_get(url, timeout=30):
            class R:
                status_code = 404
            return R()

        monkeypatch.setattr(mod.requests, "get", fake_get)
        assert mod.fetch_1433pred("NOTANID") is None

    def test_returns_none_on_network_error(self, tmp_path, monkeypatch):
        monkeypatch.setattr(mod, "CACHE_DIR", tmp_path)
        monkeypatch.setattr(
            mod.requests, "get",
            lambda *a, **k: (_ for _ in ()).throw(mod.requests.RequestException("timeout"))
        )
        assert mod.fetch_1433pred("Q06124") is None
