"""Unit tests for 14-3-3 annotation functions in scripts/4_annotate.py."""
import json
from pathlib import Path

import pandas as pd
import pytest

from conftest import import_script

mod = import_script("4_annotate.py")


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
        assert scores[59] == pytest.approx(-0.026), (
            f"position 59's consensus score should be extracted as-is, got {scores.get(59)}"
        )
        assert scores[76] == pytest.approx(1.046), (
            f"position 76's consensus score should be extracted as-is, got {scores.get(76)}"
        )

    def test_skips_malformed_entries(self):
        bad = [{"Site": "x", "Consensus": "1.0"}, {"Site": 5, "Consensus": "bad"}]
        result = mod.build_site_score_map(bad)
        assert result == {}, (
            f"a non-numeric Site or Consensus value should be skipped via the per-entry "
            f"try/except, not raise or produce a garbage entry, got {result}"
        )

    def test_empty_response(self):
        result = mod.build_site_score_map([])
        assert result == {}, f"an empty API response should map to an empty score dict, got {result}"


class TestAnnotateRow:
    def setup_method(self):
        self.scores = mod.build_site_score_map(SAMPLE_API_RESPONSE)

    def test_positive_consensus_is_yes(self):
        binding, consensus = mod.annotate_1433_row("S76", self.scores)
        assert binding == "Yes", f"a positive consensus score (1.046) should classify as binding, got {binding!r}"
        assert float(consensus) == pytest.approx(1.046), (
            f"the reported consensus score should match the source data, got {consensus}"
        )

    def test_negative_consensus_is_no(self):
        binding, consensus = mod.annotate_1433_row("T59", self.scores)
        assert binding == "No", f"a negative consensus score (-0.026) should classify as non-binding, got {binding!r}"
        assert float(consensus) == pytest.approx(-0.026), f"expected the negative score preserved, got {consensus}"

    def test_zero_consensus_is_no(self):
        binding, _ = mod.annotate_1433_row("S100", self.scores)
        assert binding == "No", (
            f"a consensus score of exactly 0 is the boundary case -- 'positive' means "
            f"strictly > 0, so 0 must classify as No, got {binding!r}"
        )

    def test_non_ser_thr_is_blank(self):
        binding, consensus = mod.annotate_1433_row("Y62", self.scores)
        assert binding == "" and consensus == "", (
            f"14-3-3 binding predictions only apply to Ser/Thr residues -- a Tyr site "
            f"should be left blank, not scored, got binding={binding!r} consensus={consensus!r}"
        )

    def test_lysine_is_blank(self):
        binding, consensus = mod.annotate_1433_row("K43", self.scores)
        assert binding == "" and consensus == "", (
            f"a non-Ser/Thr residue like Lysine should also be blank, got binding={binding!r} consensus={consensus!r}"
        )

    def test_position_not_in_response_is_blank(self):
        binding, consensus = mod.annotate_1433_row("T999", self.scores)
        assert binding == "" and consensus == "", (
            f"a position the API never returned a score for should be blank, not raise "
            f"a KeyError, got binding={binding!r} consensus={consensus!r}"
        )

    def test_malformed_ptm_site_is_blank(self):
        binding, consensus = mod.annotate_1433_row("", self.scores)
        assert binding == "" and consensus == "", (
            f"an empty/unparseable ptm_site string should be handled gracefully as blank, "
            f"not raise, got binding={binding!r} consensus={consensus!r}"
        )


class TestFetch1433pred:
    def test_caches_response_and_serves_from_cache_on_second_call(self, tmp_path, monkeypatch):
        monkeypatch.setattr(mod, "_1433_CACHE_DIR", tmp_path)
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

        assert data1 == SAMPLE_API_RESPONSE, "the first call's parsed JSON should match the mocked API response"
        assert data2 == SAMPLE_API_RESPONSE, "the second (cached) call should return the identical data"
        assert len(calls) == 1, (
            f"the second call for the same protein must be served from the on-disk cache "
            f"with no new HTTP request, got {len(calls)} call(s)"
        )

    def test_returns_none_on_non_200(self, tmp_path, monkeypatch):
        monkeypatch.setattr(mod, "_1433_CACHE_DIR", tmp_path)

        def fake_get(url, timeout=30):
            class R:
                status_code = 404
            return R()

        monkeypatch.setattr(mod.requests, "get", fake_get)
        result = mod.fetch_1433pred("NOTANID")
        assert result is None, f"a 404 response should return None, not raise or return an empty list, got {result!r}"

    def test_returns_none_on_network_error(self, tmp_path, monkeypatch):
        monkeypatch.setattr(mod, "_1433_CACHE_DIR", tmp_path)
        monkeypatch.setattr(
            mod.requests, "get",
            lambda *a, **k: (_ for _ in ()).throw(mod.requests.RequestException("timeout"))
        )
        result = mod.fetch_1433pred("Q06124")
        assert result is None, (
            f"a network-level exception (timeout, connection error) should be caught and "
            f"return None rather than propagate and crash the whole annotation phase, got {result!r}"
        )


class TestLoadConfirmedSites:
    def _write_excel(self, path: Path, rows: list[dict]) -> None:
        pd.DataFrame(rows).to_excel(path, index=False)

    def test_loads_valid_ser_thr_entries(self, tmp_path):
        xls = tmp_path / "sites.xlsx"
        self._write_excel(xls, [
            {"Uniprot Name": "X_HUMAN", "Uniprot ID": "P12345", "Site": 100, "Residue": "S", "Motif Sequence": "AAAA", "PMID": "11111111"},
            {"Uniprot Name": "X_HUMAN", "Uniprot ID": "P12345", "Site": 200, "Residue": "T", "Motif Sequence": "BBBB", "PMID": "22222222"},
        ])
        result = mod.load_confirmed_sites(xls)
        assert result[("P12345", 100)] == "11111111", "a valid Ser entry should map (uniprot, position) -> PMID"
        assert result[("P12345", 200)] == "22222222", "a valid Thr entry should also be loaded"

    def test_skips_non_ser_thr_residues(self, tmp_path):
        xls = tmp_path / "sites.xlsx"
        self._write_excel(xls, [
            {"Uniprot Name": "X_HUMAN", "Uniprot ID": "P12345", "Site": 50, "Residue": "Y", "Motif Sequence": "AAAA", "PMID": "11111111"},
            {"Uniprot Name": "X_HUMAN", "Uniprot ID": "P12345", "Site": 51, "Residue": "1", "Motif Sequence": "BBBB", "PMID": "22222222"},
        ])
        result = mod.load_confirmed_sites(xls)
        assert len(result) == 0, (
            f"14-3-3 only binds Ser/Thr sites -- a Tyr residue and a garbage '1' residue "
            f"should both be excluded, got {result}"
        )

    def test_strips_whitespace_from_residue_and_pmid(self, tmp_path):
        xls = tmp_path / "sites.xlsx"
        self._write_excel(xls, [
            {"Uniprot Name": "X_HUMAN", "Uniprot ID": "P12345", "Site": 75, "Residue": "S ", "Motif Sequence": "AAAA", "PMID": "\xa033333333"},
        ])
        result = mod.load_confirmed_sites(xls)
        assert ("P12345", 75) in result, (
            f"a trailing space on 'Residue' must not prevent the Ser/Thr check from "
            f"matching, but the entry is missing from {result}"
        )
        assert result[("P12345", 75)] == "33333333", (
            f"a leading non-breaking space (\\xa0) in PMID must be stripped, got {result[('P12345', 75)]!r}"
        )

    def test_returns_empty_dict_when_file_missing(self, tmp_path):
        result = mod.load_confirmed_sites(tmp_path / "nonexistent.xlsx")
        assert result == {}, (
            f"a missing confirmed-sites spreadsheet should degrade to an empty dict "
            f"(no confirmed sites known), not raise, got {result}"
        )


class TestAnnotateConfirmed:
    def setup_method(self):
        self.confirmed = {("P12345", 100): "11111111", ("P12345", 200): "22222222"}

    def test_confirmed_site_returns_yes_and_pmid(self):
        site, pmid = mod.annotate_confirmed("P12345", "S100", self.confirmed)
        assert site == "Yes", f"a (uniprot, position) pair present in the confirmed-sites dict should report Yes, got {site!r}"
        assert pmid == "11111111", f"the associated PMID should be returned alongside, got {pmid!r}"

    def test_unmatched_position_returns_blank(self):
        site, pmid = mod.annotate_confirmed("P12345", "S999", self.confirmed)
        assert site == "" and pmid == "", (
            f"a position not in the confirmed-sites dict should be blank, not raise, "
            f"got site={site!r} pmid={pmid!r}"
        )

    def test_wrong_uniprot_returns_blank(self):
        site, pmid = mod.annotate_confirmed("Q99999", "S100", self.confirmed)
        assert site == "" and pmid == "", (
            f"the same position (100) confirmed for a DIFFERENT protein must not match "
            f"here -- the lookup key is (uniprot, position) jointly, got site={site!r} pmid={pmid!r}"
        )

    def test_malformed_ptm_site_returns_blank(self):
        site, pmid = mod.annotate_confirmed("P12345", "", self.confirmed)
        assert site == "" and pmid == "", (
            f"an empty/unparseable ptm_site should be handled gracefully as blank, not "
            f"raise, got site={site!r} pmid={pmid!r}"
        )
