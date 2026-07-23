"""Unit tests for AIUPred disorder-score (Phase 4) functions in scripts/4_annotate.py."""
import json

import pandas as pd
import pytest

from conftest import import_script

mod = import_script("4_annotate.py")


class FakeAiupredResponse:
    def __init__(self, json_data, status_code=200):
        self._json_data = json_data
        self.status_code = status_code

    def raise_for_status(self):
        if self.status_code >= 400:
            raise RuntimeError(f"HTTP {self.status_code}")

    def json(self):
        return self._json_data


class TestAiupredCache:
    def test_round_trip(self, tmp_path, monkeypatch):
        monkeypatch.setattr(mod, "_AIUPRED_CACHE_FILE", tmp_path / "aiupred.tsv")
        cache = {("P04637", "general"): {1: 0.1, 2: 0.9}, ("P04637", "binding"): {1: 0.05}}
        mod._aiupred_save_cache(cache)
        result = mod._aiupred_load_cache()
        assert result == cache, (
            f"loading immediately after saving should reproduce the same "
            f"(uniprot_id, analysis_type) -> {{position: score}} dict, got {result}"
        )

    def test_missing_file_returns_empty_dict(self, tmp_path, monkeypatch):
        monkeypatch.setattr(mod, "_AIUPRED_CACHE_FILE", tmp_path / "does_not_exist.tsv")
        result = mod._aiupred_load_cache()
        assert result == {}, f"a never-written cache file should load as {{}}, got {result}"

    def test_malformed_scores_json_degrades_to_empty_scores(self, tmp_path, monkeypatch):
        cache_file = tmp_path / "aiupred.tsv"
        monkeypatch.setattr(mod, "_AIUPRED_CACHE_FILE", cache_file)
        pd.DataFrame([
            {"uniprot_id": "P04637", "analysis_type": "general", "scores_json": "not valid json{{{"},
        ]).to_csv(cache_file, sep="\t", index=False)

        result = mod._aiupred_load_cache()
        assert result == {("P04637", "general"): {}}, (
            f"a corrupted scores_json cell should degrade to an empty scores dict for "
            f"that entry (not skip the row or crash the whole cache load), got {result}"
        )

    def test_unreadable_file_returns_empty_dict(self, tmp_path, monkeypatch):
        cache_file = tmp_path / "aiupred.tsv"
        cache_file.write_bytes(b"\xff\xfe\x00garbage-not-a-real-tsv")
        monkeypatch.setattr(mod, "_AIUPRED_CACHE_FILE", cache_file)

        result = mod._aiupred_load_cache()
        assert result == {}, (
            f"a cache file that fails to parse entirely should degrade to {{}} via the "
            f"outer try/except, not raise, got {result}"
        )


class TestFetchAiupredAll:
    def test_binding_call_returns_general_and_binding_scores(self, monkeypatch):
        monkeypatch.setattr(
            mod._aiupred_session, "get",
            lambda url, params=None, timeout=None: FakeAiupredResponse(
                {"AIUPred": [0.1, 0.2, 0.9], "AIUPred-binding": [0.05, 0.15, 0.85]},
            ),
        )
        result = mod.fetch_aiupred_all("binding", "P04637")
        assert result == {
            "general": {1: 0.1, 2: 0.2, 3: 0.9},
            "binding": {1: 0.05, 2: 0.15, 3: 0.85},
        }, (
            f"the API's 'AIUPred'/'AIUPred-binding' keys should map to canonical "
            f"'general'/'binding' names, with 1-based positions, got {result}"
        )

    def test_unrecognized_key_is_ignored(self, monkeypatch):
        monkeypatch.setattr(
            mod._aiupred_session, "get",
            lambda url, params=None, timeout=None: FakeAiupredResponse(
                {"SomeOtherAnalysis": [0.1, 0.2]},
            ),
        )
        result = mod.fetch_aiupred_all("binding", "P04637")
        assert result == {}, (
            f"a response key not in _AIUPRED_KEY_MAP should be silently ignored, not "
            f"raise or appear under a garbage type name, got {result}"
        )

    def test_empty_score_list_is_ignored(self, monkeypatch):
        monkeypatch.setattr(
            mod._aiupred_session, "get",
            lambda url, params=None, timeout=None: FakeAiupredResponse({"AIUPred": []}),
        )
        result = mod.fetch_aiupred_all("binding", "P04637")
        assert result == {}, (
            f"an empty score list has no per-position data to report and must be "
            f"skipped, not produce an empty-but-present 'general' key, got {result}"
        )

    def test_non_numeric_score_list_is_ignored(self, monkeypatch):
        monkeypatch.setattr(
            mod._aiupred_session, "get",
            lambda url, params=None, timeout=None: FakeAiupredResponse({"AIUPred": ["not", "numbers"]}),
        )
        result = mod.fetch_aiupred_all("binding", "P04637")
        assert result == {}, (
            f"a malformed (non-numeric) score list should be rejected by the type "
            f"check, not raise a ValueError trying to float() a string, got {result}"
        )

    def test_returns_empty_dict_on_network_error(self, monkeypatch):
        monkeypatch.setattr(
            mod._aiupred_session, "get",
            lambda url, params=None, timeout=None: (_ for _ in ()).throw(mod.requests.RequestException("timeout")),
        )
        result = mod.fetch_aiupred_all("binding", "P04637")
        assert result == {}, (
            f"a network exception must be caught and return {{}}, not propagate, got {result}"
        )

    def test_returns_empty_dict_when_response_is_not_a_dict(self, monkeypatch):
        monkeypatch.setattr(
            mod._aiupred_session, "get",
            lambda url, params=None, timeout=None: FakeAiupredResponse(["unexpected", "list", "response"]),
        )
        result = mod.fetch_aiupred_all("binding", "P04637")
        assert result == {}, (
            f"a JSON response that isn't a dict (e.g. an error returned as a bare list) "
            f"should be rejected, not indexed into with .items() and crash, got {result}"
        )


class TestRunAiupredPhase:
    def test_all_cached_makes_no_fetch_calls(self, tmp_path, monkeypatch):
        monkeypatch.setattr(mod, "_AIUPRED_CACHE_FILE", tmp_path / "aiupred.tsv")
        mod._aiupred_save_cache({
            ("P04637", "general"): {1: 0.5}, ("P04637", "binding"): {1: 0.4},
        })
        calls = []
        monkeypatch.setattr(mod, "fetch_aiupred_all", lambda api_type, uid: calls.append(uid) or {})

        df = pd.DataFrame([{"UniProt": "P04637"}])
        result = mod.run_aiupred_phase(df)

        assert calls == [], (
            f"a protein already cached for both general and binding scores should not "
            f"trigger any fetch, got {len(calls)} call(s)"
        )
        assert result["general"]["P04637"] == {1: 0.5}, (
            f"the cached scores should be returned even without a fresh fetch, got {result['general']}"
        )

    def test_fetches_only_uncached_proteins(self, tmp_path, monkeypatch):
        monkeypatch.setattr(mod, "_AIUPRED_CACHE_FILE", tmp_path / "aiupred.tsv")
        mod._aiupred_save_cache({("P00001", "general"): {1: 0.1}, ("P00001", "binding"): {1: 0.1}})
        calls = []

        def fake_fetch(api_type, uid):
            calls.append(uid)
            return {"general": {1: 0.9}, "binding": {1: 0.8}}

        monkeypatch.setattr(mod, "fetch_aiupred_all", fake_fetch)

        df = pd.DataFrame([{"UniProt": "P00001"}, {"UniProt": "P00002"}])
        result = mod.run_aiupred_phase(df)

        assert calls == ["P00002"], (
            f"only the uncached protein (P00002) should be fetched -- P00001 already "
            f"has both score types cached, got fetch calls for {calls}"
        )
        assert result["general"]["P00002"] == {1: 0.9}, (
            f"the newly fetched protein's scores should appear in the result, got {result['general']}"
        )

    def test_saves_cache_after_fetching_new_data(self, tmp_path, monkeypatch):
        cache_file = tmp_path / "aiupred.tsv"
        monkeypatch.setattr(mod, "_AIUPRED_CACHE_FILE", cache_file)
        monkeypatch.setattr(mod, "fetch_aiupred_all", lambda api_type, uid: {"general": {1: 0.5}, "binding": {1: 0.5}})

        df = pd.DataFrame([{"UniProt": "P00003"}])
        mod.run_aiupred_phase(df)

        assert cache_file.exists(), (
            "after fetching new scores, the cache file must be written to disk so "
            "subsequent runs don't re-fetch the same protein"
        )
