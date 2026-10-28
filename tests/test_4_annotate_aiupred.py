"""Unit tests for AIUPred disorder-score (Phase 4) functions in scripts/4_annotate.py."""

import json
import subprocess
from types import SimpleNamespace

import pandas as pd
import pytest
from conftest import import_script

mod = import_script("4_annotate.py")
_real_run_local = mod._aiupred_run_local


@pytest.fixture(autouse=True)
def web_only_unless_testing_local(request, monkeypatch):
    """The phase tests below exercise the web path; keep them off the real local helper."""
    if "local" not in request.node.nodeid.lower():
        monkeypatch.setattr(mod, "_aiupred_run_local", lambda need, cache, on_progress: set())


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

    def test_saves_cache_periodically_while_fetching(self, tmp_path, monkeypatch):
        monkeypatch.setattr(mod, "_AIUPRED_CACHE_FILE", tmp_path / "aiupred.tsv")
        monkeypatch.setattr(mod, "_AIUPRED_SAVE_EVERY", 2)
        monkeypatch.setattr(mod, "_AIUPRED_MAX_WORKERS", 1)
        monkeypatch.setattr(mod, "fetch_aiupred_all", lambda api_type, uid: {"general": {1: 0.5}, "binding": {1: 0.5}})
        saves = []
        real_save = mod._aiupred_save_cache
        monkeypatch.setattr(mod, "_aiupred_save_cache", lambda c: (saves.append(len(c)), real_save(c)))

        df = pd.DataFrame([{"UniProt": f"P0000{i}"} for i in range(5)])
        mod.run_aiupred_phase(df)

        assert saves == [4, 8, 10], (
            "with 5 proteins (2 cache entries each) and a save every 2 fetches, the cache should be "
            f"written after proteins 2 and 4 and once more for the remainder, got saves of sizes {saves}"
        )

    def test_interrupted_run_keeps_what_was_fetched(self, tmp_path, monkeypatch):
        import threading
        cache_file = tmp_path / "aiupred.tsv"
        monkeypatch.setattr(mod, "_AIUPRED_CACHE_FILE", cache_file)
        monkeypatch.setattr(mod, "_AIUPRED_MAX_WORKERS", 1)
        two_processed = threading.Event()
        processed = []

        def count_progress(*args, **kwargs):
            processed.append(1)
            if len(processed) >= 2:
                two_processed.set()

        def fake_fetch(api_type, uid):
            if uid == "P00003":
                two_processed.wait(timeout=10)  # interrupt only after two results were handled
                raise KeyboardInterrupt
            return {"general": {1: 0.5}, "binding": {1: 0.5}}

        monkeypatch.setattr(mod, "_emit_progress", count_progress)
        monkeypatch.setattr(mod, "fetch_aiupred_all", fake_fetch)

        df = pd.DataFrame([{"UniProt": f"P0000{i}"} for i in range(5)])
        try:
            mod.run_aiupred_phase(df)
        except BaseException:
            pass

        fetched_ids = {uid for uid, _ in mod._aiupred_load_cache()}
        assert len(fetched_ids) >= 2 and "P00003" not in fetched_ids, (
            "proteins handled before the interruption must already be on disk and the interrupted one "
            f"must not be, got {sorted(fetched_ids)}"
        )

class TestAiupredFailureLogging:
    def test_failure_reasons_are_recorded(self, monkeypatch):
        import requests
        mod._aiupred_failures.clear()

        def timeout(*a, **k):
            raise requests.exceptions.ConnectTimeout("connect timed out")

        monkeypatch.setattr(mod._aiupred_session, "get", timeout)
        assert mod.fetch_aiupred_all("binding", "P00001") == {}
        assert mod._aiupred_failures["P00001"][0] == "ConnectTimeout", (
            f"a connection timeout should be recorded by its type, got {mod._aiupred_failures}"
        )

        err = requests.exceptions.HTTPError("429 Too Many Requests")
        err.response = FakeAiupredResponse({}, status_code=429)

        class Resp(FakeAiupredResponse):
            def raise_for_status(self):
                raise err

        monkeypatch.setattr(mod._aiupred_session, "get", lambda *a, **k: Resp({}, 429))
        mod.fetch_aiupred_all("binding", "P00002")
        assert mod._aiupred_failures["P00002"][0] == "HTTP 429", (
            "a rate-limit response must be distinguishable from a server that is unreachable, "
            f"got {mod._aiupred_failures['P00002']}"
        )

        monkeypatch.setattr(mod._aiupred_session, "get", lambda *a, **k: FakeAiupredResponse({"other": 1}))
        mod.fetch_aiupred_all("binding", "P00003")
        assert mod._aiupred_failures["P00003"][0] == "no scores in response"

    def test_successful_fetch_records_no_failure(self, monkeypatch):
        mod._aiupred_failures.clear()
        monkeypatch.setattr(mod._aiupred_session, "get",
                            lambda *a, **k: FakeAiupredResponse({"AIUPred": [0.1, 0.2], "AIUPred-binding": [0.3, 0.4]}))
        assert mod.fetch_aiupred_all("binding", "P00004")["general"] == {1: 0.1, 2: 0.2}
        assert mod._aiupred_failures == {}

    def test_phase_prints_summary_and_writes_log(self, tmp_path, monkeypatch, capsys):
        import requests
        monkeypatch.setattr(mod, "_AIUPRED_CACHE_FILE", tmp_path / "aiupred.tsv")

        def timeout(*a, **k):
            raise requests.exceptions.ConnectTimeout("connect timed out")

        monkeypatch.setattr(mod._aiupred_session, "get", timeout)
        log_dir = tmp_path / "logs"
        mod.run_aiupred_phase(pd.DataFrame([{"UniProt": "P00001"}, {"UniProt": "P00002"}]), log_dir=log_dir)

        out = capsys.readouterr().out
        assert "no scores for 2 of 2" in out and "ConnectTimeout" in out, (
            f"the run output should say how many proteins failed and why, got: {out}"
        )
        log = pd.read_csv(log_dir / "aiupred_failures.tsv", sep="\t")
        assert sorted(log["uniprot_id"]) == ["P00001", "P00002"] and set(log["reason"]) == {"ConnectTimeout"}

    def test_no_failures_means_no_warning_or_log(self, tmp_path, monkeypatch, capsys):
        monkeypatch.setattr(mod, "_AIUPRED_CACHE_FILE", tmp_path / "aiupred.tsv")
        monkeypatch.setattr(mod, "fetch_aiupred_all", lambda api_type, uid: {"general": {1: 0.5}, "binding": {1: 0.5}})
        mod.run_aiupred_phase(pd.DataFrame([{"UniProt": "P00001"}]), log_dir=tmp_path / "logs")
        assert "WARNING" not in capsys.readouterr().out and not (tmp_path / "logs").exists()


def completed(stdout="", returncode=0, stderr=""):
    return SimpleNamespace(stdout=stdout, returncode=returncode, stderr=stderr)


class TestLocalAiupred:
    def _setup(self, tmp_path, monkeypatch, local_scores, local_errors=None, problem=None):
        monkeypatch.setattr(mod, "_AIUPRED_CACHE_FILE", tmp_path / "aiupred.tsv")
        monkeypatch.setattr(mod, "_aiupred_run_local", _real_run_local)
        monkeypatch.setattr(mod, "_aiupred_local_sequences", lambda uids: {u: "MKV" for u in uids if u != "NOSEQ"})
        monkeypatch.setattr(mod, "check_aiupred_helper", lambda: problem)
        monkeypatch.setattr(mod, "_aiupred_local_chunk",
                            lambda seqs: ({u: local_scores for u in seqs if u not in (local_errors or {})},
                                          {u: e for u, e in (local_errors or {}).items() if u in seqs}))
        web_calls = []

        def web(api_type, uid):
            web_calls.append(uid)
            return {"general": {1: 0.9}, "binding": {1: 0.8}}

        monkeypatch.setattr(mod, "fetch_aiupred_all", web)
        return web_calls

    def test_local_scores_are_cached_and_the_web_server_is_not_called(self, tmp_path, monkeypatch):
        scores = {"general": {1: 0.1, 2: 0.2, 3: 0.3}, "binding": {1: 0.4, 2: 0.5, 3: 0.6}}
        web_calls = self._setup(tmp_path, monkeypatch, scores)
        result = mod.run_aiupred_phase(pd.DataFrame([{"UniProt": "P00001"}, {"UniProt": "P00002"}]))
        assert web_calls == [], f"proteins scored locally must not hit the web server, got {web_calls}"
        assert result["binding"]["P00002"] == {1: 0.4, 2: 0.5, 3: 0.6}
        assert ("P00001", "general") in mod._aiupred_load_cache(), "local scores must be saved to the shared cache"

    def test_setup_failure_falls_back_to_the_web_server(self, tmp_path, monkeypatch, capsys):
        web_calls = self._setup(tmp_path, monkeypatch, {}, problem="uv isn't on PATH")
        result = mod.run_aiupred_phase(pd.DataFrame([{"UniProt": "P00001"}]))
        assert web_calls == ["P00001"] and result["general"]["P00001"] == {1: 0.9}
        out = capsys.readouterr().out
        assert "falling back to the AIUPred web server" in out and "uv isn't on PATH" in out, (
            f"the user must be told why local AIUPred wasn't used, got: {out}"
        )

    def test_only_proteins_the_local_run_could_not_score_go_to_the_web(self, tmp_path, monkeypatch, capsys):
        scores = {"general": {1: 0.1, 2: 0.2, 3: 0.3}, "binding": {1: 0.4, 2: 0.5, 3: 0.6}}
        web_calls = self._setup(tmp_path, monkeypatch, scores, local_errors={"P00002": "ValueError: boom"})
        mod.run_aiupred_phase(pd.DataFrame([{"UniProt": u} for u in ("P00001", "P00002", "NOSEQ")]))
        assert sorted(web_calls) == ["NOSEQ", "P00002"], (
            f"a protein with no usable sequence and one whose local run errored both need the web, got {web_calls}"
        )
        assert "ValueError: boom" in capsys.readouterr().out

    def test_chunk_parses_helper_output_and_checks_lengths(self, monkeypatch):
        monkeypatch.setattr(mod, "_uv_script_command", lambda script, offline=False: ["helper"])
        lines = [
            json.dumps({"id": "P1", "general": [0.1, 0.2], "binding": [0.3, 0.4]}),
            json.dumps({"id": "P2", "error": "ValueError: bad"}),
            json.dumps({"id": "P3", "general": [0.1], "binding": [0.3]}),
        ]
        monkeypatch.setattr(subprocess, "run", lambda *a, **k: completed("\n".join(lines)))
        scores, errors = mod._aiupred_local_chunk({"P1": "MK", "P2": "MK", "P3": "MK", "P4": "MK"})
        assert scores == {"P1": {"general": {1: 0.1, 2: 0.2}, "binding": {1: 0.3, 2: 0.4}}}
        assert errors["P2"] == "ValueError: bad"
        assert "wrong length" in errors["P3"], "a score list shorter than the sequence must not be cached"
        assert "no result" in errors["P4"], "a protein the helper never answered for must be reported"

    def test_local_sequences_need_a_complete_chain_without_unknown_residues(self, monkeypatch):
        chains = {"GOOD": {1: "M", 2: "K"}, "GAP": {1: "M", 3: "K"}, "UNK": {1: "M", 2: "X"}}
        monkeypatch.setattr(mod, "MODELS_ROOT", mod.PROJECT_ROOT)  # any existing directory
        monkeypatch.setattr(mod.Path, "is_dir", lambda self: True)
        monkeypatch.setattr(mod, "find_canonical_cif", lambda d: d.name)
        monkeypatch.setattr(mod, "load_first_chain", lambda cif: cif if cif in chains else None)
        monkeypatch.setattr(mod, "extract_sequence", lambda chain: chains[chain])
        assert mod._aiupred_local_sequences(["GOOD", "GAP", "UNK", "NONE"]) == {"GOOD": "MK"}, (
            "only a chain numbered 1..N with standard residues is safe to score locally"
        )


class TestHelperSetupMessage:
    def _run(self, monkeypatch, capsys, offline_ok):
        calls = []

        def fake_run(command, **kwargs):
            calls.append(command)
            if "--offline" in command:
                return completed('{"ok": true}' if offline_ok else "", 0 if offline_ok else 1)
            return completed('{"ok": true}')

        monkeypatch.setattr(subprocess, "run", fake_run)
        monkeypatch.setattr(mod.shutil, "which", lambda name: "uv")
        problem = mod.check_aiupred_helper()
        return problem, calls, capsys.readouterr().out

    def test_first_run_says_dependencies_are_downloading(self, monkeypatch, capsys):
        problem, calls, out = self._run(monkeypatch, capsys, offline_ok=False)
        assert problem is None and len(calls) == 2
        assert "downloading dependencies" in out and "several minutes" in out, (
            f"users must be told the first setup downloads dependencies and can take a while, got: {out!r}"
        )

    def test_already_built_environment_is_silent_and_skips_the_download(self, monkeypatch, capsys):
        problem, calls, out = self._run(monkeypatch, capsys, offline_ok=True)
        assert problem is None and len(calls) == 1 and out == "", (
            f"a cached environment should be verified offline with nothing printed, got {len(calls)} call(s), {out!r}"
        )

    def test_failed_setup_reports_the_last_error_line(self, monkeypatch):
        monkeypatch.setattr(mod.shutil, "which", lambda name: "uv")
        monkeypatch.setattr(subprocess, "run", lambda *a, **k: completed("", 1, "boom\nNo matching distribution"))
        problem = mod.check_aiupred_helper()
        assert problem and "No matching distribution" in problem
