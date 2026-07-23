"""Unit tests for scripts/2_download_structures.py.

Every function gets at least one positive (happy-path) test and at least one
negative/edge-case test, and every assertion states what real-world behavior
is being verified. Network access is always mocked -- no test in this file
makes a real HTTP request.
"""
import pytest


@pytest.fixture
def mod(download_module):
    return download_module


class FakeHTTPResponse:
    """Stand-in for requests.Response covering both plain GET (fetch_prediction,
    which reads .status_code/.json()) and streaming GET (download, which uses
    the response as a context manager and reads .iter_content()).

    conftest.py's FakeResponse only models .text/.headers, which doesn't cover
    either of those shapes, so this is a local extension rather than a reuse.
    """

    def __init__(self, status_code=200, json_data=None, chunks=None):
        self.status_code = status_code
        self._json_data = json_data
        self._chunks = chunks if chunks is not None else [b"chunk-data"]

    def json(self):
        return self._json_data

    def iter_content(self, chunk_size=None):
        yield from self._chunks

    def raise_for_status(self):
        if self.status_code >= 400:
            raise RuntimeError(f"HTTP {self.status_code}")

    def __enter__(self):
        return self

    def __exit__(self, *exc_info):
        return False


class FakeSession:
    """queue: list of FakeHTTPResponse objects returned in order across
    successive .get() calls (models retry sequences); a single
    FakeHTTPResponse is also accepted for calls that never retry."""

    def __init__(self, queue):
        self._queue = list(queue) if isinstance(queue, list) else [queue]
        self.calls = []

    def get(self, url, timeout=None, stream=False):
        self.calls.append(url)
        if len(self._queue) > 1:
            return self._queue.pop(0)
        return self._queue[0]


# ── clean_accession ───────────────────────────────────────────────────────────

class TestCleanAccession:
    def test_returns_a_clean_accession_unchanged(self, mod):
        assert mod.clean_accession("P04637") == "P04637", (
            "a well-formed accession with no surrounding junk should pass through unchanged"
        )

    def test_picks_first_valid_accession_among_delimited_junk(self, mod):
        result = mod.clean_accession("P04637; some note, Q99999")
        assert result == "P04637", (
            f"the first valid-looking accession among delimiter-separated tokens should win, got {result}"
        )

    def test_returns_none_for_cell_with_no_valid_accession(self, mod):
        assert mod.clean_accession("not an accession at all") is None, (
            "a cell containing no token matching the UniProt accession pattern must "
            "return None, not a garbage guess"
        )

    def test_returns_none_for_none_input(self, mod):
        assert mod.clean_accession(None) is None, (
            "a None cell (e.g. a blank spreadsheet row) must return None, not raise"
        )

    def test_returns_none_for_empty_string(self, mod):
        assert mod.clean_accession("   ") is None, (
            "a whitespace-only cell must return None after stripping, not match spuriously"
        )


# ── pick_urls ──────────────────────────────────────────────────────────────────

class TestPickUrls:
    def test_prefers_cif_when_prefer_is_cif(self, mod):
        record = {
            "cifUrl": "https://example.org/AF-P1-F1-model_v4.cif",
            "pdbUrl": "https://example.org/AF-P1-F1-model_v4.pdb",
        }
        result = mod.pick_urls(record, prefer="cif")
        assert result["structure_url"].endswith(".cif"), (
            f"with prefer='cif', the .cif URL should win over .pdb, got {result['structure_url']}"
        )

    def test_prefers_pdb_when_prefer_is_pdb(self, mod):
        record = {
            "cifUrl": "https://example.org/AF-P1-F1-model_v4.cif",
            "pdbUrl": "https://example.org/AF-P1-F1-model_v4.pdb",
        }
        result = mod.pick_urls(record, prefer="pdb")
        assert result["structure_url"].endswith(".pdb"), (
            f"with prefer='pdb', the .pdb URL should win over .cif, got {result['structure_url']}"
        )

    def test_gzipped_variant_scores_higher_than_plain(self, mod):
        record = {
            "cifUrlPlain": "https://example.org/AF-P1-F1-model_v4.cif",
            "cifUrlGz": "https://example.org/AF-P1-F1-model_v4.cif.gz",
        }
        result = mod.pick_urls(record, prefer="cif")
        assert result["structure_url"].endswith(".cif.gz"), (
            "the .gz scoring bonus should make the gzipped URL win over the plain one when "
            f"both are same-format candidates, got {result['structure_url']}"
        )

    def test_empty_structure_url_when_no_matching_format_present(self, mod):
        record = {"someOtherUrl": "https://example.org/AF-P1-F1-notes.txt"}
        result = mod.pick_urls(record, prefer="cif")
        assert result["structure_url"] == "", (
            "if the record has no cif/bcif/pdb-shaped URL at all, structure_url must be "
            f"empty rather than an incorrect guess, got {result['structure_url']!r}"
        )

    def test_pae_doc_url_field_takes_precedence_over_scoring(self, mod):
        record = {
            "paeDocUrl": "https://example.org/explicit-pae.json",
            "paeUrl": "https://example.org/AF-P1-F1-predicted_aligned_error_v4.json",
        }
        result = mod.pick_urls(record)
        assert result["pae_url"] == "https://example.org/explicit-pae.json", (
            "an explicit paeDocUrl field should be used directly rather than falling back "
            "to filename-pattern scoring"
        )

    def test_non_string_and_non_url_values_are_ignored(self, mod):
        record = {"cifUrl": "https://example.org/AF-P1-F1-model_v4.cif",
                   "someNumber": 42, "someNone": None}
        result = mod.pick_urls(record, prefer="cif")
        assert result["structure_url"].endswith(".cif"), (
            "non-string / non-URL record values must be filtered out without raising, "
            f"got {result}"
        )


# ── fetch_prediction ───────────────────────────────────────────────────────────

class TestFetchPrediction:
    def test_returns_parsed_json_on_200(self, mod):
        session = FakeSession(FakeHTTPResponse(status_code=200, json_data={"a": 1}))
        result = mod.fetch_prediction("P04637", session)
        assert result == {"a": 1}, "a 200 response's JSON body should be returned as-is"

    def test_returns_none_on_404(self, mod):
        session = FakeSession(FakeHTTPResponse(status_code=404))
        result = mod.fetch_prediction("P99999", session)
        assert result is None, (
            "a 404 means AlphaFold DB has no record for this accession -- must return None, "
            "not raise, so the caller can distinguish 'no entry' from a real error"
        )

    def test_raises_after_exhausting_retries_on_persistent_error(self, mod, monkeypatch):
        monkeypatch.setattr(mod.time, "sleep", lambda *_a, **_k: None)  # skip real backoff delays
        session = FakeSession(FakeHTTPResponse(status_code=500))
        with pytest.raises(RuntimeError):
            mod.fetch_prediction("P04637", session, retries=3)

    def test_retries_on_transient_error_then_succeeds(self, mod, monkeypatch):
        monkeypatch.setattr(mod.time, "sleep", lambda *_a, **_k: None)
        session = FakeSession([
            FakeHTTPResponse(status_code=500),
            FakeHTTPResponse(status_code=200, json_data={"ok": True}),
        ])
        result = mod.fetch_prediction("P04637", session, retries=3)
        assert result == {"ok": True}, (
            "a transient error followed by a successful response should retry and return "
            "the eventual success, not give up after the first failure"
        )


# ── download ───────────────────────────────────────────────────────────────────

class TestDownload:
    def test_writes_file_via_atomic_rename(self, mod, tmp_path):
        outpath = tmp_path / "sub" / "model.cif"
        session = FakeSession(FakeHTTPResponse(status_code=200, chunks=[b"hello ", b"world"]))

        mod.download("https://example.org/model.cif", outpath, session)

        assert outpath.exists(), "download must create the output file (and its parent dir)"
        assert outpath.read_bytes() == b"hello world", (
            "all streamed chunks must be concatenated and written to the final path"
        )
        assert not outpath.with_suffix(outpath.suffix + ".part").exists(), (
            "the .part temp file must be renamed away, not left behind after a successful download"
        )

    def test_skips_download_when_already_cached_nonempty(self, mod, tmp_path):
        outpath = tmp_path / "model.cif"
        outpath.write_bytes(b"already here")
        session = FakeSession(FakeHTTPResponse(status_code=200))

        mod.download("https://example.org/model.cif", outpath, session)

        assert session.calls == [], (
            "an already-downloaded non-empty file should be skipped entirely -- no network "
            f"call should be made, but got {len(session.calls)} call(s)"
        )
        assert outpath.read_bytes() == b"already here", (
            "the existing file's content must be left untouched when the download is skipped"
        )

    def test_raises_after_exhausting_retries(self, mod, tmp_path, monkeypatch):
        monkeypatch.setattr(mod.time, "sleep", lambda *_a, **_k: None)
        outpath = tmp_path / "model.cif"
        session = FakeSession(FakeHTTPResponse(status_code=500))

        with pytest.raises(RuntimeError):
            mod.download("https://example.org/model.cif", outpath, session, retries=2)
        assert not outpath.exists(), (
            "a download that never succeeds must not leave a partial/empty output file behind"
        )


# ── _cached_version / _remove_stale_version ────────────────────────────────────

class TestCachedVersion:
    def test_finds_highest_version_across_fragments(self, mod, tmp_path):
        acc_dir = tmp_path / "P04637"
        acc_dir.mkdir()
        (acc_dir / "AF-P04637-F1-model_v4.cif").write_text("x")
        (acc_dir / "AF-P04637-F2-model_v6.cif").write_text("x")

        result = mod._cached_version(acc_dir, "P04637", "model")
        assert result == 6, (
            f"multi-fragment proteins share one version -- the max across fragments (6) "
            f"should be returned, got {result}"
        )

    def test_returns_none_when_directory_does_not_exist(self, mod, tmp_path):
        assert mod._cached_version(tmp_path / "nope", "P04637", "model") is None, (
            "a protein with no download directory at all has nothing cached -- must "
            "return None, not raise"
        )

    def test_returns_none_when_directory_empty(self, mod, tmp_path):
        acc_dir = tmp_path / "P99999"
        acc_dir.mkdir()
        assert mod._cached_version(acc_dir, "P99999", "model") is None, (
            "an empty directory has no cached version to report"
        )


class TestRemoveStaleVersion:
    def test_deletes_only_files_not_matching_keep_version(self, mod, tmp_path):
        acc_dir = tmp_path / "P04637"
        acc_dir.mkdir()
        old = acc_dir / "AF-P04637-F1-model_v4.cif"
        new = acc_dir / "AF-P04637-F1-model_v6.cif"
        old.write_text("old")
        new.write_text("new")

        mod._remove_stale_version(acc_dir, "P04637", keep_version=6)

        assert not old.exists(), "the stale (v4) file must be deleted after an update to v6"
        assert new.exists(), "the current (v6) file must NOT be deleted"

    def test_no_op_when_no_files_match_accession(self, mod, tmp_path):
        acc_dir = tmp_path / "P04637"
        acc_dir.mkdir()
        unrelated = acc_dir / "AF-Q99999-F1-model_v4.cif"
        unrelated.write_text("x")

        mod._remove_stale_version(acc_dir, "P04637", keep_version=6)

        assert unrelated.exists(), (
            "a file belonging to a different accession must never be touched, even though "
            "it lives in the same directory"
        )


# ── read_table ─────────────────────────────────────────────────────────────────

class TestReadTable:
    def test_reads_tsv_by_extension(self, mod, tmp_path):
        path = tmp_path / "data.tsv"
        path.write_text("UniProt\tGene\nP04637\tTP53\n")
        df = mod.read_table(str(path))
        assert list(df.columns) == ["UniProt", "Gene"], (
            f"a .tsv file should be parsed tab-delimited, got columns {list(df.columns)}"
        )

    def test_reads_csv_by_extension(self, mod, tmp_path):
        path = tmp_path / "data.csv"
        path.write_text("UniProt,Gene\nP04637,TP53\n")
        df = mod.read_table(str(path))
        assert list(df.columns) == ["UniProt", "Gene"], (
            f"a .csv file should be parsed comma-delimited, got columns {list(df.columns)}"
        )

    def test_reads_xlsx_by_extension(self, mod, tmp_path):
        pd = pytest.importorskip("pandas")
        path = tmp_path / "data.xlsx"
        pd.DataFrame({"UniProt": ["P04637"], "Gene": ["TP53"]}).to_excel(path, index=False)
        df = mod.read_table(str(path))
        assert list(df.columns) == ["UniProt", "Gene"], (
            f"a .xlsx file should be parsed as a spreadsheet, got columns {list(df.columns)}"
        )

    def test_unrecognized_extension_falls_back_to_tsv_parsing(self, mod, tmp_path):
        path = tmp_path / "data.txt"
        path.write_text("UniProt\tGene\nP04637\tTP53\n")
        df = mod.read_table(str(path))
        assert list(df.columns) == ["UniProt", "Gene"], (
            "an unrecognized extension (.txt) should fall back to tab-delimited parsing "
            f"per read_table's documented default, got columns {list(df.columns)}"
        )


# ── _check_and_download ────────────────────────────────────────────────────────

class TestCheckAndDownload:
    """Mocks fetch_prediction/download at the module level (both are already
    covered by their own tests above) so each branch of this orchestration
    function can be exercised in isolation."""

    def test_up_to_date_is_reported_cached_without_downloading(self, mod, tmp_path, monkeypatch):
        acc_dir = tmp_path / "P04637"
        acc_dir.mkdir()
        (acc_dir / "AF-P04637-F1-model_v6.cif").write_text("cached")

        monkeypatch.setattr(mod, "fetch_prediction", lambda acc, session: {
            "uniprotAccession": acc, "cifUrl": "https://example.org/AF-P04637-F1-model_v6.cif",
        })
        download_calls = []
        monkeypatch.setattr(mod, "download", lambda *a, **k: download_calls.append(a))

        row = mod._check_and_download("P04637", tmp_path, "cif", False, 0, session=None)

        assert row["status"] == "ALREADY_CACHED", (
            f"a cached version matching AlphaFold DB's current version must not trigger a "
            f"re-download, got status={row['status']!r}"
        )
        assert download_calls == [], "no download() call should happen when already up to date"

    def test_missing_local_cache_triggers_download(self, mod, tmp_path, monkeypatch):
        monkeypatch.setattr(mod, "fetch_prediction", lambda acc, session: {
            "uniprotAccession": acc, "cifUrl": "https://example.org/AF-P04637-F1-model_v6.cif",
        })
        download_calls = []
        monkeypatch.setattr(mod, "download", lambda *a, **k: download_calls.append(a))

        row = mod._check_and_download("P04637", tmp_path, "cif", False, 0, session=None)

        assert row["status"] == "DOWNLOADED", (
            f"an accession with nothing cached locally should be freshly downloaded, "
            f"got status={row['status']!r}"
        )
        assert len(download_calls) == 1, "exactly one file (structure only, also_pae=False) should be downloaded"

    def test_newer_version_available_triggers_update_and_removes_stale(self, mod, tmp_path, monkeypatch):
        acc_dir = tmp_path / "P04637"
        acc_dir.mkdir()
        stale = acc_dir / "AF-P04637-F1-model_v4.cif"
        stale.write_text("old")

        monkeypatch.setattr(mod, "fetch_prediction", lambda acc, session: {
            "uniprotAccession": acc, "cifUrl": "https://example.org/AF-P04637-F1-model_v6.cif",
        })
        monkeypatch.setattr(mod, "download", lambda url, path, session: path.write_text("new"))

        row = mod._check_and_download("P04637", tmp_path, "cif", False, 0, session=None)

        assert row["status"] == "UPDATED", (
            f"a newer AlphaFold version than what's cached should be reported as UPDATED, "
            f"got status={row['status']!r}"
        )
        assert not stale.exists(), "the old v4 file must be removed after updating to v6"

    def test_404_with_existing_cache_keeps_cached_copy(self, mod, tmp_path, monkeypatch):
        acc_dir = tmp_path / "P04637"
        acc_dir.mkdir()
        (acc_dir / "AF-P04637-F1-model_v6.cif").write_text("cached")

        monkeypatch.setattr(mod, "fetch_prediction", lambda acc, session: None)  # simulates 404

        row = mod._check_and_download("P04637", tmp_path, "cif", False, 0, session=None)

        assert row["status"] == "ALREADY_CACHED", (
            "if AlphaFold DB now 404s for an accession we already have cached (e.g. the "
            f"entry was removed upstream), the local copy must be kept, got {row['status']!r}"
        )

    def test_404_with_no_cache_reports_no_entry(self, mod, tmp_path, monkeypatch):
        monkeypatch.setattr(mod, "fetch_prediction", lambda acc, session: None)

        row = mod._check_and_download("Q99999", tmp_path, "cif", False, 0, session=None)

        assert row["status"] == "NO_ENTRY", (
            f"a 404 with nothing cached locally should be reported as NO_ENTRY, got {row['status']!r}"
        )

    def test_isoform_only_records_report_no_canonical_model(self, mod, tmp_path, monkeypatch):
        monkeypatch.setattr(mod, "fetch_prediction", lambda acc, session: [
            {"uniprotAccession": f"{acc}-9", "cifUrl": "https://example.org/isoform.cif"},
        ])

        row = mod._check_and_download("P11362", tmp_path, "cif", False, 0, session=None)

        assert row["status"] == "NO_CANONICAL_MODEL", (
            "records whose uniprotAccession carries an isoform suffix (e.g. 'P11362-9') "
            f"must not be mistaken for a canonical model, got {row['status']!r}"
        )

    def test_exception_during_processing_is_caught_as_error_row(self, mod, tmp_path, monkeypatch):
        def _raise(*_a, **_k):
            raise RuntimeError("simulated network failure")
        monkeypatch.setattr(mod, "fetch_prediction", _raise)

        row = mod._check_and_download("P04637", tmp_path, "cif", False, 0, session=None)

        assert row["status"] == "ERROR", (
            "an unexpected exception anywhere in this per-accession worker must be caught "
            "into an ERROR report row, not propagate and kill the whole thread pool"
        )
        assert "simulated network failure" in row["note"], (
            f"the error row's note should include the actual exception message for "
            f"debugging, got {row['note']!r}"
        )


# ── main ───────────────────────────────────────────────────────────────────────

class TestMain:
    def test_raises_when_id_column_missing(self, mod, tmp_path):
        input_path = tmp_path / "input.tsv"
        input_path.write_text("SomeOtherColumn\nvalue\n")

        with pytest.raises(ValueError):
            mod.main(
                str(input_path), id_column="UniProt", out_dir=str(tmp_path / "out"),
                prefer="cif", also_pae=False, delay=0, logs_dir=tmp_path / "logs",
            )

    def test_writes_report_for_valid_input(self, mod, tmp_path, monkeypatch):
        input_path = tmp_path / "input.tsv"
        input_path.write_text("UniProt\nP04637\nQ99999\n")

        monkeypatch.setattr(mod, "_check_and_download", lambda acc, *a, **k: {
            "UniProt": acc, "status": "DOWNLOADED", "structure_file": "x", "pae_file": "", "note": "",
        })

        out_dir = tmp_path / "out"
        logs_dir = tmp_path / "logs"
        mod.main(
            str(input_path), id_column="UniProt", out_dir=str(out_dir),
            prefer="cif", also_pae=False, delay=0, logs_dir=logs_dir,
        )

        report_path = logs_dir / "download_report.tsv"
        assert report_path.exists(), (
            "main() must write download_report.tsv summarizing every accession processed"
        )
        report_text = report_path.read_text()
        assert "P04637" in report_text and "Q99999" in report_text, (
            f"the report must include every accession from the input file, got:\n{report_text}"
        )
