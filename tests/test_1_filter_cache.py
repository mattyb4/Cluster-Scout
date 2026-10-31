"""Unit tests for the local UniProt API caches in scripts/1_filter.py.

Each fetch/compute function should:
 - serve previously-seen keys from the on-disk cache without hitting the API
 - persist newly-fetched results (including "not found" / "no restriction"
   markers) so they aren't re-queried on the next run
"""
import pytest
from conftest import FakeResponse


@pytest.fixture
def mod(filter_module, tmp_path, monkeypatch):
    # Point the cache at a per-test temp directory so tests don't touch (or
    # depend on) the real data/cache/ contents.
    monkeypatch.setattr(filter_module, "CACHE_DIR", tmp_path)
    return filter_module


class TestLoadSaveCache:
    def test_round_trip(self, mod):
        cache = {"P12345": ("GENEA",), "Q99999": ("",)}
        mod._save_cache("test.tsv", cache, ["UniProt", "gene"])
        result = mod._load_cache("test.tsv", ["UniProt", "gene"])
        assert result == cache, (
            f"loading a cache immediately after saving it should reproduce the exact same "
            f"dict, including entries with an empty-string value, got {result}"
        )

    def test_missing_file_returns_empty_dict(self, mod):
        result = mod._load_cache("does_not_exist.tsv", ["UniProt", "gene"])
        assert result == {}, (
            f"a cache file that has never been written should load as an empty dict, not "
            f"raise FileNotFoundError, got {result}"
        )


class TestFetchUniprotGeneMapping:
    def test_fetches_and_caches(self, mod, monkeypatch):
        calls = []

        def fake_get(url, params=None, **kwargs):
            calls.append(params)
            return FakeResponse(
                "Entry\tGene Names\nP04637\tTP53 p53\n",
                headers={"X-UniProt-Release": "2026_02"},
            )

        monkeypatch.setattr(mod.requests, "get", fake_get)

        df = mod.fetch_uniprot_gene_mapping(["P04637"])
        assert df.to_dict("records") == [{"UniProt": "P04637", "gene": "TP53"}], (
            "the primary gene name (first token of 'Gene Names') should be extracted, "
            f"got {df.to_dict('records')}"
        )
        assert len(calls) == 1, f"exactly one API call should be made for one uncached accession, got {len(calls)}"

    def test_second_call_is_served_from_cache(self, mod, monkeypatch):
        calls = []

        def fake_get(url, params=None, **kwargs):
            calls.append(params)
            return FakeResponse("Entry\tGene Names\nP04637\tTP53 p53\n")

        monkeypatch.setattr(mod.requests, "get", fake_get)

        mod.fetch_uniprot_gene_mapping(["P04637"])
        df2 = mod.fetch_uniprot_gene_mapping(["P04637"])

        assert df2.to_dict("records") == [{"UniProt": "P04637", "gene": "TP53"}], (
            "the second call should return the same mapping even though it's served from cache"
        )
        assert len(calls) == 1, (
            f"an accession already cached from the first call must not trigger a second "
            f"HTTP request, got {len(calls)} calls"
        )

    def test_strips_variant_suffix(self, mod, monkeypatch):
        monkeypatch.setattr(
            mod.requests, "get",
            lambda url, params=None, **kwargs: FakeResponse("Entry\tGene Names\nQ16613\tBAD1\n"),
        )

        df = mod.fetch_uniprot_gene_mapping(["Q16613_VAR_A129T"])
        assert df.to_dict("records") == [{"UniProt": "Q16613", "gene": "BAD1"}], (
            "a COSMIC-style '_VAR_...' variant suffix must be stripped before the accession "
            f"is used/reported, got {df.to_dict('records')}"
        )

    def test_not_found_accession_is_cached_and_excluded(self, mod, monkeypatch):
        monkeypatch.setattr(
            mod.requests, "get",
            lambda url, params=None, **kwargs: FakeResponse("Entry\tGene Names\n"),
        )

        df = mod.fetch_uniprot_gene_mapping(["P99999"])
        assert df.empty, (
            "an accession UniProt's API doesn't recognize should not appear in the "
            f"returned mapping at all, got {df.to_dict('records')}"
        )

        cache = mod._load_cache(mod.UNIPROT_GENE_CACHE_FILE, ["UniProt", "gene"])
        assert cache["P99999"] == ("",), (
            "the 'not found' result must still be cached (as an empty-string gene) so the "
            "next run doesn't re-query the same dead accession"
        )

    def test_empty_input_returns_empty_without_request(self, mod, monkeypatch):
        called = []
        monkeypatch.setattr(mod.requests, "get", lambda *a, **k: called.append(1))

        df = mod.fetch_uniprot_gene_mapping([])
        assert df.empty, "an empty accession list has nothing to map -- result should be an empty DataFrame"
        assert called == [], "no HTTP request should be made when there's nothing to look up"


class TestFetchGeneToUniprotMapping:
    def test_fetches_and_caches(self, mod, monkeypatch):
        calls = []

        def fake_get(url, params=None, **kwargs):
            calls.append(params)
            return FakeResponse("Entry\tGene Names\nQ06124\tPTPN11\n")

        monkeypatch.setattr(mod.requests, "get", fake_get)

        df = mod.fetch_gene_to_uniprot_mapping(["PTPN11"])
        assert df.to_dict("records") == [{"gene": "PTPN11", "UniProt": "Q06124"}], (
            f"gene->UniProt lookup should return the reviewed accession, got {df.to_dict('records')}"
        )
        assert len(calls) == 1, "the first lookup for a gene should make exactly one API call"

        df2 = mod.fetch_gene_to_uniprot_mapping(["PTPN11"])
        assert df2.to_dict("records") == [{"gene": "PTPN11", "UniProt": "Q06124"}], (
            "a second lookup of the same gene should return the same result from cache"
        )
        assert len(calls) == 1, (
            f"the second lookup must be served from cache with no new HTTP call, got {len(calls)} calls total"
        )

    def test_unmapped_gene_is_cached_and_excluded(self, mod, monkeypatch):
        monkeypatch.setattr(
            mod.requests, "get",
            lambda url, params=None, **kwargs: FakeResponse("Entry\tGene Names\n"),
        )

        df = mod.fetch_gene_to_uniprot_mapping(["NOTAGENE"])
        assert df.empty, (
            f"a gene symbol with no reviewed human UniProt match should not appear in the "
            f"result, got {df.to_dict('records')}"
        )

        cache = mod._load_cache(mod.GENE_TO_UNIPROT_CACHE_FILE, ["gene", "UniProt"])
        assert cache["NOTAGENE"] == ("",), (
            "the unmapped result must be cached (empty-string UniProt) to avoid re-querying "
            "the same unmapped gene on every future run"
        )


class TestUniprotOutage:
    """UniProt's search endpoint goes down for stretches: retry, then fail briefly."""

    @pytest.fixture(autouse=True)
    def no_sleep(self, mod, monkeypatch):
        monkeypatch.setattr(mod.time, "sleep", lambda s: None)

    def test_temporary_failure_is_retried(self, mod, monkeypatch):
        responses = [FakeResponse("", status_code=503),
                     FakeResponse("Entry\tGene Names\nP04637\tTP53 p53\n")]
        monkeypatch.setattr(mod.requests, "get", lambda url, params=None, **kw: responses.pop(0))

        df = mod.fetch_uniprot_gene_mapping(["P04637"])
        assert df.to_dict("records") == [{"UniProt": "P04637", "gene": "TP53"}], (
            f"a 503 followed by a good response should succeed on retry, got {df.to_dict('records')}"
        )

    def test_persistent_outage_raises_short_message_without_the_url(self, mod, monkeypatch):
        monkeypatch.setattr(mod.requests, "get", lambda url, params=None, **kw: FakeResponse("", status_code=503))

        with pytest.raises(mod.UniProtUnavailableError) as exc:
            mod.fetch_uniprot_gene_mapping(["P04637"])
        message = str(exc.value)
        assert "HTTP 503" in message and "try again" in message, (
            f"the message should name the cause and say to retry, got {message!r}"
        )
        assert "accession:" not in message and len(message) < 300, (
            f"the (very long) query URL must not appear in the message, got {message!r}"
        )

    def test_timeout_is_reported_as_such(self, mod, monkeypatch):
        def fake_get(url, params=None, **kw):
            raise mod.requests.Timeout()

        monkeypatch.setattr(mod.requests, "get", fake_get)
        with pytest.raises(mod.UniProtUnavailableError, match="timed out"):
            mod.fetch_gene_to_uniprot_mapping(["PTPN11"])

    def test_client_error_is_not_retried(self, mod, monkeypatch):
        calls = []

        def fake_get(url, params=None, **kw):
            calls.append(1)
            return FakeResponse("", status_code=400)

        monkeypatch.setattr(mod.requests, "get", fake_get)
        with pytest.raises(mod.UniProtUnavailableError, match="HTTP 400"):
            mod.fetch_uniprot_gene_mapping(["P04637"])
        assert len(calls) == 1, f"a 400 won't fix itself, so it shouldn't be retried; got {len(calls)} calls"

    def test_batches_fetched_before_the_outage_are_saved(self, mod, monkeypatch):
        ids = [f"P{n:05d}" for n in range(150)]  # two batches of 100
        calls = []

        def fake_get(url, params=None, **kw):
            calls.append(1)
            if len(calls) > 1:
                return FakeResponse("", status_code=503)
            batch = [q.split(":")[1] for q in params["query"].split(" OR ")]
            body = "".join(f"{u}\tG{u}\n" for u in batch)
            return FakeResponse("Entry\tGene Names\n" + body)

        monkeypatch.setattr(mod.requests, "get", fake_get)
        with pytest.raises(mod.UniProtUnavailableError):
            mod.fetch_uniprot_gene_mapping(ids)

        cache = mod._load_cache(mod.UNIPROT_GENE_CACHE_FILE, ["UniProt", "gene"])
        assert len(cache) == 100, (
            f"the first batch finished before the outage and should be on disk, found {len(cache)} entries"
        )
