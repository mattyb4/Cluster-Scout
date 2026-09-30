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

        def fake_get(url, params=None):
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

        def fake_get(url, params=None):
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
            lambda url, params=None: FakeResponse("Entry\tGene Names\nQ16613\tBAD1\n"),
        )

        df = mod.fetch_uniprot_gene_mapping(["Q16613_VAR_A129T"])
        assert df.to_dict("records") == [{"UniProt": "Q16613", "gene": "BAD1"}], (
            "a COSMIC-style '_VAR_...' variant suffix must be stripped before the accession "
            f"is used/reported, got {df.to_dict('records')}"
        )

    def test_not_found_accession_is_cached_and_excluded(self, mod, monkeypatch):
        monkeypatch.setattr(
            mod.requests, "get",
            lambda url, params=None: FakeResponse("Entry\tGene Names\n"),
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

        def fake_get(url, params=None):
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
            lambda url, params=None: FakeResponse("Entry\tGene Names\n"),
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
