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


class TestFetchUniprotSequence:
    def test_returns_sequence_string_on_200(self, mod, monkeypatch):
        monkeypatch.setattr(
            mod.requests, "get",
            lambda url: FakeResponse(">sp|P04637|P53_HUMAN Cellular tumor antigen p53\nMEEPQSDPSV\nCNTSSPQP\n"),
        )

        result = mod._fetch_uniprot_sequence("P04637")
        assert result == "MEEPQSDPSVCNTSSPQP", (
            "the FASTA header line (starting with '>') must be dropped and the remaining "
            f"sequence lines joined into one string with no newlines, got {result!r}"
        )

    def test_returns_none_on_non_200(self, mod, monkeypatch):
        monkeypatch.setattr(
            mod.requests, "get",
            lambda url: FakeResponse("Not Found", status_code=404),
        )

        result = mod._fetch_uniprot_sequence("DELETED1")
        assert result is None, (
            "a non-200 response (e.g. a withdrawn/merged UniProt entry) must return None "
            f"rather than an empty or garbage sequence string, got {result!r}"
        )


class TestComputeIsoformSafeLengths:
    def test_no_mismatch_is_cached_as_no_restriction(self, mod, monkeypatch):
        # The transcript's xref accession matches the canonical accession,
        # so there's no isoform mismatch.
        monkeypatch.setattr(
            mod.requests, "get",
            lambda url, params=None: FakeResponse(
                "Entry\txref_Ensembl\tSequence\nP04637\tENST00000269305.9 [P04637]\tMEEPQSDPSV\n"
            ),
        )

        gene_to_transcript = {"TP53": "ENST00000269305.9"}
        gene_to_uniprot = {"TP53": "P04637"}

        df = mod.compute_isoform_safe_lengths(gene_to_transcript, gene_to_uniprot)
        assert df.empty, (
            "when the transcript's isoform accession equals the canonical accession, there "
            f"is no restriction to report -- expected an empty result, got {df.to_dict('records')}"
        )

        cache = mod._load_cache(
            mod.ISOFORM_SAFE_LENGTH_CACHE_FILE, ["gene", "transcript_accession", "isoform_safe_length"]
        )
        assert cache["TP53"] == ("ENST00000269305.9", ""), (
            "the 'no restriction' result must still be cached (empty-string length) so this "
            "gene isn't rechecked every run"
        )

    def test_cached_gene_with_unchanged_transcript_skips_refetch(self, mod, monkeypatch):
        monkeypatch.setattr(
            mod.requests, "get",
            lambda url, params=None: FakeResponse(
                "Entry\txref_Ensembl\tSequence\nP04637\tENST00000269305.9 [P04637]\tMEEPQSDPSV\n"
            ),
        )
        gene_to_transcript = {"TP53": "ENST00000269305.9"}
        gene_to_uniprot = {"TP53": "P04637"}
        mod.compute_isoform_safe_lengths(gene_to_transcript, gene_to_uniprot)

        calls = []
        monkeypatch.setattr(mod.requests, "get", lambda *a, **k: calls.append(1) or FakeResponse(""))

        df = mod.compute_isoform_safe_lengths(gene_to_transcript, gene_to_uniprot)
        assert df.empty, "the second call should still report no restriction for this gene"
        assert calls == [], (
            f"a gene cached under the SAME transcript accession must not be re-fetched, "
            f"got {len(calls)} unexpected call(s)"
        )

    def test_changed_transcript_triggers_recheck(self, mod, monkeypatch):
        monkeypatch.setattr(
            mod.requests, "get",
            lambda url, params=None: FakeResponse(
                "Entry\txref_Ensembl\tSequence\nP04637\tENST00000269305.9 [P04637]\tMEEPQSDPSV\n"
            ),
        )
        gene_to_uniprot = {"TP53": "P04637"}
        mod.compute_isoform_safe_lengths({"TP53": "ENST00000269305.9"}, gene_to_uniprot)

        calls = []

        def fake_get(url, params=None):
            calls.append(params)
            return FakeResponse(
                "Entry\txref_Ensembl\tSequence\nP04637\tENST00000999999.1 [P04637]\tMEEPQSDPSV\n"
            )

        monkeypatch.setattr(mod.requests, "get", fake_get)

        # Transcript accession changed -> gene should be re-checked.
        mod.compute_isoform_safe_lengths({"TP53": "ENST00000999999.1"}, gene_to_uniprot)
        assert len(calls) == 1, (
            "if COSMIC's transcript accession for this gene changes since the last run, the "
            f"cache entry is stale and must trigger exactly one re-fetch, got {len(calls)} calls"
        )

    def test_isoform_mismatch_computes_safe_length(self, mod, monkeypatch):
        def fake_get(url, params=None):
            if url.endswith(".fasta"):
                accession = url.rsplit("/", 1)[-1].removesuffix(".fasta")
                seqs = {
                    "P04637": "MEEPQSDPSV",
                    "P04637-2": "MEEPQSDPS",
                }
                return FakeResponse(f">sp|{accession}|TEST\n{seqs[accession]}\n")
            return FakeResponse(
                "Entry\txref_Ensembl\tSequence\nP04637-2\tENST00000269305.9 [P04637-2]\tMEEPQSDPS\n"
            )

        monkeypatch.setattr(mod.requests, "get", fake_get)

        gene_to_transcript = {"TP53": "ENST00000269305.9"}
        gene_to_uniprot = {"TP53": "P04637"}

        df = mod.compute_isoform_safe_lengths(gene_to_transcript, gene_to_uniprot)
        assert df.to_dict("records") == [{"gene": "TP53", "isoform_safe_length": 9}], (
            "when the transcript maps to an isoform accession different from canonical, the "
            "safe length should be the longest common prefix between canonical ('MEEPQSDPSV') "
            f"and isoform ('MEEPQSDPS') sequences -- 9 matching leading residues -- got {df.to_dict('records')}"
        )
