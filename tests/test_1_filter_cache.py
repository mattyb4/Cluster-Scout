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
        assert mod._load_cache("test.tsv", ["UniProt", "gene"]) == cache

    def test_missing_file_returns_empty_dict(self, mod):
        assert mod._load_cache("does_not_exist.tsv", ["UniProt", "gene"]) == {}


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
        assert df.to_dict("records") == [{"UniProt": "P04637", "gene": "TP53"}]
        assert len(calls) == 1

    def test_second_call_is_served_from_cache(self, mod, monkeypatch):
        calls = []

        def fake_get(url, params=None):
            calls.append(params)
            return FakeResponse("Entry\tGene Names\nP04637\tTP53 p53\n")

        monkeypatch.setattr(mod.requests, "get", fake_get)

        mod.fetch_uniprot_gene_mapping(["P04637"])
        df2 = mod.fetch_uniprot_gene_mapping(["P04637"])

        assert df2.to_dict("records") == [{"UniProt": "P04637", "gene": "TP53"}]
        assert len(calls) == 1  # no new HTTP request on the second call

    def test_strips_variant_suffix(self, mod, monkeypatch):
        monkeypatch.setattr(
            mod.requests, "get",
            lambda url, params=None: FakeResponse("Entry\tGene Names\nQ16613\tBAD1\n"),
        )

        df = mod.fetch_uniprot_gene_mapping(["Q16613_VAR_A129T"])
        assert df.to_dict("records") == [{"UniProt": "Q16613", "gene": "BAD1"}]

    def test_not_found_accession_is_cached_and_excluded(self, mod, monkeypatch):
        monkeypatch.setattr(
            mod.requests, "get",
            lambda url, params=None: FakeResponse("Entry\tGene Names\n"),
        )

        df = mod.fetch_uniprot_gene_mapping(["P99999"])
        assert df.empty

        cache = mod._load_cache(mod.UNIPROT_GENE_CACHE_FILE, ["UniProt", "gene"])
        assert cache["P99999"] == ("",)

    def test_empty_input_returns_empty_without_request(self, mod, monkeypatch):
        called = []
        monkeypatch.setattr(mod.requests, "get", lambda *a, **k: called.append(1))

        df = mod.fetch_uniprot_gene_mapping([])
        assert df.empty
        assert called == []


class TestFetchGeneToUniprotMapping:
    def test_fetches_and_caches(self, mod, monkeypatch):
        calls = []

        def fake_get(url, params=None):
            calls.append(params)
            return FakeResponse("Entry\tGene Names\nQ06124\tPTPN11\n")

        monkeypatch.setattr(mod.requests, "get", fake_get)

        df = mod.fetch_gene_to_uniprot_mapping(["PTPN11"])
        assert df.to_dict("records") == [{"gene": "PTPN11", "UniProt": "Q06124"}]
        assert len(calls) == 1

        df2 = mod.fetch_gene_to_uniprot_mapping(["PTPN11"])
        assert df2.to_dict("records") == [{"gene": "PTPN11", "UniProt": "Q06124"}]
        assert len(calls) == 1  # served from cache

    def test_unmapped_gene_is_cached_and_excluded(self, mod, monkeypatch):
        monkeypatch.setattr(
            mod.requests, "get",
            lambda url, params=None: FakeResponse("Entry\tGene Names\n"),
        )

        df = mod.fetch_gene_to_uniprot_mapping(["NOTAGENE"])
        assert df.empty

        cache = mod._load_cache(mod.GENE_TO_UNIPROT_CACHE_FILE, ["gene", "UniProt"])
        assert cache["NOTAGENE"] == ("",)


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
        assert df.empty

        cache = mod._load_cache(
            mod.ISOFORM_SAFE_LENGTH_CACHE_FILE, ["gene", "transcript_accession", "isoform_safe_length"]
        )
        assert cache["TP53"] == ("ENST00000269305.9", "")

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
        assert df.empty
        assert calls == []  # nothing re-fetched

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
        assert len(calls) == 1

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
        assert df.to_dict("records") == [{"gene": "TP53", "isoform_safe_length": 9}]
