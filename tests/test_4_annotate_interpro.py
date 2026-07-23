"""Unit tests for InterPro functional-domain (Phase 5) functions in scripts/4_annotate.py."""
import pandas as pd
import pytest

from conftest import import_script

mod = import_script("4_annotate.py")


class FakeInterproResponse:
    def __init__(self, json_data, status_code=200):
        self._json_data = json_data
        self.status_code = status_code

    def raise_for_status(self):
        if self.status_code >= 400:
            raise RuntimeError(f"HTTP {self.status_code}")

    def json(self):
        return self._json_data


def _entry(name, etype, start, end):
    return {
        "metadata": {"name": name, "type": etype},
        "proteins": [{
            "entry_protein_locations": [{"fragments": [{"start": start, "end": end}]}],
        }],
    }


class TestInterproCache:
    def test_round_trip(self, tmp_path, monkeypatch):
        monkeypatch.setattr(mod, "_INTERPRO_CACHE_FILE", tmp_path / "interpro.tsv")
        cache = {"P04637": [{"name": "p53 domain", "type": "domain", "start": 1, "end": 50}]}
        mod._interpro_save_cache(cache)
        result = mod._interpro_load_cache()
        assert result == cache, f"loading immediately after saving should reproduce the same dict, got {result}"

    def test_missing_file_returns_empty_dict(self, tmp_path, monkeypatch):
        monkeypatch.setattr(mod, "_INTERPRO_CACHE_FILE", tmp_path / "does_not_exist.tsv")
        result = mod._interpro_load_cache()
        assert result == {}, f"a never-written cache file should load as {{}}, got {result}"

    def test_malformed_entries_json_degrades_to_empty_list(self, tmp_path, monkeypatch):
        cache_file = tmp_path / "interpro.tsv"
        monkeypatch.setattr(mod, "_INTERPRO_CACHE_FILE", cache_file)
        pd.DataFrame([
            {"uniprot_id": "P04637", "entries_json": "not valid json{{{"},
        ]).to_csv(cache_file, sep="\t", index=False)

        result = mod._interpro_load_cache()
        assert result == {"P04637": []}, (
            f"a corrupted entries_json cell should degrade to an empty entry list for "
            f"that protein (not skip the row entirely), got {result}"
        )

    def test_unreadable_file_returns_empty_dict(self, tmp_path, monkeypatch):
        cache_file = tmp_path / "interpro.tsv"
        cache_file.write_bytes(b"\xff\xfe\x00garbage-not-a-real-tsv")
        monkeypatch.setattr(mod, "_INTERPRO_CACHE_FILE", cache_file)

        result = mod._interpro_load_cache()
        assert result == {}, f"a totally unparseable cache file should degrade to {{}}, not raise, got {result}"


class TestFetchInterproDomains:
    def test_extracts_name_type_and_fragment_range(self, monkeypatch):
        monkeypatch.setattr(
            mod.requests, "get",
            lambda url, timeout=30: FakeInterproResponse({
                "results": [_entry("p53 DNA-binding domain", "domain", 94, 289)],
                "next": None,
            }),
        )
        result = mod.fetch_interpro_domains("P04637")
        assert result == [{"name": "p53 DNA-binding domain", "type": "domain", "start": 94, "end": 289}], (
            f"a single-fragment entry should produce one {{name,type,start,end}} dict, got {result}"
        )

    def test_entry_with_multiple_fragments_produces_multiple_entries(self, monkeypatch):
        multi_fragment_entry = {
            "metadata": {"name": "Split domain", "type": "domain"},
            "proteins": [{
                "entry_protein_locations": [{
                    "fragments": [{"start": 10, "end": 50}, {"start": 200, "end": 240}],
                }],
            }],
        }
        monkeypatch.setattr(
            mod.requests, "get",
            lambda url, timeout=30: FakeInterproResponse({"results": [multi_fragment_entry], "next": None}),
        )
        result = mod.fetch_interpro_domains("P04637")
        assert len(result) == 2, (
            f"a discontinuous (2-fragment) entry should produce 2 separate range "
            f"entries, each independently usable for containment checks, got {len(result)}"
        )
        assert {"name": "Split domain", "type": "domain", "start": 10, "end": 50} in result
        assert {"name": "Split domain", "type": "domain", "start": 200, "end": 240} in result

    def test_follows_pagination_next_link(self, monkeypatch):
        page1_url = "https://www.ebi.ac.uk/interpro/api/entry/interpro/protein/uniprot/P04637/"
        page2_url = "https://www.ebi.ac.uk/interpro/api/entry/interpro/protein/uniprot/P04637/?page=2"
        pages = {
            page1_url: {"results": [_entry("Domain A", "domain", 1, 50)], "next": page2_url},
            page2_url: {"results": [_entry("Domain B", "family", 100, 150)], "next": None},
        }
        monkeypatch.setattr(mod.requests, "get", lambda url, timeout=30: FakeInterproResponse(pages[url]))

        result = mod.fetch_interpro_domains("P04637")
        names = {e["name"] for e in result}
        assert names == {"Domain A", "Domain B"}, (
            f"entries from BOTH pages must be collected by following the 'next' link "
            f"until it's None, got names {names}"
        )

    def test_returns_partial_results_on_error_mid_pagination(self, monkeypatch):
        page1_url = "https://www.ebi.ac.uk/interpro/api/entry/interpro/protein/uniprot/P04637/"

        def fake_get(url, timeout=30):
            if url == page1_url:
                return FakeInterproResponse({"results": [_entry("Domain A", "domain", 1, 50)],
                                              "next": "https://broken-next-page/"})
            raise mod.requests.RequestException("network error on page 2")

        monkeypatch.setattr(mod.requests, "get", fake_get)
        result = mod.fetch_interpro_domains("P04637")
        assert result == [{"name": "Domain A", "type": "domain", "start": 1, "end": 50}], (
            f"if a later page fails, whatever was already collected from earlier pages "
            f"should still be returned (not discarded), got {result}"
        )

    def test_returns_empty_list_on_immediate_network_error(self, monkeypatch):
        monkeypatch.setattr(
            mod.requests, "get",
            lambda url, timeout=30: (_ for _ in ()).throw(mod.requests.RequestException("timeout")),
        )
        result = mod.fetch_interpro_domains("P04637")
        assert result == [], f"a network error on the very first request should return [], not raise, got {result}"

    def test_no_curated_entries_returns_empty_list(self, monkeypatch):
        monkeypatch.setattr(
            mod.requests, "get",
            lambda url, timeout=30: FakeInterproResponse({"results": [], "next": None}),
        )
        result = mod.fetch_interpro_domains("P99999")
        assert result == [], f"a protein with no curated InterPro entries should return [], got {result}"


class TestRunInterproPhase:
    def test_all_cached_makes_no_fetch_calls(self, tmp_path, monkeypatch):
        monkeypatch.setattr(mod, "_INTERPRO_CACHE_FILE", tmp_path / "interpro.tsv")
        mod._interpro_save_cache({"P04637": [{"name": "X", "type": "domain", "start": 1, "end": 10}]})
        calls = []
        monkeypatch.setattr(mod, "fetch_interpro_domains", lambda uid: calls.append(uid) or [])

        df = pd.DataFrame([{"UniProt": "P04637"}])
        result = mod.run_interpro_phase(df)

        assert calls == [], f"an already-cached protein should not be re-fetched, got {len(calls)} call(s)"
        assert result["P04637"] == [{"name": "X", "type": "domain", "start": 1, "end": 10}], (
            f"the cached entries should be returned, got {result['P04637']}"
        )

    def test_fetches_only_uncached_proteins_and_saves_cache(self, tmp_path, monkeypatch):
        cache_file = tmp_path / "interpro.tsv"
        monkeypatch.setattr(mod, "_INTERPRO_CACHE_FILE", cache_file)
        mod._interpro_save_cache({"P00001": [{"name": "X", "type": "domain", "start": 1, "end": 10}]})
        calls = []

        def fake_fetch(uid):
            calls.append(uid)
            return [{"name": "Y", "type": "family", "start": 5, "end": 20}]

        monkeypatch.setattr(mod, "fetch_interpro_domains", fake_fetch)

        df = pd.DataFrame([{"UniProt": "P00001"}, {"UniProt": "P00002"}])
        result = mod.run_interpro_phase(df)

        assert calls == ["P00002"], f"only the uncached protein should be fetched, got {calls}"
        assert result["P00002"] == [{"name": "Y", "type": "family", "start": 5, "end": 20}], (
            f"the newly fetched entries should appear in the result, got {result['P00002']}"
        )
        assert cache_file.exists(), "the cache should be persisted after fetching new data"


class TestFindDomainAtPosition:
    ENTRIES = [
        {"name": "Domain A", "type": "domain", "start": 1, "end": 100},
        {"name": "Superfamily B", "type": "homologous_superfamily", "start": 50, "end": 150},
        {"name": "Domain C", "type": "domain", "start": 200, "end": 250},
    ]

    def test_returns_single_match(self):
        result = mod.find_domain_at_position(self.ENTRIES, 220)
        assert result == "Domain C (domain, 200-250)", (
            f"a position inside exactly one entry should format as 'name (type, start-end)', got {result!r}"
        )

    def test_returns_all_overlapping_matches_joined(self):
        result = mod.find_domain_at_position(self.ENTRIES, 75)
        assert "Domain A (domain, 1-100)" in result and "Superfamily B (homologous_superfamily, 50-150)" in result, (
            f"a position inside TWO overlapping entries (a nested domain inside a "
            f"broader superfamily call) must report both, got {result!r}"
        )
        assert result.count(";") == 1, f"exactly 2 matches should be joined by exactly 1 semicolon, got {result!r}"

    def test_returns_empty_string_for_no_match(self):
        result = mod.find_domain_at_position(self.ENTRIES, 500)
        assert result == "", f"a position outside every entry's range should return '', got {result!r}"

    def test_returns_empty_string_for_none_position(self):
        result = mod.find_domain_at_position(self.ENTRIES, None)
        assert result == "", (
            f"a None position (e.g. an unparseable mutation/PTM label upstream) should "
            f"return '' rather than raise a TypeError comparing None <= int, got {result!r}"
        )

    def test_boundary_positions_are_inclusive(self):
        assert mod.find_domain_at_position(self.ENTRIES, 1) == "Domain A (domain, 1-100)", (
            "the start boundary (position == start) must be included (inclusive range)"
        )
        # Position 100 is Domain A's end boundary AND Superfamily B's interior (50-150
        # overlaps), so both legitimately match here -- check Domain A's inclusion
        # specifically rather than exact equality against the whole joined string.
        assert "Domain A (domain, 1-100)" in mod.find_domain_at_position(self.ENTRIES, 100), (
            "the end boundary (position == end) must also be included (inclusive range)"
        )
