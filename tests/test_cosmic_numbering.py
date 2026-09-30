"""Unit tests for scripts/cosmic_numbering.py (COSMIC -> canonical-sequence
mutation numbering) and step 1's apply_canonical_numbering, which uses it.

Sequences here are small synthetic proteins: an "isoform" with an extra
exon inserted, so positions after the insertion are shifted relative to the
canonical sequence -- the case that used to be tagged (isoform?) wholesale.
"""
import sys

import pandas as pd
import pytest
from conftest import SCRIPTS_DIR, FakeResponse

sys.path.insert(0, str(SCRIPTS_DIR))
import cosmic_numbering as cn  # noqa: E402

# Canonical: 40 residues. Isoform -2: an extra 6-residue exon ("WWWWWW")
# inserted after canonical position 10, so canonical position p > 10 is
# isoform position p + 6.
CANONICAL = "MSTAPKLRVE" + "QDGSHNCYFI" + "PTRKAELGVM" + "DNSQYHCFAW"
ISOFORM = CANONICAL[:10] + "WWWWWW" + CANONICAL[10:]
ISOFORMS = {"P12345": CANONICAL, "P12345-2": ISOFORM}


def label(seq, pos, alt="A"):
    """A mutation label whose reference residue is correct for *seq* at *pos*."""
    return f"{seq[pos - 1]}{pos}{alt}"


class TestChooseIsoform:
    def test_isoform_explaining_more_mutations_wins(self):
        muts = [label(ISOFORM, 20), label(ISOFORM, 30), label(ISOFORM, 35)]
        assert cn.choose_numbering_isoform("P12345", ISOFORMS, muts) == "P12345-2"

    def test_canonical_wins_ties(self):
        # Positions before the inserted exon read the same in both sequences
        muts = [label(CANONICAL, 3), label(CANONICAL, 8)]
        assert cn.choose_numbering_isoform("P12345", ISOFORMS, muts) == "P12345", (
            "another isoform must explain strictly MORE mutations than the canonical "
            "sequence to be chosen -- a tie keeps canonical numbering"
        )


class TestNumberProtein:
    def test_shifted_mutations_move_to_canonical_positions(self):
        cosmic = label(ISOFORM, 26, "V")          # isoform 26 = canonical 20
        result = cn.number_protein("P12345", ISOFORMS, [cosmic, label(ISOFORM, 30)])
        assert result.isoform == "P12345-2"
        assert result.labels[cosmic] == f"{CANONICAL[19]}20V", (
            f"a mutation after the inserted exon should move back by the exon's length, "
            f"got {result.labels[cosmic]}"
        )
        assert result.unmapped == []

    def test_mutations_before_the_difference_keep_their_label(self):
        early = label(ISOFORM, 5)
        result = cn.number_protein("P12345", ISOFORMS, [early, label(ISOFORM, 26), label(ISOFORM, 30)])
        assert result.labels[early] == early

    def test_mutation_in_isoform_only_sequence_is_unmapped(self):
        in_exon = label(ISOFORM, 13)              # inside the inserted WWWWWW exon
        result = cn.number_protein("P12345", ISOFORMS, [in_exon, label(ISOFORM, 26), label(ISOFORM, 30)])
        assert result.unmapped == [in_exon], (
            "a residue that only exists in the isoform has no canonical position -- it "
            "must be reported as unmapped, keeping COSMIC's label"
        )
        assert result.labels[in_exon] == in_exon

    def test_canonical_numbering_leaves_matching_mutations_alone(self):
        muts = [label(CANONICAL, 25), label(CANONICAL, 33)]
        result = cn.number_protein("P12345", ISOFORMS, muts)
        assert result.isoform == ""
        assert result.labels == {m: m for m in muts}
        assert result.unmapped == []

    def test_reference_residue_matching_nowhere_is_unmapped(self):
        wrong = "W3A"                              # position 3 is T in every isoform
        result = cn.number_protein("P12345", ISOFORMS, [wrong, label(CANONICAL, 25)])
        assert result.unmapped == [wrong]

    def test_no_canonical_sequence_changes_nothing(self):
        result = cn.number_protein("P99999", {}, ["R5A"])
        assert result.labels == {"R5A": "R5A"} and result.unmapped == [] and result.isoform == ""


class TestPositionMap:
    def test_gap_positions_are_absent(self):
        pmap = cn.position_map(ISOFORM, CANONICAL)
        assert pmap[5] == 5 and pmap[26] == 20
        assert all(p not in pmap for p in range(11, 17)), "the inserted exon has no canonical position"


class TestFetchIsoformSequences:
    FASTA = (">sp|P12345|TEST_HUMAN Test protein\n" + CANONICAL + "\n"
             ">sp|P12345-2|TEST_HUMAN Isoform 2 of Test protein\n" + ISOFORM + "\n")

    def test_fetches_groups_and_caches(self, tmp_path):
        calls = []

        class Session:
            def get(self, url, params, timeout):
                calls.append(params["query"])
                return FakeResponse(TestFetchIsoformSequences.FASTA)

        cache = tmp_path / "iso.tsv"
        result = cn.fetch_isoform_sequences(["P12345", "O00000"], cache, session=Session())
        assert result["P12345"] == ISOFORMS
        assert result["O00000"] == {}, "an accession UniProt returned nothing for is cached as empty"

        again = cn.fetch_isoform_sequences(["P12345", "O00000"], cache, session=Session())
        assert again == result and len(calls) == 1, "a second run must be served from the cache"

    def test_failed_batch_falls_back_to_single_fetches(self, tmp_path):
        class Session:
            def get(self, url, params, timeout):
                if " OR " in params["query"] or "O00000" in params["query"]:
                    return FakeResponse("", status_code=500)
                return FakeResponse(TestFetchIsoformSequences.FASTA)

        cache = tmp_path / "iso.tsv"
        result = cn.fetch_isoform_sequences(["P12345", "O00000"], cache, session=Session(), backoff=0)
        assert result["P12345"] == ISOFORMS, "a failed batch must be retried one accession at a time"
        assert "O00000" not in cn.load_isoform_cache(cache), (
            "an accession that still fails must NOT be cached, so the next run retries it"
        )


class TestApplyCanonicalNumbering:
    @pytest.fixture
    def mod(self, filter_module, monkeypatch):
        monkeypatch.setattr(filter_module, "fetch_isoform_sequences", lambda *a, **k: {"P12345": ISOFORMS})
        return filter_module

    def test_rewrites_labels_and_records_changes(self, mod):
        moved, early, exon = label(ISOFORM, 26, "V"), label(ISOFORM, 5), label(ISOFORM, 13)
        table = pd.DataFrame({
            "uniprot_id": ["P12345"],
            "mutations_on_protein": [f"{moved} (7); {early} (3); {exon} (2); {label(ISOFORM, 30)} (4)"],
        })
        out = mod.apply_canonical_numbering(table).iloc[0]
        canonical_moved = f"{CANONICAL[19]}20V"
        assert out["mutations_on_protein"].startswith(f"{canonical_moved} (7); {early} (3); {exon} (2)"), (
            f"moved mutations take canonical numbering, patient counts carry over, got "
            f"{out['mutations_on_protein']!r}"
        )
        assert out["cosmic_numbering_isoform"] == "P12345-2"
        assert f"{canonical_moved}={moved}" in out["cosmic_mutation_labels"], (
            "every moved mutation must record COSMIC's original label"
        )
        assert out["unmapped_mutations"] == exon
