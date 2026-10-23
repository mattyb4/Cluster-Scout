"""Tests for scripts/validate_numbering.py (UniProt VAR_SEQ ground truth) and
randomized checks that cosmic_numbering.position_map recovers known edits.
"""
import random
import sys

import pytest
from conftest import SCRIPTS_DIR

sys.path.insert(0, str(SCRIPTS_DIR))
import cosmic_numbering as cn  # noqa: E402
import validate_numbering as vn  # noqa: E402

AA = "ACDEFGHIKLMNPQRSTVWY"


def random_protein(rng, n):
    return "".join(rng.choice(AA) for _ in range(n))


def random_edits(rng, canonical):
    """Non-overlapping VAR_SEQ-style edits: pure insertions/deletions of
    exon-sized blocks (replacement is empty or a sequence of random residues)."""
    edits, pos = [], 20
    while pos < len(canonical) - 60:
        start = rng.randint(pos, pos + 40)
        kind = rng.choice(["insert", "delete"])
        if kind == "delete":
            end = start + rng.randint(5, 30)
            edits.append((start, end, ""))
        else:
            # An insertion is expressed as replacing one residue with itself + new block
            edits.append((start, start, canonical[start - 1] + random_protein(rng, rng.randint(5, 30))))
        pos = edits[-1][1] + 20
    return edits


class TestBuildTruth:
    def test_deletion(self):
        canon = "ABCDEFGHIJ"
        seq, truth, replaced = vn.build_truth(canon, [(4, 6, "")])
        assert seq == "ABCGHIJ"
        assert truth == {1: 1, 2: 2, 3: 3, 4: 7, 5: 8, 6: 9, 7: 10}
        assert replaced == set()

    def test_replacement_positions_have_no_truth(self):
        seq, truth, replaced = vn.build_truth("ABCDEFGH", [(3, 4, "XYZ")])
        assert seq == "ABXYZEFGH"
        assert replaced == {3, 4, 5} and all(p not in truth for p in replaced)
        assert truth[6] == 5

    def test_parse_var_seq(self):
        entry = {
            "features": [
                {"type": "Alternative sequence", "featureId": "VSP_1",
                 "location": {"start": {"value": 3}, "end": {"value": 5}},
                 "alternativeSequence": {"alternativeSequences": ["XY"]}},
                {"type": "Alternative sequence", "featureId": "VSP_2",
                 "location": {"start": {"value": 1}, "end": {"value": 1}},
                 "alternativeSequence": {}},   # deletion: no alternativeSequences key
            ],
            "comments": [{"commentType": "ALTERNATIVE PRODUCTS", "isoforms": [
                {"isoformIds": ["P1-1"]},
                {"isoformIds": ["P1-2"], "sequenceIds": ["VSP_2", "VSP_1"]},
            ]}],
        }
        assert vn.parse_var_seq(entry) == {"P1-2": [(1, 1, ""), (3, 5, "XY")]}


class TestPositionMapRecoversKnownEdits:
    @pytest.mark.parametrize("seed", range(200))
    def test_random_isoforms(self, seed):
        rng = random.Random(seed)
        canon = random_protein(rng, rng.randint(200, 400))
        edits = random_edits(rng, canon)
        iso, truth, _ = vn.build_truth(canon, edits)
        pmap = cn.position_map(iso, canon)
        wrong = {p: pmap.get(p) for p in truth if pmap.get(p) != truth[p]}
        # A gap can legitimately slide when the residue at one edge of the block
        # equals the residue at the other (e.g. deleting "XAB" from "..XAB.." vs
        # "..AB X.."): equal-scoring alignments, one of which differs from UniProt's.
        # Those are the only disagreements allowed, and they must be rare. The
        # residue check in number_protein cannot catch them -- see methods.md.
        for p, got in wrong.items():
            if got is not None:
                assert canon[got - 1] == canon[truth[p] - 1], "mapped to a different residue, not a gap slide"
        assert len(wrong) <= 0.05 * len(truth)

    def test_edits_with_distinct_edges_are_recovered_exactly(self):
        canon = "MSTAPKLRVEQDGSHNCYFIPTRKAELGVMDNSQYHCFAW"
        for edits in ([(12, 18, "")], [(10, 10, "EWWWW")], [(5, 8, ""), (25, 25, "LGGGG")]):
            iso, truth, _ = vn.build_truth(canon, edits)
            pmap = cn.position_map(iso, canon)
            assert {p: pmap.get(p) for p in truth} == truth

    @pytest.mark.parametrize("seed", range(50))
    def test_map_is_monotonic_and_injective(self, seed):
        rng = random.Random(1000 + seed)
        canon = random_protein(rng, 300)
        iso, _, _ = vn.build_truth(canon, random_edits(rng, canon))
        pmap = cn.position_map(iso, canon)
        targets = [pmap[p] for p in sorted(pmap)]
        assert targets == sorted(set(targets))

    def test_identical_sequences_give_identity_map(self):
        canon = random_protein(random.Random(7), 150)
        assert cn.position_map(canon, canon) == {i: i for i in range(1, 151)}


class TestCompare:
    def test_reports_disagreements_with_both_positions(self):
        canon, iso = "ACDEFGHIKL", "ACDEFGHIKL"
        truth = {i: i for i in range(1, 11)}
        truth[5] = 6                      # pretend the truth says position 5 -> 6
        res = vn.compare(canon, iso, truth, set())
        assert res["wrong"] == {5: (6, 5)} and res["agree"] == 9
