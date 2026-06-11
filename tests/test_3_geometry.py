"""Unit tests for the 3D-geometry helpers in scripts/3_find_nearby_mutations.py."""
import numpy as np
import pytest


@pytest.fixture
def mod(nearby_module):
    return nearby_module


class FakeChain:
    """Minimal stand-in for a biotite AtomArray chain.

    `atoms` is a list of (res_id, atom_name, [x, y, z]) tuples.
    """

    def __init__(self, atoms):
        self.res_id = np.array([a[0] for a in atoms])
        self.atom_name = np.array([a[1] for a in atoms], dtype=object)
        self.coord = np.array([a[2] for a in atoms], dtype=float)


@pytest.fixture
def chain():
    return FakeChain([
        (10, "N", [0.0, 0.0, 0.0]),
        (10, "CA", [0.0, 0.0, 0.0]),  # PTM site at the origin
        (10, "C", [0.0, 1.0, 0.0]),
        (15, "CA", [5.0, 0.0, 0.0]),  # 5 Å from PTM site -> within cutoff
        (30, "CA", [50.0, 0.0, 0.0]),  # 50 Å from PTM site -> beyond cutoff
    ])


class TestGetCaCoord:
    def test_returns_ca_coordinate(self, mod, chain):
        coord = mod.get_ca_coord(chain, 10)
        assert np.array_equal(coord, [0.0, 0.0, 0.0])

    def test_returns_none_for_missing_residue(self, mod, chain):
        assert mod.get_ca_coord(chain, 999) is None


class TestComputeDistance:
    def test_basic_distance(self, mod):
        assert mod.compute_distance(np.array([0.0, 0.0, 0.0]), np.array([3.0, 4.0, 0.0])) == 5.0


class TestFindNearbyMutations:
    def test_only_returns_hits_within_cutoff(self, mod, chain):
        mutation_entries = [("A15B", 15), ("C30D", 30)]
        results = mod.find_nearby_mutations(chain, ptm_pos=10, mutation_entries=mutation_entries, cutoff=10.0)

        assert len(results) == 1
        assert results[0]["mutation"] == "A15B"
        assert results[0]["mutation_pos"] == 15
        assert results[0]["distance"] == pytest.approx(5.0)
        assert results[0]["pae"] is None

    def test_returns_empty_if_ptm_position_missing(self, mod, chain):
        results = mod.find_nearby_mutations(chain, ptm_pos=999, mutation_entries=[("A15B", 15)], cutoff=10.0)
        assert results == []

    def test_skips_mutations_with_no_ca_atom(self, mod, chain):
        results = mod.find_nearby_mutations(chain, ptm_pos=10, mutation_entries=[("X20Y", 20)], cutoff=10.0)
        assert results == []

    def test_pae_is_averaged_from_both_directions(self, mod, chain):
        pae_matrix = np.zeros((30, 30))
        pae_matrix[9, 14] = 2.0  # (ptm_pos-1, mut_pos-1)
        pae_matrix[14, 9] = 4.0

        results = mod.find_nearby_mutations(
            chain, ptm_pos=10, mutation_entries=[("A15B", 15)], pae_matrix=pae_matrix, cutoff=10.0
        )
        assert results[0]["pae"] == pytest.approx(3.0)


class TestFindMutationClusters:
    def test_groups_mutations_within_cutoff_of_each_other(self, mod):
        chain = FakeChain([
            (1, "CA", [0.0, 0.0, 0.0]),
            (2, "CA", [3.0, 0.0, 0.0]),  # 3 Å from residue 1
            (3, "CA", [100.0, 0.0, 0.0]),  # far from everything
        ])
        mutation_entries = [("A1B", 1), ("C2D", 2), ("E3F", 3)]

        clusters = mod.find_mutation_clusters(chain, mutation_entries, cutoff=10.0)

        assert ("A1B", 1) in clusters
        nearby = clusters[("A1B", 1)]
        assert len(nearby) == 1
        assert nearby[0]["mutation"] == "C2D"
        assert nearby[0]["mutation_pos"] == 2

        # Residue 3 has nothing nearby, so it shouldn't appear as an anchor.
        assert ("E3F", 3) not in clusters
