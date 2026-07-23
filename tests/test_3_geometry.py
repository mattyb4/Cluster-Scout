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
        assert np.array_equal(coord, [0.0, 0.0, 0.0]), (
            "should return the CA atom's coordinate for the given residue, ignoring "
            f"the N/C atoms at the same residue -- got {coord}"
        )

    def test_returns_none_for_missing_residue(self, mod, chain):
        assert mod.get_ca_coord(chain, 999) is None, (
            "a residue number with no matching atom in the chain must return None, "
            "not raise or return a garbage coordinate"
        )


class TestComputeDistance:
    def test_basic_distance(self, mod):
        result = mod.compute_distance(np.array([0.0, 0.0, 0.0]), np.array([3.0, 4.0, 0.0]))
        assert result == 5.0, f"a 3-4-5 right triangle should give Euclidean distance 5.0, got {result}"


class TestFindNearbyMutations:
    def test_only_returns_hits_within_cutoff(self, mod, chain):
        mutation_entries = [("A15B", 15), ("C30D", 30)]
        results = mod.find_nearby_mutations(chain, ptm_pos=10, mutation_entries=mutation_entries, cutoff=10.0)

        assert len(results) == 1, (
            "residue 15 (5A away) is within the 10A cutoff and residue 30 (50A away) is "
            f"not -- exactly one hit should be returned, got {len(results)}"
        )
        assert results[0]["mutation"] == "A15B", "the surviving hit should be the in-cutoff mutation A15B"
        assert results[0]["mutation_pos"] == 15, "mutation_pos should carry through the residue position unchanged"
        assert results[0]["distance"] == pytest.approx(5.0), (
            f"distance should be the real 3D Euclidean distance (5A), got {results[0]['distance']}"
        )
        assert results[0]["pae"] is None, "no pae_matrix was supplied, so pae must be None, not 0 or omitted"

    def test_returns_empty_if_ptm_position_missing(self, mod, chain):
        results = mod.find_nearby_mutations(chain, ptm_pos=999, mutation_entries=[("A15B", 15)], cutoff=10.0)
        assert results == [], (
            "if the PTM site itself has no CA atom in the structure, there's no anchor to "
            "measure distances from -- must return an empty list, not raise"
        )

    def test_skips_mutations_with_no_ca_atom(self, mod, chain):
        results = mod.find_nearby_mutations(chain, ptm_pos=10, mutation_entries=[("X20Y", 20)], cutoff=10.0)
        assert results == [], (
            "a mutation at a residue position with no CA atom in the structure should be "
            "silently skipped, not counted as a hit or raise"
        )

    def test_pae_is_averaged_from_both_directions(self, mod, chain):
        pae_matrix = np.zeros((30, 30))
        pae_matrix[9, 14] = 2.0  # (ptm_pos-1, mut_pos-1)
        pae_matrix[14, 9] = 4.0

        results = mod.find_nearby_mutations(
            chain, ptm_pos=10, mutation_entries=[("A15B", 15)], pae_matrix=pae_matrix, cutoff=10.0
        )
        assert results[0]["pae"] == pytest.approx(3.0), (
            "PAE is directional (AlphaFold's PAE[i,j] != PAE[j,i] in general) -- the "
            f"reported value should be the average of both directions ((2+4)/2=3), got {results[0]['pae']}"
        )

    def test_max_pae_excludes_hits_above_the_threshold(self, mod, chain):
        pae_matrix = np.zeros((30, 30))
        pae_matrix[9, 14] = 8.0  # (ptm_pos-1, mut_pos-1) -- averages to 8.0
        pae_matrix[14, 9] = 8.0

        results = mod.find_nearby_mutations(
            chain, ptm_pos=10, mutation_entries=[("A15B", 15)], pae_matrix=pae_matrix,
            cutoff=10.0, max_pae=5.0,
        )
        assert results == [], (
            "a hit whose averaged PAE (8.0) exceeds max_pae (5.0) is too structurally "
            "uncertain to trust and must be excluded, not just flagged"
        )

    def test_max_pae_keeps_hits_at_or_below_the_threshold(self, mod, chain):
        pae_matrix = np.zeros((30, 30))
        pae_matrix[9, 14] = 5.0
        pae_matrix[14, 9] = 5.0

        results = mod.find_nearby_mutations(
            chain, ptm_pos=10, mutation_entries=[("A15B", 15)], pae_matrix=pae_matrix,
            cutoff=10.0, max_pae=5.0,
        )
        assert len(results) == 1, (
            "a hit exactly AT the max_pae threshold (5.0 == 5.0) should be kept, not "
            "excluded -- the filter is 'pae > max_pae', a strict inequality"
        )


class TestFindMutationClusters:
    def test_groups_mutations_within_cutoff_of_each_other(self, mod):
        chain = FakeChain([
            (1, "CA", [0.0, 0.0, 0.0]),
            (2, "CA", [3.0, 0.0, 0.0]),  # 3 Å from residue 1
            (3, "CA", [100.0, 0.0, 0.0]),  # far from everything
        ])
        mutation_entries = [("A1B", 1), ("C2D", 2), ("E3F", 3)]

        clusters = mod.find_mutation_clusters(chain, mutation_entries, cutoff=10.0)

        assert ("A1B", 1) in clusters, "residue 1 has a neighbor (residue 2) within cutoff, so it should anchor a cluster"
        nearby = clusters[("A1B", 1)]
        assert len(nearby) == 1, f"residue 1 should have exactly one neighbor within cutoff, got {len(nearby)}"
        assert nearby[0]["mutation"] == "C2D", "the neighbor found should be the in-cutoff mutation C2D"
        assert nearby[0]["mutation_pos"] == 2, "the neighbor's mutation_pos should carry through unchanged"

        assert ("E3F", 3) not in clusters, (
            "residue 3 is 100A from everything else -- with no neighbors within cutoff it "
            "must not appear as an anchor key at all (not even mapped to an empty list)"
        )

    def test_clustering_is_symmetric_both_directions_are_anchors(self, mod):
        chain = FakeChain([
            (1, "CA", [0.0, 0.0, 0.0]),
            (2, "CA", [3.0, 0.0, 0.0]),
        ])
        clusters = mod.find_mutation_clusters(chain, [("A1B", 1), ("C2D", 2)], cutoff=10.0)

        assert ("A1B", 1) in clusters and ("C2D", 2) in clusters, (
            "clustering is symmetric by design (every mutation with a neighbor becomes its "
            f"own anchor row) -- both residues should anchor a cluster, got keys {list(clusters.keys())}"
        )

    def test_mutation_at_position_missing_from_structure_is_skipped_as_anchor(self, mod):
        chain = FakeChain([(2, "CA", [0.0, 0.0, 0.0])])  # residue 1 has no CA atom at all
        clusters = mod.find_mutation_clusters(chain, [("A1B", 1), ("C2D", 2)], cutoff=10.0)

        assert ("A1B", 1) not in clusters, (
            "a mutation whose position has no CA atom in the structure has no coordinate "
            "to anchor a cluster from -- it must be skipped, not raise or appear with an "
            "empty/garbage neighbor list"
        )

    def test_neighbor_missing_from_structure_is_not_counted(self, mod):
        chain = FakeChain([(1, "CA", [0.0, 0.0, 0.0])])  # residue 2 has no CA atom
        clusters = mod.find_mutation_clusters(chain, [("A1B", 1), ("C2D", 2)], cutoff=10.0)

        assert ("A1B", 1) not in clusters, (
            "residue 2 (the only potential neighbor) has no coordinate in the structure, "
            "so residue 1 ends up with zero real neighbors and must not anchor a cluster"
        )

    def test_pae_is_averaged_from_both_directions(self, mod):
        chain = FakeChain([
            (1, "CA", [0.0, 0.0, 0.0]),
            (2, "CA", [3.0, 0.0, 0.0]),
        ])
        pae_matrix = np.zeros((2, 2))
        pae_matrix[0, 1] = 2.0  # (anchor_pos-1, other_pos-1)
        pae_matrix[1, 0] = 6.0

        clusters = mod.find_mutation_clusters(
            chain, [("A1B", 1), ("C2D", 2)], pae_matrix=pae_matrix, cutoff=10.0,
        )
        assert clusters[("A1B", 1)][0]["pae"] == pytest.approx(4.0), (
            "like find_nearby_mutations, PAE should be averaged across both directions "
            f"((2+6)/2=4), got {clusters[('A1B', 1)][0]['pae']}"
        )

    def test_max_pae_excludes_neighbors_above_the_threshold(self, mod):
        chain = FakeChain([
            (1, "CA", [0.0, 0.0, 0.0]),
            (2, "CA", [3.0, 0.0, 0.0]),
        ])
        pae_matrix = np.full((2, 2), 9.0)

        clusters = mod.find_mutation_clusters(
            chain, [("A1B", 1), ("C2D", 2)], pae_matrix=pae_matrix, cutoff=10.0, max_pae=5.0,
        )
        assert ("A1B", 1) not in clusters, (
            "residue 2 is within distance cutoff but its averaged PAE (9.0) exceeds "
            "max_pae (5.0) -- it must be filtered out, leaving residue 1 with no anchor-worthy neighbors"
        )

    def test_max_pae_keeps_neighbors_at_or_below_the_threshold(self, mod):
        chain = FakeChain([
            (1, "CA", [0.0, 0.0, 0.0]),
            (2, "CA", [3.0, 0.0, 0.0]),
        ])
        pae_matrix = np.full((2, 2), 5.0)

        clusters = mod.find_mutation_clusters(
            chain, [("A1B", 1), ("C2D", 2)], pae_matrix=pae_matrix, cutoff=10.0, max_pae=5.0,
        )
        assert ("A1B", 1) in clusters, (
            "a neighbor exactly AT the max_pae threshold (5.0 == 5.0) should be kept, "
            "not excluded -- the filter is 'pae > max_pae', a strict inequality"
        )
