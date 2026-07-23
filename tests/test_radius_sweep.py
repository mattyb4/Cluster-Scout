"""Unit tests for scripts/radius_sweep.py.

_ptm_tsv_index_cache/_cosmic_cache/_cosmic_counts_cache are module-level
globals that persist across every test in this file (mod is imported once at
collection time) -- the autouse `_reset_module_caches` fixture below clears
them before AND after each test so no test's cached data can leak into another.
"""
import numpy as np
import pandas as pd
import pytest

from conftest import import_script

mod = import_script("radius_sweep.py")


@pytest.fixture(autouse=True)
def _reset_module_caches():
    mod._ptm_tsv_index_cache = None
    mod._cosmic_cache = None
    mod._cosmic_counts_cache = None
    yield
    mod._ptm_tsv_index_cache = None
    mod._cosmic_cache = None
    mod._cosmic_counts_cache = None


class FakeChain:
    """Supports both attribute access (.res_id/.atom_name/.coord) and
    boolean-mask __getitem__/__len__, since radius_sweep.py's functions use
    both styles (get_ca_coord reads .res_id/.atom_name/.coord directly;
    get_all_ca_coords does chain[ca_mask] first).
    """

    def __init__(self, res_id, atom_name, coord):
        self.res_id = np.array(res_id)
        self.atom_name = np.array(atom_name, dtype=object)
        self.coord = np.array(coord, dtype=float)

    def __getitem__(self, mask):
        return FakeChain(self.res_id[mask], self.atom_name[mask], self.coord[mask])

    def __len__(self):
        return len(self.res_id)


class TestLooksLikeUniprotId:
    @pytest.mark.parametrize("token,expected", [
        ("P04637", True),
        ("A0A099Z4Y8", True),
        ("TP53", False),
        ("EGFR", False),
        ("", False),
    ])
    def test_matches_expected(self, token, expected):
        result = mod.looks_like_uniprot_id(token)
        assert result is expected, (
            f"looks_like_uniprot_id({token!r}) should be {expected} based on UniProt's "
            f"accession format, got {result}"
        )


class TestLoadPtmTsvIndex:
    def test_loads_gene_and_uniprot_columns(self, tmp_path, monkeypatch):
        tsv = tmp_path / "hotspots.tsv"
        pd.DataFrame([{"gene": "TP53", "uniprot_id": "P04637", "extra_col": "ignored"}]).to_csv(tsv, sep="\t", index=False)
        monkeypatch.setattr(mod, "PTM_TSV", tsv)

        df = mod._load_ptm_tsv_index()
        assert list(df.columns) == ["gene", "uniprot_id"], (
            f"only the gene/uniprot_id columns should be loaded (usecols), got {list(df.columns)}"
        )
        assert df.iloc[0]["gene"] == "TP53"

    def test_missing_file_returns_empty_dataframe_uncached(self, tmp_path, monkeypatch):
        monkeypatch.setattr(mod, "PTM_TSV", tmp_path / "does_not_exist.tsv")

        df1 = mod._load_ptm_tsv_index()
        assert df1.empty, "a missing PTM_TSV should give an empty DataFrame, not raise"
        assert mod._ptm_tsv_index_cache is None, (
            "the missing-file case must NOT populate the module cache -- so that if the "
            "file appears later (e.g. after running step 1), a subsequent call retries "
            "reading it instead of staying permanently empty"
        )

    def test_second_call_is_served_from_cache(self, tmp_path, monkeypatch):
        tsv = tmp_path / "hotspots.tsv"
        pd.DataFrame([{"gene": "TP53", "uniprot_id": "P04637"}]).to_csv(tsv, sep="\t", index=False)
        monkeypatch.setattr(mod, "PTM_TSV", tsv)

        first = mod._load_ptm_tsv_index()
        tsv.write_text("gene\tuniprot_id\nCHANGED\tQ99999\n")  # mutate the file on disk
        second = mod._load_ptm_tsv_index()

        assert second is first, "a second call should return the cached DataFrame object, not re-read the file"
        assert second.iloc[0]["gene"] == "TP53", (
            "the cached (stale) content should still reflect what was on disk at the "
            "time of the FIRST call, proving the file wasn't re-read"
        )


class TestLoadKnownGenes:
    def test_returns_gene_set(self, tmp_path, monkeypatch):
        tsv = tmp_path / "hotspots.tsv"
        pd.DataFrame([
            {"gene": "TP53", "uniprot_id": "P04637"},
            {"gene": "TP53", "uniprot_id": "P04637"},  # duplicate row
            {"gene": "EGFR", "uniprot_id": "P00533"},
        ]).to_csv(tsv, sep="\t", index=False)
        monkeypatch.setattr(mod, "PTM_TSV", tsv)

        result = mod.load_known_genes()
        assert result == {"TP53", "EGFR"}, f"duplicate rows should collapse to a set of unique genes, got {result}"

    def test_missing_file_returns_empty_set(self, tmp_path, monkeypatch):
        monkeypatch.setattr(mod, "PTM_TSV", tmp_path / "does_not_exist.tsv")
        assert mod.load_known_genes() == set(), "a missing PTM_TSV should give an empty set, not raise"


class TestResolveGeneToken:
    @pytest.fixture
    def tsv(self, tmp_path, monkeypatch):
        path = tmp_path / "hotspots.tsv"
        pd.DataFrame([{"gene": "TP53", "uniprot_id": "P04637"}]).to_csv(path, sep="\t", index=False)
        monkeypatch.setattr(mod, "PTM_TSV", path)
        return path

    def test_resolves_by_gene_symbol(self, tsv):
        result = mod.resolve_gene_token("TP53")
        assert result == ("TP53", "P04637"), f"a gene symbol should resolve to (gene, uniprot_id), got {result}"

    def test_resolves_by_uniprot_accession(self, tsv):
        result = mod.resolve_gene_token("P04637")
        assert result == ("TP53", "P04637"), (
            f"a UniProt accession should also resolve, via looks_like_uniprot_id routing "
            f"to the uniprot_id column instead of gene, got {result}"
        )

    def test_case_insensitive(self, tsv):
        result = mod.resolve_gene_token("tp53")
        assert result == ("TP53", "P04637"), f"gene symbol matching should be case-insensitive, got {result}"

    def test_unknown_token_returns_none(self, tsv):
        assert mod.resolve_gene_token("NOTAGENE") is None, "a token not present in the dataset should return None"

    def test_empty_index_returns_none(self, tmp_path, monkeypatch):
        monkeypatch.setattr(mod, "PTM_TSV", tmp_path / "does_not_exist.tsv")
        assert mod.resolve_gene_token("TP53") is None, (
            "with no PTM_TSV data available at all, nothing can be resolved -- must "
            "return None rather than raise"
        )


class TestHasCif:
    def test_true_when_canonical_cif_present(self, tmp_path, monkeypatch):
        monkeypatch.setattr(mod, "MODELS_ROOT", tmp_path)
        cif_dir = tmp_path / "P04637"
        cif_dir.mkdir()
        (cif_dir / "AF-P04637-F1-model_v4.cif").write_text("x")
        assert mod.has_cif("P04637") is True

    def test_false_when_directory_missing(self, tmp_path, monkeypatch):
        monkeypatch.setattr(mod, "MODELS_ROOT", tmp_path)
        assert mod.has_cif("Q99999") is False, "a protein with no download directory at all has no CIF"

    def test_false_when_only_isoform_present(self, tmp_path, monkeypatch):
        monkeypatch.setattr(mod, "MODELS_ROOT", tmp_path)
        cif_dir = tmp_path / "P04637"
        cif_dir.mkdir()
        (cif_dir / "AF-P04637-2-F1-model_v4.cif").write_text("x")
        assert mod.has_cif("P04637") is False, "an isoform-only download must not count as having a canonical CIF"


class TestHasMultipleFragments:
    def test_true_with_two_fragments(self, tmp_path, monkeypatch):
        monkeypatch.setattr(mod, "MODELS_ROOT", tmp_path)
        cif_dir = tmp_path / "P04637"
        cif_dir.mkdir()
        (cif_dir / "AF-P04637-F1-model_v4.cif").write_text("x")
        (cif_dir / "AF-P04637-F2-model_v4.cif").write_text("x")
        assert mod.has_multiple_fragments("P04637") is True

    def test_false_with_one_fragment(self, tmp_path, monkeypatch):
        monkeypatch.setattr(mod, "MODELS_ROOT", tmp_path)
        cif_dir = tmp_path / "P04637"
        cif_dir.mkdir()
        (cif_dir / "AF-P04637-F1-model_v4.cif").write_text("x")
        assert mod.has_multiple_fragments("P04637") is False

    def test_false_when_directory_missing(self, tmp_path, monkeypatch):
        monkeypatch.setattr(mod, "MODELS_ROOT", tmp_path)
        assert mod.has_multiple_fragments("Q99999") is False


class TestGetCaCoord:
    def test_returns_coord(self):
        chain = FakeChain(res_id=[1, 2], atom_name=["CA", "CA"], coord=[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]])
        result = mod.get_ca_coord(chain, 2)
        assert np.array_equal(result, [1.0, 0.0, 0.0]), f"should return residue 2's CA coordinate, got {result}"

    def test_returns_none_for_missing_residue(self):
        chain = FakeChain(res_id=[1], atom_name=["CA"], coord=[[0.0, 0.0, 0.0]])
        assert mod.get_ca_coord(chain, 999) is None


class TestLoadPtmPositions:
    def test_extracts_coords_for_gene(self):
        chain = FakeChain(res_id=[15, 18], atom_name=["CA", "CA"], coord=[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]])
        df = pd.DataFrame([{"gene": "TP53", "ptms_on_protein": "S15:Phosphorylation; T18:Phosphorylation"}])

        result = mod.load_ptm_positions("TP53", df, chain)
        assert [p for p, _ in result] == [15, 18], f"both PTM positions should resolve to CA coords, got {[p for p, _ in result]}"

    def test_gene_not_in_df_returns_empty_list(self):
        chain = FakeChain(res_id=[15], atom_name=["CA"], coord=[[0.0, 0.0, 0.0]])
        df = pd.DataFrame([{"gene": "TP53", "ptms_on_protein": "S15:Phosphorylation"}])
        result = mod.load_ptm_positions("EGFR", df, chain)
        assert result == [], "a gene with no row in the DataFrame should give no PTM coordinates"

    def test_position_with_no_ca_atom_is_excluded(self):
        chain = FakeChain(res_id=[999], atom_name=["CA"], coord=[[0.0, 0.0, 0.0]])
        df = pd.DataFrame([{"gene": "TP53", "ptms_on_protein": "S15:Phosphorylation"}])
        result = mod.load_ptm_positions("TP53", df, chain)
        assert result == [], (
            "a PTM site at a position with no corresponding CA atom in the structure "
            "must be excluded, not raise or produce a None coordinate"
        )


class TestGetAllCaCoords:
    def test_returns_all_ca_pairs(self):
        chain = FakeChain(
            res_id=[1, 2, 3], atom_name=["N", "CA", "CA"],
            coord=[[9, 9, 9], [0, 0, 0], [1, 0, 0]],
        )
        result = mod.get_all_ca_coords(chain)
        assert [p for p, _ in result] == [2, 3], (
            f"only CA atoms should be included (the N atom at residue 1 excluded), got {[p for p, _ in result]}"
        )


class TestSweepRadii:
    def test_counts_increase_with_radius(self):
        ptm_coords = [(10, np.array([0.0, 0.0, 0.0]))]
        mutation_coords = [(11, np.array([3.0, 0.0, 0.0])), (12, np.array([8.0, 0.0, 0.0]))]
        result = mod.sweep_radii(ptm_coords, mutation_coords, [5.0, 10.0])
        assert result[5.0] == 1.0, f"only the 3A-away mutation should be within radius 5, got {result[5.0]}"
        assert result[10.0] == 2.0, f"both mutations should be within radius 10, got {result[10.0]}"

    def test_empty_ptm_coords_returns_zeros(self):
        result = mod.sweep_radii([], [(1, np.array([0.0, 0.0, 0.0]))], [5.0, 10.0])
        assert result == {5.0: 0.0, 10.0: 0.0}, "no PTM sites means nothing to average over -- should be all zeros, not raise"

    def test_empty_mutation_coords_returns_zeros(self):
        result = mod.sweep_radii([(1, np.array([0.0, 0.0, 0.0]))], [], [5.0, 10.0])
        assert result == {5.0: 0.0, 10.0: 0.0}, "no mutations means nothing can be within range -- should be all zeros"


class TestRandomBaseline:
    def test_zero_mutations_returns_zeros(self):
        ptm_coords = [(1, np.array([0.0, 0.0, 0.0]))]
        all_ca = [(i, np.array([float(i), 0.0, 0.0])) for i in range(10)]
        result = mod.random_baseline(ptm_coords, all_ca, n_mutations=0, radii=[5.0, 10.0])
        assert result == {5.0: 0.0, 10.0: 0.0}, "with n_mutations=0 there's nothing to place randomly -- should be all zeros"

    def test_insufficient_residues_returns_zeros(self):
        ptm_coords = [(1, np.array([0.0, 0.0, 0.0]))]
        all_ca = [(i, np.array([float(i), 0.0, 0.0])) for i in range(3)]
        # Asking to place more "random mutations" than residues exist in the protein.
        result = mod.random_baseline(ptm_coords, all_ca, n_mutations=10, radii=[5.0])
        assert result == {5.0: 0.0}, (
            "sampling without replacement can't place more mutations than there are "
            "residues -- must degrade to zeros rather than raise"
        )

    def test_deterministic_with_mocked_sampling(self, monkeypatch):
        # Force np.random.choice to always pick the first n_mutations indices,
        # making the result fully deterministic and independently verifiable.
        monkeypatch.setattr(mod.np.random, "choice", lambda n, size, replace: np.arange(size))

        ptm_coords = [(1, np.array([0.0, 0.0, 0.0]))]
        all_ca = [(i, np.array([float(i), 0.0, 0.0])) for i in range(5)]
        result = mod.random_baseline(ptm_coords, all_ca, n_mutations=2, radii=[10.0], n_permutations=3)

        # Every permutation picks residues at x=0,1 (indices 0,1) -- both within
        # radius 10 of the PTM at x=0 -- so the average should be exactly 2.0.
        assert result[10.0] == pytest.approx(2.0), (
            f"with sampling forced to always pick the 2 closest residues, every "
            f"permutation should count exactly 2 within radius 10, got {result[10.0]}"
        )


class TestRunSweep:
    def test_raises_filenotfounderror_when_ptm_tsv_missing(self, tmp_path, monkeypatch):
        monkeypatch.setattr(mod, "PTM_TSV", tmp_path / "does_not_exist.tsv")
        with pytest.raises(FileNotFoundError):
            mod.run_sweep(["TP53"], [5.0, 10.0], log_cb=lambda *_: None)

    def test_raises_valueerror_when_no_data_collected(self, tmp_path, monkeypatch):
        tsv = tmp_path / "hotspots.tsv"
        pd.DataFrame([{"gene": "TP53", "uniprot_id": "P04637", "ptms_on_protein": "S15:Phosphorylation"}]).to_csv(tsv, sep="\t", index=False)
        monkeypatch.setattr(mod, "PTM_TSV", tsv)
        monkeypatch.setattr(mod, "MODELS_ROOT", tmp_path / "cif_models")  # no CIF downloaded -> every gene skipped

        with pytest.raises(ValueError):
            mod.run_sweep(["TP53"], [5.0, 10.0], log_cb=lambda *_: None)

    def test_full_run_produces_hotspot_results(self, tmp_path, monkeypatch):
        tsv = tmp_path / "hotspots.tsv"
        pd.DataFrame([{"gene": "TP53", "uniprot_id": "P04637", "ptms_on_protein": "S1:Phosphorylation"}]).to_csv(tsv, sep="\t", index=False)
        monkeypatch.setattr(mod, "PTM_TSV", tsv)
        monkeypatch.setattr(mod, "MODELS_ROOT", tmp_path / "cif_models")

        chain = FakeChain(
            res_id=[1, 2, 3], atom_name=["CA", "CA", "CA"],
            coord=[[0.0, 0.0, 0.0], [3.0, 0.0, 0.0], [8.0, 0.0, 0.0]],
        )
        monkeypatch.setattr(mod, "find_canonical_cif", lambda cif_dir: "fake.cif")
        monkeypatch.setattr(mod, "find_canonical_cifs", lambda cif_dir: ["fake.cif"])
        monkeypatch.setattr(mod, "load_first_chain", lambda cif_file: chain)
        (tmp_path / "cif_models" / "P04637").mkdir(parents=True)

        counts = pd.DataFrame([{"gene": "TP53", "aa_change": "A2B", "affected_cases": 5}])
        monkeypatch.setattr(mod, "_load_cosmic_counts_df", lambda log_cb=print: counts)

        result = mod.run_sweep(["TP53"], [5.0, 10.0], min_cases=3, log_cb=lambda *_: None)

        assert set(result.result_df["dataset"]) == {"hotspot"}, (
            f"with unfiltered=False, only the 'hotspot' dataset should be present, got {set(result.result_df['dataset'])}"
        )
        assert len(result.result_df) == 2, f"one row per radius tested (2 radii) should be produced, got {len(result.result_df)}"
        assert "TP53" in result.elbows, "an elbow-detection entry should exist for the analyzed protein"

    def test_skips_gene_with_no_cif(self, tmp_path, monkeypatch):
        tsv = tmp_path / "hotspots.tsv"
        pd.DataFrame([{"gene": "TP53", "uniprot_id": "P04637", "ptms_on_protein": "S1:Phosphorylation"}]).to_csv(tsv, sep="\t", index=False)
        monkeypatch.setattr(mod, "PTM_TSV", tsv)
        monkeypatch.setattr(mod, "MODELS_ROOT", tmp_path / "cif_models")  # empty -- no CIF dir for P04637

        with pytest.raises(ValueError):
            # Should skip TP53 (no CIF) and end up with zero results overall.
            mod.run_sweep(["TP53"], [5.0], log_cb=lambda *_: None)


class TestBuildSweepFigure:
    def _minimal_result(self):
        df = pd.DataFrame([
            {"protein": "TP53", "radius": 5.0, "dataset": "hotspot", "avg_mutation_count": 1.0,
             "random_baseline": 0.5, "avg_normalized": 2.0, "random_normalized": 1.0,
             "protein_length": 500, "n_mutations": 3},
            {"protein": "TP53", "radius": 10.0, "dataset": "hotspot", "avg_mutation_count": 3.0,
             "random_baseline": 1.5, "avg_normalized": 6.0, "random_normalized": 3.0,
             "protein_length": 500, "n_mutations": 3},
        ])
        return mod.SweepResult(result_df=df, radii=[5.0, 10.0], elbows={"TP53": None}, has_unfiltered=False)

    def test_returns_2x1_layout_without_unfiltered(self):
        import matplotlib.pyplot as plt
        fig = mod.build_sweep_figure(self._minimal_result())
        try:
            assert len(fig.axes) == 2, f"without unfiltered data, the figure should have 2 panels, got {len(fig.axes)}"
        finally:
            plt.close(fig)

    def test_uses_injected_figure(self):
        from matplotlib.figure import Figure
        injected = Figure()
        returned = mod.build_sweep_figure(self._minimal_result(), fig=injected)
        assert returned is injected, "build_sweep_figure must draw onto an injected Figure rather than creating a new one"
