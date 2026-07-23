"""Unit tests for scripts/cif_variance.py."""
import numpy as np
import pandas as pd
import pytest

from conftest import import_script

mod = import_script("cif_variance.py")


def _write_cif(path, residues, chain_id="A"):
    """Write a minimal valid mmCIF file via Biopython's own StructureBuilder +
    MMCIFIO (round-tripping through the same library cif_variance.py uses to
    parse, via Bio.PDB.MMCIFParser) rather than hand-writing mmCIF text.

    *residues* is a list of (resname, position, [x,y,z], bfactor, hetero) --
    hetero=True marks a HETATM (e.g. water), which load_ca_data must skip.
    """
    from Bio.PDB import StructureBuilder, MMCIFIO

    sb = StructureBuilder.StructureBuilder()
    sb.init_structure(path.stem)
    sb.init_model(0)
    sb.init_chain(chain_id)
    sb.init_seg(" ")
    for resname, pos, coord, bfactor, hetero in residues:
        hetfield = "H" if hetero else " "
        sb.init_residue(resname, hetfield, pos, " ")
        sb.init_atom("CA", np.array(coord, dtype=float), bfactor, 1.0, " ", "CA", pos, element="C")
    structure = sb.get_structure()

    io = MMCIFIO()
    io.set_structure(structure)
    io.save(str(path))


def _simple_residues(positions, x_offset=0.0):
    """Build a straight-line chain of Alanines at the given positions, each
    1A apart along x (optionally shifted), with a fixed pLDDT of 80.
    """
    return [
        ("ALA", p, [float(i) + x_offset, 0.0, 0.0], 80.0, False)
        for i, p in enumerate(positions)
    ]


class TestLoadCaData:
    def test_extracts_positions_coords_plddt_and_residue_names(self, tmp_path):
        cif = tmp_path / "model.cif"
        _write_cif(cif, [
            ("ALA", 1, [0.0, 0.0, 0.0], 55.5, False),
            ("SER", 2, [1.0, 0.0, 0.0], 70.0, False),
        ])
        positions, coords, plddts, names = mod.load_ca_data(cif)

        assert positions == [1, 2], f"residue positions should be extracted in order, got {positions}"
        assert coords.shape == (2, 3), f"coords should be an (N,3) array, got shape {coords.shape}"
        assert list(plddts) == pytest.approx([55.5, 70.0]), (
            f"pLDDT should come from the CA atom's B-factor column, got {list(plddts)}"
        )
        assert names == ["A", "S"], f"residue names should be converted to one-letter codes, got {names}"

    def test_hetero_residues_are_skipped(self, tmp_path):
        cif = tmp_path / "model.cif"
        _write_cif(cif, [
            ("ALA", 1, [0.0, 0.0, 0.0], 80.0, False),
            ("HOH", 2, [9.0, 9.0, 9.0], 0.0, True),  # a water molecule, hetero=True
        ])
        positions, coords, plddts, names = mod.load_ca_data(cif)

        assert positions == [1], (
            f"a HETATM record (e.g. a water molecule) must be excluded via the "
            f"residue.get_id()[0] != ' ' check, got positions {positions}"
        )


class TestAlignToReference:
    def test_aligns_translated_structure_to_near_zero_rmsd(self):
        ref = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0], [3.0, 0.0, 0.0]])
        mobile = ref + np.array([5.0, 5.0, 5.0])  # pure translation, same shape
        positions = [1, 2, 3, 4]

        aligned, rmsd = mod.align_to_reference(ref, mobile, positions, positions)
        assert rmsd == pytest.approx(0.0, abs=1e-6), (
            f"a purely translated copy of the same shape should superimpose to ~0 RMSD, got {rmsd}"
        )
        assert np.allclose(aligned, ref, atol=1e-6), (
            "after alignment, the mobile coordinates should coincide with the reference"
        )

    def test_fewer_than_three_shared_positions_returns_inf_unchanged(self):
        ref = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]])
        mobile = np.array([[9.0, 9.0, 9.0], [8.0, 8.0, 8.0]])
        # Only position 1 is shared between ref (1,2) and mobile (1,5) -- below the
        # minimum of 3 needed to define a rigid-body superposition.
        aligned, rmsd = mod.align_to_reference(ref, mobile, [1, 2], [1, 5])

        assert rmsd == float("inf"), (
            f"fewer than 3 shared positions can't define a unique rigid-body alignment "
            f"-- must signal this via rmsd=inf rather than attempting a degenerate fit, got {rmsd}"
        )
        assert np.array_equal(aligned, mobile), (
            "with too few shared positions, the mobile coordinates should be returned unchanged"
        )


class TestIterativeAverageAlignment:
    def test_identical_structures_converge_with_near_zero_shift(self):
        coords = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0], [3.0, 0.0, 0.0]])
        all_coords = [coords.copy(), coords.copy(), coords.copy()]
        all_positions = [[1, 2, 3, 4]] * 3

        result = mod.iterative_average_alignment(all_coords, all_positions, log_cb=lambda *_: None)
        assert len(result) == 3, f"one aligned coordinate array should be returned per input structure, got {len(result)}"
        for aligned in result:
            assert np.allclose(aligned, coords, atol=1e-4), (
                "three IDENTICAL structures should align to themselves with negligible "
                "movement -- the average reference is already exactly each structure"
            )

    def test_restricting_align_positions_still_transforms_all_coordinates(self):
        # Structure 1 shifted by a pure translation; align only on positions {1,2},
        # but confirm ALL 4 residues (including 3,4 outside the align set) get moved.
        base = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0], [3.0, 0.0, 0.0]])
        shifted = base + np.array([10.0, 0.0, 0.0])
        all_coords = [base.copy(), shifted.copy()]
        all_positions = [[1, 2, 3, 4], [1, 2, 3, 4]]

        result = mod.iterative_average_alignment(
            all_coords, all_positions, align_positions={1, 2}, log_cb=lambda *_: None,
        )
        # After alignment to a shared average, structure 2's residue 4 (outside the
        # align set) must have moved from its original (unaligned) coordinate.
        assert not np.allclose(result[1][3], shifted[3]), (
            "the rotation/translation computed from the align_positions subset must be "
            "applied to ALL coordinates, including ones outside that subset -- residue 4 "
            "should have moved from its pre-alignment position, not been left untouched"
        )


class TestComputePairwiseRmsd:
    def test_diagonal_is_zero_and_matrix_is_symmetric(self):
        coords_a = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]])
        coords_b = coords_a + np.array([3.0, 0.0, 0.0])
        positions = [1, 2, 3]

        df = mod.compute_pairwise_rmsd([coords_a, coords_b], [positions, positions], ["struct_a", "struct_b"])

        assert list(df.index) == ["struct_a", "struct_b"] and list(df.columns) == ["struct_a", "struct_b"], (
            f"the RMSD matrix should be indexed/labeled by the given structure names, got {df.index}/{df.columns}"
        )
        assert df.loc["struct_a", "struct_a"] == 0.0, "a structure's RMSD against itself must be exactly 0"
        assert df.loc["struct_a", "struct_b"] == df.loc["struct_b", "struct_a"], (
            "the RMSD matrix must be symmetric (dist(A,B) == dist(B,A))"
        )
        assert df.loc["struct_a", "struct_b"] == pytest.approx(0.0, abs=1e-6), (
            "a pure translation should align to ~0 RMSD, same as the align_to_reference test"
        )


class TestLoadPtmAndMutationPositions:
    def test_parses_positions_from_tsv(self, tmp_path, monkeypatch):
        tsv = tmp_path / "hotspots.tsv"
        pd.DataFrame([{
            "uniprot_id": "P04637",
            "ptms_on_protein": "S15:Phosphorylation; T18:Phosphorylation",
            "mutations_on_protein": "R175H (10); R273H (5)",
        }]).to_csv(tsv, sep="\t", index=False)
        monkeypatch.setattr(mod, "PTM_TSV", tsv)

        ptm_positions, mutation_positions = mod.load_ptm_and_mutation_positions("P04637")
        assert ptm_positions == {15, 18}, f"PTM positions should be parsed from ptms_on_protein, got {ptm_positions}"
        assert mutation_positions == {175, 273}, (
            f"mutation positions should be parsed from mutations_on_protein, got {mutation_positions}"
        )

    def test_missing_file_returns_empty_sets(self, tmp_path, monkeypatch):
        monkeypatch.setattr(mod, "PTM_TSV", tmp_path / "does_not_exist.tsv")
        ptm_positions, mutation_positions = mod.load_ptm_and_mutation_positions("P04637")
        assert ptm_positions == set() and mutation_positions == set(), (
            "with no intermediate TSV available yet, both position sets should be "
            "empty (cross-referencing simply skipped), not raise FileNotFoundError"
        )

    def test_unknown_protein_returns_empty_sets(self, tmp_path, monkeypatch):
        tsv = tmp_path / "hotspots.tsv"
        pd.DataFrame([{"uniprot_id": "P04637", "ptms_on_protein": "S15:Phosphorylation", "mutations_on_protein": ""}]).to_csv(tsv, sep="\t", index=False)
        monkeypatch.setattr(mod, "PTM_TSV", tsv)

        ptm_positions, mutation_positions = mod.load_ptm_and_mutation_positions("Q99999")
        assert ptm_positions == set() and mutation_positions == set(), (
            f"a UniProt ID not present in the TSV has no rows to cross-reference -- "
            f"expected empty sets, got ptm={ptm_positions} mut={mutation_positions}"
        )


class TestRunVarianceAnalysis:
    def test_raises_when_fewer_than_two_cif_files(self, tmp_path):
        input_dir = tmp_path / "cifs"
        input_dir.mkdir()
        _write_cif(input_dir / "only_one.cif", _simple_residues([1, 2, 3]))

        with pytest.raises(ValueError):
            mod.run_variance_analysis(input_dir, output_dir=tmp_path / "out", log_cb=lambda *_: None)

    def test_full_run_with_two_identical_structures(self, tmp_path, monkeypatch):
        monkeypatch.setattr(mod, "PTM_TSV", tmp_path / "does_not_exist.tsv")  # skip cross-referencing
        input_dir = tmp_path / "cifs"
        input_dir.mkdir()
        residues = _simple_residues([1, 2, 3, 4, 5])
        _write_cif(input_dir / "seed_a.cif", residues)
        _write_cif(input_dir / "seed_b.cif", residues)

        result = mod.run_variance_analysis(
            input_dir, output_dir=tmp_path / "out", uniprot="P00000", log_cb=lambda *_: None,
        )

        assert result.shared_positions == [1, 2, 3, 4, 5], (
            f"two identical structures share all 5 positions, got {result.shared_positions}"
        )
        assert np.allclose(result.per_residue_variance, 0.0, atol=1e-4), (
            "two IDENTICAL structures should have ~zero positional variance at every "
            f"residue, got {result.per_residue_variance}"
        )
        assert (tmp_path / "out" / "pairwise_rmsd.tsv").exists(), "the RMSD matrix should be written to disk"
        assert (tmp_path / "out" / "variance_data.tsv").exists(), "the per-residue data table should be written to disk"

    def test_report_range_filters_output_positions(self, tmp_path, monkeypatch):
        monkeypatch.setattr(mod, "PTM_TSV", tmp_path / "does_not_exist.tsv")
        input_dir = tmp_path / "cifs"
        input_dir.mkdir()
        residues = _simple_residues([1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
        _write_cif(input_dir / "seed_a.cif", residues)
        _write_cif(input_dir / "seed_b.cif", residues)

        result = mod.run_variance_analysis(
            input_dir, output_dir=tmp_path / "out", range_=(3, 6), log_cb=lambda *_: None,
        )
        assert result.shared_positions == [3, 4, 5, 6], (
            f"--range should restrict the REPORTED positions to 3-6 even though all 10 "
            f"were used/available, got {result.shared_positions}"
        )


class TestBuildVarianceFigure:
    def test_returns_figure_with_two_axes(self, tmp_path, monkeypatch):
        monkeypatch.setattr(mod, "PTM_TSV", tmp_path / "does_not_exist.tsv")
        input_dir = tmp_path / "cifs"
        input_dir.mkdir()
        residues = _simple_residues([1, 2, 3, 4, 5])
        _write_cif(input_dir / "seed_a.cif", residues)
        _write_cif(input_dir / "seed_b.cif", residues)
        result = mod.run_variance_analysis(input_dir, output_dir=tmp_path / "out", log_cb=lambda *_: None)

        fig = mod.build_variance_figure(result)
        try:
            assert len(fig.axes) == 2, (
                f"the figure should have exactly 2 panels (variance + pLDDT), got {len(fig.axes)}"
            )
        finally:
            import matplotlib.pyplot as plt
            plt.close(fig)

    def test_uses_injected_figure_instead_of_creating_a_new_one(self, tmp_path, monkeypatch):
        monkeypatch.setattr(mod, "PTM_TSV", tmp_path / "does_not_exist.tsv")
        input_dir = tmp_path / "cifs"
        input_dir.mkdir()
        residues = _simple_residues([1, 2, 3, 4, 5])
        _write_cif(input_dir / "seed_a.cif", residues)
        _write_cif(input_dir / "seed_b.cif", residues)
        result = mod.run_variance_analysis(input_dir, output_dir=tmp_path / "out", log_cb=lambda *_: None)

        from matplotlib.figure import Figure
        injected = Figure()
        returned = mod.build_variance_figure(result, fig=injected)
        assert returned is injected, (
            "when a Figure is injected (the GUI-embedding path), build_variance_figure "
            "must draw onto that SAME object rather than silently creating a new one"
        )
