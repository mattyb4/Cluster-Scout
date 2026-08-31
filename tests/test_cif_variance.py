"""Unit tests for scripts/cif_variance.py."""
import json

import numpy as np
import pandas as pd
import pytest

from conftest import FakeResponse, import_script

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


class TestResolveUniprot:
    def test_explicit_uniprot_wins(self, tmp_path):
        cif = tmp_path / "model.cif"
        cif.write_text("_ma_target_ref_db_details.db_accession   Q00000\n")
        result = mod.resolve_uniprot([cif], uniprot="P00001", gene="IGNORED")
        assert result == "P00001", (
            "an explicit uniprot arg must win outright, without even reading the CIF metadata"
        )

    def test_falls_back_to_cif_metadata(self, tmp_path):
        cif = tmp_path / "model.cif"
        cif.write_text("_ma_target_ref_db_details.db_accession   P04637\n")
        result = mod.resolve_uniprot([cif], uniprot=None, gene=None)
        assert result == "P04637", (
            "with no explicit uniprot, the accession embedded in the first CIF file's "
            "metadata should be used"
        )

    def test_falls_back_to_gene_lookup(self, tmp_path, monkeypatch):
        tsv = tmp_path / "hotspots.tsv"
        pd.DataFrame([{"gene": "TP53", "uniprot_id": "P04637"}]).to_csv(tsv, sep="\t", index=False)
        monkeypatch.setattr(mod, "PTM_TSV", tsv)

        cif = tmp_path / "model.cif"
        cif.write_text("no matching metadata line\n")
        result = mod.resolve_uniprot([cif], uniprot=None, gene="tp53")
        assert result == "P04637", (
            "with no explicit uniprot and no CIF metadata, the gene override should be "
            "looked up (case-insensitively) against the pipeline's intermediate TSV"
        )

    def test_returns_none_when_nothing_resolves(self, tmp_path, monkeypatch):
        monkeypatch.setattr(mod, "PTM_TSV", tmp_path / "does_not_exist.tsv")
        result = mod.resolve_uniprot([], uniprot=None, gene=None)
        assert result is None, (
            "with no uniprot, no CIF files, and no gene match, resolution must return "
            "None rather than an empty string or raising"
        )


class TestFetchUniprotSequence:
    def test_returns_sequence_string_on_200(self, monkeypatch):
        monkeypatch.setattr(
            mod.requests, "get",
            lambda url, timeout=None: FakeResponse(">sp|P04637|P53_HUMAN\nMEEPQSDPSV\nCNTSSPQP\n"),
        )
        result = mod.fetch_uniprot_sequence("P04637")
        assert result == "MEEPQSDPSVCNTSSPQP", (
            f"the FASTA header line must be dropped and the remaining lines joined with no "
            f"newlines, got {result!r}"
        )

    def test_returns_none_on_non_200(self, monkeypatch):
        monkeypatch.setattr(
            mod.requests, "get",
            lambda url, timeout=None: FakeResponse("Not Found", status_code=404),
        )
        result = mod.fetch_uniprot_sequence("DELETED1")
        assert result is None, (
            f"a non-200 response (e.g. a withdrawn UniProt entry) must return None, got {result!r}"
        )


class TestBuildAlphafoldSeedJson:
    def test_builds_one_separate_job_per_seed(self):
        payload = mod.build_alphafold_seed_json("MEEPQSDPSV", "TP53", [1, 2, 3])
        assert payload == [
            {
                "name": "TP53_seed1",
                "modelSeeds": [1],
                "sequences": [{"proteinChain": {"sequence": "MEEPQSDPSV", "count": 1}}],
            },
            {
                "name": "TP53_seed2",
                "modelSeeds": [2],
                "sequences": [{"proteinChain": {"sequence": "MEEPQSDPSV", "count": 1}}],
            },
            {
                "name": "TP53_seed3",
                "modelSeeds": [3],
                "sequences": [{"proteinChain": {"sequence": "MEEPQSDPSV", "count": 1}}],
            },
        ], f"unexpected AlphaFold Server batch payload shape: {payload}"

    def test_top_level_is_a_list(self):
        # AlphaFold Server's own JSON dialect is detected by the top-level
        # value being a list (as opposed to alphafold3's own dict-based format).
        payload = mod.build_alphafold_seed_json("MEEPQSDPSV", "job", [1])
        assert isinstance(payload, list), (
            "the top-level payload must be a JSON list -- alphafoldserver.com detects its "
            "own dialect by this, and a dict here would be misread as the alphafold3 dialect"
        )


class TestGenerateAlphafoldSeedJson:
    def test_writes_json_file_using_explicit_uniprot(self, tmp_path, monkeypatch):
        monkeypatch.setattr(
            mod.requests, "get",
            lambda url, timeout=None: FakeResponse(">sp|P04637|P53_HUMAN\nMEEPQSDPSV\n"),
        )
        input_dir = tmp_path / "cifs"
        input_dir.mkdir()
        out_dir = tmp_path / "out"

        path = mod.generate_alphafold_seed_json(
            input_dir=input_dir, output_dir=out_dir, uniprot="P04637", gene="TP53",
            log_cb=lambda *_: None,
        )

        assert path.parent == out_dir, f"the JSON should be written inside output_dir, got {path}"
        payload = json.loads(path.read_text())
        assert len(payload) == 10, f"default seeds 1-10 must produce 10 separate jobs, got {len(payload)}"
        assert payload[0]["sequences"][0]["proteinChain"]["sequence"] == "MEEPQSDPSV", (
            "each written job should embed the sequence fetched from UniProt"
        )
        assert [job["modelSeeds"] for job in payload] == [[s] for s in mod.DEFAULT_SEEDS] == [[s] for s in range(1, 11)], (
            f"each job must carry exactly one seed, covering 1-10 across the 10 jobs, "
            f"got {[job['modelSeeds'] for job in payload]}"
        )
        assert [job["name"] for job in payload] == [f"TP53_seed{s}" for s in range(1, 11)], (
            f"job names should be derived from the gene override with a per-seed suffix, "
            f"got {[job['name'] for job in payload]}"
        )

    def test_custom_seed_count_produces_that_many_jobs(self, tmp_path, monkeypatch):
        monkeypatch.setattr(
            mod.requests, "get",
            lambda url, timeout=None: FakeResponse(">sp|P04637|P53_HUMAN\nMEEPQSDPSV\n"),
        )
        input_dir = tmp_path / "cifs"
        input_dir.mkdir()

        path = mod.generate_alphafold_seed_json(
            input_dir=input_dir, output_dir=tmp_path / "out", uniprot="P04637", gene="TP53",
            seeds=list(range(1, 21)), log_cb=lambda *_: None,
        )

        payload = json.loads(path.read_text())
        assert len(payload) == 20, f"a custom seed count should produce exactly that many jobs, got {len(payload)}"
        assert "seeds1-20" in path.name, (
            f"the output filename should reflect the actual seed range requested, got {path.name!r}"
        )

    def test_empty_seeds_raises_value_error(self, tmp_path, monkeypatch):
        input_dir = tmp_path / "cifs"
        input_dir.mkdir()
        with pytest.raises(ValueError):
            mod.generate_alphafold_seed_json(
                input_dir=input_dir, output_dir=tmp_path / "out", uniprot="P04637",
                seeds=[], log_cb=lambda *_: None,
            )

    def test_raises_when_uniprot_cannot_be_resolved(self, tmp_path, monkeypatch):
        monkeypatch.setattr(mod, "PTM_TSV", tmp_path / "does_not_exist.tsv")
        input_dir = tmp_path / "cifs"
        input_dir.mkdir()

        with pytest.raises(ValueError):
            mod.generate_alphafold_seed_json(
                input_dir=input_dir, output_dir=tmp_path / "out", log_cb=lambda *_: None,
            )

    def test_raises_when_sequence_fetch_fails(self, tmp_path, monkeypatch):
        monkeypatch.setattr(
            mod.requests, "get",
            lambda url, timeout=None: FakeResponse("Not Found", status_code=404),
        )
        input_dir = tmp_path / "cifs"
        input_dir.mkdir()

        with pytest.raises(ValueError):
            mod.generate_alphafold_seed_json(
                input_dir=input_dir, output_dir=tmp_path / "out", uniprot="Q99999",
                log_cb=lambda *_: None,
            )
