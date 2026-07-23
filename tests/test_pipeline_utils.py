"""Unit tests for the shared helpers in scripts/pipeline_utils.py.

Every function gets at least one positive (happy-path) test and at least one
negative/edge-case test, and every assertion states what real-world behavior
is being verified -- not just the raw values pytest already shows on failure.
"""
import json

import numpy as np
import pytest


@pytest.fixture
def mod(pipeline_utils_module):
    return pipeline_utils_module


def _write_synthetic_cif(path, res_ids, res_names, atom_names, coords,
                          b_factors=None, chain_ids=None):
    """Build a minimal valid mmCIF file via biotite itself (round-trip
    through the same library pipeline_utils.py uses to parse), rather than
    hand-writing mmCIF text -- the format is strict enough that a
    hand-rolled fixture would be its own source of bugs.
    """
    import biotite.structure as struc
    import biotite.structure.io.pdbx as pdbx

    n = len(res_ids)
    arr = struc.AtomArray(n)
    arr.coord = np.asarray(coords, dtype=float)
    arr.chain_id = np.asarray(chain_ids if chain_ids is not None else ["A"] * n)
    arr.res_id = np.asarray(res_ids)
    arr.res_name = np.asarray(res_names)
    arr.atom_name = np.asarray(atom_names)
    arr.element = np.asarray(["C"] * n)
    arr.set_annotation("b_factor", np.asarray(
        b_factors if b_factors is not None else [50.0] * n, dtype=float,
    ))

    cif = pdbx.CIFFile()
    pdbx.set_structure(cif, arr, data_block="test")
    cif.write(str(path))


# ── project_root / fmt_time ──────────────────────────────────────────────────

class TestProjectRoot:
    def test_derives_project_root_from_a_script_path(self, mod, tmp_path):
        scripts_dir = tmp_path / "myproject" / "scripts"
        scripts_dir.mkdir(parents=True)
        script_path = scripts_dir / "1_filter.py"

        result = mod.project_root(str(script_path))
        assert result == tmp_path / "myproject", (
            "project_root should strip 'scripts/<file>' to reach the project root "
            f"(parent of the scripts/ dir), got {result}"
        )


class TestFmtTime:
    def test_seconds_only_below_one_minute(self, mod):
        assert mod.fmt_time(45) == "45s", (
            "durations under 60s should render as plain seconds, not minutes"
        )

    def test_zero_seconds(self, mod):
        assert mod.fmt_time(0) == "0s", "zero duration should still format without erroring"

    def test_exactly_sixty_seconds_rolls_over_to_minutes(self, mod):
        assert mod.fmt_time(60) == "1m 00s", (
            "the 60s boundary should roll over to minutes, not display as '60s'"
        )

    def test_minutes_and_seconds_are_zero_padded(self, mod):
        assert mod.fmt_time(245) == "4m 05s", (
            "245s = 4 minutes 5 seconds, and the seconds portion must be zero-padded to 2 digits"
        )

    def test_fractional_seconds_are_truncated_not_rounded(self, mod):
        assert mod.fmt_time(59.9) == "59s", (
            "fmt_time truncates via int(), so 59.9s should show as 59s, not round up to 1m 00s"
        )


# ── input_dir / resolve_input_file ───────────────────────────────────────────

class TestInputDir:
    def test_creates_and_returns_the_expected_path(self, mod, tmp_path):
        result = mod.input_dir(tmp_path, "cosmic")
        expected = tmp_path / "data" / "input" / "cosmic"
        assert result == expected, "input_dir should compose root/data/input/<subfolder>"
        assert result.is_dir(), "input_dir must create the directory, not just return the path"

    def test_idempotent_when_directory_already_exists(self, mod, tmp_path):
        mod.input_dir(tmp_path, "ptmd")
        # Calling a second time must not raise (mkdir(exist_ok=True) is the point).
        result = mod.input_dir(tmp_path, "ptmd")
        assert result.is_dir(), (
            "calling input_dir twice for the same subfolder must not raise "
            "FileExistsError -- pipeline scripts call this on every run"
        )


class TestResolveInputFile:
    def test_returns_the_single_matching_file(self, mod, tmp_path):
        folder = tmp_path / "cosmic"
        folder.mkdir()
        target = folder / "data.tsv"
        target.write_text("x")

        result = mod.resolve_input_file(folder)
        assert result == target, "the one matching file present should be returned"

    def test_raises_filenotfounderror_when_folder_missing(self, mod, tmp_path):
        missing = tmp_path / "does_not_exist"
        with pytest.raises(FileNotFoundError):
            mod.resolve_input_file(missing)

    def test_raises_filenotfounderror_when_no_matching_extension(self, mod, tmp_path):
        folder = tmp_path / "cosmic"
        folder.mkdir()
        (folder / "readme.txt").write_text("not a data file")

        with pytest.raises(FileNotFoundError):
            mod.resolve_input_file(folder)

    def test_raises_runtimeerror_when_multiple_matches(self, mod, tmp_path):
        folder = tmp_path / "cosmic"
        folder.mkdir()
        (folder / "a.tsv").write_text("x")
        (folder / "b.tsv").write_text("y")

        with pytest.raises(RuntimeError):
            mod.resolve_input_file(folder)

    def test_respects_custom_extensions_argument(self, mod, tmp_path):
        folder = tmp_path / "custom"
        folder.mkdir()
        (folder / "data.tsv").write_text("x")  # should be ignored
        target = folder / "data.json"
        target.write_text("{}")

        result = mod.resolve_input_file(folder, extensions=(".json",))
        assert result == target, (
            "passing a custom extensions tuple should restrict matching to only those "
            "extensions, ignoring files with the default extensions"
        )


# ── validate_*_file ───────────────────────────────────────────────────────────

class TestValidateColumns:
    def test_no_problems_when_all_required_columns_present(self, mod, tmp_path):
        path = tmp_path / "cosmic.tsv"
        path.write_text("GENE_SYMBOL\tMUTATION_AA\tCOSMIC_SAMPLE_ID\t"
                         "MUTATION_SOMATIC_STATUS\tTRANSCRIPT_ACCESSION\n")

        problems = mod.validate_cosmic_file(path)
        assert problems == [], (
            f"a file with every required COSMIC column should validate cleanly, got {problems}"
        )

    def test_reports_each_missing_column_by_name(self, mod, tmp_path):
        path = tmp_path / "cosmic.tsv"
        path.write_text("GENE_SYMBOL\tMUTATION_AA\n")  # missing 3 required columns

        problems = mod.validate_cosmic_file(path)
        assert len(problems) == 3, (
            f"exactly the 3 missing required columns should each produce one problem, "
            f"got {len(problems)}: {problems}"
        )
        assert any("COSMIC_SAMPLE_ID" in p for p in problems), (
            "the missing-column message should name the specific column, not just say "
            f"'columns missing' -- got {problems}"
        )

    def test_ptmd_validator_checks_ptmd_specific_columns(self, mod, tmp_path):
        path = tmp_path / "ptmd.tsv"
        path.write_text("State\tUniProt\tDisease\tMutationSite\tResidue\tPosition\tType\n")

        assert mod.validate_ptmd_file(path) == [], (
            "validate_ptmd_file must check PTMD_REQUIRED_COLUMNS, not accidentally "
            "reuse COSMIC's required-columns list"
        )

    def test_1433_validator_reads_as_excel_not_tsv(self, mod, tmp_path):
        # A real .xlsx written with pandas, containing the required columns.
        pd = pytest.importorskip("pandas")
        path = tmp_path / "interactors.xlsx"
        pd.DataFrame({"Residue": ["S123"], "PMID": ["12345"]}).to_excel(path, index=False)

        assert mod.validate_1433_file(path) == [], (
            "validate_1433_file must parse the input as an Excel spreadsheet "
            "(is_excel=True), not fall back to TSV parsing"
        )

    def test_unreadable_file_is_reported_not_raised(self, mod, tmp_path):
        path = tmp_path / "interactors.xlsx"
        path.write_bytes(b"this is not a real xlsx file, just garbage bytes")

        problems = mod.validate_1433_file(path)
        assert problems != [], (
            "a file that can't even be parsed should surface as a validation problem, "
            "not silently report 'no problems found'"
        )
        assert not any(isinstance(p, Exception) for p in problems), (
            "read errors must be caught and converted to a string message, not propagate "
            "as a raised exception (the caller expects a list of problem strings)"
        )


# ── extract_uniprot_from_cif ─────────────────────────────────────────────────

class TestExtractUniprotFromCif:
    def test_finds_accession_on_matching_line(self, mod, tmp_path):
        path = tmp_path / "model.cif"
        path.write_text(
            "data_test\n"
            "_ma_target_ref_db_details.db_accession   P04637\n"
        )
        assert mod.extract_uniprot_from_cif(path) == "P04637", (
            "the accession should be extracted from the "
            "'_ma_target_ref_db_details.db_accession' line's last token"
        )

    def test_returns_none_when_no_matching_line(self, mod, tmp_path):
        path = tmp_path / "model.cif"
        path.write_text("data_test\nsome_other_field   value\n")
        assert mod.extract_uniprot_from_cif(path) is None, (
            "a CIF file without the target-ref-db-details line has no accession to find "
            "and must return None rather than raising or returning garbage"
        )

    def test_returns_none_for_unreadable_file(self, mod, tmp_path):
        path = tmp_path / "does_not_exist.cif"
        assert mod.extract_uniprot_from_cif(path) is None, (
            "a missing/unreadable file should be caught by the except OSError clause "
            "and return None, not raise"
        )


# ── find_canonical_cif / find_canonical_cifs ─────────────────────────────────

class TestFindCanonicalCif:
    def test_picks_canonical_over_isoform_model(self, mod, tmp_path):
        uniprot_dir = tmp_path / "P12345"
        uniprot_dir.mkdir()
        (uniprot_dir / "AF-P12345-2-F1-model_v6.cif").write_text("isoform")
        (uniprot_dir / "AF-P12345-F1-model_v6.cif").write_text("canonical")

        result = mod.find_canonical_cif(uniprot_dir)
        assert result is not None and result.name == "AF-P12345-F1-model_v6.cif", (
            "the isoform-numbered file ('-2-F1-') must not be mistaken for the canonical "
            f"model; got {result}"
        )

    def test_returns_none_for_empty_directory(self, mod, tmp_path):
        uniprot_dir = tmp_path / "P99999"
        uniprot_dir.mkdir()
        assert mod.find_canonical_cif(uniprot_dir) is None, (
            "a protein directory with no CIF files at all should return None, not raise"
        )

    def test_returns_none_when_only_isoform_models_present(self, mod, tmp_path):
        uniprot_dir = tmp_path / "P55555"
        uniprot_dir.mkdir()
        (uniprot_dir / "AF-P55555-2-F1-model_v6.cif").write_text("isoform only")
        assert mod.find_canonical_cif(uniprot_dir) is None, (
            "if only isoform-numbered models are present, there is no canonical model to "
            "return -- must not fall back to an isoform file"
        )


class TestFindCanonicalCifs:
    def test_returns_all_fragments_sorted_by_fragment_number(self, mod, tmp_path):
        uniprot_dir = tmp_path / "Q00001"
        uniprot_dir.mkdir()
        # Written out of order on purpose to verify sort-by-fragment-number, not filename.
        (uniprot_dir / "AF-Q00001-F2-model_v6.cif").write_text("frag2")
        (uniprot_dir / "AF-Q00001-F10-model_v6.cif").write_text("frag10")
        (uniprot_dir / "AF-Q00001-F1-model_v6.cif").write_text("frag1")

        result = mod.find_canonical_cifs(uniprot_dir)
        names = [p.name for p in result]
        assert names == [
            "AF-Q00001-F1-model_v6.cif", "AF-Q00001-F2-model_v6.cif", "AF-Q00001-F10-model_v6.cif",
        ], (
            "fragments must sort numerically by fragment index (F1, F2, F10), not "
            f"alphabetically (which would put F10 before F2) -- got {names}"
        )

    def test_returns_empty_list_for_empty_directory(self, mod, tmp_path):
        uniprot_dir = tmp_path / "Q99999"
        uniprot_dir.mkdir()
        assert mod.find_canonical_cifs(uniprot_dir) == [], (
            "no CIF files present should give an empty list, not None or an error"
        )


# ── load_first_chain / get_plddt_map ─────────────────────────────────────────

class TestLoadFirstChain:
    def test_parses_ca_atoms_and_b_factor(self, mod, tmp_path):
        path = tmp_path / "model.cif"
        _write_synthetic_cif(
            path,
            res_ids=[1, 2, 3], res_names=["ALA", "GLY", "SER"],
            atom_names=["CA", "CA", "CA"],
            coords=[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]],
            b_factors=[50.0, 60.0, 70.0],
        )

        chain = mod.load_first_chain(path)
        assert chain is not None, "a valid, well-formed CIF should parse successfully"
        assert list(chain.res_id) == [1, 2, 3], (
            "the parsed chain's residue IDs should match what was written"
        )
        assert list(chain.b_factor) == [50.0, 60.0, 70.0], (
            "b_factor (pLDDT) must be requested/preserved as an extra field during parsing"
        )

    def test_only_returns_the_first_chain(self, mod, tmp_path):
        path = tmp_path / "model.cif"
        _write_synthetic_cif(
            path,
            res_ids=[1, 1], res_names=["ALA", "ALA"], atom_names=["CA", "CA"],
            coords=[[0.0, 0.0, 0.0], [9.0, 9.0, 9.0]],
            chain_ids=["A", "B"],
        )

        chain = mod.load_first_chain(path)
        assert set(chain.chain_id) == {"A"}, (
            "load_first_chain should only return the first chain encountered, "
            f"not mix in atoms from chain B -- got chain_ids {set(chain.chain_id)}"
        )

    def test_returns_none_for_malformed_cif(self, mod, tmp_path):
        path = tmp_path / "garbage.cif"
        path.write_text("this is not valid mmCIF syntax at all {{{")

        assert mod.load_first_chain(path) is None, (
            "an unparseable file should be caught by the except Exception clause and "
            "return None rather than propagating a parser exception"
        )

    def test_returns_none_for_missing_file(self, mod, tmp_path):
        assert mod.load_first_chain(tmp_path / "nope.cif") is None, (
            "a nonexistent file path should return None, not raise FileNotFoundError"
        )


class TestGetPlddtMap:
    def test_maps_residue_position_to_b_factor_for_ca_atoms_only(self, mod, tmp_path):
        path = tmp_path / "model.cif"
        _write_synthetic_cif(
            path,
            res_ids=[1, 1, 2], res_names=["ALA", "ALA", "GLY"],
            atom_names=["N", "CA", "CA"],  # residue 1 has an N atom too -- must be ignored
            coords=[[0.0, 0.0, 0.0], [0.1, 0.1, 0.1], [1.0, 0.0, 0.0]],
            b_factors=[10.0, 55.5, 80.0],
        )
        chain = mod.load_first_chain(path)

        result = mod.get_plddt_map(chain)
        assert result == {1: 55.5, 2: 80.0}, (
            "get_plddt_map must key by CA atoms only (using the CA b_factor, not the N "
            f"atom's), one entry per residue position -- got {result}"
        )

    def test_empty_dict_when_chain_has_no_ca_atoms(self, mod, tmp_path):
        path = tmp_path / "model.cif"
        _write_synthetic_cif(
            path, res_ids=[1], res_names=["ALA"], atom_names=["N"],
            coords=[[0.0, 0.0, 0.0]],
        )
        chain = mod.load_first_chain(path)

        assert mod.get_plddt_map(chain) == {}, (
            "a chain with no CA atoms at all should produce an empty pLDDT map, not raise"
        )


class TestGetProteinLength:
    def test_single_fragment_returns_max_residue_position(self, mod, tmp_path):
        uniprot_dir = tmp_path / "P11111"
        uniprot_dir.mkdir()
        _write_synthetic_cif(
            uniprot_dir / "AF-P11111-F1-model_v6.cif",
            res_ids=[1, 2, 3], res_names=["ALA"] * 3, atom_names=["CA"] * 3,
            coords=[[0, 0, 0], [1, 0, 0], [2, 0, 0]],
        )

        assert mod.get_protein_length(uniprot_dir) == 3, (
            "for a single-fragment protein, length should equal the highest modeled "
            "residue position"
        )

    def test_multi_fragment_positions_are_continuous_not_reset(self, mod, tmp_path):
        uniprot_dir = tmp_path / "Q22222"
        uniprot_dir.mkdir()
        _write_synthetic_cif(
            uniprot_dir / "AF-Q22222-F1-model_v6.cif",
            res_ids=[1, 2], res_names=["ALA"] * 2, atom_names=["CA"] * 2,
            coords=[[0, 0, 0], [1, 0, 0]],
        )
        _write_synthetic_cif(
            uniprot_dir / "AF-Q22222-F2-model_v6.cif",
            # AlphaFold fragment numbering continues from fragment 1, not reset to 1.
            res_ids=[1401, 1402, 1403], res_names=["ALA"] * 3, atom_names=["CA"] * 3,
            coords=[[0, 0, 0], [1, 0, 0], [2, 0, 0]],
        )

        assert mod.get_protein_length(uniprot_dir) == 1403, (
            "length must be the max across ALL fragments (continuous canonical numbering), "
            "not just the last fragment parsed or the first fragment's own max"
        )

    def test_returns_none_when_no_cifs_present(self, mod, tmp_path):
        uniprot_dir = tmp_path / "R33333"
        uniprot_dir.mkdir()
        assert mod.get_protein_length(uniprot_dir) is None, (
            "a protein with no downloaded CIF files has no known length and must return "
            "None, not 0 or raise"
        )


# ── load_pae_matrix ───────────────────────────────────────────────────────────

class TestLoadPaeMatrix:
    def test_loads_matrix_from_list_wrapped_json(self, mod, tmp_path):
        uniprot_dir = tmp_path / "P12345"
        uniprot_dir.mkdir()
        data = [{"predicted_aligned_error": [[0, 1], [1, 0]]}]
        (uniprot_dir / "AF-P12345-F1-predicted_aligned_error_v6.json").write_text(json.dumps(data))

        matrix = mod.load_pae_matrix(uniprot_dir)
        assert matrix.tolist() == [[0, 1], [1, 0]], (
            "AlphaFold's PAE JSON is sometimes a single-element list wrapping the real "
            "object -- that shape must be unwrapped, not treated as an error"
        )

    def test_loads_matrix_from_plain_dict_json(self, mod, tmp_path):
        uniprot_dir = tmp_path / "P23456"
        uniprot_dir.mkdir()
        data = {"predicted_aligned_error": [[0, 2], [2, 0]]}
        (uniprot_dir / "AF-P23456-F1-predicted_aligned_error_v6.json").write_text(json.dumps(data))

        matrix = mod.load_pae_matrix(uniprot_dir)
        assert matrix.tolist() == [[0, 2], [2, 0]], (
            "the plain-dict PAE JSON shape (not list-wrapped) must also parse correctly"
        )

    def test_returns_none_when_no_json_file_present(self, mod, tmp_path):
        uniprot_dir = tmp_path / "P99999"
        uniprot_dir.mkdir()
        assert mod.load_pae_matrix(uniprot_dir) is None, (
            "no PAE file downloaded for this protein should return None, not raise "
            "FileNotFoundError"
        )

    def test_returns_none_when_expected_key_missing(self, mod, tmp_path):
        uniprot_dir = tmp_path / "P77777"
        uniprot_dir.mkdir()
        (uniprot_dir / "AF-P77777-F1-predicted_aligned_error_v6.json").write_text(
            json.dumps({"unexpected_key": 123}),
        )
        assert mod.load_pae_matrix(uniprot_dir) is None, (
            "a PAE JSON file missing the 'predicted_aligned_error' key should return None "
            "rather than raising a KeyError"
        )


# ── permutation-test statistics (prototype primitives) ──────────────────────

class TestSamplePermutationIndices:
    def test_shape_matches_permutations_by_mutations(self, mod):
        rng = np.random.default_rng(0)
        result = mod.sample_permutation_indices(
            n_residues=20, n_mutations=5, n_permutations=100, rng=rng,
        )
        assert result.shape == (100, 5), (
            f"expected (n_permutations, n_mutations) = (100, 5), got {result.shape}"
        )

    def test_each_row_samples_without_replacement(self, mod):
        rng = np.random.default_rng(1)
        result = mod.sample_permutation_indices(
            n_residues=10, n_mutations=10, n_permutations=5, rng=rng,
        )
        for row in result:
            assert len(set(row.tolist())) == 10, (
                "sampling without replacement means every row (a simulated random "
                f"mutation set) must contain 10 distinct residue indices, got {row}"
            )

    def test_indices_stay_within_residue_range(self, mod):
        rng = np.random.default_rng(2)
        result = mod.sample_permutation_indices(
            n_residues=15, n_mutations=4, n_permutations=50, rng=rng,
        )
        assert result.min() >= 0 and result.max() < 15, (
            "sampled residue indices must fall within [0, n_residues), got range "
            f"[{result.min()}, {result.max()}]"
        )

    def test_clamps_n_mutations_to_n_residues_when_larger(self, mod):
        rng = np.random.default_rng(3)
        # More "mutations" requested than residues exist -- can't sample without
        # replacement beyond the population size, so this must clamp rather than error.
        result = mod.sample_permutation_indices(
            n_residues=3, n_mutations=10, n_permutations=2, rng=rng,
        )
        assert result.shape == (2, 3), (
            "n_mutations must be clamped down to n_residues (3) rather than raising or "
            f"producing out-of-range/duplicate indices, got shape {result.shape}"
        )


class TestPermutationPvalue:
    def test_pvalue_is_never_exactly_zero(self, mod):
        # Every null count is far below observed_count -> the "raw" fraction would be
        # 0/n, but the +1 smoothing must prevent an exact-zero p-value.
        site_coord = np.array([0.0, 0.0, 0.0])
        residue_coords = np.array([[100.0, 0.0, 0.0]] * 5)  # all far away
        sampled_idx = np.zeros((10, 1), dtype=int)  # every trial samples residue 0

        p, null_counts = mod.permutation_pvalue(
            site_coord, residue_coords, observed_count=1, cutoff=1.0, sampled_idx=sampled_idx,
        )
        assert p > 0, (
            "the +1 smoothing (Davison & Hinkley) exists specifically so a p-value is "
            f"never exactly 0 even with zero matching permutations, got p={p}"
        )
        assert p == pytest.approx(1 / 11), (
            f"with 0 of 10 null trials meeting/exceeding the observed count, "
            f"p should be (0+1)/(10+1) = 1/11, got {p}"
        )

    def test_pvalue_is_one_when_every_permutation_meets_observed(self, mod):
        site_coord = np.array([0.0, 0.0, 0.0])
        residue_coords = np.array([[0.5, 0.0, 0.0]] * 3)  # all within cutoff
        sampled_idx = np.array([[0, 1, 2]] * 4)  # every trial "hits" all 3

        p, null_counts = mod.permutation_pvalue(
            site_coord, residue_coords, observed_count=3, cutoff=1.0, sampled_idx=sampled_idx,
        )
        assert p == 1.0, (
            f"if every single permutation matches or exceeds the observed count, the "
            f"empirical p-value must be the maximum possible (1.0), got {p}"
        )


class TestBenjaminiHochberg:
    def test_identical_pvalues_produce_identical_qvalues(self, mod):
        # A textbook case: all raw BH-adjusted values come out equal, so the
        # monotonicity-enforcing cummin should leave them unchanged.
        q = mod.benjamini_hochberg([0.01, 0.02, 0.03, 0.04, 0.05])
        assert q == pytest.approx([0.05] * 5), (
            f"for p=[0.01..0.05] (m=5), each p_i*m/i equals 0.05, so all q-values should "
            f"come out equal to 0.05, got {q.tolist()}"
        )

    def test_qvalues_are_never_smaller_than_input_pvalues(self, mod):
        rng = np.random.default_rng(4)
        pvalues = rng.uniform(0, 1, size=20)
        q = mod.benjamini_hochberg(pvalues)
        assert np.all(q >= pvalues - 1e-12), (
            "BH-adjusted q-values can never be smaller than the raw p-value they came "
            "from -- that would make the correction anti-conservative"
        )

    def test_qvalues_are_monotonic_in_sorted_pvalue_order(self, mod):
        pvalues = [0.2, 0.001, 0.05, 0.5, 0.01]
        q = mod.benjamini_hochberg(pvalues)
        order = np.argsort(pvalues)
        q_sorted = q[order]
        assert np.all(np.diff(q_sorted) >= -1e-12), (
            "the step-up procedure's cumulative-min pass must guarantee q-values are "
            f"non-decreasing when read in ascending p-value order, got {q_sorted.tolist()}"
        )

    def test_qvalues_are_clipped_to_the_unit_interval(self, mod):
        q = mod.benjamini_hochberg([0.9, 0.95, 0.99])
        assert np.all((q >= 0) & (q <= 1)), (
            f"q-values must stay within [0, 1] even for large p-values close to 1, got {q.tolist()}"
        )
