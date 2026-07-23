"""Unit tests for the AlphaFold model/PAE file lookup helpers in scripts/3_find_nearby_mutations.py."""
import json

import pytest


@pytest.fixture
def mod(nearby_module):
    return nearby_module


class TestFindModelFile:
    def test_picks_canonical_cif_over_isoform(self, mod, tmp_path):
        uniprot_dir = tmp_path / "P12345"
        uniprot_dir.mkdir()
        (uniprot_dir / "AF-P12345-2-F1-model_v6.cif").write_text("isoform")
        (uniprot_dir / "AF-P12345-F1-model_v6.cif").write_text("canonical")

        result = mod.find_model_file(uniprot_dir)
        assert result.name == "AF-P12345-F1-model_v6.cif", (
            f"the isoform-numbered file ('-2-F1-') must not be mistaken for the canonical "
            f"model, got {result.name if result else None}"
        )

    def test_returns_none_if_only_isoform_models(self, mod, tmp_path):
        uniprot_dir = tmp_path / "P99999"
        uniprot_dir.mkdir()
        (uniprot_dir / "AF-P99999-2-F1-model_v6.cif").write_text("isoform")

        assert mod.find_model_file(uniprot_dir) is None, (
            "with only an isoform-numbered model present, there is no canonical model to "
            "return -- must not fall back to the isoform file"
        )

    def test_returns_none_for_empty_directory(self, mod, tmp_path):
        uniprot_dir = tmp_path / "P00000"
        uniprot_dir.mkdir()

        assert mod.find_model_file(uniprot_dir) is None, (
            "a protein directory with no CIF files at all should return None, not raise"
        )


class TestLoadPaeMatrix:
    def test_loads_matrix_from_list_wrapped_json(self, mod, tmp_path):
        uniprot_dir = tmp_path / "P12345"
        uniprot_dir.mkdir()
        data = [{"predicted_aligned_error": [[0, 1], [1, 0]]}]
        (uniprot_dir / "AF-P12345-F1-predicted_aligned_error_v6.json").write_text(json.dumps(data))

        matrix = mod.load_pae_matrix(uniprot_dir)
        assert matrix.tolist() == [[0, 1], [1, 0]], (
            f"AlphaFold's list-wrapped PAE JSON shape must be unwrapped correctly, got {matrix.tolist()}"
        )

    def test_loads_matrix_from_plain_dict_json(self, mod, tmp_path):
        uniprot_dir = tmp_path / "P23456"
        uniprot_dir.mkdir()
        data = {"predicted_aligned_error": [[0, 2], [2, 0]]}
        (uniprot_dir / "AF-P23456-F1-predicted_aligned_error_v6.json").write_text(json.dumps(data))

        matrix = mod.load_pae_matrix(uniprot_dir)
        assert matrix.tolist() == [[0, 2], [2, 0]], (
            f"the plain-dict PAE JSON shape must also parse correctly, got {matrix.tolist()}"
        )

    def test_returns_none_when_file_missing(self, mod, tmp_path):
        uniprot_dir = tmp_path / "P99999"
        uniprot_dir.mkdir()

        assert mod.load_pae_matrix(uniprot_dir) is None, (
            "no PAE file downloaded for this protein should return None, not raise FileNotFoundError"
        )
