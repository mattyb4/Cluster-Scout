"""Unit tests for Kinase Library (Phase 3) functions in scripts/4_annotate.py.

The real `kinase_library` package isn't installed in the test environment (and
even if it were, hitting its actual prediction model in a unit test would be
slow and non-deterministic) -- tests that need it inject a fake module into
sys.modules via the fake_kinase_library fixture below.
"""
import sys
import types

import numpy as np
import pandas as pd
import pytest

from conftest import import_script

mod = import_script("4_annotate.py")


class FakeChain:
    """Minimal stand-in for a biotite AtomArray chain supporting the
    boolean-mask __getitem__ that extract_sequence relies on
    (chain[ca_mask]), unlike the simpler FakeChain in test_3_geometry.py
    which only needs attribute access.
    """

    def __init__(self, atom_name, res_id, res_name):
        self.atom_name = np.array(atom_name, dtype=object)
        self.res_id = np.array(res_id)
        self.res_name = np.array(res_name, dtype=object)

    def __getitem__(self, mask):
        return FakeChain(self.atom_name[mask], self.res_id[mask], self.res_name[mask])

    def __len__(self):
        return len(self.res_id)


class TestKinCache:
    def test_round_trip(self, tmp_path, monkeypatch):
        monkeypatch.setattr(mod, "_KIN_CACHE_FILE", tmp_path / "kinase.tsv")
        cache = {"aaaaaaaSaaaaaaa": "CDK1(3.50,99.0%)"}
        mod._kin_save_cache(cache)
        result = mod._kin_load_cache()
        assert result == cache, f"loading immediately after saving should reproduce the same dict, got {result}"

    def test_missing_file_returns_empty_dict(self, tmp_path, monkeypatch):
        monkeypatch.setattr(mod, "_KIN_CACHE_FILE", tmp_path / "does_not_exist.tsv")
        result = mod._kin_load_cache()
        assert result == {}, f"a cache file that's never been written should load as {{}}, got {result}"


class TestExtractSequence:
    def test_maps_position_to_one_letter_code(self):
        chain = FakeChain(
            atom_name=["N", "CA", "CA"],
            res_id=[1, 1, 2],
            res_name=["SER", "SER", "THR"],
        )
        result = mod.extract_sequence(chain)
        assert result == {1: "S", 2: "T"}, (
            f"only CA atoms should contribute (the N atom at residue 1 must be ignored), "
            f"and three-letter codes converted to one-letter, got {result}"
        )

    def test_unknown_residue_name_maps_to_x(self):
        chain = FakeChain(atom_name=["CA"], res_id=[5], res_name=["XYZ"])
        result = mod.extract_sequence(chain)
        assert result == {5: "X"}, (
            f"a residue name not in AA3TO1 (a non-standard/modified residue) should "
            f"fall back to 'X', not raise a KeyError, got {result}"
        )

    def test_empty_chain_returns_empty_dict(self):
        chain = FakeChain(atom_name=[], res_id=[], res_name=[])
        result = mod.extract_sequence(chain)
        assert result == {}, f"a chain with no atoms at all should return an empty dict, got {result}"


class TestBuildKinaseWindow:
    def _pos_to_aa(self):
        # Positions 1-21, all Alanine except position 11 = Serine (the "site").
        d = {i: "A" for i in range(1, 22)}
        d[11] = "S"
        return d

    def test_builds_15mer_centered_on_site_with_lowercase_center(self):
        window = mod.build_kinase_window(self._pos_to_aa(), 11)
        assert window is not None and len(window) == 15, (
            f"a 15-mer window (±7 residues) should always be exactly 15 characters, got {window!r}"
        )
        assert window[7] == "s", (
            f"the center (phosphosite) residue must be lowercased to distinguish it from "
            f"context residues, got window={window!r}"
        )
        assert window == "A" * 7 + "s" + "A" * 7, f"expected all-A context around lowercase s, got {window!r}"

    def test_returns_none_for_non_phospho_residue(self):
        pos_to_aa = {10: "A"}  # Alanine can't be phosphorylated
        result = mod.build_kinase_window(pos_to_aa, 10)
        assert result is None, (
            f"a site residue that isn't S/T/Y has no phosphorylation window to build, "
            f"got {result!r}"
        )

    def test_returns_none_when_site_position_missing(self):
        result = mod.build_kinase_window({}, 100)
        assert result is None, (
            f"a site position with no residue at all in pos_to_aa should return None, "
            f"not raise a KeyError, got {result!r}"
        )

    def test_fills_missing_boundary_positions_with_underscore(self):
        # Site at position 3 -- positions -3..0 don't exist (protein starts at 1).
        pos_to_aa = {1: "A", 2: "A", 3: "S", 4: "A", 5: "A"}
        window = mod.build_kinase_window(pos_to_aa, 3)
        assert window is not None, "a site near the sequence start should still build a window"
        assert window[:4] == "____", (
            f"positions before residue 1 don't exist and must be filled with '_' "
            f"placeholders, got window={window!r}"
        )
        assert window[7] == "s", f"the center should still be the lowercased site residue, got {window!r}"


@pytest.fixture
def fake_kinase_library(monkeypatch):
    """Injects a fake kinase_library package (root module + .modules.data
    submodule + a Substrate class) into sys.modules, and returns the fake
    data-loader call counters so tests can assert on caching behavior.
    """
    calls = {"get_kinase_list": 0, "get_kinome_info": 0}

    def get_kinase_list(*a, **k):
        calls["get_kinase_list"] += 1
        return ["CDK1", "CDK2", "MAPK1"]

    def get_kinome_info(*a, **k):
        calls["get_kinome_info"] += 1
        return {}

    data_module = types.ModuleType("kinase_library.modules.data")
    data_module.get_kinase_list = get_kinase_list
    data_module.get_kinome_info = get_kinome_info

    modules_pkg = types.ModuleType("kinase_library.modules")
    modules_pkg.data = data_module

    class FakeSubstrate:
        def __init__(self, window):
            self.window = window

        def predict(self):
            return pd.DataFrame(
                {"Score": [3.5, 2.1, 1.0], "Percentile": [99.0, 95.0, 80.0]},
                index=["CDK1", "CDK2", "MAPK1"],
            )

    root_pkg = types.ModuleType("kinase_library")
    root_pkg.modules = modules_pkg
    root_pkg.Substrate = FakeSubstrate

    monkeypatch.setitem(sys.modules, "kinase_library", root_pkg)
    monkeypatch.setitem(sys.modules, "kinase_library.modules", modules_pkg)
    monkeypatch.setitem(sys.modules, "kinase_library.modules.data", data_module)
    return calls


class TestSpeedUpKinaseLibrary:
    def test_wraps_loaders_with_caching(self, fake_kinase_library):
        mod._speed_up_kinase_library()

        import kinase_library.modules.data as kl_data
        kl_data.get_kinase_list()
        kl_data.get_kinase_list()  # same args -> should hit the lru_cache

        assert fake_kinase_library["get_kinase_list"] == 1, (
            "after wrapping, calling get_kinase_list twice with the same arguments "
            f"should only invoke the real underlying loader once, got "
            f"{fake_kinase_library['get_kinase_list']} real calls"
        )

    def test_idempotent_when_called_twice(self, fake_kinase_library):
        mod._speed_up_kinase_library()
        import kinase_library.modules.data as kl_data
        first_wrapped = kl_data.get_kinase_list

        mod._speed_up_kinase_library()  # second call must be a no-op
        assert kl_data.get_kinase_list is first_wrapped, (
            "calling _speed_up_kinase_library a second time must not re-wrap an "
            "already-wrapped loader (that would silently discard the first cache "
            "and could double-wrap indefinitely) -- the function object should be unchanged"
        )


class TestPredictKinases:
    def test_formats_top_k_predictions(self, fake_kinase_library):
        result = mod.predict_kinases("a" * 7 + "s" + "a" * 7)
        assert result == "CDK1(3.50,99.0%); CDK2(2.10,95.0%); MAPK1(1.00,80.0%)", (
            f"predictions should be formatted as 'KINASE(score,percentile%)', "
            f"semicolon-joined in the model's own ranked order, got {result!r}"
        )


class TestRunKinasePhase:
    """MODELS_ROOT is monkeypatched to a tmp_path in every test here, with real
    (empty) directories created for proteins that should be found -- relying
    on find_canonical_cif's mock alone isn't enough, since run_kinase_phase
    checks `uniprot_dir.is_dir()` against the real filesystem path FIRST and
    only calls find_canonical_cif if that's True.
    """

    def test_skips_proteins_with_no_cif_directory(self, monkeypatch, tmp_path):
        monkeypatch.setattr(mod, "MODELS_ROOT", tmp_path)  # P99999 dir intentionally not created
        df = pd.DataFrame([
            {"UniProt": "P99999", "ptm_site": "S10", "ptm_type": "Phosphorylation"},
        ])

        seq_maps, kin_cache = mod.run_kinase_phase(df)

        assert seq_maps == {}, (
            f"a protein with no downloaded CIF directory has no sequence to extract -- "
            f"seq_maps should be empty, got {seq_maps}"
        )
        assert df["kinase_predictions"].tolist() == [""], (
            f"with no sequence available, the row's kinase_predictions must be blank, "
            f"not raise, got {df['kinase_predictions'].tolist()}"
        )

    def test_predicts_and_caches_new_windows(self, monkeypatch, tmp_path):
        monkeypatch.setattr(mod, "_KIN_CACHE_FILE", tmp_path / "kinase.tsv")
        monkeypatch.setattr(mod, "MODELS_ROOT", tmp_path)
        (tmp_path / "P04637").mkdir()
        monkeypatch.setattr(mod, "find_canonical_cif", lambda uniprot_dir: "fake.cif")
        pos_to_aa = {i: "A" for i in range(1, 22)}
        pos_to_aa[10] = "S"
        monkeypatch.setattr(mod, "load_first_chain", lambda cif_file: object())
        monkeypatch.setattr(mod, "extract_sequence", lambda chain: pos_to_aa)
        monkeypatch.setattr(mod, "predict_kinases", lambda window: "CDK1(3.50,99.0%)")

        df = pd.DataFrame([
            {"UniProt": "P04637", "ptm_site": "S10", "ptm_type": "Phosphorylation"},
        ])
        seq_maps, kin_cache = mod.run_kinase_phase(df)

        assert df["kinase_predictions"].tolist() == ["CDK1(3.50,99.0%)"], (
            f"a phosphorylation site with a valid S/T/Y residue and a resolvable "
            f"sequence window should get the (mocked) prediction, got "
            f"{df['kinase_predictions'].tolist()}"
        )
        assert "P04637" in seq_maps, "the protein's extracted sequence should be returned in seq_maps"

    def test_non_phosphorylation_ptm_type_is_left_blank(self, monkeypatch, tmp_path):
        monkeypatch.setattr(mod, "_KIN_CACHE_FILE", tmp_path / "kinase.tsv")
        monkeypatch.setattr(mod, "MODELS_ROOT", tmp_path)
        (tmp_path / "P04637").mkdir()
        monkeypatch.setattr(mod, "find_canonical_cif", lambda uniprot_dir: "fake.cif")
        monkeypatch.setattr(mod, "load_first_chain", lambda cif_file: object())
        monkeypatch.setattr(mod, "extract_sequence", lambda chain: {10: "S"})
        calls = []
        monkeypatch.setattr(mod, "predict_kinases", lambda window: calls.append(window) or "SHOULD_NOT_APPEAR")

        df = pd.DataFrame([
            {"UniProt": "P04637", "ptm_site": "S10", "ptm_type": "Ubiquitination"},
        ])
        mod.run_kinase_phase(df)

        assert df["kinase_predictions"].tolist() == [""], (
            f"kinase predictions only apply to phosphorylation sites -- a Ubiquitination "
            f"row must stay blank, got {df['kinase_predictions'].tolist()}"
        )
        assert calls == [], "predict_kinases should never even be called for a non-phosphorylation row"

    def test_shared_window_across_rows_is_predicted_only_once(self, monkeypatch, tmp_path):
        monkeypatch.setattr(mod, "_KIN_CACHE_FILE", tmp_path / "kinase.tsv")
        monkeypatch.setattr(mod, "MODELS_ROOT", tmp_path)
        (tmp_path / "P00001").mkdir()
        (tmp_path / "P00002").mkdir()
        monkeypatch.setattr(mod, "find_canonical_cif", lambda uniprot_dir: "fake.cif")
        pos_to_aa = {i: "A" for i in range(1, 22)}
        pos_to_aa[10] = "S"
        monkeypatch.setattr(mod, "load_first_chain", lambda cif_file: object())
        monkeypatch.setattr(mod, "extract_sequence", lambda chain: pos_to_aa)
        calls = []
        monkeypatch.setattr(mod, "predict_kinases", lambda window: calls.append(window) or "CDK1(1.0,50.0%)")

        # Two different proteins, but the SAME site position/sequence context ->
        # identical 15-mer window.
        df = pd.DataFrame([
            {"UniProt": "P00001", "ptm_site": "S10", "ptm_type": "Phosphorylation"},
            {"UniProt": "P00002", "ptm_site": "S10", "ptm_type": "Phosphorylation"},
        ])
        mod.run_kinase_phase(df)

        assert len(calls) == 1, (
            f"two rows that resolve to the identical sequence window should only be "
            f"predicted once (deduped before the thread pool), got {len(calls)} calls"
        )
