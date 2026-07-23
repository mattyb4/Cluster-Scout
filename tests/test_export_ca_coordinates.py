"""Unit tests for scripts/export_ca_coordinates.py."""
import numpy as np
import pandas as pd
import pytest

from conftest import import_script, FakeResponse

mod = import_script("export_ca_coordinates.py")


def _write_synthetic_cif(path, res_ids, res_names, atom_names, coords, chain_ids=None):
    """Same biotite-round-trip technique as test_pipeline_utils.py's fixture --
    export_ca_coordinates.py's _load_ca_from_cif uses pipeline_utils.load_first_chain
    (biotite-based), not Bio.PDB, so this mirrors that file's helper rather than
    cif_variance.py's Bio.PDB-based one.
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
    arr.set_annotation("b_factor", np.asarray([80.0] * n, dtype=float))

    cif = pdbx.CIFFile()
    pdbx.set_structure(cif, arr, data_block="test")
    cif.write(str(path))


class FakeStreamResponse:
    def __init__(self, status_code=200, chunks=None):
        self.status_code = status_code
        self._chunks = chunks or [b"cif-file-bytes"]

    def iter_content(self, chunk_size=None):
        yield from self._chunks

    def __enter__(self):
        return self

    def __exit__(self, *a):
        return False


class _JsonResp:
    def __init__(self, data, status_code=200):
        self._data = data
        self.status_code = status_code

    def raise_for_status(self):
        if self.status_code >= 400:
            raise RuntimeError(f"HTTP {self.status_code}")

    def json(self):
        return self._data


class TestDownloadCif:
    def test_downloads_canonical_record(self, tmp_path, monkeypatch):
        monkeypatch.setattr(mod, "MODELS_ROOT", tmp_path)
        record = {"uniprotAccession": "P04637", "cifUrl": "https://example.org/AF-P04637-F1-model_v4.cif"}
        monkeypatch.setattr(mod.requests, "get", lambda url, timeout=None: _JsonResp([record]))

        class FakeSession:
            def __init__(self, *a, **k): self.headers = {}
            def update(self, *a, **k): pass
            def __enter__(self): return self
            def __exit__(self, *a): return False
            def get(self, url, stream=True, timeout=None):
                return FakeStreamResponse(200, [b"data"])

        monkeypatch.setattr(mod.requests, "Session", FakeSession)

        downloaded = mod._download_cif("P04637", log_cb=lambda *_: None)
        assert len(downloaded) == 1, f"one canonical record with a cifUrl should produce one downloaded file, got {len(downloaded)}"
        assert downloaded[0].exists(), "the downloaded file should actually be written to disk"

    def test_raises_value_error_on_404(self, tmp_path, monkeypatch):
        monkeypatch.setattr(mod, "MODELS_ROOT", tmp_path)
        monkeypatch.setattr(mod.requests, "get", lambda url, timeout=None: FakeResponse("Not Found", status_code=404))

        with pytest.raises(ValueError):
            mod._download_cif("Q99999999", log_cb=lambda *_: None)

    def test_raises_value_error_when_only_isoform_records(self, tmp_path, monkeypatch):
        monkeypatch.setattr(mod, "MODELS_ROOT", tmp_path)
        record = {"uniprotAccession": "P04637-2", "cifUrl": "https://example.org/isoform.cif"}
        monkeypatch.setattr(mod.requests, "get", lambda url, timeout=None: _JsonResp([record]))

        with pytest.raises(ValueError):
            mod._download_cif("P04637", log_cb=lambda *_: None)

    def test_skips_already_downloaded_nonempty_file(self, tmp_path, monkeypatch):
        monkeypatch.setattr(mod, "MODELS_ROOT", tmp_path)
        out_dir = tmp_path / "P04637"
        out_dir.mkdir()
        existing = out_dir / "AF-P04637-F1-model_v4.cif"
        existing.write_text("already here")

        record = {"uniprotAccession": "P04637", "cifUrl": "https://example.org/AF-P04637-F1-model_v4.cif"}
        monkeypatch.setattr(mod.requests, "get", lambda url, timeout=None: _JsonResp([record]))

        session_get_calls = []

        class FakeSession:
            def __init__(self, *a, **k): self.headers = {}
            def update(self, *a, **k): pass
            def __enter__(self): return self
            def __exit__(self, *a): return False
            def get(self, url, stream=True, timeout=None):
                session_get_calls.append(url)
                return FakeStreamResponse(200)

        monkeypatch.setattr(mod.requests, "Session", FakeSession)

        downloaded = mod._download_cif("P04637", log_cb=lambda *_: None)
        assert session_get_calls == [], (
            f"an already-downloaded non-empty file should be skipped -- no download "
            f"request should be made, got {len(session_get_calls)}"
        )
        assert downloaded == [existing], "the existing file should still be reported as available"


class TestLoadCaFromCif:
    def test_extracts_ca_rows(self, tmp_path):
        cif = tmp_path / "model.cif"
        _write_synthetic_cif(
            cif, res_ids=[1, 2], res_names=["ALA", "SER"], atom_names=["CA", "CA"],
            coords=[[0.0, 0.0, 0.0], [1.5, 0.0, 0.0]],
        )
        rows = mod._load_ca_from_cif(cif)
        assert rows == [
            {"residue": "A", "position": 1, "x": 0.0, "y": 0.0, "z": 0.0},
            {"residue": "S", "position": 2, "x": 1.5, "y": 0.0, "z": 0.0},
        ], f"each CA atom should produce a {{residue,position,x,y,z}} dict with coords rounded to 3dp, got {rows}"

    def test_returns_empty_list_when_chain_unparseable(self, tmp_path):
        cif = tmp_path / "garbage.cif"
        cif.write_text("not valid mmCIF at all")
        rows = mod._load_ca_from_cif(cif)
        assert rows == [], f"an unparseable CIF should return [] (load_first_chain returns None), not raise, got {rows}"


class TestLookupGene:
    def test_returns_gene_from_cache(self, tmp_path, monkeypatch):
        cache = tmp_path / "gene_cache.tsv"
        pd.DataFrame([{"UniProt": "P04637", "gene": "TP53"}]).to_csv(cache, sep="\t", index=False)
        monkeypatch.setattr(mod, "GENE_CACHE", cache)

        result = mod._lookup_gene("P04637", log_cb=lambda *_: None)
        assert result == "TP53", f"a cache hit should return the cached gene without an API call, got {result!r}"

    def test_falls_back_to_api_when_not_in_cache(self, tmp_path, monkeypatch):
        monkeypatch.setattr(mod, "GENE_CACHE", tmp_path / "does_not_exist.tsv")
        monkeypatch.setattr(
            mod.requests, "get",
            lambda url, params=None, timeout=None: FakeResponse("Gene Names\tProtein names\nTP53 p53\tCellular tumor antigen p53\n"),
        )
        result = mod._lookup_gene("P04637", log_cb=lambda *_: None)
        assert result == "TP53", f"an uncached accession should fall back to the live UniProt API, got {result!r}"

    def test_raises_on_deleted_entry(self, tmp_path, monkeypatch):
        monkeypatch.setattr(mod, "GENE_CACHE", tmp_path / "does_not_exist.tsv")
        monkeypatch.setattr(
            mod.requests, "get",
            lambda url, params=None, timeout=None: FakeResponse("Gene Names\tProtein names\n\tdeleted\n"),
        )
        with pytest.raises(ValueError):
            mod._lookup_gene("P00000", log_cb=lambda *_: None)

    def test_raises_when_no_gene_symbol_in_response(self, tmp_path, monkeypatch):
        monkeypatch.setattr(mod, "GENE_CACHE", tmp_path / "does_not_exist.tsv")
        monkeypatch.setattr(
            mod.requests, "get",
            lambda url, params=None, timeout=None: FakeResponse("Gene Names\tProtein names\n\tSome Protein\n"),
        )
        with pytest.raises(ValueError):
            mod._lookup_gene("P00000", log_cb=lambda *_: None)

    def test_returns_none_on_network_error(self, tmp_path, monkeypatch):
        monkeypatch.setattr(mod, "GENE_CACHE", tmp_path / "does_not_exist.tsv")
        monkeypatch.setattr(
            mod.requests, "get",
            lambda url, params=None, timeout=None: (_ for _ in ()).throw(mod.requests.RequestException("timeout")),
        )
        result = mod._lookup_gene("P04637", log_cb=lambda *_: None)
        assert result is None, f"a network error should be logged and return None, not raise, got {result!r}"


class TestLookupUniprotFromGene:
    def test_returns_accession(self, monkeypatch):
        monkeypatch.setattr(
            mod.requests, "get",
            lambda url, params=None, timeout=None: FakeResponse("Entry\nP04637\n"),
        )
        result = mod._lookup_uniprot_from_gene("TP53", log_cb=lambda *_: None)
        assert result == "P04637", f"a successful gene lookup should return the accession, got {result!r}"

    def test_returns_none_when_no_results(self, monkeypatch):
        monkeypatch.setattr(
            mod.requests, "get",
            lambda url, params=None, timeout=None: FakeResponse("Entry\n"),
        )
        result = mod._lookup_uniprot_from_gene("NOTAGENE", log_cb=lambda *_: None)
        assert result is None, f"a gene with no reviewed human match should return None, got {result!r}"

    def test_returns_none_on_network_error(self, monkeypatch):
        monkeypatch.setattr(
            mod.requests, "get",
            lambda url, params=None, timeout=None: (_ for _ in ()).throw(mod.requests.RequestException("timeout")),
        )
        result = mod._lookup_uniprot_from_gene("TP53", log_cb=lambda *_: None)
        assert result is None, f"a network error should return None, not raise, got {result!r}"


class TestLoadCosmicMutations:
    def _cosmic_file(self, tmp_path, rows):
        path = tmp_path / "cosmic.tsv"
        pd.DataFrame(rows, columns=[
            "GENE_SYMBOL", "MUTATION_AA", "COSMIC_SAMPLE_ID", "MUTATION_SOMATIC_STATUS",
        ]).to_csv(path, sep="\t", index=False)
        return path

    def test_aggregates_patients_per_position(self, tmp_path):
        cosmic = self._cosmic_file(tmp_path, [
            ("TP53", "p.R175H", "S1", "Confirmed somatic variant"),
            ("TP53", "p.R175H", "S2", "Confirmed somatic variant"),
            ("TP53", "p.R273H", "S3", "Confirmed somatic variant"),
            ("PTPN11", "p.E76A", "S4", "Confirmed somatic variant"),  # different gene, excluded
        ])
        pos_mutations, pos_patients = mod._load_cosmic_mutations("TP53", cosmic, log_cb=lambda *_: None)
        assert pos_mutations == {175: ["R175H"], 273: ["R273H"]}, (
            f"only TP53 rows should contribute, grouped by position, got {pos_mutations}"
        )
        assert pos_patients == {175: 2, 273: 1}, (
            f"position 175 has 2 distinct samples (S1,S2), position 273 has 1, got {pos_patients}"
        )

    def test_excludes_non_somatic_and_non_missense(self, tmp_path):
        cosmic = self._cosmic_file(tmp_path, [
            ("TP53", "p.R175H", "S1", "Variant of unknown origin"),  # not somatic
            ("TP53", "p.E11*", "S2", "Confirmed somatic variant"),  # stop-codon, not simple substitution
        ])
        pos_mutations, pos_patients = mod._load_cosmic_mutations("TP53", cosmic, log_cb=lambda *_: None)
        assert pos_mutations == {} and pos_patients == {}, (
            f"a non-somatic row and a stop-codon row should both be excluded, got "
            f"pos_mutations={pos_mutations} pos_patients={pos_patients}"
        )


class TestComputePatientsWithinRadius:
    def _ca_df(self):
        return pd.DataFrame([
            {"position": 1, "x": 0.0, "y": 0.0, "z": 0.0},
            {"position": 2, "x": 5.0, "y": 0.0, "z": 0.0},
            {"position": 3, "x": 20.0, "y": 0.0, "z": 0.0},
        ])

    def test_sums_patients_within_radius(self):
        result = mod._compute_patients_within_radius(self._ca_df(), {1: 10}, radius=10.0)
        assert result[1] == 10, "residue 1 has a mutation at its own position (distance 0) -- must count toward itself"
        assert result[2] == 10, "residue 2 (5A from residue 1) is within the 10A radius -- should also count"
        assert result[3] == 0, "residue 3 (20A from residue 1) is outside the 10A radius -- should be 0"

    def test_empty_pos_patients_returns_all_zeros(self):
        result = mod._compute_patients_within_radius(self._ca_df(), {}, radius=10.0)
        assert result == {1: 0, 2: 0, 3: 0}, (
            f"with no mutation data at all, every residue should report 0, not raise, got {result}"
        )

    def test_mutation_position_not_in_ca_df_returns_all_zeros(self):
        # pos_patients references a position that doesn't exist in ca_df at all
        # (e.g. a mutation beyond the modeled structure's range).
        result = mod._compute_patients_within_radius(self._ca_df(), {999: 5}, radius=10.0)
        assert result == {1: 0, 2: 0, 3: 0}, (
            f"a mutation position absent from ca_df has no coordinate to measure "
            f"distance from -- every residue should report 0, not raise, got {result}"
        )

    def test_radius_boundary_is_inclusive(self):
        ca_df = pd.DataFrame([
            {"position": 1, "x": 0.0, "y": 0.0, "z": 0.0},
            {"position": 2, "x": 10.0, "y": 0.0, "z": 0.0},  # exactly 10A away
        ])
        result = mod._compute_patients_within_radius(ca_df, {1: 7}, radius=10.0)
        assert result[2] == 7, (
            f"a residue exactly AT the radius boundary (10A == 10A) must be included "
            f"(inclusive <=, matching the pipeline's convention elsewhere), got {result[2]}"
        )


class TestWriteDefattrFile:
    def test_format_matches_chimerax_spec(self, tmp_path):
        ca_df = pd.DataFrame([
            {"position": 1, "patients_within_10A": 5},
            {"position": 2, "patients_within_10A": 0},
        ])
        out_path = mod.write_defattr_file(ca_df, tmp_path / "mutations.defattr", chain_id="A")
        text = out_path.read_text(encoding="utf-8")
        lines = text.split("\n")

        assert lines[0] == "attribute: patients_within_10A", f"the first line should declare the attribute name, got {lines[0]!r}"
        assert lines[1] == "recipient: residues", f"the second line should declare the recipient type, got {lines[1]!r}"
        data_line = next(l for l in lines if "1" in l and l.startswith("\t"))
        assert data_line == "\t/A:1\t5", (
            f"each data line needs a LEADING tab and a '/' before the chain letter "
            f"(ChimeraX's atom-spec grammar requires it) -- got {data_line!r}"
        )
        assert "\r" not in text, "the file must use LF-only line endings (matching ChimeraX's own shipped files), not CRLF"


class TestWriteChimeraxScript:
    def test_includes_range_clause_when_given(self, tmp_path):
        out_path = mod.write_chimerax_script(
            tmp_path / "model.cif", tmp_path / "mutations.defattr", tmp_path / "view.cxc",
            value_range=(0, 42),
        )
        text = out_path.read_text(encoding="utf-8")
        assert "range 0,42" in text, f"an explicit value_range should appear as a 'range MIN,MAX' clause, got:\n{text}"

    def test_omits_range_clause_when_not_given(self, tmp_path):
        out_path = mod.write_chimerax_script(
            tmp_path / "model.cif", tmp_path / "mutations.defattr", tmp_path / "view.cxc",
        )
        text = out_path.read_text(encoding="utf-8")
        color_line = next(l for l in text.splitlines() if l.startswith("color byattribute"))
        assert color_line.endswith("noValueColor gray"), (
            f"with no value_range given, the color command should end right after "
            f"'noValueColor gray' with no trailing range clause (ChimeraX auto-scales), "
            f"got {color_line!r}"
        )


class TestRunExport:
    def test_raises_when_neither_uniprot_nor_gene_given(self):
        with pytest.raises(ValueError):
            mod.run_export(uniprot=None, gene=None, log_cb=lambda *_: None)

    def test_raises_filenotfounderror_when_cosmic_file_missing(self, tmp_path, monkeypatch):
        monkeypatch.setattr(mod, "MODELS_ROOT", tmp_path / "cif_models")
        uid_dir = tmp_path / "cif_models" / "P04637"
        uid_dir.mkdir(parents=True)
        _write_synthetic_cif(
            uid_dir / "AF-P04637-F1-model_v4.cif",
            res_ids=[1], res_names=["ALA"], atom_names=["CA"], coords=[[0.0, 0.0, 0.0]],
        )

        with pytest.raises(FileNotFoundError):
            mod.run_export(
                uniprot="P04637", gene="TP53",
                cosmic_file=tmp_path / "does_not_exist_cosmic.tsv",
                output_dir=tmp_path / "out", log_cb=lambda *_: None,
            )

    def test_full_run_writes_both_output_files(self, tmp_path, monkeypatch):
        monkeypatch.setattr(mod, "MODELS_ROOT", tmp_path / "cif_models")
        uid_dir = tmp_path / "cif_models" / "P04637"
        uid_dir.mkdir(parents=True)
        _write_synthetic_cif(
            uid_dir / "AF-P04637-F1-model_v4.cif",
            res_ids=[1, 2, 3], res_names=["ALA", "SER", "GLY"], atom_names=["CA", "CA", "CA"],
            coords=[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]],
        )
        cosmic = tmp_path / "cosmic.tsv"
        pd.DataFrame([
            ("TP53", "p.S2A", "S1", "Confirmed somatic variant"),
            ("TP53", "p.S2A", "S2", "Confirmed somatic variant"),
        ], columns=["GENE_SYMBOL", "MUTATION_AA", "COSMIC_SAMPLE_ID", "MUTATION_SOMATIC_STATUS"]).to_csv(
            cosmic, sep="\t", index=False,
        )

        result = mod.run_export(
            uniprot="P04637", gene="TP53", cosmic_file=cosmic,
            output_dir=tmp_path / "out", log_cb=lambda *_: None,
        )

        assert result.all_out.exists(), "all_ca.tsv should be written to disk"
        assert result.mut_out.exists(), "mutation_ca.tsv should be written to disk"
        assert len(result.all_ca_df) == 3, f"all_ca_df should have one row per residue (3), got {len(result.all_ca_df)}"
        assert len(result.mut_ca_df) == 1, (
            f"mut_ca_df should only include the one residue (position 2) with a COSMIC "
            f"mutation, got {len(result.mut_ca_df)}"
        )
        assert result.mut_ca_df.iloc[0]["total_patients"] == 2, (
            f"the mutation row's total_patients should reflect the 2 distinct COSMIC "
            f"samples, got {result.mut_ca_df.iloc[0]['total_patients']}"
        )
        assert result.defattr_out is not None and result.defattr_out.exists(), (
            "a single-fragment protein should also produce ChimeraX files"
        )

    def test_multi_fragment_protein_skips_chimerax_files(self, tmp_path, monkeypatch):
        monkeypatch.setattr(mod, "MODELS_ROOT", tmp_path / "cif_models")
        uid_dir = tmp_path / "cif_models" / "P04637"
        uid_dir.mkdir(parents=True)
        _write_synthetic_cif(
            uid_dir / "AF-P04637-F1-model_v4.cif",
            res_ids=[1, 2], res_names=["ALA", "SER"], atom_names=["CA", "CA"],
            coords=[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]],
        )
        _write_synthetic_cif(
            uid_dir / "AF-P04637-F2-model_v4.cif",
            res_ids=[1401, 1402], res_names=["ALA", "SER"], atom_names=["CA", "CA"],
            coords=[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]],
        )
        cosmic = tmp_path / "cosmic.tsv"
        pd.DataFrame(columns=["GENE_SYMBOL", "MUTATION_AA", "COSMIC_SAMPLE_ID", "MUTATION_SOMATIC_STATUS"]).to_csv(
            cosmic, sep="\t", index=False,
        )

        result = mod.run_export(
            uniprot="P04637", gene="TP53", cosmic_file=cosmic,
            output_dir=tmp_path / "out", log_cb=lambda *_: None,
        )
        assert result.defattr_out is None and result.chimerax_script_out is None, (
            "a multi-fragment protein (2 CIF files) should skip ChimeraX file "
            "generation entirely, since only fragment 1 was exported"
        )
