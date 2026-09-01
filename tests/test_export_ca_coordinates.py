"""Unit tests for scripts/export_ca_coordinates.py."""
import numpy as np
import pandas as pd
import pytest
from conftest import FakeResponse, import_script

mod = import_script("export_ca_coordinates.py")


def _write_synthetic_cif(path, res_ids, res_names, atom_names, coords,
                          b_factors=None, chain_ids=None):
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
    arr.set_annotation("b_factor", np.asarray(
        b_factors if b_factors is not None else [80.0] * n, dtype=float,
    ))

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


class TestLoadCosmicDataframe:
    def test_caches_by_resolved_path(self, tmp_path, monkeypatch):
        monkeypatch.setattr(mod, "_cosmic_df_cache", {})
        cosmic = tmp_path / "cosmic.tsv"
        pd.DataFrame(
            [("TP53", "p.R175H", "S1", "Confirmed somatic variant")],
            columns=["GENE_SYMBOL", "MUTATION_AA", "COSMIC_SAMPLE_ID", "MUTATION_SOMATIC_STATUS"],
        ).to_csv(cosmic, sep="\t", index=False)

        df1 = mod._load_cosmic_dataframe(cosmic)
        # Rewrite the same path with different content -- a second call must
        # still return the first read (proving it came from cache, not disk).
        pd.DataFrame(
            [("EGFR", "p.L858R", "S2", "Confirmed somatic variant")],
            columns=["GENE_SYMBOL", "MUTATION_AA", "COSMIC_SAMPLE_ID", "MUTATION_SOMATIC_STATUS"],
        ).to_csv(cosmic, sep="\t", index=False)
        df2 = mod._load_cosmic_dataframe(cosmic)

        assert df1["GENE_SYMBOL"].tolist() == df2["GENE_SYMBOL"].tolist() == ["TP53"], (
            "a batch scanning many genes must only read this (potentially huge) file from "
            "disk once -- the second call should return the cached frame, not the rewritten "
            "file's new content"
        )

    def test_different_paths_are_not_conflated(self, tmp_path, monkeypatch):
        monkeypatch.setattr(mod, "_cosmic_df_cache", {})
        cosmic_a = tmp_path / "a.tsv"
        cosmic_b = tmp_path / "b.tsv"
        cols = ["GENE_SYMBOL", "MUTATION_AA", "COSMIC_SAMPLE_ID", "MUTATION_SOMATIC_STATUS"]
        pd.DataFrame([("TP53", "p.R175H", "S1", "Confirmed somatic variant")], columns=cols).to_csv(
            cosmic_a, sep="\t", index=False,
        )
        pd.DataFrame([("EGFR", "p.L858R", "S2", "Confirmed somatic variant")], columns=cols).to_csv(
            cosmic_b, sep="\t", index=False,
        )

        assert mod._load_cosmic_dataframe(cosmic_a)["GENE_SYMBOL"].tolist() == ["TP53"]
        assert mod._load_cosmic_dataframe(cosmic_b)["GENE_SYMBOL"].tolist() == ["EGFR"], (
            "different COSMIC file paths must be cached independently, not share one entry"
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


class TestLoadPtmPositions:
    def test_parses_positions_and_keeps_raw_token(self, tmp_path, monkeypatch):
        tsv = tmp_path / "hotspots.tsv"
        pd.DataFrame([{
            "uniprot_id": "P04637",
            "ptms_on_protein": "S15:Phosphorylation; T18:Phosphorylation",
        }]).to_csv(tsv, sep="\t", index=False)
        monkeypatch.setattr(mod, "PTM_TSV", tsv)

        result = mod.load_ptm_positions("P04637")
        assert result == [(15, "S15:Phosphorylation"), (18, "T18:Phosphorylation")], (
            f"should parse each ';'-separated token into (position, raw token), got {result}"
        )

    def test_missing_file_returns_empty_list(self, tmp_path, monkeypatch):
        monkeypatch.setattr(mod, "PTM_TSV", tmp_path / "does_not_exist.tsv")
        assert mod.load_ptm_positions("P04637") == [], (
            "with no intermediate TSV available yet, this should return an empty list, "
            "not raise FileNotFoundError"
        )

    def test_unknown_protein_returns_empty_list(self, tmp_path, monkeypatch):
        tsv = tmp_path / "hotspots.tsv"
        pd.DataFrame([{"uniprot_id": "P04637", "ptms_on_protein": "S15:Phosphorylation"}]).to_csv(
            tsv, sep="\t", index=False,
        )
        monkeypatch.setattr(mod, "PTM_TSV", tsv)
        assert mod.load_ptm_positions("Q99999") == [], (
            "a UniProt ID not present in the TSV has no PTM row to read -- expected an "
            "empty list"
        )


class TestBuildPtmMarkerLines:
    def test_builds_one_sphere_command_per_position(self):
        ca_df = pd.DataFrame([
            {"position": 15, "x": 1.0, "y": 2.0, "z": 3.0},
            {"position": 18, "x": 4.0, "y": 5.0, "z": 6.0},
        ])
        lines = mod.build_ptm_marker_lines(
            ca_df, [(15, "S15:Phosphorylation"), (18, "T18:Phosphorylation")],
        )
        assert lines == [
            "shape sphere radius 1.2 center 1,2,3 color green name S15_Phosphorylation",
            "shape sphere radius 1.2 center 4,5,6 color green name T18_Phosphorylation",
        ], f"unexpected marker command(s): {lines}"

    def test_skips_positions_missing_from_ca_df(self):
        ca_df = pd.DataFrame([{"position": 15, "x": 1.0, "y": 2.0, "z": 3.0}])
        lines = mod.build_ptm_marker_lines(ca_df, [(15, "S15:Phosphorylation"), (999, "X999:Unknown")])
        assert len(lines) == 1, (
            f"a PTM position with no matching CA coordinate (e.g. outside the exported "
            f"fragment) should be silently skipped, not raise or produce a bad command, "
            f"got {lines}"
        )

    def test_respects_custom_radius_and_color(self):
        ca_df = pd.DataFrame([{"position": 15, "x": 0.0, "y": 0.0, "z": 0.0}])
        lines = mod.build_ptm_marker_lines(ca_df, [(15, "S15:Phosphorylation")], radius=2.5, color="magenta")
        assert lines == ["shape sphere radius 2.5 center 0,0,0 color magenta name S15_Phosphorylation"], (
            f"custom radius/color should be reflected in the command, got {lines}"
        )


class TestBuildMutationMarkerLines:
    def test_builds_show_style_color_triplet_per_position(self):
        lines = mod.build_mutation_marker_lines([592])
        assert lines == [
            "show /A:592 & sidechain atoms",
            "style /A:592 & sidechain stick",
            "color /A:592 & sidechain orange target ab",
        ], f"unexpected marker command(s): {lines}"

    def test_color_targets_atoms_bonds_only_not_cartoon(self):
        # "target ab" (atoms/bonds) is critical: "c"/"r" would be cartoons,
        # which would silently overwrite a heatmap's own color at that residue.
        lines = mod.build_mutation_marker_lines([592])
        color_line = next(l for l in lines if l.startswith("color"))
        assert "target ab" in color_line, (
            f"the color command must restrict its target to atoms/bonds only, got {color_line!r}"
        )

    def test_sorts_and_deduplicates_positions(self):
        lines = mod.build_mutation_marker_lines([18, 5, 18, 5])
        specs = [l.split()[1] for l in lines if l.startswith("show")]
        assert specs == ["/A:5", "/A:18"], (
            f"duplicate positions should be collapsed and positions sorted for stable, "
            f"predictable script output, got {specs}"
        )

    def test_respects_custom_chain_and_color(self):
        lines = mod.build_mutation_marker_lines([10], chain_id="B", color="yellow")
        assert lines == [
            "show /B:10 & sidechain atoms",
            "style /B:10 & sidechain stick",
            "color /B:10 & sidechain yellow target ab",
        ], f"custom chain_id/color should be reflected in the commands, got {lines}"

    def test_empty_positions_returns_empty_list(self):
        assert mod.build_mutation_marker_lines([]) == []


class TestBuildConfidenceDimLines:
    def test_transparency_is_100_minus_plddt(self):
        lines = mod.build_confidence_dim_lines({15: 90.0, 18: 40.0})
        assert lines == [
            "transparency /A:15 10 target c",
            "transparency /A:18 60 target c",
        ], f"unexpected dim command(s): {lines}"

    def test_target_c_restricts_to_cartoon_only(self):
        # "target c" is critical: an unscoped transparency command would also
        # fade any atom-level markers (PTM spheres, mutation sticks) layered
        # into the same script, not just the cartoon this feature targets.
        lines = mod.build_confidence_dim_lines({15: 50.0})
        assert lines == ["transparency /A:15 50 target c"]

    def test_full_confidence_is_fully_opaque(self):
        lines = mod.build_confidence_dim_lines({15: 100.0})
        assert lines == ["transparency /A:15 0 target c"], (
            "pLDDT 100 (maximum confidence) should map to 0% transparent -- fully opaque"
        )

    def test_values_are_clamped_to_0_100(self):
        # pLDDT is always 0-100 in practice, but clamp defensively so a
        # slightly-out-of-range value can never produce an invalid ChimeraX
        # transparency percentage.
        lines = mod.build_confidence_dim_lines({1: -5.0, 2: 150.0})
        assert lines == [
            "transparency /A:1 100 target c",
            "transparency /A:2 0 target c",
        ], f"unexpected dim command(s): {lines}"

    def test_positions_are_sorted(self):
        lines = mod.build_confidence_dim_lines({30: 80.0, 5: 80.0, 15: 80.0})
        specs = [l.split()[1] for l in lines]
        assert specs == ["/A:5", "/A:15", "/A:30"]

    def test_respects_custom_chain_id(self):
        lines = mod.build_confidence_dim_lines({10: 70.0}, chain_id="B")
        assert lines == ["transparency /B:10 30 target c"]

    def test_empty_map_returns_empty_list(self):
        assert mod.build_confidence_dim_lines({}) == []


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

    def test_extra_lines_inserted_before_lighting(self, tmp_path):
        out_path = mod.write_chimerax_script(
            tmp_path / "model.cif", tmp_path / "mutations.defattr", tmp_path / "view.cxc",
            extra_lines=["shape sphere radius 1.2 center 1,2,3 color purple name ptm1"],
        )
        lines = out_path.read_text(encoding="utf-8").splitlines()
        assert lines.index("shape sphere radius 1.2 center 1,2,3 color purple name ptm1") == lines.index("lighting soft") - 1, (
            f"extra_lines (e.g. PTM marker spheres) should be layered in right after the "
            f"coloring command and before the final lighting command, got {lines}"
        )

    def test_lighting_defaults_to_soft(self, tmp_path):
        out_path = mod.write_chimerax_script(
            tmp_path / "model.cif", tmp_path / "mutations.defattr", tmp_path / "view.cxc",
        )
        assert "lighting soft" in out_path.read_text(encoding="utf-8").splitlines()

    def test_lighting_can_be_overridden(self, tmp_path):
        out_path = mod.write_chimerax_script(
            tmp_path / "model.cif", tmp_path / "mutations.defattr", tmp_path / "view.cxc",
            lighting="simple",
        )
        lines = out_path.read_text(encoding="utf-8").splitlines()
        assert "lighting simple" in lines and "lighting soft" not in lines, (
            f"an explicit lighting override should replace 'soft', not add to it, got {lines}"
        )


class TestBuildMutationKeyLines:
    def test_produces_a_two_stop_zero_to_max_key(self):
        lines = mod.build_mutation_key_lines(40.0, log_scale=False, low_color="#111111", high_color="#eeeeee")
        key_line = lines[0]
        assert key_line == "key #111111:0 #eeeeee:40", (
            f"the key should always be built from an explicit 2-color list, labeled 0 and "
            f"the true max, got {key_line!r}"
        )

    def test_never_names_a_palette(self):
        # A named-palette key with blank middle labels (e.g. "key Reds :0 : : : :40")
        # was observed to render garbled overlapping digits on at least one real
        # protein -- so the key must never reference a palette by name, even when
        # the heatmap's real 3D coloring is using the default "Reds" unmodified.
        lines = mod.build_mutation_key_lines(
            40.0, log_scale=False,
            low_color=mod.MUTATION_DEFAULT_LOW_COLOR, high_color=mod.MUTATION_DEFAULT_HIGH_COLOR,
        )
        assert "Reds" not in lines[0], f"the key must not name the 'Reds' palette, got {lines[0]!r}"

    def test_includes_a_title_label(self):
        lines = mod.build_mutation_key_lines(10.0, log_scale=False, low_color="white", high_color="red")
        assert any(l.startswith("2dlabels text") and "Patients within 10" in l for l in lines), (
            f"a title label should accompany the key so a screenshot is self-explanatory, got {lines}"
        )

    def test_log_scale_high_label_shows_real_patient_count_not_log_value(self):
        # log1p(9) == log(10) ~= 2.302585 -- if the key were built from the raw
        # (pre-transform) value, the high label would be 2.302585 instead of 9.
        lines = mod.build_mutation_key_lines(np.log1p(9), log_scale=True, low_color="white", high_color="red")
        key_line = lines[0]
        assert key_line == "key white:0 red:9", (
            f"the log-space max should be converted back to a real patient count via expm1 "
            f"for display, not left as a log1p value, got {key_line!r}"
        )

    def test_log_scale_title_notes_the_log_scale(self):
        lines = mod.build_mutation_key_lines(1.0, log_scale=True, low_color="white", high_color="red")
        title_line = next(l for l in lines if l.startswith("2dlabels"))
        assert "log scale" in title_line.lower(), (
            f"the title should flag that the underlying scale is logarithmic, got {title_line!r}"
        )

    def test_zero_max_produces_a_valid_key(self):
        # e.g. a protein with zero nearby mutations everywhere.
        lines = mod.build_mutation_key_lines(0.0, log_scale=False, low_color="white", high_color="red")
        assert lines[0] == "key white:0 red:0", f"got {lines[0]!r}"


class TestBuildPlddtKeyLines:
    def test_produces_a_two_stop_0_to_100_key(self):
        lines = mod.build_plddt_key_lines(low_color="orange", high_color="blue")
        assert lines[0] == "key orange:0 blue:100", (
            f"the key should always be built from an explicit 2-color list, labeled 0 and "
            f"100 (pLDDT's fixed full range), got {lines[0]!r}"
        )

    def test_never_names_a_palette(self):
        # See TestBuildMutationKeyLines.test_never_names_a_palette -- same failure
        # mode applies to a named "alphafold" palette key.
        lines = mod.build_plddt_key_lines(
            low_color=mod.PLDDT_DEFAULT_LOW_COLOR, high_color=mod.PLDDT_DEFAULT_HIGH_COLOR,
        )
        assert "alphafold" not in lines[0], f"the key must not name the 'alphafold' palette, got {lines[0]!r}"

    def test_includes_a_title_label(self):
        lines = mod.build_plddt_key_lines(low_color="orange", high_color="blue")
        assert any("pLDDT" in l for l in lines if l.startswith("2dlabels")), (
            f"a title label should accompany the key so a screenshot is self-explanatory, got {lines}"
        )


class TestWritePlddtChimeraxScript:
    def test_colors_by_bfactor_with_alphafold_palette(self, tmp_path):
        out_path = mod.write_plddt_chimerax_script(tmp_path / "model.cif", tmp_path / "plddt_view.cxc")
        lines = out_path.read_text(encoding="utf-8").splitlines()
        assert 'open "' in lines[0] and "model.cif" in lines[0], (
            f"the script should open the given CIF path, got {lines[0]!r}"
        )
        assert "color bfactor #1 palette alphafold" in lines, (
            f"the pLDDT heatmap should use ChimeraX's own bfactor coloring with its "
            f"built-in alphafold palette -- no defattr file needed since AlphaFold CIFs "
            f"already carry pLDDT in the B-factor field. Got lines: {lines}"
        )

    def test_no_defattr_open_command(self, tmp_path):
        # Unlike write_chimerax_script (mutation heatmap), this script must not
        # reference any .defattr file -- pLDDT comes straight from the CIF's own
        # B-factor column via ChimeraX's `color bfactor`. Checking for the ".defattr"
        # extension specifically (not just "defattr") avoids a false positive from
        # pytest's own tmp_path directory name possibly containing that substring.
        out_path = mod.write_plddt_chimerax_script(tmp_path / "model.cif", tmp_path / "plddt_view.cxc")
        text = out_path.read_text(encoding="utf-8")
        assert ".defattr" not in text, f"the pLDDT script should not reference a defattr file, got:\n{text}"

    def test_extra_lines_inserted_before_lighting(self, tmp_path):
        out_path = mod.write_plddt_chimerax_script(
            tmp_path / "model.cif", tmp_path / "plddt_view.cxc",
            extra_lines=["shape sphere radius 1.2 center 1,2,3 color purple name ptm1"],
        )
        lines = out_path.read_text(encoding="utf-8").splitlines()
        assert lines.index("shape sphere radius 1.2 center 1,2,3 color purple name ptm1") == lines.index("lighting soft") - 1

    def test_palette_can_be_overridden(self, tmp_path):
        out_path = mod.write_plddt_chimerax_script(
            tmp_path / "model.cif", tmp_path / "plddt_view.cxc", palette="orange:blue",
        )
        lines = out_path.read_text(encoding="utf-8").splitlines()
        assert "color bfactor #1 palette orange:blue" in lines, (
            f"a custom palette override should replace 'alphafold' in the color command, got {lines}"
        )


class TestWritePlainChimeraxScript:
    def test_opens_cif_and_shows_plain_cartoon(self, tmp_path):
        out_path = mod.write_plain_chimerax_script(tmp_path / "model.cif", tmp_path / "ptm_markers_view.cxc")
        lines = out_path.read_text(encoding="utf-8").splitlines()
        assert 'open "' in lines[0] and "model.cif" in lines[0]
        assert "cartoon" in lines
        assert not any(l.startswith("color") for l in lines), (
            f"a plain script (no heatmap) must not include any color command, got {lines}"
        )

    def test_extra_lines_inserted_before_lighting(self, tmp_path):
        out_path = mod.write_plain_chimerax_script(
            tmp_path / "model.cif", tmp_path / "ptm_markers_view.cxc",
            extra_lines=["shape sphere radius 1.2 center 1,2,3 color purple name ptm1"],
        )
        lines = out_path.read_text(encoding="utf-8").splitlines()
        assert lines.index("shape sphere radius 1.2 center 1,2,3 color purple name ptm1") == lines.index("lighting soft") - 1


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
        assert result.all_out.parent.name == "TP53_P04637", (
            f"the per-protein output folder should be named {{gene}}_{{uid}} for "
            f"readability, not just the bare UniProt accession, got {result.all_out.parent.name!r}"
        )
        assert len(result.all_ca_df) == 3, f"all_ca_df should have one row per residue (3), got {len(result.all_ca_df)}"
        assert len(result.mut_ca_df) == 1, (
            f"mut_ca_df should only include the one residue (position 2) with a COSMIC "
            f"mutation, got {len(result.mut_ca_df)}"
        )
        assert result.mut_ca_df.iloc[0]["total_patients"] == 2, (
            f"the mutation row's total_patients should reflect the 2 distinct COSMIC "
            f"samples, got {result.mut_ca_df.iloc[0]['total_patients']}"
        )
        assert result.mutation_defattr_out is not None and result.mutation_defattr_out.exists(), (
            "a single-fragment protein should also produce ChimeraX files (mutation heatmap "
            "is on by default)"
        )
        assert result.plddt_chimerax_script_out is None, (
            "pLDDT heatmap defaults to off, so no pLDDT script should be produced unless "
            "explicitly requested"
        )
        mutation_script_text = result.mutation_chimerax_script_out.read_text()
        assert "\nkey " in mutation_script_text, (
            f"the generated mutation heatmap script should include a ChimeraX color key, "
            f"got:\n{mutation_script_text}"
        )

    def test_output_folder_name_sanitizes_unsafe_gene_characters(self, tmp_path, monkeypatch):
        # Real HGNC gene symbols never contain filesystem-unsafe characters, but
        # sanitize defensively anyway rather than trust an upstream API response.
        monkeypatch.setattr(mod, "MODELS_ROOT", tmp_path / "cif_models")
        uid_dir = tmp_path / "cif_models" / "P04637"
        uid_dir.mkdir(parents=True)
        _write_synthetic_cif(
            uid_dir / "AF-P04637-F1-model_v4.cif",
            res_ids=[1], res_names=["ALA"], atom_names=["CA"], coords=[[0.0, 0.0, 0.0]],
        )
        cosmic = tmp_path / "cosmic.tsv"
        pd.DataFrame(columns=["GENE_SYMBOL", "MUTATION_AA", "COSMIC_SAMPLE_ID", "MUTATION_SOMATIC_STATUS"]).to_csv(
            cosmic, sep="\t", index=False,
        )

        result = mod.run_export(
            uniprot="P04637", gene="MY/GENE:1", cosmic_file=cosmic,
            output_dir=tmp_path / "out", log_cb=lambda *_: None,
        )
        assert result.all_out.parent.name == "MY_GENE_1_P04637", (
            f"unsafe characters in the gene name should be replaced, not break folder "
            f"creation, got {result.all_out.parent.name!r}"
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
            output_dir=tmp_path / "out", plddt_heatmap=True, log_cb=lambda *_: None,
        )
        assert (
            result.mutation_defattr_out is None
            and result.mutation_chimerax_script_out is None
            and result.plddt_chimerax_script_out is None
        ), (
            "a multi-fragment protein (2 CIF files) should skip ALL ChimeraX heatmap "
            "generation entirely, since only fragment 1 was exported"
        )

    def _defattr_values(self, path) -> dict[int, float]:
        lines = path.read_text().splitlines()
        values = {}
        for line in lines:
            if not line.startswith("\t/"):
                continue
            _, spec, value = line.split("\t")
            values[int(spec.split(":")[1])] = float(value)
        return values

    def test_log_scale_writes_log1p_values_under_separate_attribute(self, tmp_path, monkeypatch):
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
            ("TP53", "p.S2A", "S3", "Confirmed somatic variant"),
        ], columns=["GENE_SYMBOL", "MUTATION_AA", "COSMIC_SAMPLE_ID", "MUTATION_SOMATIC_STATUS"]).to_csv(
            cosmic, sep="\t", index=False,
        )

        result = mod.run_export(
            uniprot="P04637", gene="TP53", cosmic_file=cosmic,
            output_dir=tmp_path / "out", log_scale=True, log_cb=lambda *_: None,
        )

        lines = result.mutation_defattr_out.read_text().splitlines()
        assert lines[0] == "attribute: patients_within_10A_log", (
            f"log-scaled export should write under a distinct attribute name, so it isn't "
            f"confused with the raw count, got {lines[0]!r}"
        )
        values = self._defattr_values(result.mutation_defattr_out)
        # All 3 residues are within 10A of position 2 (the only mutation, 3 patients),
        # so every residue's raw patients_within_10A is 3.
        assert values[2] == pytest.approx(np.log1p(3)), (
            f"the defattr value should be log1p of the raw patient count (log1p(3) "
            f"~= {np.log1p(3):.6f}), got {values[2]}"
        )
        assert "patients_within_10A_log" not in result.all_ca_df.columns, (
            "the log-scaled column exists only for the ChimeraX heatmap and must not leak "
            "into the returned all_ca_df (or, by extension, the already-written all_ca.tsv)"
        )

    def test_log_scale_false_uses_raw_attribute_and_values(self, tmp_path, monkeypatch):
        monkeypatch.setattr(mod, "MODELS_ROOT", tmp_path / "cif_models")
        uid_dir = tmp_path / "cif_models" / "P04637"
        uid_dir.mkdir(parents=True)
        _write_synthetic_cif(
            uid_dir / "AF-P04637-F1-model_v4.cif",
            res_ids=[1, 2], res_names=["ALA", "SER"], atom_names=["CA", "CA"],
            coords=[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]],
        )
        cosmic = tmp_path / "cosmic.tsv"
        pd.DataFrame([
            ("TP53", "p.S2A", "S1", "Confirmed somatic variant"),
        ], columns=["GENE_SYMBOL", "MUTATION_AA", "COSMIC_SAMPLE_ID", "MUTATION_SOMATIC_STATUS"]).to_csv(
            cosmic, sep="\t", index=False,
        )

        result = mod.run_export(
            uniprot="P04637", gene="TP53", cosmic_file=cosmic,
            output_dir=tmp_path / "out", log_scale=False, log_cb=lambda *_: None,
        )
        lines = result.mutation_defattr_out.read_text().splitlines()
        assert lines[0] == "attribute: patients_within_10A", (
            f"with log_scale=False (the default), the attribute name must stay the plain "
            f"'patients_within_10A', got {lines[0]!r}"
        )
        values = self._defattr_values(result.mutation_defattr_out)
        assert values[2] == 1.0, (
            f"with log_scale=False the defattr value should be the raw (unscaled) patient "
            f"count, got {values[2]}"
        )

    def _single_fragment_export_kwargs(self, tmp_path, monkeypatch):
        """Shared single-fragment CIF + empty COSMIC setup for the heatmap-toggle tests
        below, which only care about which files get produced, not their content.
        """
        monkeypatch.setattr(mod, "MODELS_ROOT", tmp_path / "cif_models")
        uid_dir = tmp_path / "cif_models" / "P04637"
        uid_dir.mkdir(parents=True)
        _write_synthetic_cif(
            uid_dir / "AF-P04637-F1-model_v4.cif",
            res_ids=[1], res_names=["ALA"], atom_names=["CA"], coords=[[0.0, 0.0, 0.0]],
        )
        cosmic = tmp_path / "cosmic.tsv"
        pd.DataFrame(columns=["GENE_SYMBOL", "MUTATION_AA", "COSMIC_SAMPLE_ID", "MUTATION_SOMATIC_STATUS"]).to_csv(
            cosmic, sep="\t", index=False,
        )
        return dict(uniprot="P04637", gene="TP53", cosmic_file=cosmic,
                    output_dir=tmp_path / "out", log_cb=lambda *_: None)

    def test_mutation_heatmap_false_skips_mutation_files_only(self, tmp_path, monkeypatch):
        kwargs = self._single_fragment_export_kwargs(tmp_path, monkeypatch)
        result = mod.run_export(mutation_heatmap=False, plddt_heatmap=True, **kwargs)
        assert result.mutation_defattr_out is None and result.mutation_chimerax_script_out is None, (
            "with mutation_heatmap=False, no mutation-heatmap files should be produced"
        )
        assert result.plddt_chimerax_script_out is not None and result.plddt_chimerax_script_out.exists(), (
            "plddt_heatmap=True should still produce the pLDDT script even with the "
            "mutation heatmap turned off -- the two are independent"
        )

    def test_plddt_heatmap_true_produces_script_with_no_mutation_files(self, tmp_path, monkeypatch):
        kwargs = self._single_fragment_export_kwargs(tmp_path, monkeypatch)
        result = mod.run_export(mutation_heatmap=False, plddt_heatmap=True, **kwargs)
        text = result.plddt_chimerax_script_out.read_text()
        assert "color bfactor #1 palette alphafold" in text, (
            "the produced pLDDT script should color by bfactor with the alphafold palette"
        )
        assert "\nkey " in text, (
            f"the produced pLDDT script should include a ChimeraX color key, got:\n{text}"
        )

    def test_both_heatmaps_true_produces_all_three_files(self, tmp_path, monkeypatch):
        kwargs = self._single_fragment_export_kwargs(tmp_path, monkeypatch)
        result = mod.run_export(mutation_heatmap=True, plddt_heatmap=True, **kwargs)
        assert result.mutation_defattr_out is not None and result.mutation_defattr_out.exists()
        assert result.mutation_chimerax_script_out is not None and result.mutation_chimerax_script_out.exists()
        assert result.plddt_chimerax_script_out is not None and result.plddt_chimerax_script_out.exists()
        assert result.mutation_chimerax_script_out != result.plddt_chimerax_script_out, (
            "the two heatmaps must be written as separate .cxc scripts, not one combined "
            "script -- ChimeraX's color command replaces the previous coloring rather than "
            "layering, so a single script could only ever show the last color applied"
        )

    def test_both_heatmaps_false_skips_all_chimerax_files(self, tmp_path, monkeypatch):
        kwargs = self._single_fragment_export_kwargs(tmp_path, monkeypatch)
        result = mod.run_export(mutation_heatmap=False, plddt_heatmap=False, **kwargs)
        assert (
            result.mutation_defattr_out is None
            and result.mutation_chimerax_script_out is None
            and result.plddt_chimerax_script_out is None
        ), "with both heatmaps off, no ChimeraX files at all should be produced"

    def test_default_kwargs_enable_only_mutation_heatmap(self, tmp_path, monkeypatch):
        # Backward-compatible default: mutation_heatmap=True, plddt_heatmap=False,
        # matching this tool's original always-on mutation-heatmap behavior.
        kwargs = self._single_fragment_export_kwargs(tmp_path, monkeypatch)
        result = mod.run_export(**kwargs)
        assert result.mutation_defattr_out is not None and result.mutation_defattr_out.exists()
        assert result.plddt_chimerax_script_out is None

    def _write_ptm_tsv(self, tmp_path, monkeypatch):
        ptm_tsv = tmp_path / "hotspots.tsv"
        pd.DataFrame([{"uniprot_id": "P04637", "ptms_on_protein": "A1:Phosphorylation"}]).to_csv(
            ptm_tsv, sep="\t", index=False,
        )
        monkeypatch.setattr(mod, "PTM_TSV", ptm_tsv)

    def test_mark_ptm_sites_adds_sphere_to_mutation_script(self, tmp_path, monkeypatch):
        kwargs = self._single_fragment_export_kwargs(tmp_path, monkeypatch)
        self._write_ptm_tsv(tmp_path, monkeypatch)

        result = mod.run_export(mutation_heatmap=True, plddt_heatmap=False, mark_ptm_sites=True, **kwargs)
        text = result.mutation_chimerax_script_out.read_text()
        assert "shape sphere" in text, (
            f"a PTM marker sphere should be layered into the mutation heatmap script, got:\n{text}"
        )
        assert result.plain_chimerax_script_out is None, (
            "a heatmap was requested, so no separate 'plain' PTM-only script should be produced"
        )

    def test_mark_ptm_sites_with_no_heatmap_produces_plain_script(self, tmp_path, monkeypatch):
        kwargs = self._single_fragment_export_kwargs(tmp_path, monkeypatch)
        self._write_ptm_tsv(tmp_path, monkeypatch)

        result = mod.run_export(mutation_heatmap=False, plddt_heatmap=False, mark_ptm_sites=True, **kwargs)
        assert result.plain_chimerax_script_out is not None and result.plain_chimerax_script_out.exists(), (
            "mark_ptm_sites=True with both heatmaps off should still produce a plain "
            "(uncolored) script carrying the PTM markers"
        )
        text = result.plain_chimerax_script_out.read_text()
        assert "shape sphere" in text
        assert not any(l.startswith("color") for l in text.splitlines()), (
            f"the plain script must have no color command, got:\n{text}"
        )

    def test_mark_ptm_sites_false_produces_no_spheres(self, tmp_path, monkeypatch):
        kwargs = self._single_fragment_export_kwargs(tmp_path, monkeypatch)
        self._write_ptm_tsv(tmp_path, monkeypatch)

        result = mod.run_export(mutation_heatmap=True, mark_ptm_sites=False, **kwargs)
        assert "shape sphere" not in result.mutation_chimerax_script_out.read_text(), (
            "with mark_ptm_sites=False (the default), no spheres should appear even though "
            "matching PTM data exists"
        )

    def test_mark_ptm_sites_with_no_ptm_data_does_not_crash(self, tmp_path, monkeypatch):
        kwargs = self._single_fragment_export_kwargs(tmp_path, monkeypatch)
        monkeypatch.setattr(mod, "PTM_TSV", tmp_path / "does_not_exist.tsv")

        result = mod.run_export(mutation_heatmap=True, mark_ptm_sites=True, **kwargs)
        assert result.mutation_chimerax_script_out is not None and result.mutation_chimerax_script_out.exists(), (
            "missing PTM data should degrade gracefully (heatmap still produced, just no "
            "markers), not raise"
        )
        assert "shape sphere" not in result.mutation_chimerax_script_out.read_text()

    def _export_kwargs_with_mutations(self, tmp_path, monkeypatch):
        """Single-fragment CIF + COSMIC setup with a real mutation position
        (592), for the mark_mutations tests below.
        """
        monkeypatch.setattr(mod, "MODELS_ROOT", tmp_path / "cif_models")
        uid_dir = tmp_path / "cif_models" / "P04637"
        uid_dir.mkdir(parents=True)
        _write_synthetic_cif(
            uid_dir / "AF-P04637-F1-model_v4.cif",
            res_ids=[1, 592], res_names=["ALA", "SER"], atom_names=["CA", "CA"],
            coords=[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]],
        )
        cosmic = tmp_path / "cosmic.tsv"
        pd.DataFrame([
            ("TP53", "p.S592A", "S1", "Confirmed somatic variant"),
        ], columns=["GENE_SYMBOL", "MUTATION_AA", "COSMIC_SAMPLE_ID", "MUTATION_SOMATIC_STATUS"]).to_csv(
            cosmic, sep="\t", index=False,
        )
        return dict(uniprot="P04637", gene="TP53", cosmic_file=cosmic,
                    output_dir=tmp_path / "out", log_cb=lambda *_: None)

    def test_mark_mutations_adds_sticks_to_mutation_script(self, tmp_path, monkeypatch):
        kwargs = self._export_kwargs_with_mutations(tmp_path, monkeypatch)
        result = mod.run_export(mutation_heatmap=True, mark_mutations=True, **kwargs)
        text = result.mutation_chimerax_script_out.read_text()
        assert "/A:592 & sidechain" in text and "orange" in text, (
            f"a mutation marker stick should be layered into the mutation heatmap script, "
            f"got:\n{text}"
        )

    def test_mark_mutations_false_produces_no_sticks(self, tmp_path, monkeypatch):
        kwargs = self._export_kwargs_with_mutations(tmp_path, monkeypatch)
        result = mod.run_export(mutation_heatmap=True, mark_mutations=False, **kwargs)
        assert "sidechain" not in result.mutation_chimerax_script_out.read_text(), (
            "with mark_mutations=False (the default), no sticks should appear even though "
            "COSMIC mutation data exists"
        )

    def test_mark_mutations_with_no_heatmap_produces_plain_script(self, tmp_path, monkeypatch):
        kwargs = self._export_kwargs_with_mutations(tmp_path, monkeypatch)
        result = mod.run_export(
            mutation_heatmap=False, plddt_heatmap=False, mark_mutations=True, **kwargs,
        )
        assert result.plain_chimerax_script_out is not None and result.plain_chimerax_script_out.exists(), (
            "mark_mutations=True with both heatmaps off should still produce a plain "
            "(uncolored) script carrying the mutation markers"
        )
        text = result.plain_chimerax_script_out.read_text()
        assert "/A:592 & sidechain" in text
        assert "color byattribute" not in text and "color bfactor" not in text, (
            f"the plain script must have no heatmap color command, got:\n{text}"
        )

    def test_mark_mutations_with_no_cosmic_data_does_not_crash(self, tmp_path, monkeypatch):
        kwargs = self._single_fragment_export_kwargs(tmp_path, monkeypatch)  # empty COSMIC
        result = mod.run_export(mutation_heatmap=True, mark_mutations=True, **kwargs)
        assert result.mutation_chimerax_script_out is not None and result.mutation_chimerax_script_out.exists(), (
            "no COSMIC mutations for this protein should degrade gracefully (heatmap still "
            "produced, just no sticks), not raise"
        )
        assert "sidechain" not in result.mutation_chimerax_script_out.read_text()

    def test_mark_ptm_sites_and_mark_mutations_coexist(self, tmp_path, monkeypatch):
        kwargs = self._export_kwargs_with_mutations(tmp_path, monkeypatch)
        self._write_ptm_tsv(tmp_path, monkeypatch)  # PTM at position 1, present in this CIF

        result = mod.run_export(
            mutation_heatmap=True, mark_ptm_sites=True, mark_mutations=True, **kwargs,
        )
        text = result.mutation_chimerax_script_out.read_text()
        assert "shape sphere" in text and "sidechain" in text, (
            f"both marker types should be able to coexist in the same script, got:\n{text}"
        )

    def _export_kwargs_with_varied_confidence(self, tmp_path, monkeypatch):
        """Single-fragment CIF with two residues at different pLDDT (b_factor)
        values, for the dim_low_confidence tests below.
        """
        monkeypatch.setattr(mod, "MODELS_ROOT", tmp_path / "cif_models")
        uid_dir = tmp_path / "cif_models" / "P04637"
        uid_dir.mkdir(parents=True)
        _write_synthetic_cif(
            uid_dir / "AF-P04637-F1-model_v4.cif",
            res_ids=[1, 2], res_names=["ALA", "SER"], atom_names=["CA", "CA"],
            coords=[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]],
            b_factors=[95.0, 30.0],
        )
        cosmic = tmp_path / "cosmic.tsv"
        pd.DataFrame(columns=["GENE_SYMBOL", "MUTATION_AA", "COSMIC_SAMPLE_ID", "MUTATION_SOMATIC_STATUS"]).to_csv(
            cosmic, sep="\t", index=False,
        )
        return dict(uniprot="P04637", gene="TP53", cosmic_file=cosmic,
                    output_dir=tmp_path / "out", log_cb=lambda *_: None)

    def test_dim_low_confidence_adds_transparency_lines(self, tmp_path, monkeypatch):
        kwargs = self._export_kwargs_with_varied_confidence(tmp_path, monkeypatch)
        result = mod.run_export(mutation_heatmap=True, dim_low_confidence=True, **kwargs)
        text = result.mutation_chimerax_script_out.read_text()
        assert "transparency /A:1 5 target c" in text, (
            f"residue 1 (pLDDT 95) should be dimmed 5% (100-95), got:\n{text}"
        )
        assert "transparency /A:2 70 target c" in text, (
            f"residue 2 (pLDDT 30) should be dimmed 70% (100-30), got:\n{text}"
        )

    def test_dim_low_confidence_false_produces_no_transparency_lines(self, tmp_path, monkeypatch):
        kwargs = self._export_kwargs_with_varied_confidence(tmp_path, monkeypatch)
        result = mod.run_export(mutation_heatmap=True, dim_low_confidence=False, **kwargs)
        assert "transparency" not in result.mutation_chimerax_script_out.read_text(), (
            "with dim_low_confidence=False (the default), no transparency commands should "
            "appear"
        )

    def test_dim_low_confidence_has_no_effect_without_mutation_heatmap(self, tmp_path, monkeypatch):
        kwargs = self._export_kwargs_with_varied_confidence(tmp_path, monkeypatch)
        result = mod.run_export(
            mutation_heatmap=False, plddt_heatmap=True, dim_low_confidence=True, **kwargs,
        )
        assert "transparency" not in result.plddt_chimerax_script_out.read_text(), (
            "dim_low_confidence is scoped to the mutation heatmap only -- it should have no "
            "effect on the pLDDT heatmap"
        )

    def test_dim_low_confidence_coexists_with_markers(self, tmp_path, monkeypatch):
        kwargs = self._export_kwargs_with_varied_confidence(tmp_path, monkeypatch)
        self._write_ptm_tsv(tmp_path, monkeypatch)  # PTM at position 1

        result = mod.run_export(
            mutation_heatmap=True, mark_ptm_sites=True, dim_low_confidence=True, **kwargs,
        )
        text = result.mutation_chimerax_script_out.read_text()
        assert "transparency /A:1 5 target c" in text and "shape sphere" in text, (
            f"dimming and PTM markers should coexist in the same script, got:\n{text}"
        )

    def test_dim_low_confidence_switches_to_simple_lighting(self, tmp_path, monkeypatch):
        # "soft" lighting's depth cues come entirely from ambient shadowing,
        # which ChimeraX doesn't compute correctly once any part of the
        # model is transparent -- opaque residues end up flatly lit at full
        # ambient intensity too, not just the dimmed ones. "simple" uses
        # real directional lights instead, which don't have that failure mode.
        kwargs = self._export_kwargs_with_varied_confidence(tmp_path, monkeypatch)
        result = mod.run_export(mutation_heatmap=True, dim_low_confidence=True, **kwargs)
        lines = result.mutation_chimerax_script_out.read_text().splitlines()
        assert "lighting simple" in lines and "lighting soft" not in lines, (
            f"dim_low_confidence=True must switch away from 'soft' lighting, got {lines}"
        )

    def test_lighting_stays_soft_without_dim_low_confidence(self, tmp_path, monkeypatch):
        kwargs = self._export_kwargs_with_varied_confidence(tmp_path, monkeypatch)
        result = mod.run_export(mutation_heatmap=True, dim_low_confidence=False, **kwargs)
        lines = result.mutation_chimerax_script_out.read_text().splitlines()
        assert "lighting soft" in lines and "lighting simple" not in lines, (
            f"without dimming there's no transparency to conflict with 'soft' lighting, so "
            f"it should stay the default, got {lines}"
        )


class TestRunExportColors:
    def _kwargs(self, tmp_path, monkeypatch):
        """Single-fragment CIF + one COSMIC mutation (position 592), matching
        the setup TestRunExport's marker tests use -- gives both heatmaps and
        both marker types something real to color/mark.
        """
        monkeypatch.setattr(mod, "MODELS_ROOT", tmp_path / "cif_models")
        uid_dir = tmp_path / "cif_models" / "P04637"
        uid_dir.mkdir(parents=True)
        _write_synthetic_cif(
            uid_dir / "AF-P04637-F1-model_v4.cif",
            res_ids=[1, 592], res_names=["ALA", "SER"], atom_names=["CA", "CA"],
            coords=[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]],
        )
        cosmic = tmp_path / "cosmic.tsv"
        pd.DataFrame([
            ("TP53", "p.S592A", "S1", "Confirmed somatic variant"),
        ], columns=["GENE_SYMBOL", "MUTATION_AA", "COSMIC_SAMPLE_ID", "MUTATION_SOMATIC_STATUS"]).to_csv(
            cosmic, sep="\t", index=False,
        )
        ptm_tsv = tmp_path / "hotspots.tsv"
        pd.DataFrame([{"uniprot_id": "P04637", "ptms_on_protein": "A1:Phosphorylation"}]).to_csv(
            ptm_tsv, sep="\t", index=False,
        )
        monkeypatch.setattr(mod, "PTM_TSV", ptm_tsv)
        return dict(uniprot="P04637", gene="TP53", cosmic_file=cosmic,
                    output_dir=tmp_path / "out", log_cb=lambda *_: None)

    def test_default_colors_still_use_named_palettes_for_the_3d_coloring(self, tmp_path, monkeypatch):
        kwargs = self._kwargs(tmp_path, monkeypatch)
        result = mod.run_export(mutation_heatmap=True, plddt_heatmap=True, **kwargs)
        mut_text = result.mutation_chimerax_script_out.read_text()
        plddt_text = result.plddt_chimerax_script_out.read_text()
        assert "palette Reds" in mut_text, (
            f"leaving colors untouched should still emit the real 'Reds' palette for the "
            f"3D coloring, unchanged from before color customization existed, got:\n{mut_text}"
        )
        assert "palette alphafold" in plddt_text, (
            f"leaving colors untouched should still emit the real 'alphafold' palette for "
            f"the 3D coloring, got:\n{plddt_text}"
        )

    def test_default_colors_key_never_names_a_palette(self, tmp_path, monkeypatch):
        # The key is always an explicit 2-color list (see build_mutation_key_lines's
        # docstring) even when the 3D coloring above uses the real named palette --
        # a named-palette key was observed to render garbled overlapping numbers.
        kwargs = self._kwargs(tmp_path, monkeypatch)
        result = mod.run_export(mutation_heatmap=True, plddt_heatmap=True, **kwargs)
        mut_text = result.mutation_chimerax_script_out.read_text()
        plddt_text = result.plddt_chimerax_script_out.read_text()
        assert f"key {mod.MUTATION_DEFAULT_LOW_COLOR}:0 {mod.MUTATION_DEFAULT_HIGH_COLOR}:" in mut_text, (
            f"got:\n{mut_text}"
        )
        assert f"key {mod.PLDDT_DEFAULT_LOW_COLOR}:0 {mod.PLDDT_DEFAULT_HIGH_COLOR}:100" in plddt_text, (
            f"got:\n{plddt_text}"
        )

    def test_customized_mutation_colors_switch_to_custom_gradient(self, tmp_path, monkeypatch):
        kwargs = self._kwargs(tmp_path, monkeypatch)
        result = mod.run_export(
            mutation_heatmap=True, mutation_low_color="#000000", mutation_high_color="#ffff00", **kwargs,
        )
        text = result.mutation_chimerax_script_out.read_text()
        assert "palette #000000:#ffff00" in text, (
            f"customized low/high colors should build a 'low:high' custom palette instead "
            f"of 'Reds', got:\n{text}"
        )
        assert "key #000000:" in text and "#ffff00:" in text, (
            f"the key should also switch to the same custom colors, got:\n{text}"
        )
        assert "palette Reds" not in text

    def test_customized_plddt_colors_switch_to_custom_gradient(self, tmp_path, monkeypatch):
        kwargs = self._kwargs(tmp_path, monkeypatch)
        result = mod.run_export(
            mutation_heatmap=False, plddt_heatmap=True,
            plddt_low_color="purple", plddt_high_color="cyan", **kwargs,
        )
        text = result.plddt_chimerax_script_out.read_text()
        assert "palette purple:cyan" in text, (
            f"customized pLDDT low/high colors should build a 'low:high' custom palette "
            f"instead of 'alphafold', got:\n{text}"
        )
        assert "key purple:0 cyan:100" in text, (
            f"the pLDDT key should switch to the same custom colors at the fixed 0/100 "
            f"range, got:\n{text}"
        )

    def test_changing_only_one_color_still_activates_the_custom_gradient(self, tmp_path, monkeypatch):
        # A swatch UI reports each color independently -- if the user has only
        # touched one of the two (the other still sitting at its default-looking
        # sentinel value), the custom gradient must still activate using BOTH
        # current values, not silently keep rendering "Reds" while the swatch
        # showing #123456 implies otherwise.
        kwargs = self._kwargs(tmp_path, monkeypatch)
        result = mod.run_export(
            mutation_heatmap=True, mutation_low_color="#123456", **kwargs,
        )
        text = result.mutation_chimerax_script_out.read_text()
        assert f"palette #123456:{mod.MUTATION_DEFAULT_HIGH_COLOR}" in text, (
            f"changing just the low color should still switch away from 'Reds', using the "
            f"still-default high color as the other end of the gradient, got:\n{text}"
        )
        assert "palette Reds" not in text

    def test_customized_marker_colors_appear_in_the_script(self, tmp_path, monkeypatch):
        kwargs = self._kwargs(tmp_path, monkeypatch)
        result = mod.run_export(
            mutation_heatmap=True, mark_ptm_sites=True, mark_mutations=True,
            ptm_marker_color="magenta", mutation_marker_color="#00ffff", **kwargs,
        )
        text = result.mutation_chimerax_script_out.read_text()
        assert "color magenta name" in text, (
            f"a customized PTM marker color should appear on the sphere command, got:\n{text}"
        )
        assert "#00ffff target ab" in text, (
            f"a customized mutation marker color should appear on the stick color command, got:\n{text}"
        )
        assert "color green name" not in text and "orange target ab" not in text, (
            f"the old default marker colors should not appear once overridden, got:\n{text}"
        )


class TestRunBatchExport:
    def _setup_two_proteins(self, tmp_path, monkeypatch):
        """Two single-fragment proteins (P04637/TP53, P00533/EGFR) with a
        gene cache so neither needs a live network call to resolve.
        """
        monkeypatch.setattr(mod, "MODELS_ROOT", tmp_path / "cif_models")
        for uid in ("P04637", "P00533"):
            uid_dir = tmp_path / "cif_models" / uid
            uid_dir.mkdir(parents=True)
            _write_synthetic_cif(
                uid_dir / f"AF-{uid}-F1-model_v4.cif",
                res_ids=[1], res_names=["ALA"], atom_names=["CA"], coords=[[0.0, 0.0, 0.0]],
            )
        gene_cache = tmp_path / "gene_cache.tsv"
        pd.DataFrame([
            {"UniProt": "P04637", "gene": "TP53"},
            {"UniProt": "P00533", "gene": "EGFR"},
        ]).to_csv(gene_cache, sep="\t", index=False)
        monkeypatch.setattr(mod, "GENE_CACHE", gene_cache)

        cosmic = tmp_path / "cosmic.tsv"
        pd.DataFrame(columns=["GENE_SYMBOL", "MUTATION_AA", "COSMIC_SAMPLE_ID", "MUTATION_SOMATIC_STATUS"]).to_csv(
            cosmic, sep="\t", index=False,
        )
        return cosmic

    def test_runs_one_export_per_token_into_separate_folders(self, tmp_path, monkeypatch):
        cosmic = self._setup_two_proteins(tmp_path, monkeypatch)
        items = mod.run_batch_export(
            ["P04637", "P00533"], cosmic_file=cosmic, output_dir=tmp_path / "out",
            log_cb=lambda *_: None,
        )
        assert [item.token for item in items] == ["P04637", "P00533"], (
            "results should be returned in the same order as the input tokens"
        )
        assert all(item.result is not None and item.error is None for item in items), (
            f"both proteins should export successfully, got "
            f"{[(i.token, i.error) for i in items]}"
        )
        assert {item.result.uid for item in items} == {"P04637", "P00533"}
        assert items[0].result.all_out.parent.name == "TP53_P04637"
        assert items[1].result.all_out.parent.name == "EGFR_P00533", (
            "each protein's outputs must land in their own {gene}_{uid} subfolder so a "
            "batch never mixes proteins' files together, and it's identifiable by gene "
            "name at a glance"
        )

    def test_classifies_uniprot_vs_gene_tokens(self, tmp_path, monkeypatch):
        cosmic = self._setup_two_proteins(tmp_path, monkeypatch)
        calls = []
        real_run_export = mod.run_export

        def spy_run_export(uniprot=None, gene=None, **kwargs):
            calls.append((uniprot, gene))
            return real_run_export(uniprot=uniprot, gene=gene, **kwargs)

        monkeypatch.setattr(mod, "run_export", spy_run_export)
        mod.run_batch_export(
            ["P04637", "EGFR"], cosmic_file=cosmic, output_dir=tmp_path / "out",
            log_cb=lambda *_: None,
        )
        assert calls == [("P04637", None), (None, "EGFR")], (
            f"a token matching the UniProt accession format should be routed to uniprot=, "
            f"and a token that doesn't should be routed to gene=, got {calls}"
        )

    def test_one_failure_does_not_abort_the_batch(self, tmp_path, monkeypatch):
        cosmic = self._setup_two_proteins(tmp_path, monkeypatch)
        # NOTFOUND1 has no AlphaFold CIF and isn't in the gene cache, so its
        # export should fail cleanly -- P00533 (after it) must still run.
        items = mod.run_batch_export(
            ["NOTFOUND1", "P00533"], cosmic_file=cosmic, output_dir=tmp_path / "out",
            log_cb=lambda *_: None,
        )
        assert items[0].result is None and items[0].error, (
            f"the unresolvable token should fail with an error message, not raise out of "
            f"run_batch_export, got {items[0]}"
        )
        assert items[1].result is not None and items[1].error is None, (
            "a failure on the first protein must not prevent the second from running"
        )

    def test_progress_cb_called_with_index_total_token(self, tmp_path, monkeypatch):
        cosmic = self._setup_two_proteins(tmp_path, monkeypatch)
        calls = []
        mod.run_batch_export(
            ["P04637", "P00533"], cosmic_file=cosmic, output_dir=tmp_path / "out",
            progress_cb=lambda i, total, token: calls.append((i, total, token)),
            log_cb=lambda *_: None,
        )
        assert calls == [(1, 2, "P04637"), (2, 2, "P00533")], (
            f"progress_cb should fire once per protein with a 1-based index, the batch "
            f"total, and that protein's own token, got {calls}"
        )

    def test_empty_token_list_returns_empty_list(self, tmp_path, monkeypatch):
        cosmic = self._setup_two_proteins(tmp_path, monkeypatch)
        assert mod.run_batch_export([], cosmic_file=cosmic, output_dir=tmp_path / "out", log_cb=lambda *_: None) == []
