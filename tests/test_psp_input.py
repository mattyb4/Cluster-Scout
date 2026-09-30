"""Unit tests for scripts/psp_input.py, the PhosphoSitePlus input reader.

Builds small synthetic PSP-format files in tmp_path rather than reading a real
download: the tests pin down the parts meant to survive a future PSP release
(files recognized by columns not names, a preamble of any length, .gz files,
unknown PTM suffixes reported rather than dropped).
"""
import gzip
import sys

import pytest
from conftest import SCRIPTS_DIR

sys.path.insert(0, str(SCRIPTS_DIR))
import psp_input  # noqa: E402

SITE_HEADER = ("GENE\tPROTEIN\tACC_ID\tHU_CHR_LOC\tMOD_RSD\tSITE_GRP_ID\tORGANISM\tMW_kD\t"
               "DOMAIN\tSITE_+/-7_AA\tLT_LIT\tMS_LIT\tMS_CST\tCST_CAT#")
DISEASE_HEADER = ("DISEASE\tALTERATION\tGENE\tPROTEIN\tACC_ID\tGENE_ID\tHU_CHR_LOC\tMW_kD\t"
                  "ORGANISM\tSITE_GRP_ID\tMOD_RSD\tDOMAIN\tSITE_+/-7_AA\tPMIDs\tLT_LIT\t"
                  "MS_LIT\tMS_CST\tCST_CAT#\tNOTES")
PREAMBLE = ["October 29 2021", "PhosphoSitePlus(R) (PSP) was created by Cell Signaling Technology Inc.", ""]


def site_row(gene, acc, mod_rsd, organism="human", ltp="", htp="1", cst=""):
    return "\t".join([gene, gene, acc, "1q1", mod_rsd, "1", organism, "50", "", "xxxxxxxSxxxxxxx",
                      ltp, htp, cst, ""])


def disease_row(disease, alteration, acc, mod_rsd, organism="human"):
    return "\t".join([disease, alteration, "G", "G", acc, "1", "1q1", "50", organism, "1",
                      mod_rsd, "", "xxxxxxxSxxxxxxx", "123", "1", "0", "0", "", ""])


def write(path, header, rows, preamble=PREAMBLE, gz=False):
    text = "\n".join([*preamble, header, *rows]) + "\n"
    path.parent.mkdir(parents=True, exist_ok=True)
    if gz:
        with gzip.open(path, "wt", encoding="latin-1") as fh:
            fh.write(text)
    else:
        path.write_text(text, encoding="latin-1")
    return path


class TestDiscovery:
    def test_files_are_recognized_by_columns_not_names(self, tmp_path):
        write(tmp_path / "renamed (1).txt", SITE_HEADER, [site_row("AKT1", "P31749", "S473-p")])
        write(tmp_path / "whatever", DISEASE_HEADER, [disease_row("breast cancer", "increased", "P31749", "S473-p")])
        (tmp_path / "notes.txt").write_text("not a PSP file\n")

        files = psp_input.discover_psp_files(tmp_path)

        assert [p.name for p in files.site] == ["renamed (1).txt"], (
            f"a browser-renamed site dataset must still be recognized by its columns, got {files.site}"
        )
        assert [p.name for p in files.disease] == ["whatever"], (
            f"the disease file has no telltale name here -- its DISEASE/ALTERATION columns "
            f"should identify it, got {files.disease}"
        )
        assert [p.name for p in files.ignored] == ["notes.txt"], (
            f"unrecognized files should be listed as ignored, got {files.ignored}"
        )

    def test_nested_folders_and_gz_files_are_found(self, tmp_path):
        write(tmp_path / "Phosphorylation_site_dataset" / "Phosphorylation_site_dataset.gz",
              SITE_HEADER, [site_row("AKT1", "P31749", "S473-p")], gz=True)
        files = psp_input.discover_psp_files(tmp_path)
        assert len(files.site) == 1, (
            f"a gzipped site file inside an extra folder level (as unzipping produces) "
            f"should be found, got {files}"
        )

    def test_header_is_found_after_a_longer_preamble(self, tmp_path):
        preamble = PREAMBLE + ["An extra license paragraph a future release might add", ""]
        write(tmp_path / "sites", SITE_HEADER, [site_row("AKT1", "P31749", "S473-p")], preamble=preamble)
        sites = psp_input.load_psp_sites(tmp_path, log=lambda *_: None)
        assert len(sites) == 1, (
            "the header line should be located by its columns, not a fixed line count, "
            f"so a longer preamble still reads correctly -- got {len(sites)} sites"
        )


class TestValidation:
    def test_no_site_files_is_an_error(self, tmp_path):
        write(tmp_path / "diseases", DISEASE_HEADER, [])
        problems = psp_input.validate_psp_folder(tmp_path)
        assert problems and "No PhosphoSitePlus site dataset" in problems[0], (
            f"a folder without any site dataset can't be used and must say so, got {problems}"
        )

    def test_two_files_with_the_same_ptm_type_is_an_error(self, tmp_path):
        write(tmp_path / "phospho_2021", SITE_HEADER, [site_row("AKT1", "P31749", "S473-p")])
        write(tmp_path / "phospho_2025", SITE_HEADER, [site_row("AKT1", "P31749", "T308-p")])
        problems = psp_input.validate_psp_folder(tmp_path)
        assert any("both contain Phosphorylation" in p for p in problems), (
            f"two phosphorylation files (e.g. an old and a new download) would double-count "
            f"sites and must be flagged, got {problems}"
        )

    def test_valid_folder_has_no_problems(self, tmp_path):
        write(tmp_path / "phospho", SITE_HEADER, [site_row("AKT1", "P31749", "S473-p")])
        write(tmp_path / "ubiq", SITE_HEADER, [site_row("AKT1", "P31749", "K8-ub")])
        assert psp_input.validate_psp_folder(tmp_path) == []

    def test_missing_required_column_is_named(self, tmp_path):
        header = SITE_HEADER.replace("\tMS_CST", "")
        write(tmp_path / "phospho", header, [])
        with pytest.raises(psp_input.PspInputError, match="MS_CST"):
            psp_input.read_psp_table(tmp_path / "phospho", psp_input.SITE_COLUMNS)


class TestLoadSites:
    def _load(self, tmp_path, rows):
        write(tmp_path / "sites", SITE_HEADER, rows)
        logged = []
        sites = psp_input.load_psp_sites(tmp_path, log=logged.append)
        return sites, "\n".join(logged)

    def test_parses_site_type_and_scores(self, tmp_path):
        sites, _ = self._load(tmp_path, [site_row("AKT1", "P31749", "S473-p", ltp="12", htp="40", cst="3")])
        row = sites.iloc[0]
        assert (row["UniProt"], row["residue"], row["position"], row["ptm_type"]) == (
            "P31749", "S", 473, "Phosphorylation"
        )
        assert (row["ltp"], row["htp"], row["cst"]) == (12, 40, 3), (
            f"LT_LIT/MS_LIT/MS_CST should become ltp/htp/cst, got {row.to_dict()}"
        )

    def test_blank_score_means_zero(self, tmp_path):
        sites, _ = self._load(tmp_path, [site_row("AKT1", "P31749", "S473-p", ltp="", htp="", cst="")])
        assert (sites.iloc[0]["ltp"], sites.iloc[0]["htp"], sites.iloc[0]["cst"]) == (0, 0, 0), (
            "PSP leaves a count blank when it's zero -- blank must read as 0, not NaN"
        )

    def test_only_human_canonical_uniprot_sites_are_kept(self, tmp_path):
        sites, log = self._load(tmp_path, [
            site_row("AKT1", "P31749", "S473-p"),
            site_row("Akt1", "P31750", "S473-p", organism="mouse"),
            site_row("AKT1", "P31749-2", "S450-p"),
            site_row("X", "AAA58698", "S10-p"),
        ])
        assert list(sites["UniProt"]) == ["P31749"], (
            f"mouse rows, isoform accessions and non-UniProt IDs must all be dropped, got {list(sites['UniProt'])}"
        )
        assert "isoform accession" in log and "not a UniProt accession" in log, (
            f"skipped rows must be counted in the log by reason, got:\n{log}"
        )

    def test_unknown_ptm_suffix_is_reported_not_silently_dropped(self, tmp_path):
        sites, log = self._load(tmp_path, [
            site_row("AKT1", "P31749", "S473-p"),
            site_row("AKT1", "P31749", "K20-zz"),
        ])
        assert len(sites) == 1
        assert "unknown PTM type suffix '-zz'" in log, (
            f"a PTM type a future release adds must show up in the log, got:\n{log}"
        )

    def test_mostly_unparseable_file_is_an_error(self, tmp_path):
        rows = [site_row("AKT1", "P31749", "Ser473")] * 10 + [site_row("AKT1", "P31749", "S473-p")]
        write(tmp_path / "sites", SITE_HEADER, rows)
        with pytest.raises(psp_input.PspInputError, match="unrecognized MOD_RSD format"):
            psp_input.load_psp_sites(tmp_path, log=lambda *_: None)

    def test_methylation_suffixes_map_to_distinct_types(self, tmp_path):
        sites, _ = self._load(tmp_path, [
            site_row("H3", "P68431", "K4-m1"),
            site_row("H3", "P68431", "K9-m3"),
        ])
        assert set(sites["ptm_type"]) == {"Monomethylation", "Trimethylation"}


class TestLoadDiseases:
    def test_alteration_is_appended_and_absent_file_is_fine(self, tmp_path):
        write(tmp_path / "sites", SITE_HEADER, [site_row("AKT1", "P31749", "S473-p")])
        empty = psp_input.load_psp_diseases(tmp_path, log=lambda *_: None)
        assert empty.empty, "the disease file is optional -- no file means an empty table, not an error"

        write(tmp_path / "diseases", DISEASE_HEADER, [
            disease_row("breast cancer", "increased", "P31749", "S473-p"),
            disease_row("Alzheimer's disease", "", "P31749", "T308-p"),
            disease_row("ignored", "increased", "P31749-2", "S473-p"),
        ])
        diseases = psp_input.load_psp_diseases(tmp_path, log=lambda *_: None)
        assert sorted(diseases["disease"]) == ["Alzheimer's disease", "breast cancer (increased)"], (
            f"alteration should be appended in parentheses when present, and isoform rows "
            f"dropped, got {sorted(diseases['disease'])}"
        )


class TestFilterSites:
    @pytest.fixture
    def loaded(self, tmp_path):
        write(tmp_path / "sites", SITE_HEADER, [
            site_row("AKT1", "P31749", "S473-p", ltp="12", htp="40"),
            site_row("AKT1", "P31749", "T308-p", ltp="0", htp="3"),
            site_row("AKT1", "P31749", "K8-ub", ltp="2", htp=""),
        ])
        write(tmp_path / "diseases", DISEASE_HEADER, [
            disease_row("breast cancer", "increased", "P31749", "S473-p"),
        ])
        quiet = {"log": lambda *_: None}
        return psp_input.load_psp_sites(tmp_path, **quiet), psp_input.load_psp_diseases(tmp_path, **quiet)

    def _sites(self, loaded, **filters):
        sites, diseases = loaded
        kept = psp_input.filter_psp_sites(sites, diseases, log=lambda *_: None, **filters)
        return sorted(kept["residue"] + kept["position"].astype(str))

    def test_defaults_keep_everything(self, loaded):
        assert self._sites(loaded) == ["K8", "S473", "T308"], (
            "default bounds (0 to no limit, disease filter off) must not drop any site"
        )

    def test_ranges_are_inclusive(self, loaded):
        assert self._sites(loaded, min_ltp=2, max_ltp=12) == ["K8", "S473"]
        assert self._sites(loaded, min_htp=3, max_htp=3) == ["T308"], (
            "both bounds are inclusive, so min=max=3 keeps exactly the HTP-3 site"
        )

    def test_blank_htp_counts_as_zero(self, loaded):
        assert self._sites(loaded, max_htp=0) == ["K8"]

    def test_disease_only(self, loaded):
        assert self._sites(loaded, disease_only=True) == ["S473"]

    def test_disease_only_without_disease_file_is_an_error(self, loaded):
        sites, diseases = loaded
        with pytest.raises(psp_input.PspInputError, match="Disease-associated_sites"):
            psp_input.filter_psp_sites(sites, diseases.iloc[0:0], disease_only=True, log=lambda *_: None)
