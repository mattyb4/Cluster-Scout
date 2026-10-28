"""Tests for ClinVar as a mutation source: scripts/clinvar_input.py (reading
variant_summary), its step-1 wiring in scripts/1_filter.py, and the ClinVar
file/step-label helpers in scripts/pipeline_utils.py. ClinVar input is a small
synthetic variant_summary file."""
import gzip
import sys
from pathlib import Path

import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "scripts"))
import clinvar_input as ci  # noqa: E402

HEADER = ["#AlleleID", "Type", "Name", "GeneID", "GeneSymbol", "ClinicalSignificance",
          "Assembly", "NumberSubmitters"]


def row(name, gene, significance, submitters, assembly="GRCh38"):
    return ["1", "single nucleotide variant", name, "0", gene, significance, assembly, str(submitters)]


def write_summary(path, rows):
    text = "\n".join("\t".join(r) for r in [HEADER, *rows]) + "\n"
    if str(path).endswith(".gz"):
        with gzip.open(path, "wt", encoding="utf-8") as handle:
            handle.write(text)
    else:
        path.write_text(text, encoding="utf-8")
    return path


ROWS = [
    row("NM_000546.6(TP53):c.524G>A (p.Arg175His)", "TP53", "Pathogenic", 6),
    # same protein change from a different codon change: submitters are summed
    row("NM_000546.6(TP53):c.525C>G (p.Arg175His)", "TP53", "Likely pathogenic", 2),
    row("NM_000546.6(TP53):c.743G>A (p.Arg248Gln)", "TP53", "Pathogenic/Likely pathogenic", 3),
    row("NM_000546.6(TP53):c.700A>G (p.Tyr234Cys)", "TP53", "Uncertain significance", 4),
    row("NM_000059.4(BRCA2):c.1A>G (p.Met1Val)", "BRCA2", "Benign", 1),
    # not single-residue substitutions
    row("NM_000546.6(TP53):c.1000C>T (p.Arg334Ter)", "TP53", "Pathogenic", 5),
    row("NM_000546.6(TP53):c.30G>A (p.Lys10=)", "TP53", "Pathogenic", 5),
    row("NM_000546.6(TP53):c.400del (p.Gly134fs)", "TP53", "Pathogenic", 5),
    row("NM_000546.6(TP53):c.4G>A", "TP53", "Pathogenic", 5),
    # same variant on the other assembly must not be double counted
    row("NM_000546.6(TP53):c.524G>A (p.Arg175His)", "TP53", "Pathogenic", 6, assembly="GRCh37"),
    # a variant spanning two genes takes the gene named in its transcript
    row("NM_000001.1(GENEA):c.10G>A (p.Gly4Asp)", "GENEA;GENEB", "Pathogenic", 0),
]


@pytest.fixture
def summary(tmp_path):
    return write_summary(tmp_path / "variant_summary.txt", ROWS)


def counts_dict(counts):
    return {(r.gene, r.aa_change): r.affected_cases for r in counts.itertuples()}


class TestParseProteinChanges:
    def test_converts_three_letter_substitutions_to_one_letter(self):
        names = pd.Series(["NM_000546.6(TP53):c.524G>A (p.Arg175His)"])
        assert ci.parse_protein_changes(names).tolist() == ["R175H"], \
            "the pipeline works in one-letter notation, so Arg175His must become R175H"

    def test_rejects_synonymous_stop_frameshift_and_missing(self):
        names = pd.Series([
            "NM_1(G):c.1A>G (p.Lys10=)", "NM_1(G):c.1A>T (p.Arg5Ter)",
            "NM_1(G):c.1del (p.Gly134fs)", "NM_1(G):c.1A>G", "NM_1(G):c.1A>G (p.Arg5Arg)",
            "NM_1(G):c.1A>G (p.Sec5Cys)",
        ])
        assert ci.parse_protein_changes(names).tolist() == [""] * 6, \
            "only a change from one standard residue to a different one is a missense mutation"


class TestLoadClinVarCounts:
    def test_default_classes_sum_submitters_per_protein_change(self, summary, tmp_path):
        counts, totals = ci.load_clinvar_counts(summary, tmp_path / "cache", ("P", "LP"))
        got = counts_dict(counts)
        assert got[("TP53", "R175H")] == 8, "6 Pathogenic + 2 Likely pathogenic submitters for one protein change"
        assert got[("TP53", "R248Q")] == 3, "a combined Pathogenic/Likely pathogenic call counts as both"
        assert ("TP53", "Y234C") not in got, "Uncertain significance is outside the default P+LP classes"
        assert ("BRCA2", "M1V") not in got, "Benign is outside the default classes"
        assert len(got) == 3, f"only three missense changes qualify here, got {sorted(got)}"
        assert totals["TP53"] == 11, "gene total = every selected-class missense submitter in the gene"

    def test_grch37_rows_are_ignored(self, summary, tmp_path):
        counts, _ = ci.load_clinvar_counts(summary, tmp_path / "cache", ("P",))
        assert counts_dict(counts)[("TP53", "R175H")] == 6, "the GRCh37 copy of the variant must not add submitters"

    def test_multi_gene_variant_uses_gene_from_its_name_and_zero_submitters_count_once(self, summary, tmp_path):
        got = counts_dict(ci.load_clinvar_counts(summary, tmp_path / "cache", ("P",))[0])
        assert got[("GENEA", "G4D")] == 1, "a variant reporting 0 submitters was still submitted once"

    def test_extra_classes_widen_the_selection(self, summary, tmp_path):
        counts, _ = ci.load_clinvar_counts(summary, tmp_path / "cache", ("P", "LP", "VUS", "B"))
        got = counts_dict(counts)
        assert got[("TP53", "Y234C")] == 4 and got[("BRCA2", "M1V")] == 1

    def test_unknown_class_is_an_error(self, summary, tmp_path):
        with pytest.raises(ci.ClinVarInputError):
            ci.load_clinvar_counts(summary, tmp_path / "cache", ("P", "NOPE"))

    def test_reads_gzipped_files(self, tmp_path):
        gz = write_summary(tmp_path / "variant_summary.txt.gz", ROWS)
        counts, _ = ci.load_clinvar_counts(gz, tmp_path / "cache", ("P",))
        assert counts_dict(counts)[("TP53", "R175H")] == 6

    def test_file_without_grch38_missense_is_rejected(self, tmp_path):
        path = write_summary(tmp_path / "variant_summary.txt", [ROWS[0][:6] + ["GRCh37", "1"]])
        with pytest.raises(ci.ClinVarInputError):
            ci.load_clinvar_counts(path, tmp_path / "cache", ("P",))


class TestCache:
    def test_second_load_comes_from_the_cache(self, summary, tmp_path, monkeypatch):
        cache = tmp_path / "cache"
        ci.load_clinvar_counts(summary, cache, ("P",))
        monkeypatch.setattr(ci, "_read_variant_summary",
                            lambda path: pytest.fail("the unchanged file should not be re-read"))
        counts, _ = ci.load_clinvar_counts(summary, cache, ("P", "LP"))
        assert counts_dict(counts)[("TP53", "R175H")] == 8, "cached table must support a different class selection"

    def test_changed_file_rebuilds_the_cache(self, summary, tmp_path):
        cache = tmp_path / "cache"
        ci.load_clinvar_counts(summary, cache, ("P",))
        write_summary(summary, [*ROWS, row("NM_9.1(NEWGENE):c.1A>G (p.Ala2Gly)", "NEWGENE", "Pathogenic", 2)])
        counts, _ = ci.load_clinvar_counts(summary, cache, ("P",))
        assert ("NEWGENE", "A2G") in counts_dict(counts), "a newer variant_summary must not be shadowed by stale cache"


class TestStepOneIntegration:
    def test_clinvar_hotspots_use_the_min_samples_threshold(self, filter_module, summary, tmp_path, monkeypatch):
        monkeypatch.setattr(filter_module, "CACHE_DIR", tmp_path / "cache")
        monkeypatch.setattr(filter_module, "HOTSPOT_MIN_AFFECTED_CASES", 3)
        hotspots, transcripts, totals = filter_module._load_and_filter_clinvar(summary, ("P", "LP"))
        assert sorted(hotspots["mutation_with_count"]) == ["R175H (8)", "R248Q (3)"], \
            "mutations with fewer than 3 submitters are not hotspots"
        assert transcripts == {} and totals["TP53"] == 11


class TestPipelineUtilsClinVar:
    def test_validate_clinvar_file(self, pipeline_utils_module, summary, tmp_path):
        assert pipeline_utils_module.validate_clinvar_file(summary) == []
        bad = tmp_path / "bad.txt"
        bad.write_text("a\tb\n1\t2\n")
        problems = pipeline_utils_module.validate_clinvar_file(bad)
        assert any("GeneSymbol" in p for p in problems), "a file missing ClinVar's columns must be reported"

    def test_clinvar_paths_never_collide_with_cosmic_paths(self, pipeline_utils_module, tmp_path):
        mod = pipeline_utils_module
        for mode in ("ptm-proximity", "mutation-clustering"):
            for ptm in ("ptmd", "psp"):
                assert mod.hotspots_tsv_path(tmp_path, mode, ptm, "clinvar") != \
                    mod.hotspots_tsv_path(tmp_path, mode, ptm, "cosmic")
        assert mod.hotspots_tsv_path(tmp_path, "ptm-proximity", "psp", "clinvar").name == \
            "PSP_ClinVar_hotspots_by_protein.tsv"
        cosmic, clinvar = mod.ptm_output_paths(tmp_path, "psp"), mod.ptm_output_paths(tmp_path, "psp", "clinvar")
        assert all(cosmic[k] != clinvar[k] for k in cosmic), "ClinVar runs must not overwrite COSMIC output"
        assert clinvar["db"].name == "clinvar_psp_ptm_mutation_proximity_db.tsv"
        assert mod.cluster_output_paths(tmp_path, "clinvar")["db"].name == "clinvar_mutation_cluster_db.tsv"

    def test_cosmic_paths_are_unchanged(self, pipeline_utils_module, tmp_path):
        mod = pipeline_utils_module
        assert mod.hotspots_tsv_path(tmp_path, "ptm-proximity").name == "PTMD_COSMIC_hotspots_by_protein.tsv"
        assert mod.ptm_output_paths(tmp_path)["db"].name == "ptm_mutation_proximity_db.tsv"
        assert mod.ptm_output_paths(tmp_path)["unmatched"].name == "ptm_genes_without_cosmic_mutations.tsv"
        assert mod.cluster_output_paths(tmp_path)["long"].name == "mutation_cluster_long.tsv"

    def test_step_labels_name_the_mutation_source(self, pipeline_utils_module):
        mod = pipeline_utils_module
        assert mod.pipeline_steps("ptm-proximity", "psp", "clinvar")[0][0] == \
            "Filter and merge PhosphoSitePlus + ClinVar data"
        assert mod.pipeline_steps("mutation-clustering", "ptmd", "clinvar")[0][0] == "Filter ClinVar hotspot mutations"
        assert mod.pipeline_steps("ptm-proximity")[0][0] == "Filter and merge PTMD + COSMIC data"


class TestResultsLabelFollowsMutationSource:
    def test_cosmic_label_headings_read_clinvar_for_clinvar_results(self):
        from types import SimpleNamespace

        from ui.results_tab import ResultsTabMixin
        tab = ResultsTabMixin.__new__(ResultsTabMixin)
        tab._results_mutation_source_var = SimpleNamespace(get=lambda: "ClinVar")
        assert tab._source_text("COSMIC label") == "ClinVar label"
        assert tab._source_text("Anchor COSMIC label") == "Anchor ClinVar label"
        assert tab._source_text("Mutation") == "Mutation", "columns that don't name the source are untouched"
        tab._results_mutation_source_var = SimpleNamespace(get=lambda: "COSMIC")
        assert tab._source_text("COSMIC label") == "COSMIC label"
