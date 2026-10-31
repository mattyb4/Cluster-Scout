"""Tests for the runtime precheck in ui/pipeline_runner.py (the estimate shown in the log and the timer).

The precheck reads caches and input files from PROJECT_ROOT, so each test builds a tiny
project under tmp_path and points the module at it.
"""
import sys
from pathlib import Path

import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

pytest.importorskip("customtkinter")
import ui.pipeline_runner as pr  # noqa: E402
from ui.common import COSMIC_INPUT_DIR, PTMD_INPUT_DIR, input_dir  # noqa: E402


class Runner(pr.PipelineRunnerMixin):
    class _Var:
        def get(self):
            return "3"

    _min_samples_var = _Var()

    def __init__(self, output_dir):
        self._out = output_dir
        self.logs = []

    @property
    def _output_dir(self):
        return self._out

    def _q(self, kind, *args):
        if kind == "log":
            self.logs.append(args[0])


def write_db(path, rows):
    """A step 3 proximity table (UTF-16 like the real one) with the columns the precheck reads."""
    df = pd.DataFrame(rows, columns=["UniProt", "gene", "mutations_within_5_positions",
                                     "mutations_more_than_5_positions"])
    df.to_csv(path, sep="\t", index=False, encoding="utf-16")


class TestUncachedPolyphenPairs:
    def test_counts_pairs_missing_from_the_cache(self, tmp_path):
        db = tmp_path / "db.tsv"
        write_db(db, [("A", "G1", "R12H(PP:D,0.9)-3.1Å, K5E-4.0Å", "T9M-12.0Å"),
                      ("B", "G2", "", "S1P-9.0Å")])
        cached = {("G1", "R12H"), ("G9", "R12H")}
        assert Runner._uncached_pp_pairs(db, cached) == 3, (
            "of the 4 distinct (gene, mutation) pairs only G1/R12H is cached; the extra cached "
            "pair for another gene must not count"
        )

    def test_missing_table_gives_none(self, tmp_path):
        assert Runner._uncached_pp_pairs(tmp_path / "nope.tsv", set()) is None, (
            "with no earlier run's table the caller needs None to fall back to an estimate"
        )


@pytest.fixture
def project(tmp_path, monkeypatch):
    monkeypatch.setattr(pr, "PROJECT_ROOT", tmp_path)
    cache = tmp_path / "data" / "cache"
    cache.mkdir(parents=True)
    (tmp_path / "Output" / "logs").mkdir(parents=True)

    # PTMD knows proteins A-D; only A and B (genes G1, G2) have a COSMIC hotspot
    ptmd_dir = input_dir(tmp_path, PTMD_INPUT_DIR)
    ptmd_dir.mkdir(parents=True, exist_ok=True)
    pd.DataFrame({"UniProt": ["A", "B", "C", "D"], "State": ["N"] * 4}).to_csv(
        ptmd_dir / "ptmd.tsv", sep="\t", index=False)
    cosmic_dir = input_dir(tmp_path, COSMIC_INPUT_DIR)
    cosmic_dir.mkdir(parents=True, exist_ok=True)
    rows = [(g, "p.R1H", f"S{i}", "Confirmed somatic variant") for g in ("G1", "G2", "G9") for i in range(3)]
    pd.DataFrame(rows, columns=["GENE_SYMBOL", "MUTATION_AA", "COSMIC_SAMPLE_ID", "MUTATION_SOMATIC_STATUS"]).to_csv(
        cosmic_dir / "cosmic.tsv", sep="\t", index=False)

    # accession -> gene is cached for every PTMD protein; G9 isn't a PTMD gene
    pd.DataFrame({"UniProt": ["A", "B", "C", "D"], "gene": ["G1", "G2", "G3", "G4"]}).to_csv(
        cache / "uniprot_gene_mapping.tsv", sep="\t", index=False)
    return tmp_path


def precheck(project):
    runner = Runner(project / "Output")
    steps, n = runner._run_precheck("ptm-proximity", "ptmd", "cosmic")
    return runner, steps, n


class TestFirstRunOfANewCombination:
    def test_protein_list_is_predicted_from_hotspot_genes_and_the_gene_cache(self, project):
        runner, _steps, n = precheck(project)
        assert n == 2, (
            "G1 and G2 are hotspot genes with PTM data (proteins A and B); G9 has no PTM "
            f"protein and C/D have no hotspot, so 2 proteins are expected, got {n}"
        )
        assert any("2 have PTM data" in line for line in runner.logs)

    def test_only_this_runs_proteins_count_as_cached(self, project):
        # A has a CIF; Z is cached from some other run and must not count
        for acc in ("A", "Z"):
            (project / "cif_models" / acc).mkdir(parents=True)
            (project / "cif_models" / acc / f"AF-{acc}-F1-model_v6.cif").write_text("")
        runner, _steps, _n = precheck(project)
        assert any("Step 2: 1 CIF files present, ~1 to download" in line for line in runner.logs), (
            "only A (of the predicted A and B) is cached; the unrelated Z must not hide B's "
            f"download, got {[line for line in runner.logs if line.startswith('Step 2')]}"
        )

    def test_unavailable_proteins_are_not_charged_as_downloads(self, project):
        report = project / "Output" / "logs" / "download_report.tsv"
        pd.DataFrame({"UniProt": ["B"], "status": ["NO_ENTRY"]}).to_csv(report, sep="\t", index=False)
        runner, _steps, _n = precheck(project)
        assert any("Step 2: 0 CIF files present, ~1 to download" in line for line in runner.logs), (
            "B has no AlphaFold model, so only A still needs a download"
        )

    def test_without_a_gene_cache_it_falls_back_to_the_gene_count(self, project):
        (project / "data" / "cache" / "uniprot_gene_mapping.tsv").unlink()
        _runner, _steps, n = precheck(project)
        assert n == 3, f"with nothing to match accessions on, the 3 hotspot genes are the best guess, got {n}"


class TestPolyphenEstimate:
    def test_exact_pairs_from_the_previous_table_set_the_cost(self, project):
        write_db(project / "Output" / "ptm_mutation_proximity_db.tsv",
                 [("A", "G1", "R1H-3.0Å, K2E-4.0Å", ""), ("B", "G2", "T3M-5.0Å", "")])
        pd.DataFrame({"gene": ["G1"], "mutation": ["R1H"], "pred": ["D"], "score": ["0.9"]}).to_csv(
            project / "data" / "cache" / "polyphen.tsv", sep="\t", index=False)
        _runner, steps, _n = precheck(project)
        # 3 distinct pairs, 1 cached -> 2 to fetch; the other step 4 phases are tiny here
        pp_seconds = 2 * pr.PipelineRunnerMixin._TIME_PER_PP_FETCH
        assert steps[3] >= pp_seconds + pr.PipelineRunnerMixin._TIME_STEP4_BASE - 1e-6, (
            f"step 4 should include the 2 uncached pairs ({pp_seconds}s) on top of the base, got {steps[3]}"
        )
