"""Tests for scripts/validate_isoform_choice.py (the check of the pipeline's
isoform choice against COSMIC's own Ensembl protein)."""
import sys

import pandas as pd
import pytest
from conftest import SCRIPTS_DIR

sys.path.insert(0, str(SCRIPTS_DIR))
import validate_isoform_choice as vic  # noqa: E402

CANONICAL = "MSTAPKLRVE" + "QDGSHNCYFI" + "PTRKAELGVM" + "DNSQYHCFAW"
ISOFORM = CANONICAL[:10] + "WWWWWW" + CANONICAL[10:]


def label(seq, pos, alt="A"):
    return f"{seq[pos - 1]}{pos}{alt}"


class TestCosmicLabels:
    def test_undoes_renumbering(self):
        row = {"mutations_on_protein": "R249Q (5); K10E (2)", "cosmic_mutation_labels": "R249Q=R248Q"}
        assert vic.cosmic_labels(row) == ["R248Q", "K10E"]

    def test_blank_labels_column(self):
        row = {"mutations_on_protein": "R249Q (5)", "cosmic_mutation_labels": ""}
        assert vic.cosmic_labels(row) == ["R249Q"]


class TestReferenceFit:
    def test_fraction_of_matching_residues(self):
        assert vic.reference_fit(CANONICAL, [label(CANONICAL, 3), "W3A"]) == 0.5

    def test_no_parsable_labels_is_nan(self):
        assert pd.isna(vic.reference_fit(CANONICAL, ["bogus"]))


class TestCompareProtein:
    ISOFORMS = {"P1": CANONICAL, "P1-2": ISOFORM}

    def test_agreeing_choice_has_no_disagreements(self):
        labels = [label(ISOFORM, 26), label(ISOFORM, 30)]
        res = vic.compare_protein("P1", self.ISOFORMS, labels, ISOFORM)
        assert res["counts"]["same_label"] == 2 and not res["disagreements"]
        assert res["chosen_equals_ensembl"] and res["pipeline_isoform"] == "P1-2"

    def test_protein_missing_from_uniprot_loses_recall(self):
        # COSMIC numbered by a protein UniProt doesn't have: the pipeline only has the
        # canonical sequence, so mutations after the insertion can't be placed
        labels = [label(ISOFORM, 26), label(ISOFORM, 30), label(ISOFORM, 35)]
        res = vic.compare_protein("P1", {"P1": CANONICAL}, labels, ISOFORM)
        assert not res["chosen_equals_ensembl"]
        assert res["counts"]["different_label"] + res["counts"]["pipeline_unmapped_only"] > 0
