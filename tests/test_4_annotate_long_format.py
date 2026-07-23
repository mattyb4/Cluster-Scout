"""Unit tests for annotate_long_format in scripts/4_annotate.py -- the function
that fills ~15 annotation columns on the long-format PTM/mutation table
in-place, aggregating results from all 5 annotation phases.
"""
import pandas as pd
import pytest

from conftest import import_script

mod = import_script("4_annotate.py")


def _kinase_window_for(site_pos: int) -> tuple[dict, str]:
    """Build a pos_to_aa dict (all Ala except a Ser at site_pos) and the
    corresponding 15-mer window string build_kinase_window would produce,
    so kin_cache can be keyed to match exactly.
    """
    pos_to_aa = {p: "A" for p in range(site_pos - 7, site_pos + 8)}
    pos_to_aa[site_pos] = "S"
    window = mod.build_kinase_window(pos_to_aa, site_pos)
    return pos_to_aa, window


class TestAnnotateLongFormat:
    def test_annotates_ser_thr_ptm_with_all_data_present(self):
        pos_to_aa, window = _kinase_window_for(100)
        df = pd.DataFrame([{
            "uniprot_id": "P04637", "ptm_position": "S100", "ptm_type": "Phosphorylation",
            "gene": "TP53", "mutation": "R175H",
        }])

        mod.annotate_long_format(
            df,
            score_maps={"P04637": {100: 0.75}},
            confirmed_sites={("P04637", 100): "12345"},
            pp_cache={("TP53", "R175H"): ("D", "0.99")},
            seq_maps={"P04637": pos_to_aa},
            kin_cache={window: "CDK1(3.50,99.0%)"},
            disorder_maps={
                "general": {"P04637": {100: 0.8, 175: 0.3}},
                "binding": {"P04637": {100: 0.2, 175: 0.9}},
            },
            domain_maps={"P04637": [{"name": "DBD", "type": "domain", "start": 90, "end": 110}]},
        )
        row = df.iloc[0]

        assert row["polyphen_class"] == "probably_damaging", (
            f"pp_cache's 'D' code should map to the 'probably_damaging' class label, "
            f"got {row['polyphen_class']!r}"
        )
        assert row["polyphen_score"] == "0.99", f"the raw PolyPhen score should pass through, got {row['polyphen_score']!r}"
        assert row["1433_predicted"] == "yes", (
            f"a positive consensus score (0.75) at a Ser site should predict binding, "
            f"got {row['1433_predicted']!r}"
        )
        assert row["1433_confirmed"] == "yes", (
            f"the (uniprot, position) pair is in confirmed_sites -- should report yes, "
            f"got {row['1433_confirmed']!r}"
        )
        assert row["kinase_predictions"] == "CDK1(3.50,99.0%)", (
            f"a phosphorylation site with a resolvable sequence window should get the "
            f"cached kinase prediction, got {row['kinase_predictions']!r}"
        )
        assert row["ptm_aiupred_general"] == "0.800", (
            f"the PTM position's general disorder score should be formatted to 3 "
            f"decimals, got {row['ptm_aiupred_general']!r}"
        )
        assert row["ptm_is_disordered"] == "yes", f"0.8 > 0.5 -- should classify as disordered, got {row['ptm_is_disordered']!r}"
        assert row["ptm_is_binding"] == "no", f"0.2 is not > 0.5 -- should not classify as binding, got {row['ptm_is_binding']!r}"
        assert row["mut_aiupred_general"] == "0.300", (
            f"the MUTATION's position (175, from R175H) should look up its own disorder "
            f"score independently of the PTM's, got {row['mut_aiupred_general']!r}"
        )
        assert row["mut_is_binding"] == "yes", f"0.9 > 0.5 at the mutation's position -- should classify as binding, got {row['mut_is_binding']!r}"
        assert row["ptm_domain"] == "DBD (domain, 90-110)", (
            f"position 100 falls inside the DBD domain (90-110), got {row['ptm_domain']!r}"
        )
        assert row["mutation_domain"] == "", (
            f"position 175 (the mutation) falls outside DBD's range -- should be blank, "
            f"got {row['mutation_domain']!r}"
        )

    def test_non_ser_thr_ptm_leaves_1433_columns_blank(self):
        df = pd.DataFrame([{
            "uniprot_id": "P04637", "ptm_position": "Y50", "ptm_type": "Phosphorylation",
            "gene": "TP53", "mutation": "R175H",
        }])
        mod.annotate_long_format(
            df, score_maps={"P04637": {50: 0.9}}, confirmed_sites={}, pp_cache={},
            seq_maps={}, kin_cache={},
        )
        row = df.iloc[0]
        assert row["1433_predicted"] == "" and row["1433_predicted_consensus"] == "" and row["1433_confirmed"] == "", (
            f"14-3-3 predictions only apply to Ser/Thr sites -- a Tyr PTM site must "
            f"leave all three 1433 columns blank, got predicted={row['1433_predicted']!r} "
            f"consensus={row['1433_predicted_consensus']!r} confirmed={row['1433_confirmed']!r}"
        )

    def test_non_phosphorylation_type_leaves_kinase_blank(self):
        pos_to_aa, window = _kinase_window_for(100)
        df = pd.DataFrame([{
            "uniprot_id": "P04637", "ptm_position": "S100", "ptm_type": "Ubiquitination",
            "gene": "TP53", "mutation": "R175H",
        }])
        mod.annotate_long_format(
            df, score_maps={}, confirmed_sites={}, pp_cache={},
            seq_maps={"P04637": pos_to_aa}, kin_cache={window: "CDK1(3.50,99.0%)"},
        )
        assert df.iloc[0]["kinase_predictions"] == "", (
            f"kinase predictions only apply to phosphorylation sites -- an "
            f"Ubiquitination row must stay blank even though a matching cached "
            f"window exists, got {df.iloc[0]['kinase_predictions']!r}"
        )

    def test_isoform_tagged_mutation_is_stripped_before_lookup(self):
        df = pd.DataFrame([{
            "uniprot_id": "P04637", "ptm_position": "S100", "ptm_type": "Phosphorylation",
            "gene": "TP53", "mutation": "R175H(isoform?)",
        }])
        mod.annotate_long_format(
            df, score_maps={}, confirmed_sites={}, pp_cache={("TP53", "R175H"): ("B", "0.02")},
            seq_maps={}, kin_cache={},
            disorder_maps={"general": {"P04637": {175: 0.9}}, "binding": {"P04637": {}}},
        )
        row = df.iloc[0]
        assert row["polyphen_class"] == "benign", (
            f"pp_cache is keyed by the CLEAN mutation label (no isoform tag) -- the "
            f"'(isoform?)' suffix must be stripped before lookup, got {row['polyphen_class']!r}"
        )
        assert row["mut_aiupred_general"] == "0.900", (
            f"the mutation's position (175) must still be correctly extracted from "
            f"the isoform-tagged label for the disorder lookup, got {row['mut_aiupred_general']!r}"
        )

    def test_disorder_and_domain_columns_blank_when_maps_not_supplied(self):
        df = pd.DataFrame([{
            "uniprot_id": "P04637", "ptm_position": "S100", "ptm_type": "Phosphorylation",
            "gene": "TP53", "mutation": "R175H",
        }])
        # disorder_maps/domain_maps intentionally omitted (default None) -- simulates
        # a run where those phases were skipped or failed upstream.
        mod.annotate_long_format(
            df, score_maps={}, confirmed_sites={}, pp_cache={}, seq_maps={}, kin_cache={},
        )
        row = df.iloc[0]
        assert row["ptm_aiupred_general"] == "" and row["mut_aiupred_binding"] == "", (
            "with disorder_maps=None, every AIUPred column must be blank rather than raise"
        )
        assert row["ptm_is_disordered"] == "no" and row["mut_is_binding"] == "no", (
            f"a blank/empty score must classify as 'no' (falsy), not crash trying to "
            f"float('') -- got ptm_is_disordered={row['ptm_is_disordered']!r} "
            f"mut_is_binding={row['mut_is_binding']!r}"
        )
        assert row["ptm_domain"] == "" and row["mutation_domain"] == "", (
            "with domain_maps=None, both domain columns must be blank rather than raise"
        )

    def test_unparseable_ptm_position_does_not_crash_kinase_or_disorder_lookup(self):
        # A malformed ptm_position (doesn't match SITE_RE) must not raise -- every
        # downstream lookup keyed on ptm_pos should just degrade to blank.
        df = pd.DataFrame([{
            "uniprot_id": "P04637", "ptm_position": "not-a-site", "ptm_type": "Phosphorylation",
            "gene": "TP53", "mutation": "R175H",
        }])
        mod.annotate_long_format(
            df, score_maps={}, confirmed_sites={}, pp_cache={}, seq_maps={}, kin_cache={},
            disorder_maps={"general": {"P04637": {}}, "binding": {"P04637": {}}},
            domain_maps={"P04637": []},
        )
        row = df.iloc[0]
        assert row["kinase_predictions"] == "", (
            f"an unparseable ptm_position has no site to build a window from -- must "
            f"stay blank, not raise, got {row['kinase_predictions']!r}"
        )
        assert row["ptm_domain"] == "", (
            f"an unparseable ptm_position also has no position for the domain lookup "
            f"-- must be blank, got {row['ptm_domain']!r}"
        )


class TestAnnotateClusterLongFormat:
    def test_fills_annotation_columns_from_mutation_position_column(self):
        df = pd.DataFrame([{
            "UniProt": "P04637", "gene": "TP53",
            "anchor_mutation": "R248Q", "mutation": "R175H", "mutation_position": "175",
        }])
        mod.annotate_cluster_long_format(
            df,
            pp_cache={("TP53", "R175H"): ("D", "0.99")},
            disorder_maps={
                "general": {"P04637": {175: 0.8}},
                "binding": {"P04637": {175: 0.2}},
            },
            domain_maps={"P04637": [{"name": "DBD", "type": "domain", "start": 100, "end": 200}]},
        )
        row = df.iloc[0]

        assert row["polyphen_class"] == "probably_damaging", (
            f"pp_cache's 'D' code should map to probably_damaging, got {row['polyphen_class']!r}"
        )
        assert row["polyphen_score"] == "0.99", f"the raw score should pass through, got {row['polyphen_score']!r}"
        assert row["mut_aiupred_general"] == "0.800", (
            f"disorder should be looked up directly from mutation_position (175) -- no "
            f"regex parsing needed, unlike annotate_long_format -- got {row['mut_aiupred_general']!r}"
        )
        assert row["mut_is_disordered"] == "yes", f"0.8 > 0.5 -- should classify as disordered, got {row['mut_is_disordered']!r}"
        assert row["mut_is_binding"] == "no", f"0.2 is not > 0.5, got {row['mut_is_binding']!r}"
        assert row["mutation_domain"] == "DBD (domain, 100-200)", (
            f"position 175 falls inside DBD (100-200), got {row['mutation_domain']!r}"
        )

    def test_isoform_tagged_mutation_is_stripped_before_polyphen_lookup(self):
        df = pd.DataFrame([{
            "UniProt": "P04637", "gene": "TP53",
            "anchor_mutation": "R248Q", "mutation": "R175H(isoform?)", "mutation_position": "175",
        }])
        mod.annotate_cluster_long_format(df, pp_cache={("TP53", "R175H"): ("B", "0.02")})
        assert df.iloc[0]["polyphen_class"] == "benign", (
            f"pp_cache is keyed by the clean mutation label -- the '(isoform?)' suffix "
            f"must be stripped before lookup, got {df.iloc[0]['polyphen_class']!r}"
        )

    def test_blank_mutation_position_does_not_crash(self):
        df = pd.DataFrame([{
            "UniProt": "P04637", "gene": "TP53",
            "anchor_mutation": "R248Q", "mutation": "R175H", "mutation_position": "",
        }])
        mod.annotate_cluster_long_format(
            df, pp_cache={},
            disorder_maps={"general": {"P04637": {}}, "binding": {"P04637": {}}},
            domain_maps={"P04637": []},
        )
        row = df.iloc[0]
        assert row["mut_aiupred_general"] == "" and row["mutation_domain"] == "", (
            "a blank mutation_position has no position to look up -- both the disorder "
            f"and domain columns must be blank rather than raise, got "
            f"aiupred={row['mut_aiupred_general']!r} domain={row['mutation_domain']!r}"
        )

    def test_maps_not_supplied_leaves_columns_blank(self):
        df = pd.DataFrame([{
            "UniProt": "P04637", "gene": "TP53",
            "anchor_mutation": "R248Q", "mutation": "R175H", "mutation_position": "175",
        }])
        mod.annotate_cluster_long_format(df, pp_cache={})
        row = df.iloc[0]
        assert row["mut_aiupred_general"] == "" and row["mutation_domain"] == "", (
            "with disorder_maps/domain_maps omitted (default None), both must be "
            "blank rather than raise"
        )
        assert row["mut_is_disordered"] == "no" and row["mut_is_binding"] == "no", (
            f"a blank score must classify as 'no', not crash on float(''), got "
            f"disordered={row['mut_is_disordered']!r} binding={row['mut_is_binding']!r}"
        )
