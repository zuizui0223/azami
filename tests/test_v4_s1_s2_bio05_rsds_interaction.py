"""Synthetic-only tests for BIO05×RSDS interaction vs elevation audit."""
from __future__ import annotations

import unittest
import numpy as np
import pandas as pd

from analysis.v4.s1_s2_extreme_lego import regression_design
from analysis.v4.s1_s2_bio05_rsds_interaction_source import (
    CLIMATE,TERRAIN,MODELS,append_weighted_interaction,fit_one_cohort,
    summarize_nested,correct_interactions
)
from analysis.v4.s1_s2_climate_only_3_vs_5 import no_private_keys


def fake_s1():
    rng=np.random.default_rng(20261031)
    taxa=40
    each=9
    keys=np.repeat([f"dummy_taxon_{i}" for i in range(taxa)],each)
    means=rng.normal(size=(taxa,4))
    x=np.repeat(means,each,axis=0)+rng.normal(0,.2,size=(taxa*each,4))
    d=pd.DataFrame(x,columns=TERRAIN)
    d["accepted_key"]=keys
    d["corolla_lab_chroma"]=(
        0.7*d.chelsa_bio18+1.3*d.chelsa_bio05*d.chelsa_rsds_mean
        +rng.normal(0,.55,size=len(d))
    )
    d["corolla_hue_sin"]=rng.normal(size=len(d))
    d["corolla_hue_cos"]=rng.normal(size=len(d))
    return d


def fake_s2():
    rng=np.random.default_rng(602)
    cells=42
    each=9
    n=cells*each
    d=pd.DataFrame(rng.normal(size=(n,4)),columns=TERRAIN)
    d["_cell"]=np.repeat(np.arange(cells),each)
    d["_month"]=rng.integers(1,13,n)
    d["_year"]=rng.integers(2015,2026,n)
    d["corolla_lab_chroma"]=(
        .5*d.chelsa_bio18+.5*d.chelsa_bio05*d.chelsa_rsds_mean
        +rng.normal(0,.5,n)
    )
    d["corolla_hue_sin"]=rng.normal(size=n)
    d["corolla_hue_cos"]=rng.normal(size=n)
    return d


class InteractionModelTests(unittest.TestCase):
    def test_four_specs_include_main_effects_and_old_baseline(self):
        self.assertEqual(list(MODELS),[
            "M3_climate","M4_original",
            "M3_plus_interaction","M4_plus_interaction"
        ])
        self.assertEqual(CLIMATE,("chelsa_bio18","chelsa_bio05","chelsa_rsds_mean"))
        self.assertEqual(TERRAIN,(*CLIMATE,"topo_elevation"))
        self.assertEqual(MODELS["M4_original"],(TERRAIN,False))
        self.assertEqual(MODELS["M4_plus_interaction"],(TERRAIN,True))

    def test_S1_between_taxon_product_constant_within_taxon(self):
        d=fake_s1()
        X,group=regression_design(d,list(TERRAIN))
        w=d.accepted_key.map(1/d.groupby("accepted_key").size()).to_numpy(float)
        extended,diag=append_weighted_interaction(X,w,mode="between_taxa")
        self.assertEqual(extended.shape[1],X.shape[1]+1)
        self.assertAlmostEqual(np.average(extended[:,-1],weights=w),0,places=10)
        self.assertAlmostEqual(np.average(extended[:,-1]**2,weights=w),1,places=10)
        for label in d.accepted_key.unique():
            within=extended[d.accepted_key.eq(label).to_numpy(),-1]
            self.assertLess(float(np.ptp(within)),1e-10)
        self.assertTrue(diag["all_main_effects_retained"])

    def test_S1_actual_fit_reports_same_cohort_AIC_and_interaction(self):
        d=fake_s1()
        base=fit_one_cohort(d,["corolla_lab_chroma"],TERRAIN,False,
                            scope="S1_between_taxa")
        plus=fit_one_cohort(d,["corolla_lab_chroma"],TERRAIN,True,
                            scope="S1_between_taxa")
        self.assertEqual(base["n_photos"],plus["n_photos"])
        self.assertEqual(len(plus["interaction_test_rows"]),1)
        self.assertEqual(plus["interaction_test_rows"][0]["scale"],
                         "between_equal_taxa")
        self.assertTrue(np.isfinite(plus["interaction_test_rows"][0]["p_uncorrected"]))
        comp=summarize_nested(base,plus)["corolla_lab_chroma"]
        self.assertGreaterEqual(comp["weighted_residual_fraction_removed"],-1e-9)

    def test_S2_within_species_range_product_and_circular_hue(self):
        d=fake_s2()
        old=fit_one_cohort(d,["corolla_hue_sin","corolla_hue_cos"],TERRAIN,False,
                           scope="S2_within_range")
        new=fit_one_cohort(d,["corolla_hue_sin","corolla_hue_cos"],TERRAIN,True,
                           scope="S2_within_range")
        self.assertEqual(old["n_photos"],new["n_photos"])
        self.assertEqual(len(new["interaction_test_rows"]),1)
        self.assertEqual(new["interaction_test_rows"][0]["kind"],
                         "joint_circular_hue_2df")
        self.assertEqual(new["interaction_test_rows"][0]["scale"],
                         "within_species_range")
        self.assertTrue(new["model"]["all_main_effects_retained"])

    def test_holm_family_exploratory_without_selective_reporting(self):
        rows=[{"p_uncorrected":0.0001,"trait":"chroma"},
              {"p_uncorrected":0.05,"trait":"angle"},
              {"p_uncorrected":0.7,"trait":"lightness"}]
        out=correct_interactions(rows,"S1_between_taxa","M4_plus_interaction")
        self.assertEqual(out["n_tests"],3)
        self.assertEqual(out["n_holm"],1)
        self.assertTrue(out["not_original_october_Holm_family"])
        self.assertTrue(out["post_result_exploratory"])
        self.assertEqual(len(out["rows"]),3)

    def test_private_provenance_does_not_leak(self):
        no_private_keys({"S1":{"models":["M3","M4"],"n_photos":168907}})
        with self.assertRaisesRegex(ValueError,"PRIVATE_SOURCE_FIELD"):
            no_private_keys({"results":{"accepted_key":"do-not-export"}})


if __name__=="__main__":
    unittest.main()
