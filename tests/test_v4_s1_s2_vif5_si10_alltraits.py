"""No-secret regression tests for V4 VIF5/main VIF10/SI all-trait update."""
import unittest

import numpy as np
import pandas as pd
import statsmodels.api as sm

from analysis.v4 import s1_s2_vif5_si10_alltraits as e
from analysis.v4.s1_s2_extreme_lego import (
    ALL_ENV, FEATURE_COLS, FEATURES, diagnostic_matrix,
)


def fake_photos(n_taxa=25, per_taxon=50, seed=11):
    rng = np.random.default_rng(seed)
    taxa = np.repeat(np.arange(n_taxa), per_taxon)
    n = len(taxa)
    data = {"obs_id":np.arange(n).astype(str), "accepted_key":taxa,
            "latitude":rng.uniform(-50,55,n),
            "longitude":rng.uniform(-170,170,n),
            "user_id":np.arange(n).astype(str),
            "observed_on":["2024-06-15"]*n,
            "taxon_name":["Cirsium arvense"]*n,
            "range_state":["native"]*n}
    for col in ALL_ENV:
        data[col] = rng.normal(size=n)+0.1*np.repeat(
            rng.normal(size=n_taxa),per_taxon)
    for col in e.ENDPOINTS:
        data[col] = rng.normal(size=n)+np.repeat(rng.normal(size=n_taxa),per_taxon)*.2
        data[col+"__eligible"] = np.ones(n,dtype=bool)
    data["orientation_image_vertical_angle"] = 55+10*rng.normal(size=n)
    return pd.DataFrame(data)


class MainSiTraitAnalysisTest(unittest.TestCase):
    def test_fixed_vif_sets_and_provenance(self):
        self.assertEqual(len(e.MAIN),5)
        self.assertEqual(len(e.SI),7)
        self.assertTrue(set(e.MAIN) <= set(e.SI))
        self.assertEqual(set(e.SI)-set(e.MAIN),
                         {"chelsa_rsds_mean","chelsa_bio13"})
        self.assertEqual(len(e.ENDPOINTS),11)
        self.assertEqual(len(FEATURE_COLS),9)

    def test_vifstep_removes_redundant_exposure_without_outcomes(self):
        rng=np.random.default_rng(2026)
        source=pd.DataFrame(rng.normal(size=(3000,len(ALL_ENV))),columns=ALL_ENV)
        source["chelsa_bio17"] = source["chelsa_bio14"] + rng.normal(size=3000)*.04
        diagnostic=diagnostic_matrix(source[ALL_ENV],ALL_ENV)
        self.assertEqual(diagnostic["status"],"OK")
        answers=e.vifstep_all_contexts({"one":diagnostic},10.)
        self.assertGreaterEqual(len(answers["drops"]),1)
        self.assertLessEqual(answers["max_vif"],10.)
        self.assertEqual(answers["status"],"VIFSTEP_CONVERGED")

    def test_source_qc_no_imputed_bracts(self):
        original=fake_photos()
        cohort,reason=e.selected_source(original,
                ["orientation_image_vertical_angle","bract_projection_p95"])
        self.assertEqual(reason["status"],"SOURCE_LINKED")
        self.assertEqual(len(cohort),len(original))
        original.drop(columns=["bract_projection_p95__eligible"],inplace=True)
        cohort,reason=e.selected_source(original,["bract_projection_p95"])
        self.assertIsNone(cohort)
        self.assertEqual(reason["status"],"ORIGINAL_S1_TRAIT_OR_QC_UNAVAILABLE")

    def test_joint_circular_hue_cluster_score_covariance(self):
        rng=np.random.default_rng(77)
        G, n_each=40,10
        n=G*n_each
        groups=np.repeat(np.arange(G),n_each)
        x=np.column_stack([np.ones(n),rng.normal(size=n),rng.normal(size=n)])
        cov_by_taxon=rng.normal(size=G)
        y1=x[:,1]*.5+np.repeat(cov_by_taxon,n_each)+rng.normal(size=n)*.5
        y2=x[:,1]*.35+np.repeat(cov_by_taxon,n_each)*.8+rng.normal(size=n)*.5
        weights=np.ones(n)
        f1=sm.WLS(y1,x,weights=weights).fit(
            cov_type="cluster",cov_kwds={"groups":groups},use_t=False)
        f2=sm.WLS(y2,x,weights=weights).fit(
            cov_type="cluster",cov_kwds={"groups":groups},use_t=False)
        cross=e._cluster_cross_cov(f1,f2,x,weights,groups)
        self.assertEqual(cross.shape,(3,3))
        test=e.hue_joint([f1,f2],x,weights,groups,
             ["chelsa_bio01","chelsa_bio05"],[("between_equal_taxa",1)])
        self.assertEqual(len(test),2)
        self.assertEqual(test[0]["status"],"JOINT_CIRCULAR_HUE_WALD_2DF")
        self.assertLess(test[0]["p_uncorrected"],.05)

    def test_all_28_v4_edges_under_five_and_seven(self):
        source=fake_photos(n_taxa=25,per_taxon=50)
        result=e.s1_joint_28(source)
        self.assertEqual(result["n_pairs"],28)
        self.assertEqual(result["taxa"],25)
        for case in e.SPEC:
            self.assertEqual(result["models"][case]["status"],"ESTIMATED")
            self.assertEqual(len(result["models"][case]["pairs"]),28)
        self.assertEqual(result["n_visible_constructs"],8)

    def test_bract_links_to_all_original_constructs(self):
        source=fake_photos(n_taxa=25,per_taxon=7)
        result=e.s1_bract_edges(source)
        self.assertEqual(set(result["bridges_to_all_eight_original_constructs"]),
                         set(FEATURES))
        for key,bridge in result["bridges_to_all_eight_original_constructs"].items():
            self.assertEqual(bridge["status"],"OBSERVED_PHOTO_BRACT_BRIDGE")
            self.assertTrue(all(bridge["models"][spec]["status"]=="ESTIMATED"
                                for spec in e.SPEC))

    def test_family_corresponds_to_all_coordinates_not_only_angle(self):
        effects=[{"environment":"chelsa_bio01","scale":"within_taxon",
                  "status":"ESTIMATED","beta_std":.3,"p_uncorrected":.001,
                  "ci95_wald":[.20,.40]},
                 {"environment":"chelsa_bio05","scale":"between_equal_taxa",
                  "status":"ESTIMATED","beta_std":.1,"p_uncorrected":.3,
                  "ci95_wald":[-.09,.29]}]
        fitted={"status":"ESTIMATED","trait_cohort_vif_pass":True,
                "n":123,"component_effects":{"bract_projection_p95":effects}}
        s1={"bract_projection":{"status":"SOURCE_LINKED",
              "models":{"main_vif5":fitted,"si_vif10":fitted}}}
        s2={"Cirsium arvense/native":{
            "source_measured_visible_traits":{"chroma":{
                "status":"SOURCE_LINKED","balanced_n":150,
                "models":{"main_vif5":fitted,"si_vif10":fitted}}}}}
        corrected=e.family_corrections(s1,s2)
        self.assertEqual(corrected["S1/main_vif5"]["n_tests"],2)
        self.assertEqual(corrected["S2/si_vif10"]["n_tests"],2)
        self.assertLess(corrected["S1/main_vif5"]["rows"][0]["p_Holm"],.05)


if __name__=="__main__":
    unittest.main()
