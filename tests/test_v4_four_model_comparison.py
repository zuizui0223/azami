"""Tests for matched four-model comparison and honest hue accounting."""
import unittest
from analysis.v4 import s1_s2_vif5_si10_alltraits as fit
from analysis.v4 import s1_s2_four_model_compare as compare


class FourClimateModelTests(unittest.TestCase):
    def test_four_explicit_models_and_climate_scopes(self):
        self.assertEqual(len(fit.SPEC),4)
        self.assertEqual(len(fit.FOUR_ANNUAL),4)
        self.assertEqual(len(fit.FOUR_WARM),4)
        self.assertEqual(
            set(fit.FOUR_ANNUAL)-set(fit.FOUR_WARM), {"chelsa_bio12"})
        self.assertEqual(
            set(fit.FOUR_WARM)-set(fit.FOUR_ANNUAL), {"chelsa_bio18"})
        self.assertEqual(set(fit.VIFSTEP_NAMES),{"main_vif5","si_vif10"})
        self.assertEqual(set(fit.SPEC),{
            "four_bio12","four_bio18","main_vif5","si_vif10"})
        self.assertEqual(set(fit.SI)-set(fit.MAIN),
                         {"chelsa_bio13","chelsa_rsds_mean"})
        self.assertEqual(compare.EXPECTED["S1"],{
            "four_bio12":80,"four_bio18":80,"main_vif5":100,"si_vif10":140})
        self.assertEqual(compare.EXPECTED["S2"],{
            "four_bio12":128,"four_bio18":128,"main_vif5":160,"si_vif10":224})

    def test_biological_hue_joint_only_one_test_per_environment(self):
        source={}
        for name in fit.SPEC:
            source[name]={
              "status":"ESTIMATED","n":400,
              "trait_cohort_vif_pass":True,
              "component_effects":{
                  "corolla_hue_sin":[{
                     "status":"ESTIMATED","environment":"chelsa_bio05",
                     "scale":"between_equal_taxa",
                     "beta_std":1.,"p_uncorrected":.00001,
                     "ci95_wald":[.5,1.5]}],
                  "corolla_hue_cos":[{
                     "status":"ESTIMATED","environment":"chelsa_bio05",
                     "scale":"between_equal_taxa",
                     "beta_std":.2,"p_uncorrected":.00003,
                     "ci95_wald":[.1,.3]}]},
              "joint_circular_hue_tests":[{
                "environment":"chelsa_bio05","scale":"between_equal_taxa",
                "status":"JOINT_CIRCULAR_HUE_WALD_2DF",
                "p_uncorrected":.002}]
            }
        out=fit.family_corrections({"hue":{"models":source}}, {})
        for model in fit.SPEC:
            family=out["S1/"+model]
            self.assertEqual(family["n_tests"],1)
            self.assertEqual(family["n_unadjusted"],1)
            self.assertEqual(family["n_BH"],1)
            self.assertEqual(family["n_Holm"],1)
            self.assertEqual(family["rows"][0]["coordinate"],
                             "hue_joint_sin_cos")
            self.assertIsNone(family["rows"][0]["beta_std"])
            self.assertIsNone(family["rows"][0]["ci95_wald"])

    def test_unadjusted_and_corrected_raw_p_are_both_preserved(self):
        source={}
        for model in fit.SPEC:
            source[model]={
              "status":"ESTIMATED","n":200,
              "trait_cohort_vif_pass":True,
              "component_effects":{
                 "bract_projection_p95":[
                   {"status":"ESTIMATED","environment":"chelsa_bio05",
                    "scale":"between_equal_taxa","beta_std":-.4,
                    "ci95_wald":[-.5,-.3],"p_uncorrected":.001},
                   {"status":"ESTIMATED","environment":"topo_elevation",
                    "scale":"within_taxon","beta_std":.03,
                    "ci95_wald":[-.01,.07],"p_uncorrected":.2}]},
              "joint_circular_hue_tests":[]}
        out=fit.family_corrections({"bract_projection":{"models":source}}, {})
        for model in fit.SPEC:
            rows=out["S1/"+model]["rows"]
            self.assertEqual(len(rows),2)
            self.assertAlmostEqual(rows[0]["p_original"],.001)
            self.assertAlmostEqual(rows[0]["p_BH"],.002)
            self.assertAlmostEqual(rows[0]["p_Holm"],.002)
            self.assertEqual(rows[0]["ci95_wald"],[-.5,-.3])
            self.assertFalse(rows[1]["raw_lt_0_05"])
            self.assertEqual(out["S1/"+model]["n_unadjusted"],1)

    def test_source_authentication_gate(self):
        with self.assertRaises(ValueError):
            compare.analyze({})
        with self.assertRaises(ValueError):
            compare.analyze({"status":
                "SOURCE_VERIFIED_S1_S2_MAIN_SI_ALL_AVAILABLE_PHENOTYPES"})
        self.assertEqual(compare.domain("bract_projection_p95"),
                         "bract_photo_proxy")
        self.assertEqual(compare.domain("hue_joint_sin_cos"),"colour")
        with self.assertRaises(ValueError):
            compare.domain("unknown")

if __name__=="__main__":
    unittest.main()
