"""Guard new S1/S2 exploratory family against post hoc outcome selection."""
import unittest
from analysis.v4 import s1_s2_multiplicity_audit as audit


class MultitestAudit(unittest.TestCase):
    def test_all_estimable_tests_adjusted(self):
        rows=[{"p_original":.001,"scenario":"A"},
              {"p_original":.02,"scenario":"B"},
              {"p_original":.9,"scenario":"C"}]
        answer=audit.correct_family(rows)
        self.assertEqual(answer["n_tests"],3)
        self.assertEqual(answer["n_BH_005"],2)
        self.assertEqual(answer["n_Holm_005"],2)
        self.assertGreaterEqual(answer["rows"][0]["p_BH_exploratory_family"],.001)
        self.assertGreaterEqual(answer["rows"][1]["p_BH_exploratory_family"],.02)

    def test_cannot_promote_invalid_probability(self):
        for invalid in (-1,1.1,float("nan")):
            with self.assertRaises(ValueError):
                audit.validated_probability(invalid)

    def test_s1_nested_scenario_exposure_family(self):
        m={"S1":{"scenario_results":{"annual":{
            "status":"ESTIMATED_POSTHOC_EXPLORATORY",
            "sample_n":120,"taxa_n":25,
            "trait_environment_component_gradients":{
                "chroma":{"components":[
                    {"trait_coordinate":"corolla_lab_chroma",
                     "components":[{"environment":"chelsa_bio12",
                                    "scale":"between","p_component":.02,
                                    "beta_std":.1}]}]}}}}}}
        rows=audit.collect_s1(m)
        self.assertEqual(len(rows),1)
        self.assertEqual(rows[0]["coordinate"],"corolla_lab_chroma")
        self.assertEqual(rows[0]["study"],"S1")

    def test_nonmatching_scientific_family_refused(self):
        with self.assertRaises(ValueError):
            audit.audit({"status":"SOURCE_VERIFIED_V4_EXTENSION_EXPLORATORY",
                         "not_V3_reanalysis":False})


if __name__=="__main__":
    unittest.main()
