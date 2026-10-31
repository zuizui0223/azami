"""Unit tests for elevation-free matched S1/S2 three-vs-five climate models.

Synthetic data only: original protected observations stay in GitHub Actions.
"""
import unittest

import numpy as np
import pandas as pd

from analysis.v4.s1_s2_climate_only_3_vs_5 import (
    MODELS, PRIMARY, HOT, FIVE, CLIMATE_REQUIRED, ANGLE,
    source_s1, apply_multiplicity, no_private_keys,
)


class ClimateOnlyThreeFiveTests(unittest.TestCase):
    def test_three_and_five_omit_elevation_and_are_nested(self):
        self.assertEqual(tuple(MODELS[PRIMARY]),
            ("chelsa_bio18","chelsa_bio01","chelsa_rsds_mean"))
        self.assertEqual(tuple(MODELS[HOT]),
            ("chelsa_bio18","chelsa_bio05","chelsa_rsds_mean"))
        self.assertEqual(tuple(MODELS[FIVE]),
            ("chelsa_bio18","chelsa_bio01","chelsa_rsds_mean",
             "chelsa_bio04","chelsa_bio15"))
        self.assertTrue(set(MODELS[PRIMARY])<set(MODELS[FIVE]))
        self.assertNotIn("topo_elevation",CLIMATE_REQUIRED)
        self.assertNotIn("topo_elevation",set().union(*MODELS.values()))

    def test_s1_climate_gate_does_not_use_elevation(self):
        rng=np.random.default_rng(31)
        n=160
        x=rng.normal(size=(n,len(CLIMATE_REQUIRED)))
        data=pd.DataFrame(x,columns=CLIMATE_REQUIRED)
        data["accepted_key"]=np.repeat([f"Taxon_{i}" for i in range(20)],8)
        data[ANGLE]=np.linspace(1,160,n)
        data[ANGLE+"__eligible"]=True
        data["topo_elevation"]=np.nan
        cohort,summary=source_s1(data,[ANGLE])
        self.assertEqual(len(cohort),n)
        self.assertEqual(summary["taxa"],20)
        self.assertEqual(summary["status"],"CLIMATE_ONLY_SOURCE_LINKED")

    def test_s1_paired_models_same_climate_complete_cohort(self):
        rng=np.random.default_rng(42)
        n=160
        data=pd.DataFrame(rng.normal(size=(n,len(CLIMATE_REQUIRED))),
                          columns=CLIMATE_REQUIRED)
        data["accepted_key"]=np.repeat([f"Taxon_{i}" for i in range(20)],8)
        data[ANGLE]=np.arange(n)*0.5+1
        data[ANGLE+"__eligible"]=True
        data["chelsa_bio04"]=np.nan
        first,_=source_s1(data,[ANGLE])
        self.assertIsNone(first)
        data["chelsa_bio04"]=rng.normal(size=n)
        first,_=source_s1(data,[ANGLE])
        self.assertEqual(len(first),n)

    def test_holm_refits_separate_biological_families(self):
        def row(p,trait):
            return {"p_uncorrected":p,"trait":trait,"environment":"chelsa_bio18",
                    "scale":"between_equal_taxa","n_photos":160}
        results=apply_multiplicity({
            PRIMARY:[row(.0001,"angle"),row(.06,"circular_hue_joint")],
            FIVE:[row(.01,"angle"),row(.04,"circular_hue_joint"),row(.9,"chroma")],
        },"S1")
        self.assertEqual(results[PRIMARY]["n_tests"],2)
        self.assertEqual(results[FIVE]["n_tests"],3)
        self.assertNotEqual(results[PRIMARY]["results"][0]["p_Holm"],.0001)
        self.assertTrue(results[FIVE]["not_original_20261024_Holm_family"])
        self.assertEqual(len(results[FIVE]["results"]),3)

    def test_private_export_fields_are_rejected(self):
        no_private_keys({"S1":{"sample_count":160,"models":list(MODELS)}})
        with self.assertRaisesRegex(ValueError,"PRIVATE_SOURCE_FIELD"):
            no_private_keys({"S1":{"nested":[{"accepted_key":"private"}]}})
        with self.assertRaisesRegex(ValueError,"PRIVATE_SOURCE_FIELD"):
            no_private_keys({"S2":{"longitude":140.5}})


if __name__=="__main__":
    unittest.main()
