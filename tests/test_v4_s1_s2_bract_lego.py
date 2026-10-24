"""Source-gated candidate bract LEGO proxy tests: never botanical spine proof."""
import unittest
import numpy as np
import pandas as pd
from analysis.v4 import s1_s2_bract_lego as b
from analysis.v4.s1_s2_extreme_lego import ALL_ENV


class BractVisibleModules(unittest.TestCase):
    def test_missing_protective_envelope_kept_unknown(self):
        frame=pd.DataFrame({"accepted_key":[1],"orientation_image_vertical_angle":[1.]})
        result=b.evaluate(frame,ALL_ENV,{
            "annual_BIO12_elevation":{"status":"UNRESOLVED"},
            "warmest_quarter_BIO18_elevation":{"status":"UNRESOLVED"}})
        self.assertEqual(result["status"],"SOURCE_ENDPOINTS_NOT_MEASURED")

    def test_same_observation_bract_three_component_when_sufficient(self):
        rng=np.random.default_rng(323)
        n=24; each=7
        names=np.repeat(np.arange(n),each)
        data={"accepted_key":names}
        for v in b.COLUMNS:
            data[v]=rng.normal(size=n*each)+.5*np.repeat(
                rng.normal(size=n),each)
            data[v+"__eligible"]=np.ones(n*each,dtype=bool)
        for v in ALL_ENV:
            data[v]=rng.normal(size=n*each)
        frame=pd.DataFrame(data)
        models={
            "annual_BIO12_elevation":{"status":"UNRESOLVED_PROTECTED_CONFOUNDING"},
            "warmest_quarter_BIO18_elevation":{"status":"UNRESOLVED_PROTECTED_CONFOUNDING"}}
        out=b.evaluate(frame,ALL_ENV,models)
        self.assertEqual(out["status"],"SOURCE_LINKED_BRACHT_PROXY_EXPLORATORY")
        self.assertEqual(out["taxa_after_min5"],n)
        self.assertEqual(len(out["raw_bract_pair_relations"]),3)
        self.assertTrue(out["raw_visible_triad_beyond_pairwise"][
                        "triple_is_visible_phyllary_projection_NOT_genetic_G123"])
        self.assertEqual(out["adjusted_scenarios"]["annual_BIO12_elevation"][
                         "status"],"UNRESOLVED_PROTECTED_CONFOUNDING")

    def test_single_taxon_too_few_not_fabricated(self):
        rng=np.random.default_rng(23)
        d={"accepted_key":np.ones(15)}
        for v in b.COLUMNS:
            d[v]=rng.normal(size=15);d[v+"__eligible"]=[True]*15
        for v in ALL_ENV:
            d[v]=rng.normal(size=15)
        cohort,report=b.valid_source(pd.DataFrame(d),ALL_ENV)
        self.assertIsNone(cohort)
        self.assertIn("INSUFFICIENT",report["status"])


if __name__=="__main__":
    unittest.main()
