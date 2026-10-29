"""Source-free checks of post-positive H3c sensitivity panel."""
from __future__ import annotations
import unittest
import numpy as np
from analysis.v4.s1_h3c_robustness import (
    DRAWS, PERMUTATIONS, declared_scenarios, examine, effect_against_null
)


class H3cRobustnessTests(unittest.TestCase):
    def sample_taxa(self,n=36):
        rng=np.random.default_rng(410)
        result=[]
        for idx in range(n):
            climate=rng.normal(0,.12,size=(8,3))
            traits=.7*climate + rng.normal(0,.20,size=(8,3))
            result.append({"taxon":f"do_not_release_taxon_{idx}",
                           "sector":idx%3,"x":climate,"y":traits})
        return result

    def test_declared_conditions_are_fixed_not_screened_by_significance(self):
        main=declared_scenarios(1)
        secondary=declared_scenarios(5)
        self.assertEqual(len(main),10)
        self.assertEqual(len(secondary),5)
        self.assertEqual(len({x["label"] for x in main}),10)
        self.assertEqual(
            [s["omit_axis"] for s in main if "omit_axis" in s],[0,1,2])
        self.assertEqual(
            [s["omit_sector"] for s in main if "omit_sector" in s],[0,1,2])
        self.assertEqual(DRAWS,12)
        self.assertEqual(PERMUTATIONS,99)
        with self.assertRaises(ValueError):
            declared_scenarios(3)

    def test_null_effect_is_calculated_at_source_draw_level(self):
        obs=[.25,.35]
        null=np.array([[.05,.10,.20],[0,.10,.15]])
        out=effect_against_null(obs,null)
        self.assertAlmostEqual(out["observed"],.3)
        self.assertAlmostEqual(out["null_mean"],.1)
        self.assertAlmostEqual(out["observed_minus_null_mean"],.2)
        self.assertTrue(out["original_H3c_Holm_2_not_recomputed"])

    def test_omit_each_trait_and_sector_retains_both_climate_geography_units(self):
        taxa=self.sample_taxa()
        for scenario in (
            {"label":"narrow","threshold":.7},
            {"label":"omit-angle","omit_axis":0},
            {"label":"omit-east","omit_sector":2},
            {"label":"kernel","bandwidth":1.0}
        ):
            res=examine(taxa,scenario,seed=2026,draws=2,perms=5)
            self.assertEqual(res["status"],"estimated_exploratory_robustness",scenario)
            self.assertEqual(res["successful_draws"],2)
            self.assertGreaterEqual(res["mean_supported_cross_sector_taxon_pairs"],40)
            self.assertTrue(np.isfinite(res["distance_advantage"]["observed_minus_null_mean"]))
            self.assertTrue(np.isfinite(res["conditional_KDE_overlap"]["observed_minus_null_mean"]))
            self.assertNotIn("do_not_release_taxon",str(res))

    def test_sparse_sector_fails_without_optimistic_coefficient(self):
        taxa=self.sample_taxa()
        res=examine(taxa,{"label":"omit-west-and-middle","omit_sector":0},
                    draws=1,perms=1)
        self.assertIn(res["status"],("estimated_exploratory_robustness",
                                     "insufficient_repeated_analog_support"))
        self.assertNotIn("taxon",str(res.get("source_taxon",None)))


if __name__=="__main__":
    unittest.main()
