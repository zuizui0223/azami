"""Source-free tests for three/elevation-free vs five-dimensional H3c climate matching."""
import unittest

import numpy as np
import pandas as pd

from analysis.v4.h3_climate_axis_comparison import (
    CLIMATE_SETS, CLIMATE_UNION, PHENOTYPE, analyse_paired,
    aggregate_cells, assert_public, build_common_pair, pack_pairs,
    source_blind_pairs, scaled_taxa,
)


class H3ClimateAxisComparisonTests(unittest.TestCase):
    def fixture(self,n_taxa=36,n_cells=8):
        rng=np.random.default_rng(19821106)
        result=[]
        for j in range(n_taxa):
            sect=j%3
            for k in range(n_cells):
                z=rng.normal(0,.33,size=3)+np.array([sect*.09,0,0])
                five=rng.normal(0,.22,size=2)
                y=.65*z+rng.normal(0,.14,size=3)
                result.append({
                    "accepted_key":f"private_taxon_{j}",
                    "longitude":[-130,35,123][sect]+k*.002,
                    "n_photos":10,
                    **dict(zip(CLIMATE_UNION,[z[0],z[1],z[2],z[0]*.8+five[0],z[1]*.85+five[1],five[0]+five[1],z[2]*.8+five[1]])),
                    **dict(zip(PHENOTYPE,y)),
                })
        return pd.DataFrame(result)

    def test_original_source_reader_includes_original_seasonality_columns(self):
        import inspect
        from analysis.v4 import h3_climate_axis_comparison as mod
        source = inspect.getsource(mod.load_common)
        self.assertIn('usecols=[*_required_cols(), "chelsa_bio04", "chelsa_bio15"]', source)

    def test_no_elevation_and_exact_variable_sets(self):
        self.assertEqual(CLIMATE_SETS["annual3_original"],("chelsa_bio12","chelsa_bio01","chelsa_rsds_mean"))
        self.assertEqual(CLIMATE_SETS["warm3_bio18_bio05"],("chelsa_bio18","chelsa_bio05","chelsa_rsds_mean"))
        self.assertEqual(CLIMATE_SETS["seasonal5_bio18_bio01"],("chelsa_bio18","chelsa_bio01","chelsa_rsds_mean","chelsa_bio04","chelsa_bio15"))
        self.assertFalse(any("elevation" in k.lower() for k in CLIMATE_UNION))

    def test_new_distances_can_use_five_axes_with_same_three_flower_axes(self):
        rng=np.random.default_rng(82)
        x=rng.normal(size=(4,8,5))
        y=rng.normal(size=(4,8,3))
        sec=np.array([0,1,2,1])
        rows,distance=source_blind_pairs(x,y,sec,0.9)
        self.assertEqual(distance.shape,(5*64,))
        self.assertEqual(len(rows),5)
        self.assertEqual(rows[0][2].shape,(8,8))
        packed=pack_pairs(rows,y,float(distance.max())+1)
        self.assertEqual(len(packed),5)
        pair=build_common_pair(packed,sorted(packed))
        self.assertEqual(pair["d"].shape,(5,8,8))
        self.assertEqual(pair["mask"].shape,(5,8,8))

    def test_synthetic_three_and_five_models_use_shared_photo_cells_and_pairs(self):
        result=analyse_paired(self.fixture(),draws=3,permutations=19,seed=37)
        self.assertEqual(result["status"],"PAIRED_MATCHED_CLIMATE_AXIS_COMPARISON_EXPLORATORY")
        self.assertEqual(result["successful_common_rarefaction_draws"],3)
        self.assertGreaterEqual(result["mean_common_eligible_taxon_pairs"],40)
        self.assertEqual(set(result["models"]),set(CLIMATE_SETS))
        for name,row in result["models"].items():
            self.assertEqual(row["n_climate_dimensions"],len(CLIMATE_SETS[name]))
            self.assertGreaterEqual(row["p_holm3_across_posthoc_climate_definitions_exploratory"],row["p_exploratory_unadjusted"])
            self.assertTrue(np.isfinite(row["observed_minus_null"]))
        self.assertNotIn("private_taxon_",str(result))

    def test_tiny_source_support_reports_no_analogue_claim(self):
        result=analyse_paired(self.fixture(n_taxa=3),draws=2,permutations=9,seed=3)
        self.assertEqual(result["status"],"INSUFFICIENT_COMMON_CLIMATE_ANALOGUE_SUPPORT")
        self.assertTrue(result["comparison_no_conclusions"])

    def test_aggregate_source_observations_keep_common_taxa(self):
        fixture=self.fixture(n_taxa=30,n_cells=8).copy()
        fixture["latitude"]=20
        # Construct ungrouped 240 photo rows in distinct 1-degree bins.
        fixture["longitude"]=np.tile(np.arange(8)*3.0,30)
        fixture["n_photos"]=1
        # Only data with 40 original photo records per taxon pass; this should fail.
        with self.assertRaisesRegex(ValueError,"INSUFFICIENT_COMMON_CLIMATE_TAXON_CELL_SUPPORT"):
            aggregate_cells(fixture,1)

    def test_identifiers_forbidden_in_aggregate_output(self):
        valid={"source":{"original_H3c_qualifying_photographs":157822},"results":{}}
        assert_public(valid)
        valid["results"]={"foo":[{"latitude":21.2}]}
        with self.assertRaisesRegex(ValueError,"PRIVATE_SOURCE_KEY"):
            assert_public(valid)


if __name__=="__main__":
    unittest.main()
