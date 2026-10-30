"""Source-free tests for the S2 original-range geographic block sensitivities."""
import unittest
import numpy as np
import pandas as pd
from analysis.v4.s2_original_colour_spatial_block_sensitivity import (
    EXPECTED, BOOTSTRAPS, spatial_block_labels, bootstrap_blocks,
)


class SpatialS2ColourTests(unittest.TestCase):
    def test_original_S2_three_holm_positives_and_native_six_panel_lock(self):
        self.assertEqual(len(EXPECTED),6)
        self.assertEqual([a[4] for a in EXPECTED],
                         [22936,26171,17335,27200,17335,27200])
        self.assertEqual(
            [round(row[5],6) for row in EXPECTED if row[1]=="introduced"],
            [-.173732,.134504,-.142801])
        self.assertEqual(BOOTSTRAPS,399)

    def test_no_latitude_or_longitude_regression_and_block_dimensions(self):
        frame=pd.DataFrame({
            "latitude":[-15,-5,0,19.9,20,39.9,40,79],
            "longitude":[-179,-161,-141,-39,-20,0,20,130],
        })
        clusters=spatial_block_labels(frame)
        self.assertEqual(len(clusters),len(frame))
        self.assertNotEqual(clusters[3],clusters[4])
        with self.assertRaises(ValueError):
            spatial_block_labels(frame,degrees=5)
        with self.assertRaises(ValueError):
            spatial_block_labels(pd.DataFrame({"latitude":[90.1],"longitude":[0]}))

    def test_block_bootstrap_spatial_reweighting_maintains_sign(self):
        rng=np.random.default_rng(29)
        nblocks=24
        nperblock=18
        block=np.repeat(np.arange(nblocks),nperblock)
        n=len(block)
        env=rng.normal(size=n)
        eta=rng.normal(size=n)
        fixed=rng.normal(scale=.10,size=nblocks)
        noise=rng.normal(scale=.12,size=n)
        X=np.column_stack([np.ones(n),env,eta])
        y=.65*env+.11*eta+fixed[block]+noise
        weight=np.ones(n)
        out=bootstrap_blocks(X,y,weight,1,block,reps=99,seed=3)
        self.assertEqual(out["status"],"estimated")
        self.assertEqual(out["n_original_source_geo_blocks"],nblocks)
        self.assertEqual(out["n_largest_source_weight_blocks_heldout"],5)
        self.assertEqual(out["n_largest_holdouts_retaining_sign"],5)
        self.assertGreater(out["source_block_bootstrap_95_percentile_interval"][0],0)
        self.assertTrue(out["block_partition_used_for_resampling_not_a_regression_predictor"])
        self.assertNotIn("latitude",out)
        self.assertNotIn("longitude",out)

    def test_too_few_spatial_blocks_does_not_issue_fake_ci(self):
        rng=np.random.default_rng(2)
        x=np.column_stack([np.ones(60),rng.normal(size=60)])
        y=rng.normal(size=60)
        r=bootstrap_blocks(x,y,np.ones(60),1,np.repeat(np.arange(6),10),reps=50)
        self.assertEqual(r["status"],"insufficient_spatial_blocks")
        self.assertNotIn("source_block_bootstrap_95_percentile_interval",r)


if __name__=="__main__":
    unittest.main()
