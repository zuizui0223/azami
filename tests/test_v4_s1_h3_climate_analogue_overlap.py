"""Source-free ecological alignment and 3D KDE-overlap tests for H3c."""
from __future__ import annotations
import unittest

import numpy as np
import pandas as pd

from analysis.v4.s1_h3_climate_analogue_alignment import (
    analogue_analysis, kde_overlap_and_permutations, matched_pairs,
)


class ClimateAnalogueOverlapTests(unittest.TestCase):
    def synthetic_original_cells(self,n_taxa=36):
        rng=np.random.default_rng(20261029)
        cells=[]
        for j in range(n_taxa):
            sec=j%3
            for k in range(8):
                climate=rng.normal(0,.35,size=3) + np.array([sec*.03,0,0])
                phenotype=.8*climate+rng.normal(0,.15,size=3)
                cells.append({
                    "accepted_key":f"private_taxon_{j}",
                    "longitude":[-130,30,125][sec]+k*.002,
                    **{"_z"+str(q):float(v)
                       for q,v in enumerate(np.r_[climate,phenotype])},
                })
        return pd.DataFrame(cells)

    def test_KDE_normalised_L2_overlap_is_bounded_and_null_finite(self):
        rng=np.random.default_rng(76)
        y=rng.normal(size=(3,8,3))
        y[1]=y[0]
        d=np.linalg.norm(y[0,:,None,:]-y[1,None,:,:],axis=2)
        mask=np.zeros((1,8,8),bool)
        mask[0,:4,:4]=True
        pair={"i":np.array([0]),"j":np.array([1]),
              "d":d[None,...],"mask":mask,
              "analog_n":np.array([16.])}
        ans,null=kde_overlap_and_permutations(pair,y,rng,reps=17)
        self.assertAlmostEqual(ans["unconditional_mean"],1.0,places=10)
        self.assertGreaterEqual(ans["conditional_mean"],0)
        self.assertLessEqual(ans["conditional_mean"],1.0+1e-8)
        self.assertEqual(len(null),17)
        self.assertTrue(np.isfinite(null).all())

    def test_qualified_cross_region_climate_analogues_and_kde_overlap(self):
        result=analogue_analysis(self.synthetic_original_cells(),draws=3,
                                permutations=19,seed=49)
        self.assertEqual(result["status"],"source_matched_exploratory_analogue_alignment")
        self.assertGreaterEqual(result["mean_eligible_cross_sector_taxon_pairs"],40)
        self.assertEqual(result["successful_rarefaction_draws"],3)
        index=result["climate_analogue_conditional_KDE_L2_overlap"]
        self.assertEqual(index["status"],"estimated")
        self.assertGreaterEqual(index["mean_conditional_analogue_support"],0)
        self.assertLessEqual(index["mean_conditional_analogue_support"],1)
        self.assertLessEqual(index["mean_unconditional_eight_cells"],1)
        self.assertGreaterEqual(index["holm2_distance_vs_overlap_p"],
                                index["two_sided_taxon_block_permutation_p"])
        self.assertNotIn("private_taxon_",str(result))

    def test_no_supported_climate_pairs_does_not_generate_false_overlap(self):
        x=np.stack([np.zeros((8,3)),np.full((8,3),20)])
        y=np.stack([np.zeros((8,3)),np.ones((8,3))])
        pair,meta=matched_pairs(x,y,np.array([0,1]))
        self.assertIsNone(pair)
        self.assertEqual(meta["status"],"no_analogous_pair_support")


if __name__=="__main__":
    unittest.main()
