"""Source-free tests for H3d conditional trait covariance and breadth."""
from __future__ import annotations
import unittest

import numpy as np
import pandas as pd

from analysis.v4.s1_h3_environment_moderated_trait_covariance import (
    within_taxon_covariance_moderation, pooled_additive_residuals,
    taxon_axis_breadth
)
from analysis.v4.s1_climate_trait_hypervolume import draw_taxon_scores


class H3ConditionalTraitTests(unittest.TestCase):
    def sample(self, n_taxa=45, n_cells=10):
        rng=np.random.default_rng(303)
        records=[]
        for taxon in range(n_taxa):
            sector=taxon%3
            for cell in range(n_cells):
                rain=(cell-4.5)/3.0 + rng.normal(scale=.012)
                temp=rng.normal()
                rsds=rng.normal()
                angle=rng.normal()
                chroma=angle*(.9 if rain>0 else -.9)+rng.normal(scale=.25)
                aspect=rng.normal()
                records.append({
                    "accepted_key":f"private_taxon_{taxon}",
                    "_z0":rain,"_z1":temp,"_z2":rsds,
                    "_z3":angle,"_z4":chroma,"_z5":aspect,
                    "longitude":-140 if sector==0 else 16 if sector==1 else 120,
                    "latitude":20+cell/10,
                    "n_photos":5,
                })
        return pd.DataFrame(records)

    def test_additive_climate_removed_before_cross_product_comparison(self):
        cells=self.sample()
        data,residual,beta=pooled_additive_residuals(cells)
        self.assertEqual(residual.shape,(450,3))
        self.assertEqual(beta.shape,(3,3))
        self.assertTrue(np.isfinite(residual).all())
        self.assertEqual(len(data),450)

    def test_moderation_detects_climate_dependent_cross_trait_covariance(self):
        cells=self.sample()
        outcome=within_taxon_covariance_moderation(
            cells,permutations=199,bootstraps=199)
        self.assertEqual(outcome["status"],"estimated")
        self.assertEqual(outcome["n_estimable_tests"],6)
        self.assertEqual(outcome["n_taxa_source"],45)
        target=next(x for x in outcome["results"]
                    if x["climate"]=="chelsa_bio12"
                    and x["trait_a"]=="orientation_image_vertical_angle"
                    and x["trait_b"]=="corolla_lab_chroma")
        self.assertGreater(target["higher_minus_lower_exposure_residual_covariance"],0)
        self.assertLessEqual(target["taxon_sign_flip_permutation_p"],.05)
        self.assertTrue(all("exploratory_Holm_6_moderations" in v
                            for v in outcome["results"]))
        self.assertNotIn("private_taxon",str(outcome))

    def test_trait_specific_diversity_has_no_new_taxon_identifiers(self):
        cells=self.sample()
        scores=draw_taxon_scores(cells,draws=3)
        result=taxon_axis_breadth(cells,scores,draws=3,permutations=39)
        self.assertEqual(result["status"],"estimated")
        self.assertEqual(result["n_original_predeclared_traits"],3)
        self.assertEqual(result["taxon_units"],45)
        self.assertEqual(len(result["results"]),3)
        self.assertTrue(all("exploratory_Holm_3_traits" in v
                            for v in result["results"]))
        self.assertNotIn("private_taxon",str(result))


if __name__=="__main__":
    unittest.main()
