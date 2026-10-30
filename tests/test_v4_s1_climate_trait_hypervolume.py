"""Source-free tests for the ecological hypervolume / climate breadth analysis."""
from __future__ import annotations
import unittest
import numpy as np
import pandas as pd

from analysis.v4.s1_climate_trait_hypervolume import (
    AXES, CLIMATE, PHENOTYPE, ellipsoid_log_volume, geo_aggregate,
    standardise_cells, draw_taxon_scores, breadth_regression,
    cross_region_convergence,
)


class HypervolumeEcologyTests(unittest.TestCase):
    def synthetic_photo_cells(self):
        rng=np.random.default_rng(20261029)
        records=[]
        for taxon in range(45):
            sector=taxon%3
            base_lon=(-135 if sector==0 else 20 if sector==1 else 115)
            wide=0.5+taxon/24
            for cell in range(10):
                for photo in range(5):
                    climate=rng.normal(size=3)*wide + np.array([taxon/70,sector*.13,0])
                    traits=climate*0.65+rng.normal(size=3)*.33
                    records.append({
                        "accepted_key":f"source_taxon_{taxon}",
                        "latitude":-10+cell*2+sector*18+photo*.002,
                        "longitude":base_lon+cell*2+photo*.003,
                        **{key:float(climate[j]) for j,key in enumerate(CLIMATE)},
                        **{key:float(traits[j]) for j,key in enumerate(PHENOTYPE)}
                    })
        return pd.DataFrame(records)

    def test_ellipsoid_volume_identifies_spread_without_convex_hull(self):
        rng=np.random.default_rng(32)
        x=rng.normal(size=(18,3))
        self.assertGreater(ellipsoid_log_volume(x*2),
                           ellipsoid_log_volume(x*.25))
        with self.assertRaises(ValueError):
            ellipsoid_log_volume(x[:, :2])
        with self.assertRaises(ValueError):
            ellipsoid_log_volume(x[:4])

    def test_geographic_support_and_equal_cell_rarefaction(self):
        frame=self.synthetic_photo_cells()
        cells,coverage=geo_aggregate(frame,1)
        self.assertEqual(coverage["source_taxa_with_cell_and_photo_support"],45)
        self.assertEqual(len(cells),450)
        cells,scales=standardise_cells(cells)
        self.assertEqual(len(scales["global_sd"]),len(AXES))
        scores=draw_taxon_scores(cells,draws=7)
        self.assertEqual(len(scores),45)
        self.assertTrue(np.isfinite(scores.climate_log_volume).all())
        self.assertTrue(np.isfinite(scores.trait_log_volume).all())
        self.assertTrue((scores.n_cells==10).all())
        self.assertEqual(set(scores.sector),{0,1,2})
        self.assertTrue((scores.geo_spread_km>0).all())

    def test_ecological_tests_use_taxa_as_permutation_units(self):
        cells,_=geo_aggregate(self.synthetic_photo_cells(),1)
        cells,_=standardise_cells(cells)
        scores=draw_taxon_scores(cells,draws=5)
        b=breadth_regression(scores,permutations=11)
        c=cross_region_convergence(scores,permutations=11)
        self.assertEqual(b["status"],"estimated")
        self.assertEqual(c["status"],"estimated")
        self.assertEqual(b["n_taxa"],45)
        self.assertEqual(b["n_permutations"],11)
        self.assertGreater(c["cross_sector_taxon_pairs"],100)
        self.assertTrue(c["phylogenetic_evolutionary_convergence_NOT_tested"])
        self.assertTrue(0<b["sector_conditioned_residual_permutation_p_exploratory"]<=1)
        self.assertTrue(0<c["within_sector_taxon_label_permutation_p_exploratory"]<=1)
        self.assertNotIn("source_taxon",str(b)+str(c))

    def test_sparse_coverage_fails_without_inventing_results(self):
        frame=self.synthetic_photo_cells()
        tiny=frame.loc[frame.accepted_key=="source_taxon_0"].copy()
        cells,coverage=geo_aggregate(tiny,5)
        self.assertLess(coverage["source_taxa_with_cell_and_photo_support"],25)
        invalid=smaller=pd.DataFrame({
            "sector":[0,1],
            "climate_log_volume":[0,0],
            "trait_log_volume":[0,0],
        })
        self.assertEqual(breadth_regression(smaller)["status"],
                         "insufficient_taxonomic_support")


if __name__ == "__main__":
    unittest.main()
