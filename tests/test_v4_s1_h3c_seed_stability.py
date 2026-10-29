"""Source-free test for post-positive H3c Monte Carlo seed stability."""
from __future__ import annotations
import unittest
import numpy as np
import pandas as pd

from analysis.v4.s1_h3c_seed_stability import SEEDS, one_panel


class H3cSeedStabilityTests(unittest.TestCase):
    def synthetic_cells(self,n=36):
        rng=np.random.default_rng(19)
        rows=[]
        for i in range(n):
            s=i%3
            for j in range(8):
                climate=rng.normal(scale=.12,size=3)
                appearance=.85*climate+rng.normal(scale=.16,size=3)
                row={
                    "accepted_key":f"synthetic_not_exported_{i}",
                    "longitude":(-130 if s==0 else 10 if s==1 else 120) + j*.002,
                }
                row.update({f"_z{k}":float(value)
                           for k,value in enumerate(np.r_[climate,appearance])})
                rows.append(row)
        return pd.DataFrame(rows)

    def test_exact_six_seed_series_are_fixed(self):
        self.assertEqual(len(SEEDS[1]),6)
        self.assertEqual(len(SEEDS[5]),6)
        self.assertEqual(len(set(SEEDS[1]+SEEDS[5])),12)

    def test_all_sequences_retained_and_stats_are_only_diagnostics(self):
        cells=self.synthetic_cells()
        result=one_panel(cells,1,draws=2,permutations=7,seeds=[101,202])
        self.assertEqual(result["status"],
                         "estimated_exploratory_random_sequence_stability")
        self.assertEqual(result["n_seeds_registered"],2)
        self.assertEqual(len(result["all_seed_results"]),2)
        self.assertEqual(result["n_cell_draws_per_seed"],2)
        self.assertTrue(result["random_sequences_are_monte_carlo_sensitivity_not_phylogenetic_replications"])
        for metric in ("summary_distance","summary_KDE_overlap"):
            s=result[metric]
            self.assertIn("min",s)
            self.assertIn("max",s)
            self.assertLessEqual(s["min"],s["median"])
            self.assertLessEqual(s["median"],s["max"])
            self.assertNotIn("synthetic_not_exported",str(s))
        self.assertNotIn("synthetic_not_exported",str(result))

    def test_unestimable_source_does_not_claim_effect(self):
        cells=self.synthetic_cells(n=12)
        result=one_panel(cells,1,draws=2,permutations=7,seeds=[101,202])
        self.assertEqual(result["status"],
                         "incomplete_not_all_seed_sequences_estimable")
        self.assertTrue(result["not_a_negative_or_equivalence_test"])


if __name__=="__main__":
    unittest.main()
