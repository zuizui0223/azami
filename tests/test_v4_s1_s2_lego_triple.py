"""Synthetic algebraic falsifiers for V4 phenotypic LEGO triad."""
import unittest
import numpy as np
from analysis.v4.s1_s2_lego_triple import (
    TRIAD,maxent_from_observed,frequencies_from_pairwise,analyze_triad)
from analysis.v4.s1_s2_extreme_lego import FEATURE_COLS


class PhenotypicLegoTripleTests(unittest.TestCase):
    def test_same_marginals_same_pairs_different_third_order(self):
        iid=np.ones(8)/8
        even_parity=np.array([.25,0,0,.25,0,.25,.25,0])
        independent=maxent_from_observed(iid)
        parity=maxent_from_observed(even_parity)
        self.assertEqual(independent["single_module_positive_probabilities"],
                         parity["single_module_positive_probabilities"])
        self.assertEqual(independent["pairwise_joint_positive_probabilities"],
                         parity["pairwise_joint_positive_probabilities"])
        self.assertAlmostEqual(independent["p111_maximum_entropy_pairwise"],.125,places=8)
        self.assertAlmostEqual(parity["p111_maximum_entropy_pairwise"],.125,places=8)
        self.assertEqual(parity["p111_observed"],0)
        self.assertAlmostEqual(parity["observed_minus_pairwise_maxent"],-.125,places=8)
        self.assertEqual(parity["n_nonempty_triple_configurations"],4)

    def test_exact_pairwise_constraint_algebra(self):
        singles=[.55,.45,.6]
        pairs=[.2,.35,.25]
        p=frequencies_from_pairwise(.1,singles,pairs)
        self.assertAlmostEqual(p.sum(),1.)
        self.assertTrue(np.min(p)>=0)
        reference=maxent_from_observed(p)
        self.assertLessEqual(reference["p111_pairwise_feasible_low"],.1)
        self.assertGreaterEqual(reference["p111_pairwise_feasible_high"],.1)
        self.assertAlmostEqual(sum(reference["all_8_maxent_reference_state_probabilities"]),1.)

    def test_observation_scale_and_taxon_scale_separate(self):
        rng=np.random.default_rng(43)
        n=24
        labels=np.repeat(np.arange(n),12)
        features=rng.normal(size=(12*n,len(FEATURE_COLS)))
        features += np.repeat(rng.normal(size=(n,1)),12,axis=0)
        answer=analyze_triad(features,labels,TRIAD,replicates=100,seed=91)
        self.assertEqual(answer["triad"],list(TRIAD))
        self.assertEqual(set(answer["within_and_among"]),{"within","among"})
        self.assertEqual(answer["within_and_among"]["among"]["n_source_taxa"],24)
        for scale in ("within","among"):
            self.assertEqual(len(answer["within_and_among"][scale][
                             "all_8_observed_state_probabilities"]),8)
            self.assertGreaterEqual(answer["within_and_among"][scale][
                                    "p111_pairwise_feasible_high"],0)
            self.assertTrue(answer["within_and_among"][scale][
                            "p111_not_genetic_joint_resource_frequency"])

    def test_invalid_distribution_fails(self):
        with self.assertRaises(ValueError):
            maxent_from_observed(np.ones(9)/9)
        with self.assertRaises(ValueError):
            maxent_from_observed(np.array([-1,2,0,0,0,0,0,0],float))


if __name__=="__main__":
    unittest.main()
