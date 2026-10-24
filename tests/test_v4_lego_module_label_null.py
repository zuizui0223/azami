"""Exact 280-label partition test has no genotype or evolutionary null."""
import unittest
from itertools import combinations
from analysis.v4 import lego_module_label_null as module


class ModuleLabelNullTests(unittest.TestCase):
    def test_all_fixed_partition_space(self):
        labels=list(module.all_partitions())
        self.assertEqual(len(labels),280)
        for p in labels:
            self.assertEqual([list(p.values()).count(v) for v in (0,1,2)],[1,3,4])

    def test_equal_relationships_give_no_fictitious_enrichment(self):
        rows=[{"left":left,"right":right,
               "within_rv":.1,"among_rv":.2}
              for left,right in combinations(module.CONSTRUCTS,2)]
        result=module.evaluate(rows)
        self.assertEqual(result["within_rv"]["exact_one_sided_module_label_reference_p"],1.)
        self.assertEqual(result["among_rv"]["observed_within_module_minus_cross_module_median_RV"],0.)

    def test_predeclared_modular_structure_can_be_detected(self):
        rows=[]
        for left,right in combinations(module.CONSTRUCTS,2):
            in_block=module.GROUPS[left]==module.GROUPS[right]
            rows.append({"left":left,"right":right,
                         "within_rv":.70 if in_block else .01,
                         "among_rv":.80 if in_block else .02})
        result=module.evaluate(rows)
        self.assertLess(result["within_rv"]["exact_one_sided_module_label_reference_p"],.02)
        self.assertGreater(result["among_rv"]["observed_within_module_minus_cross_module_median_RV"],.5)

    def test_missing_relation_fails(self):
        with self.assertRaises(ValueError):
            module.evaluate([])


if __name__=="__main__":
    unittest.main()
