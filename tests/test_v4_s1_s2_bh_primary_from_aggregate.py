"""Source-free tests for V4 trait-level BH-FDR analysis."""
import unittest

from analysis.v4.s1_s2_bh_primary_from_aggregate import (
    bh, handle, EXPECTED_NEW, EXPECTED_OLD,
)


class BhPrimaryTests(unittest.TestCase):
    def test_bh_reference_vector(self):
        q=bh([.01,.03,.04,.2])
        for a,b in zip(q,[.04,.05333333333333334,.05333333333333334,.2]):
            self.assertAlmostEqual(a,b)
        self.assertEqual(bh([]),[])
        self.assertEqual(bh([0,0,1]),[0.,0.,1.])
        self.assertEqual(bh([.03]*5),[.03]*5)

    def test_invalid_probability_stops(self):
        for p in [float("nan"),float("inf"),-1.e-5,1.0001]:
            with self.assertRaises(ValueError):
                bh([.1,p])

    def test_bh_adjusted_q_monotone_in_p(self):
        p=[.1,.001,.049,.4,.02,.0001,.01]
        q=bh(p)
        sorted_q=[x for _,x in sorted(zip(p,q))]
        self.assertEqual(sorted_q,sorted(sorted_q))
        self.assertTrue(all(0<=a<=1 for a in q))
        self.assertTrue(all(a+1.e-15>=b for a,b in zip(q,p)))

    def test_original_global_V4_scope_and_hue_dedup(self):
        self.assertEqual(EXPECTED_OLD, {
            "S1/main_vif5":120,"S1/si_vif10":168,
            "S2/main_vif5":200,"S2/si_vif10":280})
        self.assertEqual(EXPECTED_NEW,{
            "S1/main_vif5":100,"S1/si_vif10":140,
            "S2/main_vif5":160,"S2/si_vif10":224})

    def test_nonverified_aggregate_must_stop(self):
        with self.assertRaises(ValueError):
            handle({"status":"SOME_OTHER_ANALYSIS","V4_NOT_V3":True})
        with self.assertRaises(ValueError):
            handle({"status":"SOURCE_VERIFIED_S1_S2_MAIN_SI_ALL_AVAILABLE_PHENOTYPES",
                    "V4_NOT_V3":False})


if __name__=="__main__":
    unittest.main()
