"""Synthetic-only tests: no simulated outputs appear as measured ecological results."""
from __future__ import annotations
import unittest
import numpy as np
import pandas as pd
import analysis.v4.s2_four_axis_joint_niche as mod

def synthetic(n=1000):
    rng=np.random.default_rng(37)
    frames=[];j=0
    for sp in mod.SPECIES:
        for state in ("native","introduced"):
            x=rng.normal(size=(n,4))
            x[:,0]+=0.2 if state=="introduced" else 0
            y=np.arange(n)%70
            lat=-40+(y%15)*4
            lon=-150+(y//15)*25
            frame=pd.DataFrame({
                "obs_id":["synthetic_"+str(q) for q in range(j,j+n)],
                "taxon_name":sp,"range_state":state,
                "user_id":["tester_"+str(q%65) for q in range(n)],
                "latitude":lat,"longitude":lon,
                "_cell":mod.cell_id(lat,lon),
                "observed_on":["2020-06-01"]*n,
                "chelsa_bio18":220+35*x[:,0],
                "chelsa_bio05":25+2*x[:,1],
                "chelsa_rsds_mean":14+2*x[:,2],
                "topo_elevation":500+75*x[:,3],
                "_angle":68+3*x[:,0]+rng.normal(0,5,n)})
            frames.append(frame);j+=n
    return pd.concat(frames,ignore_index=True)

class JointV4SourceContract(unittest.TestCase):
    def test_exact_spatial_projection(self):
        self.assertEqual(mod.cell_id([0],[0]).tolist(),[1476])
        with self.assertRaises(ValueError):
            mod.cell_id([91],[0])

    def test_phenotype_support_synthetic(self):
        sample=synthetic()
        result=mod.phenotype_common_support(sample)
        self.assertEqual(set(result),set(mod.SPECIES))
        for row in result.values():
            self.assertIn(row["status"],{
                "OBSERVATIONAL_FOUR_AXIS_ANALOG_SUPPORT",
                "INSUFFICIENT_PHOTOGRAPHIC_COMMON_SUPPORT"})
            if row["status"]=="OBSERVATIONAL_FOUR_AXIS_ANALOG_SUPPORT":
                self.assertTrue(row["not_an_equivalence_test"])

    def test_niche_observations_no_ids(self):
        sample=synthetic()
        original=dict(mod.BIO18_MEDIANS)
        try:
            # Synthetic medians are unrelated to the original evidence. Monkeypatch
            # only this explicit source identity assertion for the test.
            for sp in mod.SPECIES:
                for state in ("native","introduced"):
                    vals=sample.loc[(sample.taxon_name==sp)&(sample.range_state==state),"chelsa_bio18"]
                    mod.BIO18_MEDIANS[(sp,state)]=float(vals.median())
            result=mod.empirical_niche(sample,boot=30)
        finally:
            mod.BIO18_MEDIANS.clear();mod.BIO18_MEDIANS.update(original)
        for row in result.values():
            self.assertTrue(row["not_environmental_availability_corrected_niche_expansion"])
            self.assertGreaterEqual(row["joint_4D_regularized_mahalanobis_centroid_distance"],0)
        import json
        payload=json.dumps(result)
        self.assertNotIn("synthetic_",payload)
        self.assertNotIn("tester_",payload)

if __name__=="__main__":
    unittest.main()
