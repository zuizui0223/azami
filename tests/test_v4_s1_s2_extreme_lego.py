"""Synthetic numerical and negative-control tests for V4 S1/S2 LEGO extension."""
import unittest
import numpy as np
import pandas as pd
from analysis.v4 import s1_s2_extreme_environment as source
from analysis.v4 import s1_s2_extreme_lego as lego
from analysis.v4 import s1_s2_lego_three_way as triple
from analysis.v4 import s1_s2_extreme_range_models as ranges


class ExtremeAndLego(unittest.TestCase):
    def test_source_registration_and_coordinates(self):
        self.assertEqual(len(source.EXTREMES),7)
        self.assertTrue(source.urls("chelsa_bio18")[0].endswith(
            "/CHELSA_bio18_1981-2010_V.2.1.tif"))
        self.assertIn("elev.zip",source.urls("topo_elevation")[0])
        with self.assertRaises(ValueError):
            source.urls("https://evil.invalid")
        native=pd.DataFrame({"obs_id":["x"],"latitude":[1.],"longitude":[2.]})
        focal=pd.DataFrame({"obs_id":["x"],"latitude":[8.],"longitude":[2.],
                             "range_state":["native"],
                             "taxon_name":["Cirsium vulgare"]})
        with self.assertRaises(ValueError):
            source.validate_coordinates(native,focal)  # wrong frozen cohorts

    def test_extreme_covariate_diagnostics(self):
        rng=np.random.default_rng(83)
        d=pd.DataFrame(rng.normal(size=(220,12)),columns=lego.ALL_ENV)
        d["accepted_key"]=np.repeat(np.arange(22),10)
        checks=lego.all_diagnostics(d,lego.ALL_ENV)
        self.assertEqual(set(checks),{"raw","among_equal_taxa",
                                      "within_taxon_demeaned"})
        self.assertEqual(checks["raw"]["status"],"OK")
        self.assertLess(lego.max_vif(checks,lego.ALL_ENV),5)
        selection=lego.outcome_blind_set(d,lego.SCENARIOS["warmest_quarter_BIO18"])
        self.assertEqual(selection["status"],"SELECTED")
        self.assertIn("chelsa_bio18",selection["selected"])
        reduced=lego.outcome_blind_set(
            d,lego.SCENARIOS["warmest_quarter_BIO18_hot_BIO05_elevation_reduced"])
        self.assertEqual(reduced["status"],"SELECTED")
        self.assertEqual(reduced["selected"],
                         ["chelsa_bio18","chelsa_bio05",
                          "chelsa_rsds_mean","topo_elevation"])
        self.assertNotIn("chelsa_bio12",selection["selected"])
        d["chelsa_bio18"]=d["chelsa_bio12"].copy()
        checks=lego.all_diagnostics(d,lego.ALL_ENV)
        self.assertEqual(checks["raw"]["status"],"NEAR_SINGULAR")

    def test_eight_v4_constructs_and_equal_n(self):
        rng=np.random.default_rng(12)
        n_taxa=24
        n_each=5
        group=np.repeat(np.arange(n_taxa),n_each)
        y=rng.normal(size=(n_taxa*n_each,len(lego.FEATURE_COLS)))
        y+=np.repeat(rng.normal(size=(n_taxa,1)),n_each,axis=0)*.3
        w,a,*_=lego.split_covariances(y,group)
        pairs=lego.rv_rows(w,a)
        self.assertEqual(len(lego.FEATURES),8)
        self.assertEqual(len(lego.PAIRS),28)
        self.assertEqual(len(pairs),28)
        self.assertEqual(sum(1 for p in pairs if p["left_module"]=="colour"
                             and p["right_module"]=="outline"),12)
        summary=lego.lego_module_summary(pairs)
        self.assertTrue(summary["within"]["n_within_module_pairs"]>0)
        x=lego.equal_n_stats(y,y,group,replicates=100)
        self.assertEqual(x["raw"],x["adjusted"])

    def test_environment_residual_removes_shared_slopes(self):
        rng=np.random.default_rng(9)
        n_taxa=24
        group=np.repeat(np.arange(n_taxa),6)
        x=pd.DataFrame({"accepted_key":group})
        x["chelsa_bio12"]=rng.normal(size=len(x))+group*.03
        x["chelsa_rsds_mean"]=rng.normal(size=len(x))
        vals=np.column_stack([4+2*x.chelsa_bio12 + rng.normal(size=len(x))*.01,
                              x.chelsa_rsds_mean])
        corrected,beta,_,design=lego.residualize(
            x,vals,["chelsa_bio12","chelsa_rsds_mean"])
        self.assertEqual(corrected.shape,vals.shape)
        self.assertLess(np.abs(corrected[:,0]).max(),.045)
        self.assertEqual(design.shape[1],5)

    def test_pairwise_maximum_entropy_holds_fixed_margins(self):
        observed=np.array([[[25.,5.],[4.,12.]],[[7.,6.],[3.,18.]]])
        expected=triple.pairwise_max_entropy(observed)
        self.assertAlmostEqual(float(expected.sum()),1.0,places=10)
        for axis in (0,1,2):
            np.testing.assert_allclose(expected.sum(axis=axis),
                                       observed.sum(axis=axis)/observed.sum(),
                                       atol=1e-9)
        parity,dev,fit=triple.three_way_excess(observed)
        self.assertTrue(np.isfinite(parity))
        self.assertGreaterEqual(dev,0.)

    def test_three_way_phenotypes_do_not_imply_genotypes(self):
        rng=np.random.default_rng(17)
        groups=np.repeat(np.arange(20),36)
        arr=rng.normal(size=(len(groups),len(lego.FEATURE_COLS)))
        out=triple.test_three_way(arr,groups,reps=100,seed=17)
        self.assertEqual(out['status'],'OBSERVED_VISIBLE_THREE_WAY_CONFIGURATION')
        self.assertEqual(out['taxa_eligible'],20)
        self.assertTrue(out['not_a_test_of_genetic_resource_co_retention'])
        self.assertEqual(len(out['observed_taxon_equal_state_fraction']),8)

    def test_real_s2_model_balancing_and_collinearity_stop(self):
        rng=np.random.default_rng(444)
        n=360
        dataframe=pd.DataFrame({
            'obs_id':[str(i) for i in range(n)],
            'user_id':[str(i) for i in range(n)],
            'observed_on':[f'{2023+(i//12)%3}-{1+i%12:02d}-15' for i in range(n)],
            'latitude':np.linspace(-65,65,n),
            'longitude':np.linspace(-175,175,n),
        })
        for name in lego.ALL_ENV:
            dataframe[name]=rng.normal(size=n)
        dataframe['orientation_image_vertical_angle']=60+10*dataframe['chelsa_bio12']+rng.normal(size=n)
        frame=ranges._get_dated_balanced(dataframe)
        self.assertGreaterEqual(frame._cell.nunique(),15)
        result=ranges.model_one(frame,'orientation_image_vertical_angle',
                      ranges.S2_MODELS['annual_BIO12'])
        self.assertEqual(result['status'],'EXPLORATORY_OBSERVED_S2_CLIMATE_GRADIENT')
        self.assertGreater(result['hydric_coefficient_per_sample_SD'],0)
        frame['chelsa_bio18']=frame['chelsa_bio12']
        stopped=ranges.model_one(frame,'orientation_image_vertical_angle',
                      ['chelsa_bio12','chelsa_bio18','chelsa_rsds_mean'])
        self.assertEqual(stopped['status'],'UNRESOLVED_S2_CLIMATE_COLLINEARITY')

    def test_s2_nonindependent_range_climate_only(self):
        rows=[]
        for species in sorted(lego.S2_NAMES):
            for state in ("native","introduced"):
                for i in range(4):
                    d={"taxon_name":species,"range_state":state,
                       "user_id":f"{state}_{i}"}
                    for name in lego.ALL_ENV:
                        d[name]=float(i+(1 if state=="introduced" else 0))
                    rows.append(d)
        frame=pd.DataFrame(rows)
        report=lego.aggregate_s2(frame)
        self.assertEqual(len(report["taxa"]),2)
        self.assertEqual(report["taxa"]["Cirsium arvense"][
                          "introduced_minus_native_climate_medians"][
                          "chelsa_bio18"],1)
        self.assertEqual(report["taxa"]["Cirsium vulgare"]["native"][
                          "available_pairwise_visible_phenotypes"],[])


if __name__=="__main__":
    unittest.main()
