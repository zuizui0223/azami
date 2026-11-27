import numpy as np
from analysis.v5.phenotype_water_interaction import (
    common_scale,difference_features,RAIN,SEASON,P,fit_water,in_hull,
    ENDPOINTS,SPECIES
)

def test_interaction_is_difference_of_products_not_product_of_differences():
    a=np.zeros((1,5));b=a.copy()
    a[:,RAIN]=2;a[:,SEASON]=3;b[:,RAIN]=4;b[:,SEASON]=5
    z=difference_features(a,b,np.zeros(5),np.ones(5))
    assert z[0,-1]==14.
    assert z[0,-1]!=4.

def test_common_physical_scale():
    rng=np.random.default_rng(3)
    a=rng.normal(size=(30,5));b=rng.normal(size=(40,5))+2
    c,s=common_scale(a,b)
    assert s[RAIN]==500 and s[SEASON]==20
    np.testing.assert_allclose(c,(a.mean(0)+b.mean(0))/2)

def test_species_specific_interaction_recovers_known_coefficients():
    rng=np.random.default_rng(42)
    prepared=[]
    truth=[np.array([.01,.005,.01,.02,.03,.01,-.05]),
           np.array([.01,.005,.01,.04,.01,.01,.06])]
    for b in truth:
        n=80
        z=np.column_stack((np.ones(n),rng.uniform(-1,1,(n,6))))
        # Add nonzero noise orthogonal to design, so the coefficient is
        # identifiable and the cluster SE is nondegenerate.
        r=rng.normal(0,.01,n);r-=z@np.linalg.lstsq(z,r,rcond=None)[0]
        a=np.full((n,2),.3);v=a.copy();v[:,0]+=z@b+r
        prepared.append(dict(design=z,blocks=np.repeat(np.arange(10),8),
                             tail_left=a,tail_right=v))
    rows=fit_water(prepared,ENDPOINTS[0])
    rr=[r for r in rows if r["term"]=="rain_x_seasonality"]
    assert len(rr)==3 and all(r["status"]=="estimated" for r in rr)
    np.testing.assert_allclose([r["coefficient"] for r in rr],[-.05,.06,.11],atol=1e-9)

def test_hull_support():
    a=np.array([[0.,0.],[0.,1.],[1.,0.],[1.,1.]])
    np.testing.assert_array_equal(in_hull(a,np.array([[.5,.5],[1.5,.5]])),[True,False])
