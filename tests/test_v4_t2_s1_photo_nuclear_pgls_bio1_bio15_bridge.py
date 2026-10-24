"""True multilocus phylogenetic GLS != nearest-tip clustering or floristic U/D."""
import numpy as np
import pytest
from analysis.v4 import t2_s1_photo_nuclear_pgls_bio1_bio15_bridge as m

def test_two_climate_axes_frozen_from_prior_and_five_axes_original_S1():
    assert m.TWO==("chelsa_bio01","chelsa_bio15")
    assert m.FIVE==("chelsa_bio01","chelsa_bio04","chelsa_bio12",
                    "chelsa_bio15","chelsa_rsds_mean")
    assert "BIO15_minus_BIO1" in m.DESIGNS["EAzami_vector_BIO15up_BIO1down"]
    assert len(m.RESPONSES)==2

def mock_rows(n=20):
    rows=[]
    for i in range(n):
        rows.append({
          "chelsa_bio01":str(10+3*np.sin(i*.7)),
          "chelsa_bio04":str(550+100*np.sin(i*.3)+.4*i),
          "chelsa_bio12":str(700+50*i+33*np.sin(i*.8)),
          "chelsa_bio15":str(30+9*np.cos(i*.37)+i*.1),
          "chelsa_rsds_mean":str(14+2*np.sin(i*.61)+i*.07),
        })
    return rows

def test_EAzami_direction_H1_is_higher_BIO15_lower_BIO1_continuous_translation():
    x,info=m.design_matrix(mock_rows(),("BIO15_minus_BIO1",))
    assert x.shape==(20,2)
    assert np.isclose(x[:,1].mean(),0,atol=1e-9)
    assert np.isclose(x[:,1].std(),1)
    assert info["EAzami_composite_is_posthoc_translation_not_CTMC_U_to_D"]
    assert m.covariance(np.eye(4),"IID").shape==(4,4)

def test_GLS_recovers_known_beta_with_nonzero_residual_and_boundary_lambda():
    x,_=m.design_matrix(mock_rows(),m.TWO)
    n=len(x)
    K=np.eye(n)*1.1+np.ones((n,n))*.12
    rng=np.random.default_rng(55)
    y=68+1.6*x[:,1]-2*x[:,2]+rng.normal(size=n)
    a=m.fit_gls(y,x,K,"IID")
    b=m.fit_gls(y,x,K,1.)
    lam,ll=m.profile_ml_lambda(y,x,K)
    assert a["n_taxa"]==20 and a["residual_degrees_freedom"]==17
    assert a["coefficients_degrees_per_SD"][1]>0
    assert a["coefficients_degrees_per_SD"][2]<0
    assert len(a["conditional_t_95CI_degrees_per_SD"])==3
    assert np.isfinite(ll) and 0<=lam<=1
    assert np.isfinite(b["maximized_gaussian_loglikelihood_fixed_covariance_shape"])

def test_partial_multiplicity_and_sufficient_species_gate():
    assert m.bh_adjust([.01,.2,.03,.9])==pytest.approx([.04,.2666666666667,.06,.9])
    badx=np.ones((7,3));bady=np.arange(7,dtype=float)
    with pytest.raises(ValueError,match="Insufficient matched taxa"):
        m.fit_gls(bady,badx,np.eye(7),"IID")
    with pytest.raises(ValueError,match="lambda"):
        m.covariance(np.eye(5),1.5)

def test_EAzami_U_D_not_implicitly_called_S1_image_angle():
    assert m.TIP_NOISE==.05
    assert "EAzami_vector_BIO15up_BIO1down" in m.DESIGNS
    assert m.MIN_TIPS==8
