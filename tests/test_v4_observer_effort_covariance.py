"""Observer/taxon/region weighting gates for the v4 native original-source audit."""
import numpy as np
import pandas as pd
import pytest
from analysis.v4.observer_effort_covariance import (
    covariance_parts,observation_weights,equal_area_cells,SCHEMES
)


def frame():
    return pd.DataFrame({
      "obs_id":[str(i) for i in range(9)],
      "taxon":["A"]*6+["B"]*3,
      "cell":[1,1,1,1,2,2,1,1,1],
      "user_id":["ua","ua","ua","ub","ua","ub","uc","uc","ud"],
      "rad":[0.,1.,2.,3.,4.,5.,6.,7.,8.],
      "chroma":[8.,7.,6.,5.,4.,3.,2.,2.,1.] })


def test_nested_observer_weights_balance_species_region_and_users():
    f=frame()
    w=observation_weights(f,"equal_taxon_region_observer")
    assert np.isclose(w.sum(),1)
    assert np.isclose(w[:6].sum(),.5)
    assert np.isclose(w[6:].sum(),.5)
    assert np.allclose(w[:3],1/24)
    assert np.isclose(w[3],1/8)
    assert np.isclose(w[4],1/8) and np.isclose(w[5],1/8)
    assert np.allclose(w[6:8],1/8) and np.isclose(w[8],1/4)
    for scheme in SCHEMES:
        assert np.isclose(observation_weights(f,scheme).sum(),1)


def test_pooled_standardized_covariance_exactly_equals_between_plus_within():
    f=frame()
    for scheme in SCHEMES:
        out=covariance_parts(f,scheme,"chroma","rad")
        assert out["identity_residual"]<1e-12
        assert np.isclose(out["pooled_standardized_coefficient"],
                          out["between_contribution"]+out["within_contribution"])
        assert -1.0000001 <= out["pooled_standardized_coefficient"] <= 1.0000001
        assert out["n_taxa"]==2


def test_unknown_observer_does_not_create_a_new_independent_sampling_unit():
    f=frame();f.loc[0,"user_id"]=None
    with pytest.raises(ValueError,match="Missing observer"):
        covariance_parts(f,"equal_taxon_region_observer","chroma","rad")
    with pytest.raises(ValueError,match="Missing observer ID"):
        observation_weights(f,"equal_taxon_region_observer")


def test_coordinate_grid_near_limits_has_no_invalid_cells():
    cells=equal_area_cells([40.,-20.,89.999],[-100.,100.,179.999])
    assert len(cells)==3
    assert (cells>=0).all()
    assert (cells<90*180).all()
    with pytest.raises(ValueError,match="Missing or invalid"):
        equal_area_cells([np.nan],[10.])


def test_omitting_dominant_species_changes_observation_but_not_method_definition():
    f=frame()
    full=covariance_parts(f,"observation","chroma","rad")
    only=f.loc[f.taxon.eq("A")].copy()
    with pytest.raises(ValueError,match="Constant predictor/response"):
        covariance_parts(only.assign(rad=1.),"observation","chroma","rad")
    assert full["n_observations"]==9
