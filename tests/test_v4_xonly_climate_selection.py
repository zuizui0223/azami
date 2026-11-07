"""Tests enforce selection only on X, and deterministic X-only coverage."""
import numpy as np
import pytest
from analysis.v4.xonly_climate_selection import (
    select_climate_variables,greedy_rank,retained_covariance_fraction,
    weighted_covariance)

def test_redundant_climate_axes_removed():
    r=np.random.default_rng(11)
    n=400;a=r.normal(size=n);b=r.normal(size=n)
    x=np.column_stack([a,a+.003*r.normal(size=n),b,
                       .4*a+.8*b+.04*r.normal(size=n)])
    o=select_climate_variables(x,["a","copy_a","b","mixed"],
                               bootstraps=6,seed=3,target=.9)
    assert len(o["selected_variables"])<=3
    assert ("a" in o["selected_variables"]) != ("copy_a" in o["selected_variables"])
    assert o["selection_never_uses_phenotype_or_previous_results"]
    assert o["n_valid_X_only_bootstraps"]==6

def test_reconstruction_independent_of_any_other_table():
    rng=np.random.default_rng(9)
    a=rng.normal(size=(200,5))
    r1=select_climate_variables(a,list("abcde"),bootstraps=4,seed=99)
    r2=select_climate_variables(a.copy(),list("abcde"),bootstraps=4,seed=99)
    assert r1==r2

def test_missing_rows_are_explicitly_reported():
    x=np.random.default_rng(3).normal(size=(60,4))
    x[0,1]=np.nan
    a=select_climate_variables(x,list("abcd"),bootstraps=2)
    assert a["n_complete_X_rows"]==59 and a["n_incomplete_X_rows"]==1

def test_constant_variable_is_not_silently_selected():
    x=np.random.default_rng(3).normal(size=(40,3))
    x[:,0]=1.
    with pytest.raises(ValueError,match="CONSTANT"):
        select_climate_variables(x,list("abc"))

def test_covariance_fraction_in_unit_interval():
    c,_,_,_=weighted_covariance(np.random.default_rng(3).normal(size=(80,6)))
    assert retained_covariance_fraction(c,[])==0
    idx,curve=greedy_rank(c,list("abcdef"))
    assert 0.<=curve[-1]["climate_reconstruction_fraction"]<=1.
    assert curve[0]["n_selected"]==1
