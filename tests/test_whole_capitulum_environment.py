import numpy as np
from analysis.v3.run_whole_capitulum_environment import statistic,test_gradient as run_gradient_test

def test_construct_weight_not_component_count():
    assert np.isclose(statistic(np.array([1.,2.,2.]),[np.array([0]),np.array([1,2])]),2.5)

def test_omnibus_deterministic_and_signal():
    x=np.linspace(-1,1,30);y=np.column_stack([x,-x]);groups=[np.array([0]),np.array([1])]
    a=run_gradient_test(y,x,groups,[np.arange(30)],99,np.random.default_rng(1))
    b=run_gradient_test(y,x,groups,[np.arange(30)],99,np.random.default_rng(1))
    assert np.allclose(a[0],[1,-1]) and a[1]==b[1] and a[2]==b[2]==.01
