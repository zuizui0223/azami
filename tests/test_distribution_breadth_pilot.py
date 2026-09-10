import importlib.util
from pathlib import Path
import numpy as np
import pytest
import pandas as pd

spec = importlib.util.spec_from_file_location('breadth', Path(__file__).resolve().parents[1]/'analysis/v3/run_distribution_breadth_pilot.py')
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


def test_density_and_geography_invariants():
    x = np.array([[-1., 0.], [0., 1.], [1., -1.]])
    points = np.array([[0.,0.], [2.,2.]])
    assert np.allclose(m.log_kde(points,x,.5),m.log_kde(points+7,x+7,.5))
    a = m.kde_volume(x,.5,m.rng_for(1),1024)
    b = m.kde_volume(2*x,1.,m.rng_for(1),1024)
    assert np.isclose(b/a,4.)
    assert m.geographic_rms(np.array([0.,0.]),np.array([179.,-179.])) < 225
    assert m.geographic_rms(np.array([0.,0.]),np.array([0.,0.])) == 0
    assert len(m.GROUPS['whole'])==18 and len(set(m.GROUPS['whole']))==18
    assert set(m.GROUPS)=={'orientation','colour','outline','involucre','whole'}
    with pytest.raises(ValueError,match='Non-finite'):
        m.weighted_scale(np.array([[1.],[np.nan]]),np.ones(2))


def test_spatial_holdout_prediction_and_support_gate():
    r=np.random.default_rng(91)
    f=pd.DataFrame({'latitude':r.uniform(-60,60,250),'longitude':r.uniform(-170,170,250)})
    env=r.normal(size=(250,4))
    y=(3*env[:,0]+r.normal(scale=.05,size=250))[:,None]
    result=m.predictive_gain(f,y,env,'test','synthetic',1)
    assert result['heldout_gain']>.95
    assert result['min_training_rank']==result['design_columns']
    assert m.predictive_gain(f.iloc[:50],y[:50],env[:50],'test','small',1) is None
