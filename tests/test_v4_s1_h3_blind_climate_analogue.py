"""Source-independent H3 X-only climate PCA and matching checks."""
import numpy as np
import pandas as pd
from analysis.v4.s1_h3_blind_climate_analogue import climate_pcs,standardized,calibrate_threshold

def test_three_pc_projection_uses_no_y():
    rng=np.random.default_rng(5)
    a=rng.normal(size=(90,7))
    frame=pd.DataFrame(a,columns=["v"+str(i) for i in range(7)])
    weights=np.repeat(1/90,90)
    pcs,details=climate_pcs(frame,list(frame.columns),weights)
    assert pcs.shape==(90,3)
    assert details["PCA_outcome_blind"]
    assert .4<details["first_three_PC_variance_fraction"]<=1.
    assert np.allclose(pcs.mean(axis=0),0,atol=1e-8)

def test_calibrate_distance_from_no_flower_values():
    rng=np.random.default_rng(7)
    data=[]
    for t in range(9):
        for i in range(8):
            data.append(dict(accepted_key=str(t),longitude=-110 if t<3 else (10 if t<6 else 120)))
    cell=pd.DataFrame(data)
    x=rng.normal(size=(len(cell),3))
    x2=x*2
    threshold,meta=calibrate_threshold(cell,x,x2,20261107)
    assert threshold>0
    assert meta["outcome_blind_X_distance_quantile_calibration"]
