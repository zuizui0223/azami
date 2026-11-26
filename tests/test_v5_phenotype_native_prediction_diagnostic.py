import numpy as np
import pandas as pd
from analysis.v5.phenotype_native_prediction_diagnostic import predictions, aggregate

def test_fixed_models_constant_target():
    x=np.column_stack([np.arange(120.),np.arange(120.)**2])
    y=np.repeat([[.25,.75]],120,axis=0)
    r=predictions(x,y,x[[2,47,89]])
    assert len(r)==4
    for pred in r.values():
        np.testing.assert_allclose(pred,np.repeat([[.25,.75]],3,axis=0))

def test_block_equal_vs_weighted():
    s=aggregate(pd.DataFrame(dict(sse=[1.,12.],baseline_sse=[2.,10.],test_cells=[8,20],train_cells=[50,50])))
    np.testing.assert_allclose(s["q2_cell_weighted"],-1/12)
    np.testing.assert_allclose(s["q2_block_equal"],.15)
    assert s["positive_blocks"]==1
