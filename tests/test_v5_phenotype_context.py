import numpy as np
from analysis.v5.phenotype_context import coarsen, weights_for

def test_coarsening_mass():
    m=np.arange(16).reshape(4,4)/120
    np.testing.assert_allclose(coarsen(m,2),[[10/120,18/120],[42/120,50/120]])
    assert np.isclose(coarsen(m,2).sum(),1)
    w=weights_for(m,np.array([0,1,2]),np.array([0,1,3]),2)
    np.testing.assert_allclose(w,[5/120,5/120,50/120])
