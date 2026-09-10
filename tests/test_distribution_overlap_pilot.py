import importlib.util
from pathlib import Path
import sys
import numpy as np
import pytest

root=Path(__file__).resolve().parents[1]/'analysis/v3'
sys.path.insert(0,str(root))
spec=importlib.util.spec_from_file_location('overlap_pilot',root/'run_distribution_overlap_pilot.py')
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)


def test_wasserstein_exact_decomposition():
    x=np.array([[-1.,0.],[0.,1.],[1.,-1.]])
    assert np.allclose(m.wasserstein_parts(x,x),[0,0,0])
    a=m.wasserstein_parts(x,x+[3,4])
    assert np.allclose(a,[5,5,0])
    y=x*np.array([2.,.5])
    total,centre,shape=m.wasserstein_parts(x,y)
    assert np.isclose(total**2,centre**2+shape**2)
    assert np.allclose(m.wasserstein_parts(x,y),m.wasserstein_parts(y,x))
    with pytest.raises(ValueError,match='equal nonempty'):
        m.wasserstein_parts(x,x[:2])


def test_overlap_identity_symmetry_and_separation():
    x=np.array([[-1.],[0.],[1.]])
    r=np.random.default_rng(8);draw=x[r.integers(3,size=256)]+r.normal(size=(256,1))*.5
    own=m.log_kde(draw,x,.5)
    assert m.overlap(x,x,.5,draw,draw,own,own)==1
    y=x+20;dy=draw+20;py=m.log_kde(dy,y,.5)
    a=m.overlap(x,y,.5,draw,dy,own,py)
    assert 0<=a<1e-10
    assert np.isclose(a,m.overlap(y,x,.5,dy,draw,py,own))
