import numpy as np
from analysis.v5.phenotype_native_calibration import native_splits, distances
def test_disjoint_native_split():
    blocks=np.repeat(np.arange(12),50)
    splits=native_splits(blocks,nref=100,ncal=100,draws=4)
    for a,b in splits:
        assert len(a)==len(b)==100
        assert not set(blocks[a])&set(blocks[b])
def test_distance():
    a=np.column_stack([np.arange(10),np.arange(10)**2])
    v=np.vstack([a[0],a[1]+1000])
    d=distances(a,v,k=1)
    assert d[0]==0 and d[1]>1
