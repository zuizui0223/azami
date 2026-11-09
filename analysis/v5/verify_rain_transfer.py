from pathlib import Path
import sys
import numpy as np
import pandas as pd

def verify(root):
    t=pd.read_csv(root/'transfer.csv');n=pd.read_csv(root/'niche.csv')
    a=pd.read_csv(root/'rain_displacement.csv');coefs=pd.read_csv(root/'native_coefficients.csv')
    assert len(t)==6 and len(coefs)==6 and len(a)==12 and set(t.group)=={'native_spatial_cv','introduced_analogous','introduced_novel'}
    assert (t.n_cells>=t.n_blocks).all() and (t.n_photos>=t.n_cells).all()
    assert np.allclose(t.observed-t.predicted,t.bias)
    assert (t.rmse>=abs(t.bias)-1e-9).all()
    prev=pd.read_csv(root.parent/'S2_full_range/niche.csv')
    assert np.allclose(n.D,prev.D) and np.allclose(n.outside_density,prev.outside_density)
    for sp in t.species.unique():
        for group in ['introduced_novel','introduced_analogous']:
            sub=a[(a.species==sp)&(a.group==group)]
            measured=t[(t.species==sp)&(t.group==group)].iloc[0]
            assert np.isfinite(sub[['angle_estimate','lower','upper']]).all().all()
            assert sub[sub.comparison.isin(['wetter','drier'])].n.sum()==measured.n_cells
            assert sub.loc[sub.comparison=='wetter','rain_difference'].item()>0
            assert sub.loc[sub.comparison=='drier','rain_difference'].item()<0
    print('Verified support, both species, transfer and displacement tables against original climate classification.')
if __name__=='__main__':verify(Path(sys.argv[1]))
