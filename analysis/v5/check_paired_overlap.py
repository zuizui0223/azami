from pathlib import Path
import sys
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

def verify(path):
    t=pd.read_csv(path/'regional.csv')
    p=pd.read_csv(path/'pairs.csv')
    w=pd.read_csv(path/'whole.csv')
    assert len(t)==len(w)==12 and ((w.climate_D>=0)&(w.climate_D<=1)).all()
    assert ((p[['climate_D','head_D']]>=0)&(p[['climate_D','head_D']]<=1)).all().all()
    for (sp,shape),d in p.groupby(['species','shape']):
        v=t[(t.species==sp)&(t['shape']==shape)&(t.h==.5)].iloc[0]
        assert np.isclose(spearmanr(d.climate_D,d.head_D).statistic,v.rho)
        assert len(d)==v.n_pairs==v.n_native_blocks*v.n_intro_blocks and not d.duplicated(['native_block','intro_block']).any()
    for _,g in p.groupby('species'):
        a=g[g['shape']];b=g[~g['shape']]
        assert np.allclose(a.climate_D,b.climate_D)
    print('Both species, all regional pairs and shared climate scores verified.')
if __name__=='__main__':verify(Path(sys.argv[1]))
