"""Allocate change in reference overlap among traits (descriptive, not causal)."""
import itertools
from math import factorial
from pathlib import Path
import sys
import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests

def run(root):
    d=pd.read_csv(root/'novel_head_overlap.csv')
    d=d[(d.design=='geographic_pairs')&(d.native95_mass==.95)]
    rows=[]
    for (sp,h),g in d.groupby(['species','bandwidth']):
        val={frozenset(s.split('+')):v for s,v in zip(g.traits,g.drop_in_overlap)}
        val[frozenset()]=0.
        axes=['angle','chroma','outline']
        phi={}
        for a in axes:
            others=[b for b in axes if b!=a]
            phi[a]=0.
            for size in range(3):
                for u in itertools.combinations(others,size):
                    phi[a]+=factorial(size)*factorial(2-size)/6*(val[frozenset([*u,a])]-val[frozenset(u)])
        assert np.isclose(sum(phi.values()),val[frozenset(axes)])
        for a in axes:
            v=g[g.traits==a].iloc[0]
            rows.append(dict(species=sp,axis=a,bandwidth=h,full_delta=val[frozenset(axes)],component=phi[a],
                single_D_sim=v.D_native_analogous,single_D_new=v.D_native_novel,single_delta=v.drop_in_overlap,p=v.p_swap))
    table=pd.DataFrame(rows)
    primary=table.bandwidth==.5
    table.loc[primary,'holm_singleton6']=multipletests(table.loc[primary,'p'],method='holm')[1]
    table.to_csv(root/'trait_overlap_components.csv',index=False)
    print(table[primary].round(4).to_string(index=False))

if __name__=='__main__':run(Path(sys.argv[1]))
