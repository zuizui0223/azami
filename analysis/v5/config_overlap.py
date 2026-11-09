"""Overlap attributable to the co-occurrence of marginal trait coordinates."""
from itertools import product
import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter
from statsmodels.stats.multitest import multipletests
from analysis.v5.novel_head_overlap import make_cells,match_cells,T

def run(sp,native,intro):
    n=make_cells(native,True);c=make_cells(intro)
    c=c[c.novel!=.5].reset_index(drop=True)
    pairs=match_cells(c)
    ns=n.block.to_numpy()
    a=c.iloc[np.concatenate([l for l,r in pairs])]
    b=c.iloc[np.concatenate([r for l,r in pairs])]
    assert np.array_equal(a.block.to_numpy(),b.block.to_numpy())
    rows=[]
    for cols in [T,T[:2]]:
        fullnative=n[cols].to_numpy()
        mean=fullnative.mean(0);sd=fullnative.std(0)
        bounds=np.quantile((np.vstack([fullnative,c[cols].to_numpy()])-mean)/sd,[.005,.995],axis=0)
        mats=[np.clip((g[cols].to_numpy()-mean)/sd,*bounds) for g in (n,a,b)]
        blocks=[g.block.to_numpy() for g in (n,a,b)]
        edges=[np.linspace(lo-1.5,hi+1.5,41) for lo,hi in bounds.T]
        def density(m):
            hist=np.histogramdd(m,bins=edges)[0]
            v=gaussian_filter(hist,[.5/np.diff(e)[0] for e in edges],mode='constant').ravel()
            return v/v.sum()
        def D(u,v):return np.minimum(u,v).sum()
        def score(matrices):
            ref,sim,new=map(density,matrices)
            return np.array([D(ref,sim),D(ref,new),D(ref,sim)-D(ref,new)])
        ob=score(mats)
        for restricted in [False,True]:
            rng=np.random.default_rng(20261109)
            indices=[[np.flatnonzero(g==block) for block in np.unique(g)] if restricted else [np.arange(len(g))] for g in blocks]
            null=[]
            for _ in range(499):
                versions=[]
                for m,ixs in zip(mats,indices):
                    u=m.copy()
                    for ix in ixs:
                        for j in range(1,m.shape[1]):
                            u[ix,j]=m[rng.permutation(ix),j]
                    versions.append(u)
                null.append(score(versions))
            null=np.array(null);avg=null.mean(0);lo,hi=np.quantile(null[:,2],[.025,.975])
            p=(1+np.count_nonzero(abs(null[:,2]-avg[2])>=abs(ob[2]-avg[2])-1e-12))/500
            rows.append(dict(species=sp,dimensions=len(cols),within_block=restricted,pairs=len(a),blocks=a.block.nunique(),
                observed_D_sim=ob[0],observed_D_new=ob[1],observed_drop=ob[2],recombined_D_sim=avg[0],recombined_D_new=avg[1],recombined_drop=avg[2],
                configuration_excess=ob[2]-avg[2],null_low=lo,null_high=hi,p=p))
    # The centroid vector uses the same paired cells as the overlap comparison.
    difference = b[T].to_numpy() - a[T].to_numpy()
    unique_blocks = np.unique(blocks[1])
    block_sums = np.stack([difference[blocks[1] == block].sum(0) for block in unique_blocks])
    flips = np.array(list(product([-1, 1], repeat=len(unique_blocks))))
    null = flips @ block_sums / len(a)
    means = difference.mean(0)
    p_values = (abs(null) >= abs(means) - 1e-12).mean(0)
    vectors = []
    for j, trait in enumerate(["angle", "chroma", "outline"]):
        vectors.append(dict(species=sp, trait=trait, pairs=len(a), blocks=len(unique_blocks),
                            native=n[T[j]].mean(), analogous=a[T[j]].mean(), novel=b[T[j]].mean(),
                            novel_minus_analogous=means[j], difference_native_sd=means[j]/n[T[j]].std(ddof=0),
                            p_block_flip=p_values[j]))
    return rows, vectors, [g[T].to_numpy() for g in (n,a,b)]
