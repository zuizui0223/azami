"""Regional support and out-of-region incremental prediction for paired displacement."""
from itertools import product
import numpy as np
from scipy.spatial import Delaunay


def audit(species,ref,left,right,axes,x,dx,dy,prediction,sd_y):
    blocks=left.block.to_numpy()
    unique,counts=np.unique(blocks,return_counts=True)
    weights=1/counts[np.searchsorted(unique,blocks)]
    mu=ref[axes].to_numpy().mean(0); scale=ref[axes].to_numpy().std(0)
    hull=Delaunay(x)
    support=dict(species=species,native_cells=len(ref),native_blocks=ref.block.nunique(),pairs=len(left),
                 native_records=ref.photos.sum(),similar_records=left.photos.sum(),novel_records=right.photos.sum(),
                 blocks=len(unique),min_pairs=counts.min(),median_pairs=np.median(counts),max_pairs=counts.max(),
                 largest_share=counts.max()/len(left),effective_blocks=len(left)**2/(counts@counts),
                 analogous_outside_hull=sum(hull.find_simplex((left[axes]-mu)/scale)<0),
                 novel_outside_hull=sum(hull.find_simplex((right[axes]-mu)/scale)<0))
    dy=dy[:,:2]/sd_y
    pred=prediction[:,:2]/sd_y
    baseline=np.zeros_like(dy)
    calibrated=np.zeros_like(dy)
    for b in unique:
        te=blocks==b
        baseline[te]=dy[~te].mean(0)
        calibrated[te]=pred[te]+(dy[~te]-pred[~te]).mean(0)
    rows=[]
    for name,ids in [("angle",[0]),("chroma",[1]),("joint",[0,1])]:
        obs=dy[:,ids];p=pred[:,ids]
        null=(obs**2).sum(1)
        err=((obs-p)**2).sum(1)
        constant=((obs-baseline[:,ids])**2).sum(1)
        adjusted=((obs-calibrated[:,ids])**2).sum(1)
        def q(e,base,w):
            return 1-np.sum(e*w)/np.sum(base*w)
        loo=[q(err[blocks!=b],null[blocks!=b],1) for b in unique]
        # Compare site-specific deviations after removing all block-level shifts.
        cov=[]; vary=0
        for b in unique:
            ix=blocks==b
            a=obs[ix]-obs[ix].mean(0); c=p[ix]-p[ix].mean(0)
            if ix.sum()>1:
                cov.append((a*c).sum())
                vary+=ix.sum()
        signs=np.array(list(product([-1,1],repeat=len(cov))))
        vals=signs@np.array(cov)
        pv=np.mean(np.abs(vals)>=abs(sum(cov))-1e-10) if cov else np.nan
        rows.append(dict(species=species,endpoint=name,q_zero=q(err,null,1),equal_block_q=q(err,null,weights),
                         loo_min=min(loo),loo_max=max(loo),loo_positive=sum(v>0 for v in loo),
                         positive_blocks=sum(err[blocks==b].sum()<null[blocks==b].sum() for b in unique),
                         q_native_vs_constant=q(err,constant,1),q_incremental=q(adjusted,constant,1),
                         q_incremental_equal_block=q(adjusted,constant,weights),
                         within_cov=sum(cov),within_signflip_p=pv,blocks_in_within_test=len(cov),pairs_in_within_test=vary))
    return [support],rows
