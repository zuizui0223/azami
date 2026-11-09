"""Within-introduced-range contrasts with region-specific intercepts."""
from itertools import product
from pathlib import Path
import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests
from analysis.v5.novel_head_overlap import xyz


def estimate(x,y,blocks):
    # Frisch-Waugh centering removes every additive 20-degree regional average.
    groups=np.unique(blocks)
    xc=x.copy();yc=y.copy()
    for group in groups:
        take=blocks==group
        xc[take]-=x[take].mean(0);yc[take]-=y[take].mean()
    if np.linalg.matrix_rank(xc)!=x.shape[1]:
        raise ValueError("Unidentified within-region direction")
    return np.linalg.lstsq(xc,yc,rcond=None)[0],xc,yc


def run(sp,left,right,axes,dx,dy,sd_y,out):
    blocks=left.block.to_numpy()
    unique,counts=np.unique(blocks,return_counts=True)
    km=6371*2*np.arcsin(np.clip(np.linalg.norm(xyz(left)-xyz(right),axis=1)/2,0,1))
    extent=right.niche_distance.to_numpy()-left.niche_distance.to_numpy()
    head_dist=np.linalg.norm(dy[:,:2]/sd_y,axis=1)
    rows=[]
    for trait,y,labels,x in [
        ("angle",dy[:,0],axes,dx),("chroma",dy[:,1],axes,dx),
        ("joint_distance",head_dist,["niche_distance"],extent[:,None])]:
        # Effect units for x: one SD among introduced cell-pair contrasts.
        sx=x.std(0)
        base=x/sx
        for geo in [False,True]:
            design=np.column_stack([base,km/km.std()]) if geo else base
            beta,xc,yc=estimate(design,y,blocks)
            for j,key in enumerate(labels):
                rest=np.delete(xc,j,axis=1)
                focal=xc[:,j]
                if rest.shape[1]:
                    focal=focal-rest@np.linalg.lstsq(rest,focal,rcond=None)[0]
                scores=np.array([np.dot(focal[blocks==b],yc[blocks==b]) for b in unique[counts>1]])
                signs=np.array(list(product([-1,1],repeat=len(scores))))
                p=np.mean(abs(signs@scores)>=abs(sum(scores))-1e-10)
                loo=[]
                for b in unique:
                    idx=blocks!=b
                    if np.linalg.matrix_rank(xc[idx])==design.shape[1]:
                        loo.append(np.linalg.lstsq(xc[idx],yc[idx],rcond=None)[0][j])
                rows.append(dict(species=sp,trait=trait,axis=key,adjust_pair_distance=geo,pairs=len(y),
                    informative_blocks=len(scores),pairs_in_informative_blocks=counts[counts>1].sum(),
                    median_pair_km=np.median(km),slope=beta[j],loo_min=min(loo),loo_max=max(loo),
                    conditional_signflip_p=p))
    path=out/"introduced_direct_contrasts.csv"
    pd.DataFrame(rows).to_csv(path,index=False,mode="a",header=not path.exists())


def finalize(out):
    path=out/"introduced_direct_contrasts.csv"
    data=pd.read_csv(path)
    data["holm"]=data.groupby("adjust_pair_distance").conditional_signflip_p.transform(lambda p:multipletests(p,method="holm")[1])
    data.to_csv(path,index=False)
