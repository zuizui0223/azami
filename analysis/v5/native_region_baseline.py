"""Native-native spatial variation as a baseline for an introduced target region."""
import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter1d
from scipy.spatial.distance import cdist
from analysis.v4.native_cached_climate import VARIABLES


def run(c,sp,traits,project,grid,density,seed):
    x,_=project(c,list(VARIABLES));y,_=project(c,traits)
    w=1/c.groupby("range_state").range_state.transform("size")
    raw=c[traits].to_numpy()
    mean=np.average(raw,axis=0,weights=w);sd=np.sqrt(np.average((raw-mean)**2,axis=0,weights=w))
    z=(raw-mean)/sd
    gx,gy=grid(x),grid(y)
    single=[np.linspace(v.min()-2.5,v.max()+2.5,65) for v in z.T]
    rs=[[g.index.to_numpy() for _,g in c[c.range_state==side].groupby("block") if len(g)>=8]
        for side in ["native","introduced"]]
    n,m=map(len,rs)
    total=np.zeros((n+m,5))
    rng=np.random.default_rng(seed)
    for draw in range(32):
        values=[[] for _ in range(5)]
        for cells in rs[0]+rs[1]:
            idx=rng.choice(cells,8,replace=False)
            ds=[density(x[idx],gx,.5),density(y[idx],gy,.5)]
            for j,e in enumerate(single):
                hist=np.histogram(z[idx,j],e)[0].astype(float)
                v=gaussian_filter1d(hist,.5/np.diff(e)[0],mode="constant");ds.append(v/v.sum())
            for dest,v in zip(values,ds):dest.append(v)
        for j,vs in enumerate(values):
            arr=np.array(vs)
            refs=(arr[:n].sum(0)-arr[:n])/(n-1)
            total[:n,j]+=(1-np.minimum(arr[:n],refs).sum(1))/32
            # Each introduced target uses the same set of n leave-one-region-out references.
            total[n:,j]+=cdist(arr[n:],refs,metric="cityblock").mean(1)/64
    rows=[];summary=[]
    for j,axis in enumerate(["climate","joint_head_pc12","angle","chroma","outline"]):
        nn,ni=total[:n,j],total[n:,j]
        lo,hi=np.quantile(nn,[.1,.9])
        summary.append(dict(species=sp,axis=axis,n_native=n,n_intro=m,native_median=np.median(nn),
                            intro_median=np.median(ni),median_difference=np.median(ni)-np.median(nn),
                            native_p10=lo,native_p90=hi,intro_above_native_p90=sum(ni>hi),
                            intro_above_native_max=sum(ni>nn.max())))
        for k,val in enumerate(total[:,j]):
            rows.append(dict(species=sp,axis=axis,range="native" if k<n else "introduced",
                             region_index=k+1 if k<n else k-n+1,dissimilarity=val))
    return summary,rows
