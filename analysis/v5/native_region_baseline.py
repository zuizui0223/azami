"""Native-native spatial variation as a baseline for an introduced target region."""
import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter1d
from scipy.spatial.distance import cdist
from analysis.v4.native_cached_climate import VARIABLES


def run(c,sp,traits,project,grid,density,seed):
    x,_=project(c,list(VARIABLES));y,_=project(c,traits)
    w=1/c.groupby("range_state").range_state.transform("size")
    source_axes=[*traits,*VARIABLES]
    raw=c[source_axes].to_numpy()
    mean=np.average(raw,axis=0,weights=w);sd=np.sqrt(np.average((raw-mean)**2,axis=0,weights=w))
    z=(raw-mean)/sd
    gx,gy=grid(x),grid(y)
    single=[np.linspace(v.min()-2.5,v.max()+2.5,65) for v in z.T]
    rs=[[g.index.to_numpy() for _,g in c[c.range_state==side].groupby("block") if len(g)>=8]
        for side in ["native","introduced"]]
    n,m=map(len,rs)
    labels=["climate","joint_head_pc12","angle","chroma","outline",*VARIABLES]
    total=np.zeros((n+m,len(labels)))
    shifts=np.zeros((n+m,len(source_axes)))
    rng=np.random.default_rng(seed)
    for draw in range(32):
        values=[[] for _ in labels]
        region_means=[]
        for cells in rs[0]+rs[1]:
            idx=rng.choice(cells,8,replace=False)
            region_means.append(raw[idx].mean(0))
            ds=[density(x[idx],gx,.5),density(y[idx],gy,.5)]
            for j,e in enumerate(single):
                hist=np.histogram(z[idx,j],e)[0].astype(float)
                v=gaussian_filter1d(hist,.5/np.diff(e)[0],mode="constant");ds.append(v/v.sum())
            for dest,v in zip(values,ds):dest.append(v)
        means=np.stack(region_means)
        reference=(means[:n].sum(0)-means[:n])/(n-1)
        shifts[:n]+=(means[:n]-reference)/32
        shifts[n:]+=(means[n:]-reference.mean(0))/32
        for j,vs in enumerate(values):
            arr=np.array(vs)
            refs=(arr[:n].sum(0)-arr[:n])/(n-1)
            total[:n,j]+=(1-np.minimum(arr[:n],refs).sum(1))/32
            # Each introduced target uses the same set of n leave-one-region-out references.
            total[n:,j]+=cdist(arr[n:],refs,metric="cityblock").mean(1)/64
    rows=[];summary=[]
    for j,axis in enumerate(labels):
        nn,ni=total[:n,j],total[n:,j]
        lo,hi=np.quantile(nn,[.1,.9])
        summary.append(dict(species=sp,axis=axis,n_native=n,n_intro=m,native_median=np.median(nn),
                            intro_median=np.median(ni),median_difference=np.median(ni)-np.median(nn),
                            native_p10=lo,native_p90=hi,intro_above_native_p90=sum(ni>hi),
                            intro_above_native_max=sum(ni>nn.max())))
        for k,val in enumerate(total[:,j]):
            rows.append(dict(species=sp,axis=axis,range="native" if k<n else "introduced",
                             region_index=k+1 if k<n else k-n+1,dissimilarity=val,
                             signed_shift=shifts[k,j-2] if j>=2 else np.nan,
                             signed_shift_sd=shifts[k,j-2]/sd[j-2] if j>=2 else np.nan))
    return summary,rows
