"""All occupied climates: couple regional climate displacement to individual traits."""
import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter1d
from scipy.stats import rankdata
from statsmodels.stats.multitest import multipletests

from analysis.v4.native_cached_climate import VARIABLES


def rank_vector(a):
    v=rankdata(a).reshape(a.shape)
    return (v-v.mean())/np.linalg.norm(v-v.mean())


def compare(x,y,rng):
    a=rank_vector(x);b=rank_vector(y)
    rho=(a*b).sum()
    null=np.array([(a*b[rng.permutation(len(b))][:,rng.permutation(b.shape[1])]).sum() for _ in range(4999)])
    p=(1+(abs(null)>=abs(rho)-1e-12).sum())/5000
    leave=[]
    for dim,n in enumerate(a.shape):
        if n>1:
            for k in range(n):
                xx=np.delete(x,k,axis=dim);yy=np.delete(y,k,axis=dim)
                leave.append((rank_vector(xx)*rank_vector(yy)).sum())
    return rho,p,min(leave),max(leave)


def analyze(c,sp,traits,project,grid,density,overlap,seed):
    from scipy.spatial.distance import cdist
    x,_=project(c,list(VARIABLES));edges=grid(x)
    native=c.range_state.eq("native").to_numpy()
    raw=c[[*VARIABLES,*traits]].to_numpy()
    weight=1/c.groupby("range_state").range_state.transform("size").to_numpy()
    avg=np.average(raw,axis=0,weights=weight)
    z=(raw-avg)/np.sqrt(np.average((raw-avg)**2,axis=0,weights=weight))
    regions=[[g.index.to_numpy() for _,g in c[c.range_state==kind].groupby("block") if len(g)>=8] for kind in ("native","introduced")]
    n,m=map(len,regions)
    sums=np.zeros((n+1,m,4))
    delta=np.zeros((n+1,m,z.shape[1]))
    cut=[np.linspace(v.min()-2.5,v.max()+2.5,65) for v in z[:,-3:].T]
    def surfaces(ids):
        values=[density(x[ids],edges,.5)]
        for j,e in enumerate(cut):
            hist=np.histogram(z[ids,j-3],bins=e)[0].astype(float)
            s=gaussian_filter1d(hist,.5/np.diff(e)[0],mode="constant")
            values.append(s/s.sum())
        return values
    reference=surfaces(native);ref_mean=z[native].mean(0)
    rng=np.random.default_rng(seed)
    for _ in range(32):
        picks=[[rng.choice(ids,8,replace=False) for ids in side] for side in regions]
        ns=[surfaces(ids) for ids in picks[0]]+[reference]
        ins=[surfaces(ids) for ids in picks[1]]
        for j in range(4):
            sums[:,:,j]+=overlap([s[j] for s in ns],[s[j] for s in ins])/32
        left=np.stack([z[ids].mean(0) for ids in picks[0]]+[ref_mean])
        right=np.stack([z[ids].mean(0) for ids in picks[1]])
        delta+=(right[None,:,:]-left[:,None,:])/32
    records=[];pairs=[]
    for design,sl in [("regional",slice(0,n)),("all_native_reference",slice(n,n+1))]:
        climate=1-sums[sl,:,0]
        for j,trait in enumerate(["angle","chroma","outline"]):
            tests=[("overlap","climate_D",climate,1-sums[sl,:,j+1]),
                   ("novelty_direction","climate_D",climate,delta[sl,:,5+j])]
            tests.extend(("signed",axis,delta[sl,:,k],delta[sl,:,5+j]) for k,axis in enumerate(VARIABLES))
            for test,axis,a,b in tests:
                rho,p,lo,hi=compare(a,b,np.random.default_rng(seed+9))
                records.append(dict(species=sp,design=design,test=test,axis=axis,trait=trait,
                                    n_native_regions=n if design=="regional" else 1,n_intro_regions=m,
                                    rho=rho,p=p,leave_min=lo,leave_max=hi))
        for (i,j),ds in np.ndenumerate(climate):
            row=dict(species=sp,design=design,native_region=i+1,intro_region=j+1,climate_D=1-ds)
            for k,trait in enumerate(["angle","chroma","outline"]):
                row[trait+"_D"]=sums[sl,:,k+1][i,j]
                row[trait+"_shift"]=delta[sl,:,5+k][i,j]
            pairs.append(row)
    return records,pairs


def save(rows,pairs,out):
    df=pd.DataFrame(rows)
    df["holm"]=df.groupby(["design","test"]).p.transform(lambda v:multipletests(v,method="holm")[1])
    df.to_csv(out/"trait_coupling.csv",index=False)
    pd.DataFrame(pairs).to_csv(out/"trait_regional_pairs.csv",index=False)
