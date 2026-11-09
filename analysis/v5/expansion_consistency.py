"""Connect climate occupancy expansion, signed trait shifts and cross-range slopes."""
import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.ndimage import gaussian_filter
from statsmodels.stats.multitest import multipletests
from analysis.v4.native_cached_climate import VARIABLES
from analysis.v5.novel_head_overlap import T


def mass_rank(p):
    values, inverse=np.unique(p,return_inverse=True)
    mass=np.bincount(inverse,weights=p.ravel(),minlength=len(values))
    rank=1-mass.cumsum()+mass/2
    return rank[inverse.ravel()].reshape(p.shape)


def attach_cells(frame,cells):
    f=frame[frame[T[0]].between(0,180)].dropna(subset=T).copy()
    f["i"]=np.minimum(179,np.floor(f.latitude+90)).astype(int)
    f["j"]=np.floor((f.longitude+180)%360).astype(int)
    climate=f.groupby(["i","j"])[list(VARIABLES)].median()
    return cells.merge(climate,on=["i","j"],validate="one_to_one")


def compare(sp,n,i,native,intro,nat,common,out):
    n,i=[attach_cells(f,c) for f,c in [(native,n),(intro,i)]]
    r=mass_rank(nat)
    np.testing.assert_allclose(np.sum(r*nat),.5,atol=1e-9)
    # Fixed native phenotype density used only to score introduced cells.
    y0=n[T].to_numpy();mu=y0.mean(0);sd=y0.std(0)
    y=(i[T].to_numpy()-mu)/sd
    lo,hi=np.quantile(np.vstack([(y0-mu)/sd,y]),[.005,.995],axis=0)
    edges=[np.linspace(a-1.5,b+1.5,41) for a,b in zip(lo,hi)]
    hist=np.histogramdd(np.clip((y0-mu)/sd,lo,hi),bins=edges)[0]
    p=gaussian_filter(hist,[.5/np.diff(e)[0] for e in edges],mode="constant")
    p/=p.sum();rank=mass_rank(p)
    idx=tuple(np.searchsorted(e,np.clip(y[:,j],lo[j],hi[j]),side="right")-1 for j,e in enumerate(edges))
    i["head_rarity"]=rank[idx]
    for c in [n,i]:
        c["climate_rarity"]=r[c.u,c.v]
        c["common"]=common[c.u,c.v]
        c["block"]=c.i.floordiv(20)*18+c.j.floordiv(20)
    slopes=[];differences=[];support=[]
    g=i[i.common].copy()
    rarity=g.climate_rarity.to_numpy()
    for col,name in zip([*T,"head_rarity"],["angle","chroma","outline","joint_rarity"]):
        # Per 0.1 increment of native climate density rank (no split at .95).
        fit=sm.OLS(g[col].to_numpy(),sm.add_constant(rarity*10)).fit(cov_type="cluster",cov_kwds={"groups":g.block},use_t=True)
        l,h=fit.conf_int()[1]
        slopes.append(dict(species=sp,trait=name,n_cells=len(g),blocks=g.block.nunique(),beta=fit.params[1],lo=l,hi=h,p=fit.pvalues[1]))
    for binsize in [4,8]:
        groups=[c[c.common].copy() for c in [n,i]]
        bins=[((c.u//binsize)*80+c.v//binsize).to_numpy() for c in groups]
        freq=[np.bincount(ids,minlength=6400)/len(ids) for ids in bins]
        target=np.minimum(*freq);shared=target.sum()
        for k,(c,ids) in enumerate(zip(groups,bins)):
            c["w"]=target[ids]/freq[k][ids]/(len(ids)*shared)
            c["introduced"]=k
        d=pd.concat(groups,ignore_index=True).query("w>0").copy()
        for k,c in d.groupby("introduced"):
            support.append(dict(species=sp,coarsen=binsize,introduced=k,n_cells=len(c),blocks=c.block.nunique(),ess=c.w.sum()**2/(c.w@c.w),shared_mass=shared))
        z=d.introduced.to_numpy()
        w=d.w.to_numpy()
        for model,axes in [("occupancy_rank",["climate_rarity"]),("annual",list(VARIABLES))]:
            x=d[axes].to_numpy();cen=np.average(x,axis=0,weights=w)
            sx=np.sqrt(np.average((x-cen)**2,axis=0,weights=w));x=(x-cen)/sx
            X=np.column_stack([np.ones(len(d)),z,x,x*z[:,None]])
            assert np.linalg.matrix_rank(X)==X.shape[1]
            for col,name,sy in zip(T,["angle","chroma","outline"],sd):
                fit=sm.WLS(d[col].to_numpy()/sy,X,weights=w).fit(cov_type="cluster",cov_kwds={"groups":d.block},use_t=True)
                for j,axis in enumerate(axes):
                    for term,ids in [("native",[2+j]),("introduced",[2+j,2+len(axes)+j]),("difference",[2+len(axes)+j])]:
                        v=np.zeros(X.shape[1]);v[ids]=1
                        t=fit.t_test(v);l,h=t.conf_int().ravel()
                        differences.append(dict(species=sp,coarsen=binsize,model=model,trait=name,axis=axis,term=term,beta=t.effect.item(),lo=l,hi=h,p=t.pvalue.item(),axis_sd=sx[j],native_trait_sd=sy))
    for file,rows in [("rarity_gradient",slopes),("range_gradient",differences),("gradient_support",support)]:
        path=out/(file+".csv");pd.DataFrame(rows).to_csv(path,index=False,mode="a",header=not path.exists())


def finalize(out):
    for file,keys in [("rarity_gradient",[]),("range_gradient",["coarsen","model","term"])]:
        path=out/(file+".csv");d=pd.read_csv(path)
        adjust=lambda x:multipletests(x,method="holm")[1]
        d["holm"]=d.groupby(keys).p.transform(adjust) if keys else adjust(d.p)
        d.to_csv(path,index=False)
