"""Signed differences on the exact S1 analogue draws and a within-taxon analogue test."""
import numpy as np
import pandas as pd
from analysis.v5.trait_convergence import root
from analysis.v4.s1_climate_trait_hypervolume import load_source,geo_aggregate,standardise_cells
from analysis.v4 import s1_h3_climate_analogue_alignment as h
from statsmodels.stats.multitest import multipletests


def moments(dx,dy,weights):
    return np.einsum("nka,nkb,nk->nab",dx,dx,weights),np.einsum("nka,nkb,nk->nab",dx,dy,weights)


def main():
    decoded=root/"private/decoded"
    df,_=load_source(decoded/"observations.csv.gz",decoded/"measurement_receipt.json")
    out=[];sim=[]
    for deg in (1,5):
        cells,_=geo_aggregate(df,deg)
        z,_=standardise_cells(cells)
        taxa=h.eligible_taxa(z)
        rng=np.random.default_rng(20261109)
        perm_rng=np.random.default_rng(20261111)
        # Fixed taxon bootstrap weights across rarefactions: pairs sharing taxa remain linked.
        mult=np.random.default_rng(20261112).multinomial(len(taxa),np.full(len(taxa),1/len(taxa)),size=499)
        mult=np.vstack([np.ones(len(taxa)),mult])
        acc={};within=[]
        for rep in range(24):
            x,y,sector=h._sample_taxa(taxa,rng)
            p,support=h.matched_pairs(x,y,sector)
            if support["status"]!="qualified":continue
            i,j=p["i"],p["j"]
            dx=(x[i,:,None,:]-x[j,None,:,:]).reshape(len(i),64,3)
            dy=(y[i,:,None,:]-y[j,None,:,:]).reshape(len(i),64,3)
            w=mult[:,i]*mult[:,j]
            w=w/w.sum(1,keepdims=True)
            for label,mask in [("all_cells",np.ones(dx.shape[:2])),("analogues",p["mask"].reshape(len(i),64))]:
                weights=mask/mask.sum(1,keepdims=True)
                mx=(dx*weights[:,:,None]).sum(1,keepdims=True)
                my=(dy*weights[:,:,None]).sum(1,keepdims=True)
                for scale,parts in [("pair_mean",(mx,my,np.ones((len(i),1)))),
                                    ("within_pair",(dx-mx,dy-my,weights))]:
                    xx,xy=moments(*parts)
                    key=(label,scale)
                    if key not in acc: acc[key]=[]
                    acc[key].append((w@xx.reshape(-1,9),w@xy.reshape(-1,9)))
            u,v=np.triu_indices(8,1)
            near=np.linalg.norm(x[:,u]-x[:,v],axis=2)<=h.CLIMATE_ANALOG_Z
            ok=(near.sum(1)>=2)&(near.sum(1)<len(u))
            if not ok.any():continue
            yy=y[ok];mm=near[ok]
            def advantage(vectors):
                dist=np.linalg.norm(vectors[:,u]-vectors[:,v],axis=2)
                return np.mean(dist.mean(1)-(dist*mm).sum(1)/mm.sum(1))
            scores=[advantage(yy)]
            for _ in range(499):
                order=np.argsort(perm_rng.random(yy.shape[:2]),axis=1)
                scores.append(advantage(np.take_along_axis(yy,order[:,:,None],axis=1)))
            within.append([ok.sum(),*scores])
        for (sample,scale),reps in acc.items():
            xx,xy=np.mean(reps,axis=0).reshape(2,-1,3,3)
            coef=np.linalg.solve(xx,xy)
            lo,hi=np.quantile(coef[1:],[.025,.975],axis=0)
            for a,axis in enumerate(["annual_rain","mean_temperature","radiation"]):
                for t,trait in enumerate(["angle","chroma","outline"]):
                    out.append(dict(degrees=deg,sample=sample,scale=scale,axis=axis,trait=trait,
                                    beta=coef[0,a,t],lo=lo[a,t],hi=hi[a,t],draws=len(reps)))
        k,obs,*null=np.mean(within,axis=0)
        null=np.array(null);center=null.mean()
        sim.append(dict(degrees=deg,mean_supported_taxa=k,draws=len(within),observed=obs,
                        null_mean=center,excess=obs-center,p=(1+(abs(null-center)>=abs(obs-center)).sum())/500))
    pd.DataFrame(out).to_csv(root/"public/v5_signed_pairs.csv",index=False)
    tab=pd.DataFrame(sim);tab["holm"]=multipletests(tab.p,method="holm")[1]
    tab.to_csv(root/"public/v5_within_analogue.csv",index=False)

if __name__=="__main__":
    main()
