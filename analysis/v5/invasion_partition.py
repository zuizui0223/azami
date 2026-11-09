"""Availability-aware invasion contrasts, using the V5 fixed environmental grid."""
import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter
from analysis.v5.native_pca_env import _probability_content, KDE_SIGMA_BINS
from analysis.v5.novel_head_overlap import T
from analysis.v5.record_density import density


def background(env,transform,edges):
    hist=np.histogramdd(transform(env),bins=edges,weights=env._area_weight)[0]
    h=gaussian_filter(hist,KDE_SIGMA_BINS,mode="constant")
    return h/h.sum()


def head_cells(frame,transform,edges):
    d=frame[frame[T[0]].between(0,180)].dropna(subset=T).copy()
    pc=transform(d)
    d["p1"],d["p2"]=pc.T
    d["i"]=np.minimum(179,np.floor(d.latitude+90)).astype(int)
    d["j"]=np.floor((d.longitude+180)%360).astype(int)
    c=d.groupby(["i","j"])[[*T,"p1","p2"]].median().reset_index()
    bins=np.column_stack([np.searchsorted(e,c[k],side="right")-1 for e,k in zip(edges,["p1","p2"])])
    ok=((bins>=0)&(bins<len(edges[0])-1)).all(1)
    c=c[ok].copy()
    c["u"],c["v"]=bins[ok].T
    return c


def run(sp,env_n,env_i,native,intro,transform,edges,nat,inv,out):
    bn,bi=[background(g,transform,edges) for g in [env_n,env_i]]
    records=[];traits=[];balanced=[]
    n,i=[head_cells(d,transform,edges) for d in [native,intro]]
    for bg in [.95,.99,1.]:
        an,ai=[(_probability_content(v,bg) if bg<1 else v>=v.max()*.0001) for v in (bn,bi)]
        for occ in [.9,.95,.99]:
            on,oi=[_probability_content(v,occ) for v in [nat,inv]]
            for side,mass,own,other,occupied in [("introduced",inv,ai,an,on),("native",nat,an,ai,oi)]:
                common=own&other
                partition=[mass[common&occupied].sum(),mass[common&~occupied].sum(),
                           mass[own&~other].sum(),mass[~own].sum()]
                np.testing.assert_allclose(sum(partition),1)
                records.append(dict(species=sp,background=bg,occ_core=occ,side=side,
                    shared_core=partition[0],shared_outside_other_core=partition[1],
                    outside_other_background=partition[2],outside_own_background=partition[3],
                    outside_core_conditional_common=partition[1]/sum(partition[:2])))
            if not (bg==.99 and occ==.95):continue
            from analysis.v5.expansion_consistency import compare
            compare(sp,n,i,native,intro,nat,an&ai,out)
            u,v=i.u,i.v
            groups={"occupied_native_core":(ai&an&on)[u,v],
                    "expansion_in_common_climate":(ai&an&~on)[u,v],
                    "outside_native_background":(ai&~an)[u,v],
                    "outside_intro_background":(~ai)[u,v]}
            assert np.all(np.sum(list(groups.values()),axis=0)==1)
            for cols,name in [(T,"joint"),(T[:2],"angle_chroma"),([T[0]],"angle"),([T[1]],"chroma"),([T[2]],"outline")]:
                mu=n[cols].mean().to_numpy();sd=n[cols].std(ddof=0).to_numpy()
                a=(n[cols].to_numpy()-mu)/sd;b=(i[cols].to_numpy()-mu)/sd
                bounds=np.quantile(np.vstack([a,b]),[.005,.995],axis=0)
                ref=density(a,None,bounds,.5)
                for key,mask in groups.items():
                    row=dict(species=sp,trait=name,group=key,cells=sum(mask),native_cells=len(n))
                    if sum(mask)>=8:
                        row.update(D=np.minimum(ref,density(b[mask],None,bounds,.5)).sum(),
                            shift_raw=float((i.loc[mask,cols].mean().to_numpy()-mu)[0]) if len(cols)==1 else np.nan)
                    traits.append(row)
                # Shared-support standardisation: give both ranges identical PC-bin frequencies.
                for factor in [4,8]:
                    weights=[];hist=[]
                    for g in [n,i]:
                        ids=(g.u.to_numpy()//factor)*80+(g.v.to_numpy()//factor)
                        weights.append(ids);hist.append(np.bincount(ids,minlength=6400)/len(ids))
                    h=np.minimum(*hist);overlap=h.sum()
                    for k,ids in enumerate(weights):
                        weights[k]=h[ids]/hist[k][ids]
                    wn,wi=weights
                    q=dict(species=sp,trait=name,bin_width=factor,shared_mass=overlap,
                        native_support=int((wn>0).sum()),intro_support=int((wi>0).sum()),
                        native_ess=wn.sum()**2/(wn@wn),intro_ess=wi.sum()**2/(wi@wi),
                        D_raw=np.minimum(ref,density(b,None,bounds,.5)).sum(),
                        D_balanced=np.minimum(density(a,wn,bounds,.5),density(b,wi,bounds,.5)).sum())
                    if len(cols)==1:q.update(shift_raw=float(i[cols[0]].mean()-mu[0]),
                        shift_balanced=np.average(i[cols[0]],weights=wi)-np.average(n[cols[0]],weights=wn))
                    balanced.append(q)
    for file,data in [("availability_partition",records),("availability_heads",traits),("climate_balanced_heads",balanced)]:
        path=out/(file+".csv")
        pd.DataFrame(data).to_csv(path,index=False,mode="a",header=not path.exists())
