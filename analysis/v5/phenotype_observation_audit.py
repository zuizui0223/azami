"""Audit the previously selected climate-rarity/chroma associations.

Freeze original V5 climate bins and cell cohort. Reproduce the published
coefficients before any sensitivity. Aggregate outputs only. Observation
coverage and image lightness are diagnostics, not colour calibration.
"""
from pathlib import Path
from itertools import product
import json
import numpy as np
import pandas as pd
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests
from analysis.v5.phenotype_redistribution import ROOT, SEED, eligible
from analysis.v5.phenotype_endpoint_scan import load_measures
from analysis.v5.novel_head_overlap import T
from analysis.v5.invasion_partition import head_cells, background
from analysis.v5.expansion_consistency import mass_rank
from analysis.v5.native_pca_env import background_mesh, _background_species_audit, _pca_env_2d, _probability_content
from analysis.v5.native_pca_env_fullcoverage import apply_original_five_CHELSA
from analysis.v4.native_introduced_opportunity import authority_source, angular_cells
from analysis.v4.native_cached_climate import VARIABLES

OUT=ROOT/"phenotype_observation_audit"
CHROMA=T[1]
HEADS="corolla_lab_chroma_n_usable_heads"
LIGHT="scan_corolla_lab_lightness"
PROXIES=["sin_month","cos_month","year","research_fraction","log_records","log_usable_heads"]
EXPECTED={
 "Cirsium arvense": (850, .0238707832093976, -.0753041100747241),
 "Cirsium vulgare": (1023, .7612325886661637, .1750536045403736),
}

def enrich_cells(f,c):
    q=eligible(f)
    dt=pd.to_datetime(q.observed_on,format="mixed",errors="coerce")
    q["sin_month"]=np.sin(2*np.pi*dt.dt.month/12)
    q["cos_month"]=np.cos(2*np.pi*dt.dt.month/12)
    q["year"]=dt.dt.year-2020
    q["research_fraction"]=q.quality_grade.eq("research").astype(float)
    q["log_usable_heads"]=np.log1p(q[HEADS].where(q[HEADS]>0))
    g=q.groupby(["i","j"])
    m=g[["sin_month","cos_month","year","research_fraction","log_usable_heads",LIGHT]].mean()
    m["records"]=g.size()
    m["observers"]=g.user_id.nunique()
    m["log_records"]=np.log(m.records)
    c=c.merge(m,on=["i","j"],validate="one_to_one")
    assert (c.records>=1).all()
    c["block"]=c.i.floordiv(20)*18+c.j.floordiv(20)
    return c

def center_within(a,w,g):
    a=np.asarray(a,float)
    result=np.empty_like(a)
    for k in np.unique(g):
        take=g==k
        result[take]=a[take]-np.average(a[take],axis=0,weights=w[take])
    return result

def estimate(y,X,w,block,target):
    y,X,w,block=map(np.asarray,(y,X,w,block))
    groups=np.unique(block)
    if len(y)<max(25,3*X.shape[1]) or len(groups)<8 or np.linalg.matrix_rank(X)<X.shape[1]:
        return dict(status="insufficient_support",p=1.,p_sign=1.,cells=len(y),blocks=len(groups))
    fit=sm.WLS(y,X,weights=w).fit(cov_type="cluster",cov_kwds={"groups":block},use_t=True)
    lo,hi=np.asarray(fit.conf_int())[target]
    # Approximate independent-block score sign symmetry under a restricted model.
    # Weighted Frisch-Waugh residuals, keeping all within-group transformations.
    sw=np.sqrt(w)
    a=X*sw[:,None];b=y*sw
    z=np.delete(a,target,axis=1)
    if z.shape[1]:
        xr=a[:,target]-z@np.linalg.lstsq(z,a[:,target],rcond=None)[0]
        yr=b-z@np.linalg.lstsq(z,b,rcond=None)[0]
    else:xr,yr=a[:,target],b
    contributions=np.array([xr[block==g]@yr[block==g] for g in groups])
    k=len(groups);exact=k<=15
    signs=np.array(list(product([-1,1],repeat=k))) if exact else np.random.default_rng(SEED).choice([-1,1],(9999,k))
    null=abs(signs@contributions)
    p_sign=(np.count_nonzero(null>=abs(contributions.sum())-1e-12)+int(not exact))/(len(signs)+int(not exact))
    loo=[]
    for g in groups:
        keep=block!=g
        if np.linalg.matrix_rank(a[keep])==a.shape[1]:
            loo.append(np.linalg.lstsq(a[keep],b[keep],rcond=None)[0][target])
    return dict(status="estimated",cells=len(y),blocks=k,beta=float(fit.params[target]),
                lo=float(lo),hi=float(hi),p=float(fit.pvalues[target]),p_sign=float(p_sign),
                loo_min=min(loo,default=np.nan),loo_max=max(loo,default=np.nan),
                effective_cells=float(w.sum()**2/(w@w)))

def design(d,estimand,within,extra):
    x=d.x.to_numpy()
    if estimand=="introduced":
        X=x[:,None] if within else np.column_stack([np.ones(len(d)),x])
    else:
        z=d.introduced.to_numpy()
        X=np.column_stack([x,x*z]) if within else np.column_stack([np.ones(len(d)),z,x,x*z])
    target=X.shape[1]-1
    if extra:X=np.column_stack([X,d[extra].to_numpy(float)])
    y=d.y.to_numpy();w=d.w.to_numpy();block=d.block.to_numpy()
    if within:
        groups=block*2+d.introduced.to_numpy() if estimand=="interaction" else block
        X=center_within(X,w,groups);y=center_within(y,w,groups)
    return estimate(y,X,w,block,target)

def run_sensitivities(sp,estimand,d,out):
    base=design(d,estimand,False,[])
    target=EXPECTED[sp][1 if estimand=="introduced" else 2]
    assert base["status"]=="estimated" and abs(base["beta"]-target)<1e-8,(sp,estimand,base.get("beta"),target)
    if estimand=="introduced":assert len(d)==EXPECTED[sp][0]
    common=d.dropna(subset=PROXIES+[LIGHT]).copy()
    qualified=d.loc[(d.records>=5)&(d.observers>=3)].copy()
    specs=[("original",d,False,[]),
           ("region_within",d,True,[]),
           ("complete_proxy_baseline",common,False,[]),
           ("complete_proxy_region",common,True,[]),
           ("region_observation_proxies",common,True,PROXIES),
           ("region_observation_proxies_plus_Lstar",common,True,PROXIES+[LIGHT]),
           ("higher_record_support_baseline",qualified,False,[]),
           ("higher_record_support_region",qualified,True,[])]
    for name,g,within,extra in specs:
        result=design(g,estimand,within,extra) if len(g) else dict(status="insufficient_support",p=1.,p_sign=1.,cells=0,blocks=0)
        out.append(dict(species=sp,estimand=estimand,specification=name,
            reference_beta=target,n_native=int((g.introduced==0).sum()),n_intro=int((g.introduced==1).sum()),**result))

def shared_weights(n,i,factor=4):
    groups=[g[g.common].copy() for g in [n,i]]
    ids=[(g.u.to_numpy()//factor)*80+g.v.to_numpy()//factor for g in groups]
    freq=[np.bincount(a,minlength=6400)/len(a) for a in ids]
    target=np.minimum(*freq);mass=target.sum()
    for k,(g,a) in enumerate(zip(groups,ids)):
        g["w"]=target[a]/freq[k][a]/(len(a)*mass)
        g["introduced"]=k
    d=pd.concat(groups,ignore_index=True).query("w>0").copy()
    np.testing.assert_allclose(d.groupby("introduced").w.sum().to_numpy(),[1,1])
    return d

def observation_coverage(sp,side,f,transform,edges,nat,common,out):
    q=f.copy()
    pc=transform(q);q["p1"],q["p2"]=pc.T
    q["i"]=np.minimum(179,np.floor(q.latitude+90)).astype(int)
    q["j"]=np.floor((q.longitude+180)%360).astype(int)
    q["measured_row"]=q.has_measured_row.fillna(False).astype(float)
    q["colour_usable"]=(q[HEADS].gt(0)&np.isfinite(q[CHROMA])).astype(float)
    q["three_complete"]=(q[T[0]].between(0,180)&np.isfinite(q[T].to_numpy(float)).all(1)).astype(float)
    q["research_grade"]=q.quality_grade.eq("research").astype(float)
    q["date_available"]=pd.to_datetime(q.observed_on,format="mixed",errors="coerce").notna().astype(float)
    metrics=["measured_row","colour_usable","three_complete","research_grade","date_available"]
    g=q.groupby(["i","j"])
    c=g[metrics].mean().join(g[["p1","p2"]].median()).reset_index()
    c["source_records"]=g.size().to_numpy()
    uv=np.column_stack([np.searchsorted(e,c[key],side="right")-1 for e,key in zip(edges,["p1","p2"])])
    for k,e in enumerate(edges):uv[c[["p1","p2"][k]].to_numpy()==e[-1],k]=len(e)-2
    ok=np.all((uv>=0)&(uv<np.array([len(e)-1 for e in edges])),axis=1)
    c=c.loc[ok].copy();u,v=uv[ok].T
    c["common"]=common[u,v]
    c["outside"]=~_probability_content(nat,.95)[u,v]
    c["rarity"]=mass_rank(nat)[u,v]
    c["block"]=c.i.floordiv(20)*18+c.j.floordiv(20)
    for scope,part in [("all_classified",c),("common_background",c[c.common])]:
        for outside,b in part.groupby("outside"):
            for m in metrics:
                out["coverage"].append(dict(species=sp,range_state=side,scope=scope,
                    native_core_outside=bool(outside),metric=m,cells=len(b),blocks=b.block.nunique(),
                    source_records=int(b.source_records.sum()),cell_equal_fraction=b[m].mean()))
    # Primary missingness diagnostic: within-region change per .1 rarity,
    # in all source cells in common background, including unmeasured cells.
    c=c[c.common].copy();w=np.ones(len(c));b=c.block.to_numpy()
    X=center_within(c.rarity.to_numpy()[:,None]*10,w,b)
    for m in ["measured_row","colour_usable","three_complete"]:
        y=center_within(c[m].to_numpy(),w,b)
        if np.std(y)<1e-12:
            r=dict(status="constant_response",cells=len(c),blocks=len(np.unique(b)),p=1.,p_sign=1.)
        else:r=estimate(y,X,w,b,0)
        out["coverage_slopes"].append(dict(species=sp,range_state=side,metric=m,**r))

def run_species(sp,sn,si,env,transform,edges,nat,inv,out):
    bn,bi=[background(env[env._range==k],transform,edges) for k in ["native","introduced"]]
    common=_probability_content(bn,.99)&_probability_content(bi,.99)
    n,i=[enrich_cells(f,head_cells(f,transform,edges)) for f in [sn,si]]
    r=mass_rank(nat)
    for side,c in enumerate([n,i]):
        c["rarity"]=r[c.u,c.v];c["common"]=common[c.u,c.v];c["introduced"]=side
    d=i[i.common].copy();d["x"]=d.rarity*10;d["y"]=d[CHROMA];d["w"]=1.
    run_sensitivities(sp,"introduced",d,out["sensitivity"])
    d=shared_weights(n,i)
    x=d.rarity.to_numpy();w=d.w.to_numpy()
    mu=np.average(x,weights=w);sd=np.sqrt(np.average((x-mu)**2,weights=w))
    d["x"]=(x-mu)/sd;d["y"]=d[CHROMA]/n[CHROMA].std(ddof=0)
    run_sensitivities(sp,"interaction",d,out["sensitivity"])
    for side,f in [("native",sn),("introduced",si)]:
        observation_coverage(sp,side,f,transform,edges,nat,common,out)

def main():
    OUT.mkdir(parents=True,exist_ok=True)
    auth=authority_source((ROOT/"private/authority350.zip").read_bytes())
    wild=pd.read_csv(ROOT/"private/focal_wild.csv.gz",dtype={"obs_id":str},low_memory=False)
    assert len(wild)==297622
    wild["_cell"]=angular_cells(wild.latitude,wild.longitude)
    wild=wild.sort_values("obs_id").drop_duplicates(["taxon_name","range_state","_cell","user_id"]).copy()
    land=background_mesh(auth)
    wild,land,_=apply_original_five_CHELSA(wild,land,source_meta=Path("analysis/v4/environment_source_metadata.json"),work=ROOT/"private/rasters")
    path=ROOT/"private/original-measured.csv.gz"
    measured,mapping=load_measures(path)
    counts=pd.read_csv(path,usecols=["obs_id",HEADS],dtype={"obs_id":str})
    counts[HEADS]=pd.to_numeric(counts[HEADS],errors="coerce")
    measured=measured.merge(counts,on="obs_id",validate="one_to_one")
    measured["has_measured_row"]=True
    assert mapping["corolla_lab_chroma"][1]==HEADS
    wild=wild.drop(columns=[c for c in measured if c!="obs_id"],errors="ignore").merge(measured,on="obs_id",how="left",validate="one_to_one")
    frozen=json.loads(Path("analysis/v4/results/20261021_unified_abcde/S2_TWO_SPECIES_FULL_WCVP_NICHE_VERIFIED_20261021.json").read_text())
    out={k:[] for k in ["sensitivity","coverage","coverage_slopes"]}
    for sp,rows in wild.groupby("taxon_name"):
        env,_=_background_species_audit(auth,land.dropna(subset=list(VARIABLES)),sp,land)
        complete=rows.dropna(subset=list(VARIABLES))
        sn=complete.loc[complete.range_state=="native"];si=complete.loc[complete.range_state=="introduced"]
        transform,edges,nat,inv=_pca_env_2d(env[env._range=="native"],env[env._range=="introduced"],sn,si,return_grid=True)
        target=frozen["species"][sp]["PCA_env_background_corrected_Schoener_D"]
        assert abs(1-abs(nat-inv).sum()/2-target)<1e-5
        run_species(sp,sn,si,env,transform,edges,nat,inv,out)
    for name,rows in out.items():
        d=pd.DataFrame(rows)
        if name=="sensitivity":
            assert len(d)==32
            for p in ["p","p_sign"]:
                d["holm_"+p]=d.groupby("specification")[p].transform(lambda s:multipletests(s.fillna(1),method="holm")[1])
        if name=="coverage_slopes":
            assert len(d)==12
            for p in ["p","p_sign"]:
                d["holm_"+p]=multipletests(d[p].fillna(1),method="holm")[1]
        d.to_csv(OUT/(name+".csv"),index=False)
    print("Complete; aggregate outputs only.",flush=True)

if __name__=="__main__":
    main()
