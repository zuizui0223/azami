"""Fixed three-axis dispersion: geographic differentiation vs within-cell spread.

Continuous moments, not histogram occupancy or geometric hypervolume.
All axes scaled once per species by full native equal-cell record moments.
The same V5 climate definitions, source cohort, and geographic pairs are used.
No lightness, time axis, genetic inference or new prediction model.
"""
from pathlib import Path
from itertools import product
import json
import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests
from analysis.v5.phenotype_three_axis import (
    ROOT, SEED, eligible, load_three_axis, AXES, joint_records, T,
    make_cells, match_cells, head_cells, background_mesh,
    _background_species_audit, _pca_env_2d, _probability_content,
    apply_original_five_CHELSA, authority_source, angular_cells, VARIABLES,
)
OUT=ROOT/"phenotype_three_axis_dispersion"
NAMES=list(AXES)
COLS=["scan_"+key for key,_ in AXES.values()]
METRICS=("total_trace","between_cell_trace")
EXPECTED={"Cirsium arvense":(1366,860,30,183,18505),
          "Cirsium vulgare":(1113,1139,146,316,13115)}

def cell_moments(v,ix,n):
    """One unit per geographic cell, retain all eligible observations within it."""
    v=np.asarray(v,float)
    ix=np.asarray(ix,int)
    if v.ndim!=2 or v.shape[1]!=3 or not np.isfinite(v).all():
        raise ValueError("Expected finite joint three-axis records")
    if len(ix)!=len(v) or np.any((ix<0)|(ix>=n)):
        raise ValueError("Invalid cell indices")
    counts=np.bincount(ix,minlength=n)
    means=np.stack([np.bincount(ix,weights=v[:,j],minlength=n) for j in range(3)],axis=1)
    seconds=np.stack([np.bincount(ix,weights=v[:,j]**2,minlength=n) for j in range(3)],axis=1)
    good=counts>0
    means[good]/=counts[good,None]
    seconds[good]/=counts[good,None]
    means[~good]=np.nan
    seconds[~good]=np.nan
    # Per-cell sufficient moments for exact mixture and its decomposition:
    # [first moment, record second moment, squared cell centroid].
    m=np.concatenate([means,seconds,means**2],axis=1)
    return m,counts

def summarize(m):
    """Vectorized group mean of sufficient moments -> variance decomposition."""
    a,q,b=np.split(np.asarray(m,float),3,axis=-1)
    total=q-a*a
    between=b-a*a
    within=q-b
    if min(np.min(total),np.min(between),np.min(within))< -1e-8:
        raise ValueError("Negative variance beyond numerical tolerance")
    total=np.maximum(total,0.)
    between=np.maximum(between,0.)
    within=np.maximum(within,0.)
    np.testing.assert_allclose(total,between+within,atol=1e-8,rtol=1e-8)
    return a,total,between,within

def group_metrics(m):
    _,total,between,_=summarize(m)
    return np.stack([total.sum(-1),between.sum(-1)],axis=-1)

def reference_scale(nv,ni,n):
    m,count=cell_moments(nv,ni,n)
    good=count>0
    mu,var,_,_=summarize(m[good].mean(0))
    if np.any(var<=1e-10):
        raise ValueError("Zero native variance")
    return mu,np.sqrt(var)

def paired_block_sums(m,pairs):
    return (np.stack([m[l].sum(0) for l,r in pairs]),
            np.stack([m[r].sum(0) for l,r in pairs]),
            np.array([len(l) for l,r in pairs]))

def delta_from_sums(l,r,ns):
    ml=l.sum(0)/ns.sum()
    mr=r.sum(0)/ns.sum()
    return group_metrics(mr)-group_metrics(ml)

def block_test(m,pairs):
    """Recompute both centroids for every whole-region label swap and bootstrap."""
    l,r,ns=paired_block_sums(m,pairs)
    b=len(ns)
    assert 8<=b<=15
    signs=np.array(list(product([1.,-1.],repeat=b)))
    centre=(l+r).sum(0)/2
    offset=signs@(r-l)/2
    nl=(centre-offset)/ns.sum()
    nr=(centre+offset)/ns.sum()
    delta=group_metrics(nr)-group_metrics(nl)
    obs=delta_from_sums(l,r,ns)
    np.testing.assert_allclose(delta[0],obs,atol=1e-10)
    p=np.mean(np.abs(delta)>=np.abs(obs)[None,:]-1e-12,axis=0)
    rng=np.random.default_rng(SEED)
    boots=rng.multinomial(b,np.full(b,1/b),size=3999)
    denom=boots@ns
    dl=boots@l/denom[:,None]
    dr=boots@r/denom[:,None]
    boot=group_metrics(dr)-group_metrics(dl)
    lo,hi=np.quantile(boot,[.025,.975],axis=0)
    first,second=group_metrics(l.sum(0)/ns.sum()),group_metrics(r.sum(0)/ns.sum())
    null_lo,null_hi=np.quantile(delta,[.025,.975],axis=0)
    return [dict(metric=metric,analogue=float(first[j]),novel=float(second[j]),
                 difference=float(obs[j]),ratio=float(second[j]/first[j]) if first[j]>0 else np.nan,
                 ci_low=float(lo[j]),ci_high=float(hi[j]),p=float(p[j]),
                 null_q025=float(null_lo[j]),null_q975=float(null_hi[j]),
                 permutations=len(signs),pairs=int(ns.sum()),blocks=b)
            for j,metric in enumerate(METRICS)]

def report_group(out,sp,definition,group,m,counts,indices,mu,sd):
    indices=np.asarray(indices,int)
    assert len(indices)>0 and len(np.unique(indices))==len(indices)
    assert np.all(counts[indices]>0) and np.isfinite(m[indices]).all()
    a,total,between,within=summarize(m[indices].mean(0))
    d=dict(species=sp,definition=definition,group=group,
           cells=len(indices),records=int(counts[indices].sum()),
           total_trace=float(total.sum()),between_cell_trace=float(between.sum()),
           within_cell_trace=float(within.sum()),
           centroid_distance_native=float(np.linalg.norm(a)))
    for j,name in enumerate(NAMES):
        d.update({name+"_mean":float(mu[j]+sd[j]*a[j]),name+"_mean_z":float(a[j]),
                  name+"_sd_z":float(np.sqrt(total[j])),
                  name+"_between_variance_z":float(between[j]),
                  name+"_within_variance_z":float(within[j])})
    out["groups"].append(d)

def analyze(sp,definition,n,c,pairs,fn,fi,out,mu,sd):
    iv,ii,_=joint_records(fi,c,COLS,False)
    m,counts=cell_moments((iv-mu)/sd,ii,len(c))
    good=counts>0
    left=np.concatenate([l for l,r in pairs])
    right=np.concatenate([r for l,r in pairs])
    assert good[left].all() and good[right].all()
    for group,idx in [
        ("analogous_all",np.flatnonzero(good & (c.novel.to_numpy()==0))),
        ("novel_all",np.flatnonzero(good & (c.novel.to_numpy()==1))),
        ("analogous_matched",left),
        ("novel_matched",right),
    ]:
        report_group(out,sp,definition,group,m,counts,idx,mu,sd)
    for d in block_test(m,pairs):
        out["tests"].append(dict(species=sp,definition=definition,**d))

def run_species(sp,sn,si,transform,edges,nat,inv,out):
    fn,fi=eligible(sn),eligible(si)
    n=make_cells(sn,True)
    full=make_cells(si,True)
    ex=EXPECTED[sp]
    assert (len(n),len(full))==ex[:2]
    nv,ni,_=joint_records(fn,n,COLS,False)
    assert len(nv)==ex[4]
    mu,sd=reference_scale(nv,ni,len(n))
    m,count=cell_moments((nv-mu)/sd,ni,len(n))
    assert np.all(count>0)
    report_group(out,sp,"all_ranges","native",m,count,np.flatnonzero(count>0),mu,sd)
    np.testing.assert_allclose(out["groups"][-1]["total_trace"],3.,atol=1e-8)
    iv,ii,_=joint_records(fi,full,COLS,False)
    m,count=cell_moments((iv-mu)/sd,ii,len(full))
    report_group(out,sp,"all_ranges","introduced",m,count,np.flatnonzero(count>0),mu,sd)
    # Reuse exact original climate labels, including exclusions of mixed cells.
    coords=transform(si)
    uv=np.column_stack([np.searchsorted(e,coords[:,k],side="right")-1 for k,e in enumerate(edges)])
    for k,e in enumerate(edges):uv[coords[:,k]==e[-1],k]=len(e)-2
    inside=np.all((uv>=0)&(uv<np.array([len(e)-1 for e in edges])),axis=1)
    labeled=si.loc[inside].copy()
    labeled["novel"]=(~_probability_content(nat,.95)[uv[inside,0],uv[inside,1]]).astype(int)
    c=make_cells(labeled);c=c.loc[c.novel!=.5].reset_index(drop=True)
    pairs=match_cells(c)
    assert sum(len(l) for l,r in pairs)==ex[2]
    analyze(sp,"native_core95",n,c,pairs,fn,fi,out,mu,sd)
    hc=head_cells(si,transform,edges)[["i","j","u","v"]]
    c=full.merge(hc,on=["i","j"],validate="one_to_one")
    u,v=c.u.to_numpy(),c.v.to_numpy()
    c=c.loc[inv[u,v]>0].reset_index(drop=True);u,v=c.u.to_numpy(),c.v.to_numpy()
    c["novel"]=(inv[u,v]>nat[u,v]).astype(int)
    pairs=match_cells(c)
    assert sum(len(l) for l,r in pairs)==ex[3]
    analyze(sp,"density_surplus",n,c,pairs,fn,fi,out,mu,sd)

def main():
    OUT.mkdir(parents=True,exist_ok=True)
    auth=authority_source((ROOT/"private/authority350.zip").read_bytes())
    wild=pd.read_csv(ROOT/"private/focal_wild.csv.gz",dtype={"obs_id":str},low_memory=False)
    assert len(wild)==297622
    wild["_cell"]=angular_cells(wild.latitude,wild.longitude)
    wild=wild.sort_values("obs_id").drop_duplicates(["taxon_name","range_state","_cell","user_id"]).copy()
    land=background_mesh(auth)
    wild,land,_=apply_original_five_CHELSA(wild,land,source_meta=Path("analysis/v4/environment_source_metadata.json"),work=ROOT/"private/rasters")
    measured=load_three_axis(ROOT/"private/original-measured.csv.gz")
    wild=wild.drop(columns=[c for c in measured if c!="obs_id"],errors="ignore").merge(measured,on="obs_id",how="left",validate="one_to_one")
    frozen=json.loads(Path("analysis/v4/results/20261021_unified_abcde/S2_TWO_SPECIES_FULL_WCVP_NICHE_VERIFIED_20261021.json").read_text())
    out={"groups":[],"tests":[]}
    for sp,rows in wild.groupby("taxon_name"):
        env,_=_background_species_audit(auth,land.dropna(subset=list(VARIABLES)),sp,land)
        complete=rows.dropna(subset=list(VARIABLES))
        sn=complete.loc[complete.range_state=="native"];si=complete.loc[complete.range_state=="introduced"]
        transform,edges,nat,inv=_pca_env_2d(env[env._range=="native"],env[env._range=="introduced"],sn,si,return_grid=True)
        assert abs(1-abs(nat-inv).sum()/2-frozen["species"][sp]["PCA_env_background_corrected_Schoener_D"])<1e-5
        run_species(sp,sn,si,transform,edges,nat,inv,out)
    d=pd.DataFrame(out["tests"])
    assert len(d)==8
    d["holm_all8"]=multipletests(d.p.fillna(1),method="holm")[1]
    d.to_csv(OUT/"dispersion_tests.csv",index=False)
    g=pd.DataFrame(out["groups"])
    assert len(g)==20
    g.to_csv(OUT/"dispersion_groups.csv",index=False)
    print("Complete; aggregate outputs only.",flush=True)

if __name__=="__main__":
    main()
