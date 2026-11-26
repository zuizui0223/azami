"""Final fixed-three-axis block-label null diagnostic.

Do existing-bin positive gains, unseen-bin gains, total redistribution or
existing share exceed the amounts expected under whole-block label exchange?
This is a conditional cross-sectional diagnostic; no selection or ancestry
inference. The source, bins, splits and matched pairs reproduce the preceding
three-axis analysis exactly. No lightness or time variables.
"""
from pathlib import Path
from itertools import product
import json
import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests
from analysis.v5.phenotype_three_axis import (
    ROOT, SEED, eligible, load_three_axis, AXES, REFS, DEFS,
    N_REF, N_CAL, N_DRAWS, joint_records, codes, hist_by_cell,
    native_splits, one_per_cell, T, make_cells, match_cells, head_cells,
    background_mesh, _background_species_audit, _pca_env_2d, _probability_content,
    apply_original_five_CHELSA, authority_source, angular_cells, VARIABLES,
)
OUT=ROOT/"phenotype_three_axis_null"
SPACES={"joint_3axis":[key for key,_ in AXES.values()]}
METRICS=["existing_gain","unobserved_gain","total_redistribution","existing_share"]
# Previous verified aggregate, rounded to 8 decimal places. Reproduce both
# gain components before computing or interpreting the null comparison.
EXPECTED={
("Cirsium arvense","native_core95","full_native"):(.58458509,0.),
("Cirsium arvense","native_core95","rarefied_native"):(.53522809,.07631142),
("Cirsium arvense","density_surplus","full_native"):(.27506629,.00110666),
("Cirsium arvense","density_surplus","rarefied_native"):(.24590334,.04254531),
("Cirsium vulgare","native_core95","full_native"):(.36294241,0.),
("Cirsium vulgare","native_core95","rarefied_native"):(.30432533,.06263947),
("Cirsium vulgare","density_surplus","full_native"):(.19751446,.00033368),
("Cirsium vulgare","density_surplus","rarefied_native"):(.15732524,.04262595),
}

def partition_null(block_differences,support,signs,count):
    """A whole block flips all its paired cells; retain actual pair counts."""
    block_differences=np.asarray(block_differences,float)
    support=np.asarray(support,bool)
    np.testing.assert_allclose(block_differences.sum(1),0,atol=1e-9)
    keep=np.any(block_differences!=0,axis=0)
    x=block_differences[:,keep]/count
    s=support[keep]
    out=np.zeros((len(signs),2))
    for i in range(0,len(signs),512):
        a=signs[i:i+512]@x
        pos=np.maximum(a,0)
        e=pos[:,s].sum(1)
        u=pos[:,~s].sum(1)
        np.testing.assert_allclose(e+u,.5*np.abs(a).sum(1),atol=1e-9)
        out[i:i+len(a)]=np.column_stack([e,u])
    return out

def evaluate_draw(ref,rw,iv,ii,n,matched,signs,count):
    rc,k=codes(ref,rw,ref,False)
    nc,_=codes(ref,rw,iv,False)
    support=np.bincount(rc,weights=rw,minlength=k)>0
    h=hist_by_cell(nc,ii,n,k)
    block_differences=np.stack([(h[r]-h[l]).sum(0) for l,r in matched])
    assert np.isfinite(block_differences).all()
    return partition_null(block_differences,support,signs,count)

def summarize_null(arr):
    # Average E/U over identical native reference draws before the ratio.
    total=arr.sum(1)
    ratio=np.divide(arr[:,0],total,out=np.full(len(arr),np.nan),where=total>1e-12)
    matrix=np.column_stack([arr,total,ratio])
    result=[]
    for j,label in enumerate(METRICS):
        observed=matrix[0,j]
        v=matrix[:,j]
        v=v[np.isfinite(v)]
        if not np.isfinite(observed) or not len(v):
            result.append(dict(metric=label,status="undefined",p=1.))
            continue
        # First pattern is the observed all-positive orientation; enumerate
        # every pattern. One-sided upper-tail test fixed for every metric.
        p=float(np.mean(v>=observed-1e-12))
        lo,hi=np.quantile(v,[.025,.975])
        result.append(dict(metric=label,status="estimated",observed=float(observed),
          null_mean=float(v.mean()),null_q025=float(lo),null_q975=float(hi),
          excess_over_null=float(observed-v.mean()),p=p,valid_patterns=len(v),
          patterns=len(matrix)))
    return result

def analyze(sp,definition,space,n,c,pairs,fn,fi,out):
    cols=["scan_"+x for x in SPACES[space]]
    nv,ni,nw=joint_records(fn,n,cols,False)
    iv,ii,iw=joint_records(fi,c,cols,False)
    good_native=np.unique(ni)
    good_intro=np.bincount(ii,minlength=len(c))>0
    total=sum(len(l) for l,r in pairs)
    matched=[(l[good_intro[l]&good_intro[r]],r[good_intro[l]&good_intro[r]]) for l,r in pairs]
    matched=[(l,r) for l,r in matched if len(l)]
    count=sum(len(l) for l,r in matched)
    k=len(matched)
    assert count==total and len(good_native)>=N_REF+N_CAL and 8<=k<=15
    signs=np.array(list(product([1.,-1.],repeat=k)))
    np.testing.assert_allclose(signs[0],1)
    splits=native_splits(n.block.to_numpy()[good_native],nref=N_REF,ncal=N_CAL,draws=N_DRAWS)
    lookup={int(a):np.flatnonzero(ni==a) for a in good_native}
    rng=np.random.default_rng(SEED+1)
    for mode in REFS:
        if mode=="full_native":
            arr=evaluate_draw(nv,nw,iv,ii,len(c),matched,signs,count)
            draws=1
        else:
            arr=np.zeros((len(signs),2))
            for tr,ca in splits:
                ref=one_per_cell(nv,ni,good_native[tr],rng,lookup)
                # Consume the identical calibration draws to preserve reference RNG.
                _=one_per_cell(nv,ni,good_native[ca],rng,lookup)
                arr+=evaluate_draw(ref,np.ones(len(ref)),iv,ii,len(c),matched,signs,count)
            draws=len(splits)
            arr/=draws
        np.testing.assert_allclose(arr[0],EXPECTED[(sp,definition,mode)],atol=1e-7,rtol=0)
        for r in summarize_null(arr):
            out.append(dict(species=sp,definition=definition,reference=mode,space=space,
                     pairs=count,blocks=k,draws=draws,**r))

def run_species(sp,sn,si,transform,edges,nat,inv,out):
    fn,fi=eligible(sn),eligible(si)
    n=make_cells(sn,True);full=make_cells(si,True)
    expected={"Cirsium arvense":(1366,860,30,183),"Cirsium vulgare":(1113,1139,146,316)}[sp]
    assert (len(n),len(full))==expected[:2]
    coords=transform(si)
    uv=np.column_stack([np.searchsorted(e,coords[:,k],side="right")-1 for k,e in enumerate(edges)])
    for k,e in enumerate(edges):uv[coords[:,k]==e[-1],k]=len(e)-2
    inside=np.all((uv>=0)&(uv<np.array([len(e)-1 for e in edges])),axis=1)
    labeled=si.loc[inside].copy()
    labeled["novel"]=(~_probability_content(nat,.95)[uv[inside,0],uv[inside,1]]).astype(int)
    c=make_cells(labeled);c=c.loc[c.novel!=.5].reset_index(drop=True)
    pairs=match_cells(c);assert sum(len(l) for l,r in pairs)==expected[2]
    for space in SPACES:analyze(sp,"native_core95",space,n,c,pairs,fn,fi,out)
    hc=head_cells(si,transform,edges)[["i","j","u","v"]]
    c=full.merge(hc,on=["i","j"],validate="one_to_one")
    u,v=c.u.to_numpy(),c.v.to_numpy()
    c=c.loc[inv[u,v]>0].reset_index(drop=True);u,v=c.u.to_numpy(),c.v.to_numpy()
    c["novel"]=(inv[u,v]>nat[u,v]).astype(int)
    pairs=match_cells(c);assert sum(len(l) for l,r in pairs)==expected[3]
    for space in SPACES:analyze(sp,"density_surplus",space,n,c,pairs,fn,fi,out)

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
    out=[]
    for sp,rows in wild.groupby("taxon_name"):
        env,_=_background_species_audit(auth,land.dropna(subset=list(VARIABLES)),sp,land)
        complete=rows.dropna(subset=list(VARIABLES))
        sn=complete.loc[complete.range_state=="native"];si=complete.loc[complete.range_state=="introduced"]
        transform,edges,nat,inv=_pca_env_2d(env[env._range=="native"],env[env._range=="introduced"],sn,si,return_grid=True)
        assert abs(1-abs(nat-inv).sum()/2-frozen["species"][sp]["PCA_env_background_corrected_Schoener_D"])<1e-5
        run_species(sp,sn,si,transform,edges,nat,inv,out)
    d=pd.DataFrame(out)
    assert len(d)==2*len(REFS)*len(DEFS)*len(METRICS)==32
    d["holm_all32"]=multipletests(d.p.fillna(1),method="holm")[1]
    d.to_csv(OUT/"three_axis_partition_null.csv",index=False)
    print("Complete; aggregate outputs only.",flush=True)

if __name__=="__main__":
    main()
