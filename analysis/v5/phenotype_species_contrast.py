"""Exploratory difference of fixed, within-species phenotype enrichment curves.
Species-relative quintiles and species-specific climate definitions; not
identical environmental treatments, selection coefficients or fitness.
"""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy.stats import t as student_t
from statsmodels.stats.multitest import multipletests
from analysis.v5.phenotype_enrichment_transfer import (
    ROOT,NAMES,COLS,EXPECTED,EXPECTED_BETA,joint_records,T,eligible,
    load_three_axis,make_cells,match_cells,head_cells,background_mesh,
    _background_species_audit,_pca_env_2d,_probability_content,
    apply_original_five_CHELSA,authority_source,angular_cells,VARIABLES,
    X,K,native_cuts,rank_codes,cell_hist,fit,
)
OUT=ROOT/"phenotype_species_contrast"

def fit_and_influence(a,b,blocks):
    beta,score,h=fit(a,b,X,np.arange(6))
    ids=np.unique(blocks)
    summed=np.stack([score[blocks==g].sum(0) for g in ids])
    influence=summed@np.linalg.inv(h)*np.sqrt(len(ids)/(len(ids)-1))
    return beta,dict(zip(ids,influence))

def direct_contrast(ca,cv):
    ba,ia=fit_and_influence(*ca)
    bv,iv=fit_and_influence(*cv)
    ids=sorted(set(ia)|set(iv))
    zero=np.zeros(6)
    dif=np.stack([iv.get(g,zero)-ia.get(g,zero) for g in ids])
    se=np.sqrt(np.square(dif).sum(0))
    df=min(len(ia),len(iv))-1
    delta=bv-ba
    crit=student_t.ppf(.975,df)
    p=np.ones(6)
    use=se>1e-12
    p[use]=2*student_t.sf(abs(delta[use]/se[use]),df)
    return [dict(trait=NAMES[j%3],term="directional" if j<3 else "quadratic",
                 ca_coefficient=float(ba[j]),cv_coefficient=float(bv[j]),
                 cv_minus_ca=float(delta[j]),cluster_se=float(se[j]),
                 lo=float(delta[j]-crit*se[j]),hi=float(delta[j]+crit*se[j]),
                 p=float(p[j]),df=df,ca_blocks=len(ia),cv_blocks=len(iv),
                 shared_blocks=len(set(ia)&set(iv)),union_blocks=len(ids),
                 status="estimated" if use[j] else "degenerate_se")
            for j in range(6)]

def common_support(ca,cv):
    ba,bv=ca[2],cv[2]
    ids=np.intersect1d(ba,bv)
    return ids,np.array([min(int((ba==g).sum()),int((bv==g).sum())) for g in ids],int)

def balanced_diagnostic(ca,cv):
    ids,counts=common_support(ca,cv)
    info=dict(ca_pairs=len(ca[0]),cv_pairs=len(cv[0]),ca_blocks=len(np.unique(ca[2])),
              cv_blocks=len(np.unique(cv[2])),shared_blocks=len(ids),
              balanced_pairs_per_species=int(counts.sum()),planned_draws=32,
              valid_draws=0,status="insufficient_support")
    def unavailable(status):
        return [dict(trait=NAMES[j%3],term="directional" if j<3 else "quadratic",status=status) for j in range(6)]
    if len(ids)<8 or counts.sum()<20:
        return info,unavailable("insufficient_support")
    rng=np.random.default_rng(20261126)
    allb=[]
    for _ in range(32):
        fit_draw=[]
        try:
            for a,b,block in (ca,cv):
                ix=np.concatenate([rng.choice(np.flatnonzero(block==g),k,replace=False)
                                   for g,k in zip(ids,counts)])
                beta,_,_=fit(a[ix],b[ix],X,np.arange(6))
                fit_draw.append(beta)
        except (ValueError,np.linalg.LinAlgError):
            continue
        allb.append(fit_draw)
    info["valid_draws"]=len(allb)
    info["status"]="estimated" if len(allb)==32 else "fit_failed"
    if len(allb)!=32:
        return info,unavailable("fit_failed")
    b=np.asarray(allb)
    delta=b[:,1]-b[:,0]
    return info,[dict(trait=NAMES[j%3],term="directional" if j<3 else "quadratic",
                      status="estimated",ca_median=float(np.median(b[:,0,j])),
                      cv_median=float(np.median(b[:,1,j])),
                      difference_median=float(np.median(delta[:,j])),
                      difference_min=float(delta[:,j].min()),
                      difference_max=float(delta[:,j].max()),
                      cv_gt_ca_draws=int((delta[:,j]>0).sum())) for j in range(6)]

def analyze(sp,definition,n,c,pairs,fn,fi,out):
    nv,ni,nw=joint_records(fn,n,COLS,False)
    iv,ii,_=joint_records(fi,c,COLS,False)
    assert len(nv)==EXPECTED[sp][4]
    cuts,lo,hi=native_cuts(nv,nw)
    rc,valid=rank_codes(nv,cuts,lo,hi)
    assert valid.all()
    support=np.bincount(rc,minlength=K)>0
    assert support.sum()==125
    h,tot,kept=cell_hist(iv,ii,len(c),cuts,lo,hi,support)
    left=np.concatenate([l for l,r in pairs]);right=np.concatenate([r for l,r in pairs])
    blocks=c.iloc[left].block.to_numpy()
    np.testing.assert_array_equal(blocks,c.iloc[right].block.to_numpy())
    assert np.all(kept[np.r_[left,right]]>0)
    a,b=h[left],h[right]
    beta,_,_=fit(a,b,X,np.arange(6))
    np.testing.assert_allclose(beta,EXPECTED_BETA[(sp,definition)],atol=1e-6,rtol=0)
    out[(sp,definition)]=(a,b,blocks)

def run_species(sp,sn,si,transform,edges,nat,inv,out):
    fn,fi=eligible(sn),eligible(si)
    n=make_cells(sn,True);full=make_cells(si,True)
    ex=EXPECTED[sp]
    assert (len(n),len(full))==ex[:2]
    coords=transform(si)
    uv=np.column_stack([np.searchsorted(e,coords[:,k],side="right")-1 for k,e in enumerate(edges)])
    for k,e in enumerate(edges):uv[coords[:,k]==e[-1],k]=len(e)-2
    inside=np.all((uv>=0)&(uv<np.array([len(e)-1 for e in edges])),axis=1)
    labeled=si.loc[inside].copy()
    labeled["novel"]=(~_probability_content(nat,.95)[uv[inside,0],uv[inside,1]]).astype(int)
    c=make_cells(labeled);c=c.loc[c.novel!=.5].reset_index(drop=True)
    pairs=match_cells(c)
    assert sum(len(l) for l,r in pairs)==ex[2]
    analyze(sp,"native_core95",n,c,pairs,fn,fi,out)
    hc=head_cells(si,transform,edges)[["i","j","u","v"]]
    c=full.merge(hc,on=["i","j"],validate="one_to_one")
    u,v=c.u.to_numpy(),c.v.to_numpy()
    c=c.loc[inv[u,v]>0].reset_index(drop=True);u,v=c.u.to_numpy(),c.v.to_numpy()
    c["novel"]=(inv[u,v]>nat[u,v]).astype(int)
    pairs=match_cells(c)
    assert sum(len(l) for l,r in pairs)==ex[3]
    analyze(sp,"density_surplus",n,c,pairs,fn,fi,out)

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
    out={}
    for sp,rows in wild.groupby("taxon_name"):
        env,_=_background_species_audit(auth,land.dropna(subset=list(VARIABLES)),sp,land)
        complete=rows.dropna(subset=list(VARIABLES))
        sn=complete.loc[complete.range_state=="native"];si=complete.loc[complete.range_state=="introduced"]
        transform,edges,nat,inv=_pca_env_2d(env[env._range=="native"],env[env._range=="introduced"],sn,si,return_grid=True)
        assert abs(1-abs(nat-inv).sum()/2-frozen["species"][sp]["PCA_env_background_corrected_Schoener_D"])<1e-5
        run_species(sp,sn,si,transform,edges,nat,inv,out)
    contrasts=[];balanced=[];coverage=[]
    for definition in ["native_core95","density_surplus"]:
        ca=out[("Cirsium arvense",definition)]
        cv=out[("Cirsium vulgare",definition)]
        contrasts.extend(dict(definition=definition,**r) for r in direct_contrast(ca,cv))
        info,rows=balanced_diagnostic(ca,cv)
        coverage.append(dict(definition=definition,**info))
        balanced.extend(dict(definition=definition,**r) for r in rows)
    d=pd.DataFrame(contrasts)
    assert len(d)==12
    d["holm12"]=multipletests(d.p,method="holm")[1]
    d["holm6"]=d.groupby("definition").p.transform(lambda p:multipletests(p,method="holm")[1])
    d.to_csv(OUT/"species_contrast.csv",index=False)
    pd.DataFrame(balanced).to_csv(OUT/"balanced_coefficients.csv",index=False)
    pd.DataFrame(coverage).to_csv(OUT/"shared_region_coverage.csv",index=False)
    print("Complete; aggregate outputs only.",flush=True)

if __name__=="__main__":
    main()
