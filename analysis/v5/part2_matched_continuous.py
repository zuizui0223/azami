"""Continuous phenotype density ratios in matched native/introduced climates.

Same 5 annual climate axes, no novelty or trait bins, no truncation to native
trait extrema. Geographic regions may differ across ranges; the estimand is
a descriptive conditional association, not selection or temporal evolution.
"""
from pathlib import Path
import sys
import numpy as np
import pandas as pd
from scipy.stats import t
from statsmodels.stats.multitest import multipletests
from analysis.v5.part2_annual_continuous import (
    ROOT, SPECIES, AXES, NAMES, RANGES, cell_data, support, authority_source,
    angular_cells, background_mesh, apply_original_five_CHELSA, load_three_axis,
    MIN_PAIRS, MIN_REGIONS, CALIPER, MAX_AXIS,
)
from analysis.v5.phenotype_environment_explanation import match_one
from analysis.v5.phenotype_continuous_enrichment import raw_pairs, fit_raw, native_reference, curve_grid
from analysis.v5.phenotype_redistribution import weighted_quantile
OUT=ROOT/"part2_matched_continuous"

def grouped_meat(score, groups):
    _,ix=np.unique(groups,return_inverse=True,axis=0)
    sums=np.zeros((ix.max()+1,score.shape[1]))
    np.add.at(sums,ix,score)
    return sums.T@sums

def dyadic_cov(score,h,gn,gi):
    # Allow dependence whenever two pairs share either geographical block,
    # regardless of whether the block occurs on the native or introduced side.
    # Inclusion-exclusion: subtract the double-counted two-distinct-block pairs.
    # Self-loop pairs are counted once. This covariance can be indefinite.
    distinct=gn!=gi
    blocks=np.r_[gn,gi[distinct]]
    s=np.vstack([score,score[distinct]])
    meat=grouped_meat(s,blocks[:,None])
    if distinct.any():
        both=np.sort(np.column_stack([gn[distinct],gi[distinct]]),axis=1)
        meat-=grouped_meat(score[distinct],both)
    inv=np.linalg.inv(h)
    clusters=min(len(np.unique(gn)),len(np.unique(gi)))
    cov=inv@meat@inv*clusters/(clusters-1)
    return (cov+cov.T)/2,clusters-1

def run_species(sp,frame,out):
    data=[cell_data(frame.loc[frame.range_state==r],1) for r in RANGES]
    (a,va,ia),(b,vb,ib)=data
    ca=a[list(AXES)].to_numpy();cb=b[list(AXES)].to_numpy()
    mu=ca.mean(0);sd=ca.std(0)
    assert np.isfinite(sd).all() and (sd>0).all()
    xa=(ca-mu)/sd; xb=(cb-mu)/sd
    shared=support(xa,xb)
    l,r=match_one(xa,xb,tuple(range(5)))
    assert shared[0][l].all() and shared[1][r].all()
    gn=a.block.to_numpy()[l]; gi=b.block.to_numpy()[r]
    n,m=len(np.unique(gn)),len(np.unique(gi))
    meta=dict(species=sp,pairs=len(l),native_regions=n,introduced_regions=m)
    for j,(c,v,ix) in enumerate(data):
        ids=(l,r)[j]
        mask=np.zeros(len(c),bool);mask[ids]=True
        for status,keep in (("all",np.ones(len(c),bool)),("continuous_shared_support",shared[j]),
                            ("matched",mask),("no_neighbour_within_caliper",~shared[j])):
            out["cohort"].append(dict(**meta,range=RANGES[j],group=status,cells=int(keep.sum()),
                records=int(c.loc[keep,"count"].sum()),regions=c.loc[keep,"block"].nunique()))
            if not keep.any(): continue
            for k,axis in enumerate(AXES):
                vals=(ca,cb)[j][keep,k]
                qs=np.quantile(vals,[.1,.5,.9])
                out["climate"].append(dict(**meta,range=RANGES[j],group=status,axis=axis,
                      mean=vals.mean(),q10=qs[0],median=qs[1],q90=qs[2]))
    if len(l):
        delta=xb[r]-xa[l]
        assert np.max(np.sqrt(np.mean(delta**2,axis=1)))<=CALIPER+1e-10
        assert np.max(abs(delta))<=MAX_AXIS+1e-10
        for j,axis in enumerate(AXES):
            out["balance"].append(dict(**meta,axis=axis,signed_mean_native_sd=delta[:,j].mean(),
                 mean_abs_native_sd=abs(delta[:,j]).mean(),max_abs_native_sd=abs(delta[:,j]).max()))
    counts_a=a["count"].to_numpy();counts_b=b["count"].to_numpy()
    center,scale,qq=native_reference(va,ia,1/counts_a[ia],len(a))
    for j,(c,v,ix) in enumerate(data):
        ids=(l,r)[j]
        k=np.isin(ix,ids)
        if not k.any(): continue
        w=1/c["count"].to_numpy()[ix[k]];w/=w.sum()
        for j2,trait in enumerate(NAMES):
            q=weighted_quantile(v[k,j2],w,[.05,.1,.5,.9,.95])
            mean=w@v[k,j2];var=w@((v[k,j2]-mean)**2)
            out["phenotypes"].append(dict(**meta,range=RANGES[j],trait=trait,mean=mean,
                 sd=np.sqrt(var),q05=q[0],q10=q[1],median=q[2],q90=q[3],q95=q[4]))
    reason="limited_support"
    beta=cov=None;df=None
    if len(l)>=MIN_PAIRS and min(n,m)>=MIN_REGIONS:
        v=np.vstack([va,vb]);ix=np.r_[ia,ib+len(a)];kept=np.r_[counts_a,counts_b]
        d=raw_pairs(v,ix,np.ones(len(v),bool),l,r+len(a),kept,center,scale,gn)
        try:
            beta,s,h=fit_raw(d,np.arange(6))
            cov,df=dyadic_cov(s,h,gn,gi)
            reason="estimated"
        except (ValueError,AssertionError,np.linalg.LinAlgError):
            reason="fit_or_covariance_failure"
            beta=None
    for j,trait in enumerate(NAMES):
        for k,term in ((j,"linear"),(j+3,"quadratic")):
            row=dict(**meta,trait=trait,term=term,status=reason,native_center=center[j],native_sd=scale[j])
            if beta is not None:
                row["beta"]=float(beta[k])
                if cov is not None and np.isfinite(cov[k,k]) and cov[k,k]>0:
                    se=np.sqrt(cov[k,k]);crit=t.ppf(.975,df)
                    row.update(se=se,lower=beta[k]-crit*se,upper=beta[k]+crit*se,p=2*t.sf(abs(beta[k]/se),df),df=df)
                else:row["status"]="descriptive_invalid_covariance"
            out["coefficients"].append(row)
        if beta is not None:
            for value in curve_grid(qq,j):
                z=(value-center[j])/scale[j];f=np.zeros(6);f[j]=z;f[j+3]=z*z
                log=float(f@beta)
                row=dict(**meta,trait=trait,value=value,native_z=z,log_ratio=log,relative_density=np.exp(log))
                var=float(f@cov@f) if cov is not None else np.nan
                if np.isfinite(var) and var>=0:
                    delta=t.ppf(.975,df)*np.sqrt(var)
                    row.update(lower=np.exp(log-delta),upper=np.exp(log+delta))
                out["curves"].append(row)

def main():
    OUT.mkdir(parents=True,exist_ok=True)
    auth=authority_source((ROOT/"private/authority350.zip").read_bytes())
    wild=pd.read_csv(ROOT/"private/focal_wild.csv.gz",dtype={"obs_id":str},low_memory=False)
    assert len(wild)==297622
    wild["_cell"]=angular_cells(wild.latitude,wild.longitude)
    wild=wild.sort_values("obs_id").drop_duplicates(["taxon_name","range_state","_cell","user_id"]).copy()
    land=background_mesh(auth)
    wild,_,_=apply_original_five_CHELSA(wild,land,source_meta=Path("analysis/v4/environment_source_metadata.json"),work=ROOT/"private/rasters")
    measured=load_three_axis(ROOT/"private/original-measured.csv.gz")
    wild=wild.drop(columns=[c for c in measured if c!="obs_id"],errors="ignore").merge(measured,on="obs_id",how="left",validate="one_to_one")
    assert set(wild.taxon_name.unique())==set(SPECIES)
    out={k:[] for k in ("cohort","climate","balance","phenotypes","coefficients","curves")}
    for sp,frame in wild.groupby("taxon_name"):
        run_species(sp,frame.dropna(subset=list(AXES)),out)
    d=pd.DataFrame(out["coefficients"])
    assert len(d)==12
    if "p" not in d:d["p"]=np.nan
    d["holm12"]=multipletests(d.p.fillna(1),method="holm")[1]
    d.loc[d.status!="estimated","holm12"]=np.nan
    d.to_csv(OUT/"coefficients.csv",index=False)
    for k in out:
        if k!="coefficients":pd.DataFrame(out[k]).to_csv(OUT/(k+".csv"),index=False)
    print({"coefficients":len(d),"estimated":int((d.status=="estimated").sum()),"cohort_rows":len(out["cohort"])},flush=True)

def self_test():
    score=np.array([[1.,2.],[3.,4.],[5.,6.],[2.,1.],[4.,2.]])
    a=np.array([1,1,3,4,7]);b=np.array([2,1,4,3,8]);h=np.eye(2)
    # Compare the inclusion-exclusion formula against direct pairwise shared-block calculation.
    expected=np.zeros((2,2))
    for i in range(len(a)):
        for j in range(len(a)):
            if {a[i],b[i]} & {a[j],b[j]}:
                expected+=np.outer(score[i],score[j])
    clusters=min(len(np.unique(a)),len(np.unique(b)))
    cov,df=dyadic_cov(score,h,a,b)
    np.testing.assert_allclose(cov,expected*clusters/(clusters-1))
    assert df==clusters-1
    print("dyadic covariance test passed",flush=True)

if __name__=="__main__":
    self_test() if "--self-test" in sys.argv else main()
