"""Threshold-independent climate occupancy excess and native phenotype redistribution."""
from pathlib import Path
from itertools import product
import json
import numpy as np
import pandas as pd
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests
from analysis.v5.phenotype_redistribution import (
    ROOT, Q, TRAITS, LABELS, SEED, eligible, cellkeys, rows_for, weighted_quantile,
    encode, cell_average, distributions, climate_profile
)
from analysis.v5.novel_head_overlap import T, make_cells, match_cells
from analysis.v5.invasion_partition import head_cells
from analysis.v5.native_pca_env import background_mesh, _background_species_audit, _pca_env_2d
from analysis.v5.native_pca_env_fullcoverage import apply_original_five_CHELSA
from analysis.v4.native_introduced_opportunity import authority_source, angular_cells
from analysis.v4.native_cached_climate import VARIABLES

OUT=ROOT/"surplus"
def targets(nat,inv,u,v):
    """Exact common and directional excess mass; only observed grid bins are estimable."""
    common=np.minimum(nat,inv)
    excess=np.maximum(inv-nat,0)
    np.testing.assert_allclose(common.sum()+excess.sum(),1,atol=1e-9)
    np.testing.assert_allclose(excess.sum(),1-common.sum(),atol=1e-9)
    keys=u*nat.shape[1]+v
    counts=np.bincount(keys,minlength=nat.size)
    support=counts.reshape(nat.shape)>0
    masses=np.array([common[support].sum(),excess[support].sum()])
    weights=np.array([m[u,v]/counts[keys] for m in [common,excess]])
    np.testing.assert_allclose(weights.sum(1),masses)
    return common,excess,masses,weights


def slope_row(sp,rep,trait,category,values,score,blocks):
    """Within-region slope with clustered uncertainty, conditional on estimated densities."""
    d=pd.DataFrame(dict(y=values,x=score,block=blocks))
    d[["y","x"]]=d[["y","x"]]-d.groupby("block")[["y","x"]].transform("mean")
    informative=d.groupby("block").x.apply(lambda v:float(v@v)>1e-12)
    d=d.loc[d.block.isin(informative[informative].index)]
    if d.block.nunique()<8 or len(d)<20 or (d.x@d.x)<1e-12:
        return dict(species=sp,representation=rep,trait=trait,category=category,
                    cells=len(d),blocks=d.block.nunique(),status="insufficient")
    x=d.x.to_numpy(); y=d.y.to_numpy()
    fit=sm.OLS(y,x[:,None]).fit(cov_type="cluster",cov_kwds={"groups":d.block},use_t=True)
    beta=float(fit.params[0])
    loo=[]
    for b in d.block.unique():
        keep=(d.block.to_numpy()!=b)
        if x[keep]@x[keep]>1e-12:loo.append(float((x[keep]@y[keep])/(x[keep]@x[keep])))
    # Independent-block multiplier sign sensitivity for the zero-slope score.
    # Whole blocks, not individual cells; conditional approximation, not random assignment.
    contributions=np.array([x[d.block==b]@y[d.block==b] for b in d.block.unique()])
    k=len(contributions); exact=k<=15
    rng=np.random.default_rng(SEED)
    signs=np.array(list(product([-1,1],repeat=k))) if exact else rng.choice([-1,1],(9999,k))
    null=np.abs(signs@contributions)
    p=(np.count_nonzero(null>=abs(contributions.sum())-1e-12)+int(not exact))/(len(signs)+int(not exact))
    lo,hi=np.asarray(fit.conf_int()).ravel()
    return dict(species=sp,representation=rep,trait=trait,category=category,
                cells=len(d),blocks=k,status="estimated",beta=beta,lo=lo,hi=hi,
                p=float(fit.pvalues[0]),p_block_sign=p,loo_min=min(loo),loo_max=max(loo))


def run_species(sp,sn,si,env,transform,edges,nat,inv,outputs):
    fn,fi=eligible(sn),eligible(si)
    n=make_cells(sn,True)
    full_i=make_cells(si,True)
    expected={"Cirsium arvense":(1366,860),"Cirsium vulgare":(1113,1139)}[sp]
    assert (len(n),len(full_i))==expected
    hc=head_cells(si,transform,edges)[["i","j","u","v"]]
    c=full_i.merge(hc,on=["i","j"],validate="one_to_one")
    u,v=c.u.to_numpy(),c.v.to_numpy()
    valid=inv[u,v]>0
    c=c.loc[valid].reset_index(drop=True);u,v=c.u.to_numpy(),c.v.to_numpy()
    c["score"]=(inv[u,v]-nat[u,v])/(inv[u,v]+nat[u,v])
    c["novel"]=(c.score>0).astype(int)
    common,excess,mass,weights=targets(nat,inv,u,v)
    pairs=match_cells(c)
    outputs["cohort"].append(dict(species=sp,n_native=len(n),n_intro=len(full_i),assigned=len(c),
        positive=int(c.novel.sum()),nonpositive=int((c.novel==0).sum()),
        pairs=sum(len(a) for a,b in pairs),blocks=len(pairs),climate_D=common.sum(),
        introduced_excess_mass=excess.sum(),common_observed_mass=mass[0],
        excess_observed_mass=mass[1],common_coverage=mass[0]/common.sum(),
        excess_coverage=mass[1]/excess.sum()))
    climate_profile(sp,"density_surplus",n,c,pairs,fn,fi,outputs["climate_profiles"])
    for rep in ["median","record_equalcell"]:
        distributions(sp,"density_surplus",rep,fn,fi,n,c,pairs,{},full_i,
                      outputs["frequencies"],outputs["paired_categories"],
                      outputs["omnibus"],outputs["thresholds"],outputs["climate_profiles"])
        vn,ni,nw=rows_for(fn,n,rep)
        vi,ii,iw=rows_for(fi,c,rep)
        score=c.score.to_numpy();blocks=c.block.to_numpy()
        def summarize(feature,labels,h,ref,kind):
            # Mass targets are fixed; unsupported climate bins are not imputed.
            if np.any(mass<=0):return
            v0,v1=weights@h/mass[:,None]
            for j,label in enumerate(labels):
                outputs["mass_weighted"].append(dict(species=sp,representation=rep,feature=feature,
                    category=label,kind=kind,native_reference=ref[j],overlap_component=v0[j],
                    excess_component=v1[j],difference=v1[j]-v0[j]))
                if kind in ("mean","tails"):
                    outputs["slopes"].append(slope_row(sp,rep,feature,label,h[:,j],score,blocks))
        for j,name in enumerate(TRAITS):
            qs=np.quantile(vn[:,j],Q) if rep=="median" else weighted_quantile(vn[:,j],nw,Q)
            hn=encode(vn,ni,len(n),[qs],[j]);h=encode(vi,ii,len(c),[qs],[j])
            summarize(name,LABELS,h,hn.mean(0),"distribution")
            low,high=vn[:,j].min(),vn[:,j].max()
            def tail(x,ix,k):
                z=x[:,j]
                return np.column_stack([cell_average(a,ix,k) for a in
                    [(z>=low)&(z<qs[1]),(z<=high)&(z>=qs[-2]),z<low,z>high]])
            summarize(name,["existing_low20","existing_high20","below_observed_min","above_observed_max"],
                      tail(vi,ii,len(c)),tail(vn,ni,len(n)).mean(0),"tails")
            summarize(name,["mean_raw"],cell_average(vi[:,j],ii,len(c))[:,None],
                      np.array([np.average(vn[:,j],weights=nw)]),"mean")
            # Dispersion measured around the fixed native median; no assumption of direction.
            center=np.quantile(vn[:,j],.5) if rep=="median" else weighted_quantile(vn[:,j],nw,[.5])[0]
            summarize(name,["abs_deviation_native_median"],cell_average(abs(vi[:,j]-center),ii,len(c))[:,None],
                      np.array([np.average(abs(vn[:,j]-center),weights=nw)]),"mean")
        for cols,name in [([0,1],"angle+chroma"),([0,1,2],"angle+chroma+outline")]:
            qs=[np.quantile(vn[:,j],[1/3,2/3]) if rep=="median" else weighted_quantile(vn[:,j],nw,[1/3,2/3]) for j in cols]
            a=encode(vn,ni,len(n),qs,cols);h=encode(vi,ii,len(c),qs,cols)
            ref=a.mean(0);labs=["/".join(x) for x in product(("low","mid","high"),repeat=len(cols))]
            summarize(name,labs,h,ref,"joint")
            cls=np.where(ref==0,0,np.where(ref<=.02,1,2))
            hn=np.column_stack([a[:,cls==k].sum(1) for k in range(3)])
            hi=np.column_stack([h[:,cls==k].sum(1) for k in range(3)])
            summarize(name,["unobserved_native_combination","rare_existing_combination","common_existing_combination"],
                      hi,hn.mean(0),"joint_support")


def main():
    OUT.mkdir(parents=True,exist_ok=True)
    auth=authority_source((ROOT/"private/authority350.zip").read_bytes())
    wild=pd.read_csv(ROOT/"private/focal_wild.csv.gz",dtype={"obs_id":str},low_memory=False)
    assert len(wild)==297622
    wild["_cell"]=angular_cells(wild.latitude,wild.longitude)
    wild=wild.sort_values("obs_id").drop_duplicates(["taxon_name","range_state","_cell","user_id"]).copy()
    land=background_mesh(auth)
    wild,land,_=apply_original_five_CHELSA(wild,land,
        source_meta=Path("analysis/v4/environment_source_metadata.json"),work=ROOT/"private/rasters")
    measured=pd.read_csv(ROOT/"private/original-measured.csv.gz",usecols=["obs_id",*T],dtype={"obs_id":str})
    wild=wild.drop(columns=T,errors="ignore").merge(measured,on="obs_id",how="left",validate="one_to_one")
    frozen=json.loads(Path("analysis/v4/results/20261021_unified_abcde/S2_TWO_SPECIES_FULL_WCVP_NICHE_VERIFIED_20261021.json").read_text())
    outputs={k:[] for k in ["frequencies","paired_categories","omnibus","thresholds","climate_profiles","cohort","mass_weighted","slopes"]}
    for sp,rows in wild.groupby("taxon_name"):
        env,_=_background_species_audit(auth,land.dropna(subset=list(VARIABLES)),sp,land)
        complete=rows.dropna(subset=list(VARIABLES))
        sn=complete.loc[complete.range_state=="native"];si=complete.loc[complete.range_state=="introduced"]
        transform,edges,nat,inv=_pca_env_2d(env[env._range=="native"],env[env._range=="introduced"],sn,si,return_grid=True)
        target=frozen["species"][sp]["PCA_env_background_corrected_Schoener_D"]
        assert abs(1-abs(nat-inv).sum()/2-target)<1e-5
        run_species(sp,sn,si,env,transform,edges,nat,inv,outputs)
    for name,rows in outputs.items():
        d=pd.DataFrame(rows)
        if name in ("paired_categories","omnibus") and len(d):
            keys=["representation","family" if name=="paired_categories" else "test"]
            d["holm"]=d.groupby(keys).p.transform(lambda p:multipletests(p,method="holm")[1])
        if name=="slopes" and len(d):
            for p in ["p","p_block_sign"]:
                d["holm_"+p]=np.nan
                for _,idx in d[d.status=="estimated"].groupby("representation").groups.items():
                    d.loc[idx,"holm_"+p]=multipletests(d.loc[idx,p],method="holm")[1]
        d.to_csv(OUT/(name+".csv"),index=False)
    print("Complete; aggregate outputs only.",flush=True)

if __name__=="__main__":
    main()
