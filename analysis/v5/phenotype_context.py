"""Coarsened occupancy mass support and exploratory matched climate contrasts."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests
from analysis.v5.phenotype_redistribution import (
    ROOT, Q, TRAITS, eligible, cellkeys, rows_for, weighted_quantile, cell_average
)
from analysis.v5.novel_head_overlap import T, make_cells, match_cells
from analysis.v5.invasion_partition import head_cells
from analysis.v5.native_pca_env import (
    background_mesh, _background_species_audit, _pca_env_2d, _probability_content
)
from analysis.v5.native_pca_env_fullcoverage import apply_original_five_CHELSA
from analysis.v4.native_introduced_opportunity import authority_source, angular_cells
from analysis.v4.native_cached_climate import VARIABLES

OUT=ROOT/"phenotype_context"


def coarsen(m,k):
    if m.shape[0]%k or m.shape[1]%k: raise ValueError("Incompatible grid")
    return m.reshape(m.shape[0]//k,k,m.shape[1]//k,k).sum((1,3))


def weights_for(m,u,v,k):
    # Preserve overlap/excess mass decomposition on original fine grid.
    a=coarsen(m,k); keys=(u//k)*a.shape[1]+v//k
    counts=np.bincount(keys,minlength=a.size)
    w=a.ravel()[keys]/counts[keys]
    assert np.isclose(w.sum(),a.ravel()[counts>0].sum())
    return w


def outcomes(v,ids,n,ref,w,rep):
    result={}
    for j,name in enumerate(TRAITS):
        x=v[:,j];r=ref[:,j]
        q=np.quantile(r,[.2,.5,.8]) if rep=="median" else weighted_quantile(r,w,[.2,.5,.8])
        lo,hi=r.min(),r.max()
        for label,y in [("mean_raw",x),("existing_low20",(x>=lo)&(x<q[0])),
                        ("existing_high20",(x>=q[2])&(x<=hi)),
                        ("abs_deviation_native_median",abs(x-q[1])),
                        ("outside_observed_range",(x<lo)|(x>hi))]:
            result[(name,label)]=cell_average(y,ids,n)
    return result


def climate_rows(frame,cells):
    d=frame.copy()
    dates=pd.to_datetime(d.observed_on,format="mixed",errors="coerce")
    d["sin_month"]=np.sin(2*np.pi*dates.dt.month/12)
    d["cos_month"]=np.cos(2*np.pi*dates.dt.month/12)
    d["year"]=dates.dt.year-2020
    # No imputation for missing dates.
    cols=list(VARIABLES)+["sin_month","cos_month","year"]
    g=d.groupby("cellkey")[cols].median()
    g["log_records"]=np.log(d.groupby("cellkey").size())
    return g.reindex(cellkeys(cells))


def contrasts(sp,definition,c,pairs,frame,n,fn,output):
    if not pairs:return
    left=np.concatenate([a for a,b in pairs]);right=np.concatenate([b for a,b in pairs])
    blocks=c.iloc[left].block.to_numpy()
    assert np.array_equal(blocks,c.iloc[right].block.to_numpy())
    cv=climate_rows(frame,c)
    native_clim=climate_rows(fn,n)[list(VARIABLES)].to_numpy()
    scales=np.std(native_clim,axis=0,ddof=0)
    raw=cv.to_numpy()[right]-cv.to_numpy()[left]
    x=raw[:,:5]/scales
    for k,axis in enumerate(VARIABLES):
        output["paired_climate"].append(dict(species=sp,definition=definition,axis=axis,
          pairs=len(left),blocks=len(np.unique(blocks)),left_median=np.median(cv.iloc[left,k]),
          right_median=np.median(cv.iloc[right,k]),mean_delta=np.mean(raw[:,k]),
          sd_native=scales[k],standardized_mean_delta=np.mean(x[:,k])))
    for rep in ["median","record_equalcell"]:
        vn,ni,nw=rows_for(fn,n,rep);vi,ii,iw=rows_for(frame,c,rep)
        ys=outcomes(vi,ii,len(c),vn,nw,rep)
        for adjustment in ["climate_only","date_count_sensitivity"]:
            X=x if adjustment=="climate_only" else np.column_stack([x,raw[:,5:]])
            valid=np.isfinite(X).all(1)
            XX=sm.add_constant(X[valid],has_constant="add")
            b=blocks[valid];groups=np.unique(b)
            if len(groups)<8 or len(b)<max(25,3*XX.shape[1]) or np.linalg.matrix_rank(XX)<XX.shape[1]:continue
            for (trait,kind),v in ys.items():
                if kind not in ("mean_raw","existing_low20","existing_high20"):continue
                y=(v[right]-v[left])[valid]
                if np.std(y)<1e-12:continue
                fit=sm.OLS(y,XX).fit(cov_type="cluster",cov_kwds={"groups":b},use_t=True)
                loo=[]
                for g in groups:
                    k=b!=g
                    if np.linalg.matrix_rank(XX[k])==XX.shape[1]:
                        loo.append(np.linalg.lstsq(XX[k],y[k],rcond=None)[0])
                loo=np.asarray(loo)
                for j,axis in enumerate(VARIABLES,1):
                    # FWL score test with a whole-block multiplier sign sensitivity.
                    nuisance=np.delete(XX,j,axis=1)
                    xr=XX[:,j]-nuisance@np.linalg.lstsq(nuisance,XX[:,j],rcond=None)[0]
                    yr=y-nuisance@np.linalg.lstsq(nuisance,y,rcond=None)[0]
                    s=np.array([xr[b==g]@yr[b==g] for g in groups])
                    from itertools import product
                    exact=len(groups)<=15
                    flips=np.array(list(product([-1,1],repeat=len(groups)))) if exact else np.random.default_rng(20261126).choice([-1,1],(9999,len(groups)))
                    p_sign=(np.count_nonzero(abs(flips@s)>=abs(s.sum())-1e-12)+int(not exact))/(len(flips)+int(not exact))
                    lo,hi=np.asarray(fit.conf_int())[j]
                    output["paired_climate_slopes"].append(dict(species=sp,definition=definition,representation=rep,
                       adjustment=adjustment,trait=trait,endpoint=kind,axis=axis,pairs=len(b),blocks=len(groups),
                       beta=fit.params[j],lo=lo,hi=hi,p=fit.pvalues[j],p_block_sign=p_sign,
                       loo_min=loo[:,j].min() if loo.size else np.nan,loo_max=loo[:,j].max() if loo.size else np.nan))


def run_species(sp,sn,si,env,transform,edges,nat,inv,output):
    fn,fi=eligible(sn),eligible(si)
    n=make_cells(sn,True);full=make_cells(si,True)
    expected={"Cirsium arvense":(1366,860,30,183),"Cirsium vulgare":(1113,1139,146,316)}[sp]
    assert (len(n),len(full))==expected[:2]
    hc=head_cells(si,transform,edges)[["i","j","u","v"]]
    c=full.merge(hc,on=["i","j"],validate="one_to_one")
    u,v=c.u.to_numpy(),c.v.to_numpy()
    c=c.loc[inv[u,v]>0].reset_index(drop=True);u,v=c.u.to_numpy(),c.v.to_numpy()
    common=np.minimum(nat,inv);excess=np.maximum(inv-nat,0)
    c["novel"]=(inv[u,v]>nat[u,v]).astype(int)
    pairs=match_cells(c)
    assert sum(len(a) for a,b in pairs)==expected[3]
    contrasts(sp,"density_surplus",c,pairs,fi,n,fn,output)
    for rep in ["median","record_equalcell"]:
        vn,ni,nw=rows_for(fn,n,rep);vi,ii,iw=rows_for(fi,c,rep)
        ys=outcomes(vi,ii,len(c),vn,nw,rep)
        for k in [1,2,4,8]:
            w=np.stack([weights_for(m,u,v,k) for m in (common,excess)])
            masses=w.sum(1)
            ess=masses**2/(w*w).sum(1)
            for (trait,kind),y in ys.items():
                means=w@y/masses
                output["coarsened_mass"].append(dict(species=sp,representation=rep,grid_factor=k,
                   trait=trait,endpoint=kind,overlap_coverage=masses[0]/common.sum(),
                   excess_coverage=masses[1]/excess.sum(),overlap_ess=ess[0],excess_ess=ess[1],
                   overlap=means[0],excess=means[1],difference=means[1]-means[0]))
    coords=transform(si)
    uv=np.column_stack([np.searchsorted(e,coords[:,k],side="right")-1 for k,e in enumerate(edges)])
    for k,e in enumerate(edges):uv[coords[:,k]==e[-1],k]=len(e)-2
    inside=np.all((uv>=0)&(uv<np.array([len(e)-1 for e in edges])),axis=1)
    labeled=si.loc[inside].copy()
    on=_probability_content(nat,.95)
    labeled["novel"]=(~on[uv[inside,0],uv[inside,1]]).astype(int)
    a=make_cells(labeled);a=a.loc[a.novel!=.5].reset_index(drop=True)
    pairs=match_cells(a)
    assert sum(len(i) for i,j in pairs)==expected[2]
    contrasts(sp,"native_core95",a,pairs,fi,n,fn,output)


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
    output={k:[] for k in ["coarsened_mass","paired_climate","paired_climate_slopes"]}
    for sp,rows in wild.groupby("taxon_name"):
        env,_=_background_species_audit(auth,land.dropna(subset=list(VARIABLES)),sp,land)
        complete=rows.dropna(subset=list(VARIABLES))
        sn=complete.loc[complete.range_state=="native"];si=complete.loc[complete.range_state=="introduced"]
        transform,edges,nat,inv=_pca_env_2d(env[env._range=="native"],env[env._range=="introduced"],sn,si,return_grid=True)
        target=frozen["species"][sp]["PCA_env_background_corrected_Schoener_D"]
        assert abs(1-abs(nat-inv).sum()/2-target)<1e-5
        run_species(sp,sn,si,env,transform,edges,nat,inv,output)
    for name,rows in output.items():
        d=pd.DataFrame(rows)
        if name=="paired_climate_slopes" and len(d):
            for p in ["p","p_block_sign"]:
                d["holm_"+p]=d.groupby(["definition","representation","adjustment"])[p].transform(lambda s:multipletests(s.fillna(1),method="holm")[1])
        d.to_csv(OUT/(name+".csv"),index=False)
    print("Complete; aggregate outputs only.",flush=True)

if __name__=="__main__":
    main()
