"""Temporal support audit before testing any trait changes.

No trait effect estimates, individual records, dates by location, or coordinates
are exported. Fixed equal-length time windows, local cell+month repeated
sampling, and observer/year support. Current climate assignment is not a date
of colonization. Source phenotype and climate filters match original V5.
"""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from analysis.v5.phenotype_redistribution import ROOT
from analysis.v5.phenotype_endpoint_scan import load_measures
from analysis.v5.native_pca_env import background_mesh, _background_species_audit, _pca_env_2d, _probability_content
from analysis.v5.native_pca_env_fullcoverage import apply_original_five_CHELSA
from analysis.v4.native_introduced_opportunity import authority_source, angular_cells
from analysis.v4.native_cached_climate import VARIABLES
from analysis.v4.native_integration import FEATURES

OUT=ROOT/"phenotype_temporal_support"
WINDOWS=[("early",2011,2015),("middle",2016,2020),("late",2021,2025)]
GRID=[.1,.25,1.]
MIN_COUNTS=[3,5]
NOVEL_GROUPS=["analogue","novel","mixed"]
ENDPOINTS={"all_source":[],**FEATURES}

def annotate_time(f):
    f=f.copy()
    dt=pd.to_datetime(f.observed_on,format="mixed",errors="coerce",utc=True)
    f["year"]=dt.dt.year;f["month"]=dt.dt.month
    return f

def grid_keys(f,size):
    i=np.minimum(round(180/size)-1,np.floor((f.latitude.to_numpy()+90)/size)).astype(int)
    j=np.floor(((f.longitude.to_numpy()+180)%360)/size).astype(int)
    block=np.floor(i*size/20).astype(int)*18+np.floor(j*size/20).astype(int)
    return i,j,block

def temporal_groups(f,size):
    d=f.copy()
    d["i"],d["j"],d["block"]=grid_keys(d,size)
    d["period"]=pd.Series(index=d.index,dtype="object")
    for key,lo,hi in WINDOWS:d.loc[d.year.between(lo,hi),"period"]=key
    d=d.dropna(subset=["period","month","year"])
    # Deduplicate per period within a local cell/month/year, not across years.
    d=d.sort_values("obs_id").drop_duplicates(["i","j","year","month","user_id"])
    keys=["i","j","month"]
    # One climate assignment for each location/month stratum, across periods,
    # based on all eligible source records, not on any phenotype.
    g=d.groupby(keys).novel.agg(["min","max"])
    g["climate_group"]=np.where(g["min"]==g["max"],np.where(g["min"]==1,"novel","analogue"),"mixed")
    d=d.merge(g[["climate_group"]],left_on=keys,right_index=True,validate="many_to_one")
    return d

def supported_strata(d,period_a,period_b,min_count):
    keys=["i","j","month","block","climate_group"]
    t=d[d.period.isin([period_a,period_b])]
    stats=t.groupby(keys+["period"]).agg(n=("obs_id","size"),users=("user_id","nunique"),years=("year","nunique"))
    qualified=stats[(stats.n>=min_count)&(stats.users>=3)&(stats.years>=2)].reset_index()
    if qualified.empty:return qualified.drop(columns=["period","n","users","years"])
    paired=qualified.groupby(keys).period.nunique()
    return paired[paired==2].reset_index()[keys]

def count_groups(sp,side,f,out):
    f=annotate_time(f)
    for name,cols in ENDPOINTS.items():
        scan=["scan_"+c for c in cols]
        if scan:
            valid=np.isfinite(f[scan].to_numpy(float)).all(1)
            if name=="hue":valid &= np.linalg.norm(f[scan].to_numpy(float),axis=1)>1e-8
        else:valid=np.ones(len(f),bool)
        g=f.loc[valid]
        for period,lo,hi in WINDOWS:
            p=g.loc[g.year.between(lo,hi)]
            out["date_coverage"].append(dict(species=sp,range_state=side,trait=name,period=period,
                       observations=len(p),users=p.user_id.nunique(),observed_years=p.year.nunique()))
    for size in GRID:
        d=temporal_groups(f,size)
        for name,cols in ENDPOINTS.items():
            scan=["scan_"+c for c in cols]
            if scan:
                good=np.isfinite(d[scan].to_numpy(float)).all(1)
                if name=="hue":good &= np.linalg.norm(d[scan].to_numpy(float),axis=1)>1e-8
                g=d.loc[good]
            else:g=d
            for pa,pb in [("middle","late"),("early","late")]:
                for min_count in MIN_COUNTS:
                    t=supported_strata(g,pa,pb,min_count)
                    for climate_group in [*NOVEL_GROUPS,"all"]:
                        q=t if climate_group=="all" else t[t.climate_group==climate_group]
                        nc=len(q.drop_duplicates(["i","j"]));nb=q.block.nunique()
                        out["repeated_support"].append(dict(species=sp,range_state=side,trait=name,
                             cell_degrees=size,period_a=pa,period_b=pb,min_records_each_period=min_count,
                             climate_group=climate_group,cell_month_strata=len(q),local_cells=nc,blocks=nb,
                             eligible=(nc>=30 and nb>=8),
                             primary=(size==.25 and pa=="middle" and min_count==3)))

def assign_climate(frame,transform,edges,nat):
    pc=transform(frame)
    uv=np.column_stack([np.searchsorted(e,pc[:,k],side="right")-1 for k,e in enumerate(edges)])
    for k,e in enumerate(edges):uv[pc[:,k]==e[-1],k]=len(e)-2
    inside=((uv>=0)&(uv<np.array([len(e)-1 for e in edges]))).all(1)
    d=frame.loc[inside].copy()
    d["novel"]=(~_probability_content(nat,.95)[uv[inside,0],uv[inside,1]]).astype(int)
    return d

def main():
    OUT.mkdir(parents=True,exist_ok=True)
    auth=authority_source((ROOT/"private/authority350.zip").read_bytes())
    wild=pd.read_csv(ROOT/"private/focal_wild.csv.gz",dtype={"obs_id":str},low_memory=False)
    assert len(wild)==297622 and wild.obs_id.is_unique
    land=background_mesh(auth)
    wild,land,_=apply_original_five_CHELSA(wild,land,source_meta=Path("analysis/v4/environment_source_metadata.json"),work=ROOT/"private/rasters")
    measured,_=load_measures(ROOT/"private/original-measured.csv.gz")
    wild=wild.drop(columns=[c for c in measured if c!="obs_id"],errors="ignore").merge(measured,on="obs_id",how="left",validate="one_to_one")
    wild["_cell"]=angular_cells(wild.latitude,wild.longitude)
    frozen=json.loads(Path("analysis/v4/results/20261021_unified_abcde/S2_TWO_SPECIES_FULL_WCVP_NICHE_VERIFIED_20261021.json").read_text())
    out={k:[] for k in ["date_coverage","repeated_support"]}
    for sp,rows in wild.groupby("taxon_name"):
        env,_=_background_species_audit(auth,land.dropna(subset=list(VARIABLES)),sp,land)
        complete=rows.dropna(subset=list(VARIABLES))
        climate=rows.sort_values("obs_id").drop_duplicates(["taxon_name","range_state","_cell","user_id"]).dropna(subset=list(VARIABLES))
        sn=climate.loc[climate.range_state=="native"];si=climate.loc[climate.range_state=="introduced"]
        transform,edges,nat,inv=_pca_env_2d(env[env._range=="native"],env[env._range=="introduced"],sn,si,return_grid=True)
        assert abs(1-abs(nat-inv).sum()/2-frozen["species"][sp]["PCA_env_background_corrected_Schoener_D"])<1e-5
        complete=assign_climate(complete,transform,edges,nat)
        for side,f in complete.groupby("range_state"):
            count_groups(sp,side,f,out)
    for name,rows in out.items():
        d=pd.DataFrame(rows)
        d.to_csv(OUT/(name+".csv"),index=False)
    print("Complete; temporal support counts only.",flush=True)

if __name__=="__main__":
    main()
