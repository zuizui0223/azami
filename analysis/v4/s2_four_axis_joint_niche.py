"""Recheck S2 in the frozen BIO18/BIO05/RSDS/elevation climate space.

PRIVATE source-row input, PUBLIC aggregate-only output. Source itself is
restored by the pinned original V4 workflow from WCVP-qualified records and
source-pixel CHELSA/WorldClim layers. The estimates are observational:
neither niche evolution nor image-angle equivalence can be inferred.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import numpy as np
import pandas as pd
import statsmodels.api as sm

ENV=("chelsa_bio18","chelsa_bio05","chelsa_rsds_mean","topo_elevation")
SPECIES=("Cirsium arvense","Cirsium vulgare")
COUNTS={("Cirsium arvense","native"):78646,
        ("Cirsium arvense","introduced"):70964,
        ("Cirsium vulgare","native"):62368,
        ("Cirsium vulgare","introduced"):85644}
BIO18_MEDIANS={("Cirsium arvense","native"):222.2,
               ("Cirsium arvense","introduced"):270.4,
               ("Cirsium vulgare","native"):218.0,
               ("Cirsium vulgare","introduced"):267.7}

def cell_id(latitude,longitude):
    a=np.asarray(latitude,float);b=np.asarray(longitude,float)
    if (not np.isfinite(a).all() or not np.isfinite(b).all()
        or (abs(a)>90).any() or (abs(b)>180).any()):
        raise ValueError("INVALID_PRIVATE_GEOREFERENCES")
    return (np.clip(np.floor((np.sin(np.deg2rad(a))+1)/.05).astype(int),0,39)*72
            +np.floor(((b+180)%360)/5).astype(int))

def restore(focal:Path, external:Path, measured:Path)->pd.DataFrame:
    p=pd.read_csv(focal,dtype={"obs_id":str},low_memory=False)
    e=pd.read_csv(external,dtype={"obs_id":str},low_memory=False)
    m=pd.read_csv(measured,dtype={"obs_id":str},low_memory=False)
    needed={"obs_id","taxon_name","range_state","user_id","latitude",
            "longitude","observed_on","chelsa_rsds_mean"}
    if not needed<=set(p) or not {"obs_id",*ENV[0:2],"topo_elevation"}<=set(e):
        raise ValueError("VERIFIED_SOURCE_COLUMNS_MISSING")
    if not(p.obs_id.is_unique and e.obs_id.is_unique and m.obs_id.is_unique):
        raise ValueError("PRIVATE_OBSERVATION_ID_DUPLICATE")
    if p.groupby(["taxon_name","range_state"]).size().to_dict()!=COUNTS:
        raise ValueError("WCVP_ORIGINAL_297622_STRATA_CHANGED")
    cols=["obs_id"]+[c for c in ENV if c not in p]
    p=p.merge(e[cols],on="obs_id",how="left",validate="one_to_one")
    angle=next((k for k in ("orientation_image_vertical_angle",
                           "orientation_angle_degrees_median") if k in m),None)
    if not angle:
        raise ValueError("ORIGINAL_CAMERA_ANGLE_NOT_PRESENT")
    p=p.merge(m[["obs_id",angle]].rename(columns={angle:"_angle"}),
              on="obs_id",how="left",validate="one_to_one")
    if len(p)!=297622 or not set(ENV)<=set(p):
        raise ValueError("S2_COARSE_SOURCE_DROPPED")
    p["_cell"]=cell_id(p.latitude,p.longitude)
    for c in (*ENV,"_angle"):
        p[c]=pd.to_numeric(p[c],errors="coerce")
    return p

def range_cell_points(t:pd.DataFrame)->pd.DataFrame:
    t=(t.sort_values("obs_id",kind="stable")
       .drop_duplicates(["user_id","_cell"]))
    if t._cell.nunique()<15:
        raise ValueError("SOURCE_CELLS_INSUFFICIENT")
    return t.groupby("_cell")[list(ENV)].mean()

def empirical_niche(t:pd.DataFrame,boot:int=399)->dict:
    answer={}; rng=np.random.default_rng(20261025)
    for sp in SPECIES:
        ranges={st:t.loc[(t.taxon_name==sp)&(t.range_state==st)] for st in
                ("native","introduced")}
        complete={st:z.loc[np.isfinite(z[list(ENV)].to_numpy(float)).all(axis=1)]
                  for st,z in ranges.items()}
        points={st:range_cell_points(z).to_numpy(float) for st,z in complete.items()}
        n,i=points["native"],points["introduced"]
        mu,mi=n.mean(0),i.mean(0)
        sd=n.std(0,ddof=1)
        if (sd<1e-8).any():
            raise ValueError("NATIVE_CLIMATE_SD_ZERO")
        diff=mi-mu; relative=diff/sd
        # IMPORTANT: raw units mix mm/month, degrees C, radiation and metres.
        # A ridge proportional to trace(raw covariance) is unit dependent.
        # Standardise the four axes by native-cell SD BEFORE covariance
        # regularisation, distance and novelty calculations.
        zn=(n-mu)/sd
        zi=(i-mu)/sd
        cov=np.cov(zn,rowvar=False)
        ridge=.05*np.trace(cov)/len(ENV)
        inv=np.linalg.inv(cov+ridge*np.eye(len(ENV)))
        distance=float(np.sqrt(relative@inv@relative))
        native_md=np.sqrt(np.einsum("ij,jk,ik->i",zn,inv,zn))
        intro_md=np.sqrt(np.einsum("ij,jk,ik->i",zi,inv,zi))
        threshold=float(np.quantile(native_md,.95))
        draws=np.empty((boot,4))
        for it in range(boot):
            x=n[rng.integers(0,len(n),len(n))].mean(0)
            y=i[rng.integers(0,len(i),len(i))].mean(0)
            draws[it]=(y-x)/sd
        axes={}
        for j,key in enumerate(ENV):
            native=ranges["native"][key].dropna().to_numpy(float)
            introduced=ranges["introduced"][key].dropna().to_numpy(float)
            if key=="chelsa_bio18":
                if (abs(np.median(native)-BIO18_MEDIANS[(sp,"native")])>.011
                    or abs(np.median(introduced)-BIO18_MEDIANS[(sp,"introduced")])>.011):
                    raise ValueError("FROZEN_BIO18_MARGINS_NOT_REPRODUCED")
            axes[key]={
                "n_native":int(len(native)),"n_introduced":int(len(introduced)),
                "native_source_median":float(np.median(native)),
                "introduced_source_median":float(np.median(introduced)),
                "introduced_minus_native_median":float(np.median(introduced)-np.median(native)),
                "native_cell_sd":float(sd[j]),
                "cell_mean_standardized_difference":float(relative[j]),
                "cell_bootstrap_95_percentile_interval":[float(v) for v in np.quantile(draws[:,j],[.025,.975])]
            }
        answer[sp]={
            "source_rows":{st:int(len(z)) for st,z in ranges.items()},
            "four_environment_complete":{st:int(len(z)) for st,z in complete.items()},
            "original_observer_deduplicated_spatial_cell_n":{st:int(len(z)) for st,z in points.items()},
            "climate_axes":axes,"joint_4D_regularized_mahalanobis_centroid_distance":distance,
            "ridge_trace_fraction":.05,
            "distance_unit_invariant_native_sd_standardized":True,
            "introduced_cell_fraction_outside_native_cell_95pct_ellipsoid":float(np.mean(intro_md>threshold)),
            "native_cell_fraction_outside_native_cell_95pct_ellipsoid":float(np.mean(native_md>threshold)),
            "not_environmental_availability_corrected_niche_expansion":True,
            "cell_resampling_intervals_not_independent_region_replicates":True,
        }
    return answer

def phenotype_common_support(t:pd.DataFrame)->dict:
    result={}
    for sp in SPECIES:
        z=t.loc[(t.taxon_name==sp)&t._angle.between(0,180)].copy()
        z=z.loc[np.isfinite(z[list(ENV)].to_numpy(float)).all(axis=1)]
        dts=pd.to_datetime(z.observed_on,errors="coerce",format="mixed")
        z=z.loc[dts.notna()].copy()
        z["_month"]=dts.loc[z.index].dt.month
        z["_year"]=dts.loc[z.index].dt.year.astype(float)
        z=z.sort_values("obs_id",kind="stable").drop_duplicates(
            ["range_state","user_id","_cell"])
        support={}
        for c in ENV:
            low=max(z.loc[z.range_state==st,c].quantile(.05) for st in ("native","introduced"))
            high=min(z.loc[z.range_state==st,c].quantile(.95) for st in ("native","introduced"))
            support[c]=[float(low),float(high)]
        if any(low>=high for low,high in support.values()):
            result[sp]={"status":"NO_FOUR_AXIS_ANALOG_SUPPORT","support":support}
            continue
        for c,(low,high) in support.items():
            z=z.loc[z[c].between(low,high)]
        n=z.groupby("range_state").size().to_dict()
        cells=z.groupby("range_state")._cell.nunique().to_dict()
        if min(n.get(st,0) for st in ("native","introduced"))<120 or min(cells.get(st,0) for st in ("native","introduced"))<15:
            result[sp]={"status":"INSUFFICIENT_PHOTOGRAPHIC_COMMON_SUPPORT",
                        "n":n,"cells":cells,"four_axis_bounds":support}
            continue
        x=z[list(ENV)].to_numpy(float)
        sd=x.std(0)
        if (sd<1e-9).any():
            raise ValueError("CLIMATE_COMMON_SUPPORT_CONSTANT")
        x=(x-x.mean(0))/sd
        intro=(z.range_state=="introduced").to_numpy(float)
        month=z._month.to_numpy(float);year=z._year.to_numpy(float)
        ys=year.std()
        X=np.column_stack([np.ones(len(z)),intro,x,intro*x[:,0],
                           np.sin(month*2*np.pi/12),np.cos(month*2*np.pi/12),
                           (year-year.mean())/ys if ys>0 else np.zeros(len(z))])
        group=z.groupby(["range_state","_cell"]).obs_id.transform("size")
        fit=sm.WLS(z._angle.to_numpy(float),X,weights=1/group.to_numpy(float)).fit(
            cov_type="cluster",cov_kwds={"groups":z._cell.to_numpy()},use_t=False)
        ci=fit.conf_int()
        result[sp]={
            "status":"OBSERVATIONAL_FOUR_AXIS_ANALOG_SUPPORT",
            "four_axis_analogue_rectangular_bounds":support,
            "observer_cell_balanced_n":{k:int(v) for k,v in n.items()},
            "occupied_cells":{k:int(v) for k,v in cells.items()},
            "median_camera_angle":{st:float(z.loc[z.range_state==st,"_angle"].median())
                                   for st in ("native","introduced")},
            "range_status_effect_at_pooled_climate_mean":float(fit.params[1]),
            "range_status_effect_95CI":[float(v) for v in ci[1]],
            "BIO18_by_introduction_interaction_per_pooled_SD":float(fit.params[6]),
            "BIO18_introduction_interaction_95CI":[float(v) for v in ci[6]],
            "not_an_equivalence_test":True,"not_genetic_evolution":True,
            "camera_roll_unmeasured":True,
        }
    return result

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--focal",type=Path,required=True)
    p.add_argument("--environment-ext",type=Path,required=True)
    p.add_argument("--s2-measurements",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    args=p.parse_args()
    t=restore(args.focal,args.environment_ext,args.s2_measurements)
    result={
        "study":"V4_S2_FOUR_AXIS_NATURAL_CONTRAST_20261025",
        "status":"SOURCE_297622_FOUR_ENVIRONMENT_OBSERVATIONAL_RECHECK",
        "climatic_occupancy":empirical_niche(t),
        "image_angle_matched_climate_support":phenotype_common_support(t),
        "source_private_identifiers_never_exported":True,
        "natural_introduction_not_randomized":True,
        "true_niche_evolution_not_identified_without_available_climate_background":True,
        "source_photo_angle_not_gravity_validated":True}
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(result,ensure_ascii=False,indent=2,allow_nan=False)+"\n")
    print({k:v["joint_4D_regularized_mahalanobis_centroid_distance"]
           for k,v in result["climatic_occupancy"].items()})
if __name__=="__main__":
    main()
