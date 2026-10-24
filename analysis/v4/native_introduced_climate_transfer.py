"""Exploratory native→introduced within-species rainfall-cline transfer test.

Only the two focal source species and WCVP EXPLICIT introduced territories are
compared; frozen source/WCVP/tile identities are checked. Five global
climatologies are held fixed (BIO1, BIO4, BIO12, BIO15, RSDS). Cluster-robust
WLS on distinct original observer×5° site strata estimates the within-range
conditional precipitation coefficient. No evolutionary time, plasticity,
genotype, true head tilt or rain-protection fitness can be identified.

Outputs have aggregate species-range values only; no photographer/location
point, original photo or observation ID may leave the Actions private runner.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd

from analysis.v4.native_cached_climate import (
    VARIABLES, physical_climate_for_focal, digest_file
)
from analysis.v4.native_spatial import cells

SOURCE_REPORT_SHA256="bbfc9a459a1a5c230c31a405671633e3032491f92f0123634e5cbb8f4db33562"
MIN_ORIENTATION=30
MIN_CELLS=15
MIN_OBSERVERS=10
LIMIT_CONDITION=1e10


def load_verified_focal(path:Path,report_path:Path)->tuple[pd.DataFrame,dict]:
    raw=report_path.read_bytes()
    if hashlib.sha256(raw).hexdigest()!=SOURCE_REPORT_SHA256:
        raise ValueError("Cannot mix original source census with different focal source")
    report=json.loads(raw)
    if (report["overall"]["native_primary_wild_observations"],
        report["overall"]["introduced_primary_wild_observations"])!=(313350,156750):
        raise ValueError("WCVP native/introduced status partition changed")
    data=pd.read_csv(path,dtype={"obs_id":str,"user_id":str,"accepted_key":str},
                     low_memory=False)
    if not data.obs_id.is_unique or data.obs_id.isna().any():
        raise ValueError("Original focal observation identity changed")
    if set(data.taxon_name)!={"Cirsium arvense","Cirsium vulgare"}:
        raise ValueError("Only the two predeclared original source taxa are permitted")
    if set(data.range_state)!={"native","introduced"}:
        raise ValueError("Only accepted WCVP wild native/introduced source states are permitted")
    if data.user_id.isna().any():
        raise ValueError("Original observer ID absent in focal sample")
    if data.orientation_angle_degrees_median.isna().any() or (
       ~data.orientation_angle_degrees_median.between(0,180)).any():
        raise ValueError("Original focal orientation response not validated")
    expected={
        ("Cirsium arvense","native"):49072,
        ("Cirsium arvense","introduced"):43195,
        ("Cirsium vulgare","native"):33637,
        ("Cirsium vulgare","introduced"):47833
    }
    actual={tuple(k):len(g) for k,g in data.groupby(["taxon_name","range_state"])}
    if actual!=expected:
        raise ValueError("Focal source-range orientation counts disagree with first independent audit")
    return data,report


def dated_and_balanced(frame:pd.DataFrame)->pd.DataFrame:
    s=frame.copy()
    date=pd.to_datetime(s.observed_on,errors="coerce",format="mixed")
    s["_month"]=date.dt.month
    s["_year"]=date.dt.year
    s=s.loc[s._month.between(1,12)&s._year.between(1980,2026)].copy()
    if not len(s):
        return s
    if s.user_id.isna().any():
        raise ValueError("Missing source observer after date filter")
    s["_cell"]=cells(s.latitude.to_numpy(float),s.longitude.to_numpy(float))
    s=s.sort_values("obs_id").drop_duplicates(
        ["taxon_name","range_state","_cell","user_id"],keep="first"
    ).copy()
    return s


def standardization(frame:pd.DataFrame)->dict:
    result={}
    for v in VARIABLES:
        values=frame[v].to_numpy(float)
        mean=float(values.mean())
        sd=float(values.std(ddof=0))
        if not np.isfinite(mean) or not np.isfinite(sd) or sd<=1e-9:
            raise ValueError("Species pooled climatic support was constant: "+v)
        result[v]={"mean":mean,"sd":sd}
    return result


def design(sample:pd.DataFrame,scales:dict)->np.ndarray:
    climate=[(sample[v].to_numpy(float)-scales[v]["mean"])/scales[v]["sd"]
             for v in VARIABLES]
    mo=sample._month.to_numpy(float)
    year=sample._year.to_numpy(float)
    return np.column_stack([
        np.ones(len(sample)),*climate,
        np.sin(2*np.pi*mo/12),np.cos(2*np.pi*mo/12),
        (year-2015)/10
    ])


def cell_weighted_fit(sample:pd.DataFrame,scales:dict)->dict:
    """OLS with one per cell×original observer, equal total cell weight.

    Sandwich meat resamples cells, so this is not two-way crossed camera
    random effects. Original same-head replicate validation remains open.
    """
    if len(sample)<MIN_ORIENTATION or sample._cell.nunique()<MIN_CELLS:
        return {"status":"insufficient_original_pixels_cells_or_focal_orientation"}
    if sample.user_id.nunique()<MIN_OBSERVERS:
        return {"status":"insufficient_original_observers"}
    x=design(sample,scales)
    y=sample.orientation_angle_degrees_median.to_numpy(float)
    n,p=x.shape
    groups=sample._cell.to_numpy(int)
    cell_n=sample.groupby("_cell").obs_id.transform("size").to_numpy(float)
    weight=1/cell_n
    xw=x*weight[:,None]
    bread=x.T@xw
    rank=int(np.linalg.matrix_rank(bread))
    cond=float(np.linalg.cond(bread))
    if rank<p or cond>LIMIT_CONDITION or not np.isfinite(cond):
        return {"status":"five_climate_and_month_design_not_identifiable",
                "n":n,"cells":int(np.unique(groups).size),
                "design_rank":rank,"design_columns":p,
                "weighted_design_condition_number":cond if np.isfinite(cond) else None}
    beta=np.linalg.solve(bread,x.T@(weight*y))
    predicted=x@beta
    residual=y-predicted
    group_ids,inv=np.unique(groups,return_inverse=True)
    sums=np.zeros((len(group_ids),p),float)
    np.add.at(sums,inv,x*(weight*residual)[:,None])
    meat=sums.T@sums
    inverse=np.linalg.inv(bread)
    correction=len(group_ids)/(len(group_ids)-1)*(n-1)/(n-p)
    cov=inverse@meat@inverse*correction
    ses=np.sqrt(np.maximum(np.diag(cov),0))
    d=VARIABLES.index("chelsa_bio12")+1
    unit=100/scales["chelsa_bio12"]["sd"]
    rain=float(beta[d]*unit)
    se=float(ses[d]*unit)
    result={
        "status":"estimated_five_exposure_cell_cluster_WLS",
        "n_source_observer_cell_balanced":n,
        "source_distinct_cells":int(len(group_ids)),
        "source_distinct_original_observers":int(sample.user_id.nunique()),
        "design_rank":rank,"design_columns":p,
        "weighted_design_condition_number":cond,
        "annual_precip_degrees_per_100mm":rain,
        "annual_precip_robust_se_per_100mm":se,
        "annual_precip_95pct_cell_cluster_interval":[rain-1.96*se,rain+1.96*se],
        "weighted_rmse_degrees":float(np.sqrt(np.average(residual**2,weights=weight))),
        "original_camera_reference":"EXIF image vertical; NOT gravity angle",
        "report":"Five climate original pixel exposures with observation-month sine/cosine and year; one original observation per 5 degree equal-sin cell×human observer, all cells weighted equally.",
    }
    return result


def climate_common_support(native:pd.DataFrame,introduced:pd.DataFrame)->dict:
    overlaps={}
    for name in VARIABLES:
        an=native[name].to_numpy(float)
        ai=introduced[name].to_numpy(float)
        na=np.quantile(an,[.05,.95])
        ia=np.quantile(ai,[.05,.95])
        lower=max(na[0],ia[0])
        upper=min(na[1],ia[1])
        frac_nat=float(np.mean((an>=lower)&(an<=upper))) if upper>lower else 0.
        frac_intro=float(np.mean((ai>=lower)&(ai<=upper))) if upper>lower else 0.
        overlaps[name]={
            "native_q05_q95":[float(na[0]),float(na[1])],
            "introduced_q05_q95":[float(ia[0]),float(ia[1])],
            "joint_90pct_central_overlap":[float(lower),float(upper)] if upper>lower else None,
            "native_fraction_in_overlap":frac_nat,
            "introduced_fraction_in_overlap":frac_intro,
        }
    return overlaps


def transfer_against_native_baseline(native:pd.DataFrame,introduced:pd.DataFrame,
                                     scales:dict)->dict:
    # Prediction is made on exactly the introduced sample (not its expected
    # phenotype) from a model trained on the WCVP-native range only.
    # Train with identical fixed climatic normalization and monthly context.
    if min(native._cell.nunique(),introduced._cell.nunique())<MIN_CELLS:
        return {"status":"not_evaluable_insufficient_cells"}
    x=design(native,scales)
    y=native.orientation_angle_degrees_median.to_numpy(float)
    wn=1/native.groupby("_cell").obs_id.transform("size").to_numpy(float)
    bread=x.T@(x*wn[:,None])
    if np.linalg.matrix_rank(bread)!=x.shape[1] or np.linalg.cond(bread)>LIMIT_CONDITION:
        return {"status":"not_evaluable_native_climate_design_singular"}
    beta=np.linalg.solve(bread,x.T@(wn*y))
    xi=design(introduced,scales)
    yi=introduced.orientation_angle_degrees_median.to_numpy(float)
    wi=1/introduced.groupby("_cell").obs_id.transform("size").to_numpy(float)
    intercept=float(np.average(y,weights=wn))
    pred=xi@beta
    baseline=np.full(len(xi),intercept,float)
    squared=float(np.average((yi-pred)**2,weights=wi))
    raw=float(np.average((yi-baseline)**2,weights=wi))
    local_constant=float(np.average(yi,weights=wi))
    local=float(np.average((yi-local_constant)**2,weights=wi))
    return {
        "status":"cross_range_native_trained_direct_prediction",
        "native_trained_model_to_introduced_skill_vs_native_constant":float(1-squared/raw) if raw>0 else None,
        "native_trained_model_to_introduced_skill_vs_introduced_local_mean":float(1-squared/local) if local>0 else None,
        "introduced_weighted_rmse":float(np.sqrt(squared)),
        "native_constant_weighted_rmse":float(np.sqrt(raw)),
        "introduced_local_mean_weighted_rmse":float(np.sqrt(local)),
        "climatic_overlap_gated":False,
        "scope":"Exploratory geographic transfer on same original photo-derived response, source/territory observer-balanced and region-cell weighted; errors include intercept and ecological/community/camera shifts, not evolutionary rates."
    }


def evaluate_focal(frame:pd.DataFrame,climate_audit:dict,
                   *,required_cells:int=MIN_CELLS)->dict:
    if required_cells!=MIN_CELLS:
        raise ValueError("Registered minimum 15 source cells changed")
    if set(frame.range_state)!={"native","introduced"}:
        raise ValueError("Original WCVP states are not exactly two eligible classes")
    full=frame.copy()
    complete=np.isfinite(full[list(VARIABLES)]).all(axis=1)
    full["_complete"]=complete
    values={}
    for name,group in full.groupby("taxon_name",sort=True):
        group=group.copy()
        fresh={}
        for state in ("native","introduced"):
            g=group.loc[group.range_state.eq(state)]
            c=g.loc[g._complete]
            fresh[state]={
                "source_orientation_observations":int(len(g)),
                "five_CHELSA_original_pixel_complete":int(len(c)),
                "five_CHELSA_coverage_fraction":float(len(c)/len(g)) if len(g) else None,
            }
        enough=all(fresh[state]["five_CHELSA_original_pixel_complete"]>=MIN_ORIENTATION
                   for state in ("native","introduced"))
        if not enough:
            values[name]={"status":"no_original_pixel_climate_support", "range_support":fresh}
            continue
        eligible=group.loc[group._complete].copy()
        scaled=standardization(eligible)
        sample=dated_and_balanced(eligible)
        splits={state:sample.loc[sample.range_state.eq(state)]
                for state in ("native","introduced")}
        for state in ("native","introduced"):
            fresh[state]["source_original_observer_cell_balanced_rows"]=len(splits[state])
            fresh[state]["distinct_original_observers"]=int(splits[state].user_id.nunique())
            fresh[state]["distinct_equal_sin_cells"]=int(splits[state]._cell.nunique())
        fitted={state:cell_weighted_fit(splits[state],scaled)
                for state in ("native","introduced")}
        overlap=climate_common_support(splits["native"],splits["introduced"])
        transferable=transfer_against_native_baseline(
            splits["native"],splits["introduced"],scaled)
        status="estimated_provisional_S2_native_introduced_climate_contrast"
        if any(fitted[state]["status"]!="estimated_five_exposure_cell_cluster_WLS"
               for state in ("native","introduced")):
            status="five_climate_design_not_identifiable_under_frozen_support"
        if status.startswith("estimated"):
            sn=fitted["native"]["annual_precip_robust_se_per_100mm"]
            si=fitted["introduced"]["annual_precip_robust_se_per_100mm"]
            delta=(fitted["introduced"]["annual_precip_degrees_per_100mm"]
                   -fitted["native"]["annual_precip_degrees_per_100mm"])
            sd=float(np.sqrt(sn**2+si**2))
            contrast={
                "introduced_minus_native_within_species_BIO12_degrees_per_100mm":float(delta),
                "descriptive_independent_range_cell_cluster_95pct_interval":[
                    float(delta-1.96*sd),float(delta+1.96*sd)],
                "contrasts_are_not_equivalence_tests":True,
            }
        else:
            contrast={"status":"not_estimated_due_to_nonidentifiable_range"}
        values[name]={
            "status":status,
            "range_support":fresh,
            "model_covariates":list(VARIABLES)+[
                "observation_month_sine","observation_month_cosine","observation_year"],
            "same_species_original_global_SD_for_both_ranges":{
                k:float(v["sd"]) for k,v in scaled.items()},
            "estimated_within_range_slope":fitted,
            "range_BIO12_difference":contrast,
            "joint_climatic_support":overlap,
            "original_native_trained_introduced_transfer":transferable,
            "source_ids_and_photographers_retained_only_on_private_runner":True,
        }
    return {
        "schema_version":"v4_S2_source_original_CHELSA_five_environment_slope_and_transfer_v1",
        "status":"PROVISIONAL_SOURCE_MEASUREMENT_CLIMATE_TRANSFER_NO_CAUSAL_EVOLUTION",
        "source_original_WCVP_territory_state":"native vs explicitly introduced, wild and un-obscured",
        "single_image_relative_response":"orientation_angle_degrees_median",
        "focal_species":values,
        "source_original_pixel_climate_audit":climate_audit,
        "source_calendar_window":"1981-2010 CHELSA; observation years recorded separately",
        "assertions":{
            "native_intro_are_separate_taxon_specific_fits":True,
            "same_5_original_pixel_predictors_without_cherry_picking":True,
            "observer_cell_source_deduplicated":True,
            "each_coarse_cell_equally_weighted":True,
            "spatial_cluster_uncertainty_not_phylogenetic":True,
            "same_head_angle_photographic_error_externally_verified":False,
            "evolution_rate_or_plasticity_causally_identified":False,
            "introduced_region_is_not_proof_each_plant_introduced":True,
        },
        "predeclared_failure_conditions":[
            "Missing verified original CHELSA pixel cache prevents full-area transfer",
            "Incomplete phenotypes/observers/native or introduced cohort",
            "Less than 30 observations, 15 independent coarse site cells or 10 photographers after balancing",
            "Ill-conditioned five-climate plus season design",
            "Little native-introduced climatic overlap implies extrapolation rather than common-support inference",
        ]
    }


def execute(focal_path:Path,source_receipt:Path,
            climate_root:Path,metadata:Path,output:Path):
    frame,report=load_verified_focal(focal_path,source_receipt)
    values,audit=physical_climate_for_focal(frame,climate_root,metadata)
    for k,x in values.items():
        frame[k]=x
    result=evaluate_focal(frame,audit)
    result["source_original_intro_status_json_sha256"]=SOURCE_REPORT_SHA256
    result["climate_metadata_sha256"]=digest_file(metadata)
    output.mkdir(parents=True,exist_ok=True)
    (output/"S2_native_introduced_climate_transfer.json").write_text(
        json.dumps(result,ensure_ascii=False,indent=2,allow_nan=False)+"\n",encoding="utf-8")
    print(json.dumps({
        "status":result["status"],
        "focal":{
            species:{
                "status":r["status"],
                "native_rain":r.get("estimated_within_range_slope",{}).get("native",{}).get("annual_precip_degrees_per_100mm"),
                "introduced_rain":r.get("estimated_within_range_slope",{}).get("introduced",{}).get("annual_precip_degrees_per_100mm")
            } for species,r in result["focal_species"].items()
        },
    },ensure_ascii=False),flush=True)
    return result


def main()->None:
    p=argparse.ArgumentParser(description=__doc__)
    for x in ("focal","source_receipt","climate_root","metadata","output"):
        p.add_argument("--"+x.replace("_","-"),type=Path,required=True)
    a=p.parse_args()
    execute(a.focal,a.source_receipt,a.climate_root,a.metadata,a.output)


if __name__=="__main__":
    main()
