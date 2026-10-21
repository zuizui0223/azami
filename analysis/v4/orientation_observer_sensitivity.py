"""S1 orientation rainfall: author-aware robustness + nested observer variance.

A) Compare the existing fixed five-CHELSA-covariate species-REML fit with
source-name-collision-excluded, one-per-taxon×cell×observer, and original
observer-effort top-1% deletion. B) Estimate residual variance attributable
to repeated observers *nested within species×5° cell* by unbalanced ANOVA
and cell-cluster resampling. The nested variance includes repeated individuals,
microhabitats, and camera style: it is NOT pure photographer causality,
same-head photo measurement error, or a crossed random-effect REML.

Private decoded native observations plus the one-to-one authenticated original
human-user join never leave Actions. Aggregate outcomes alone are published.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd

from analysis.v4.native_gates_ae import authority_mapping, source_and_observer
from analysis.v4.native_spatial import cells
from analysis.v4.native_trait_models import prepare
from analysis.v4.observation_models import fit_module
from analysis.v4.run_native_abcde import EXPECTED_SHA, FIXED_PREDICTORS

ENDPOINTS={
    "orientation":"orientation_image_vertical_angle",
    "colour_chroma":"corolla_lab_chroma",
    "outline_aspect":"capitulum_outline_aspect_ratio",
}
SOURCE_N=313350
SOURCE_TAXA=306
PRIMARY=("orientation_image_vertical_angle","chelsa_bio12__between")
BOOTSTRAP_SEED=20261021


def original_observer_cohort(frame:pd.DataFrame,observers:pd.DataFrame,mapping:list[dict])->pd.DataFrame:
    frame, joined=source_and_observer(frame,observers)
    if not joined:
        raise ValueError("Original observer input required")
    if len(frame)!=SOURCE_N or frame.accepted_key.nunique()!=SOURCE_TAXA:
        raise ValueError("Historical native population changed")
    mixed={str(r["accepted_key"]) for r in mapping if len(set(r["source_names"]))>1}
    if len(mixed)!=7:
        raise ValueError("Exact seven original taxon concepts changed")
    frame["_mixed_accepted_key"]=frame.accepted_key.astype(str).isin(mixed)
    frame["_cell"]=cells(frame.latitude,frame.longitude)
    return frame


def finite_endpoint(frame:pd.DataFrame,name:str)->pd.Series:
    pred=np.isfinite(frame[FIXED_PREDICTORS]).all(axis=1)
    response=pd.to_numeric(frame[name],errors="coerce")
    return pred &frame[name+"__eligible"].eq(True)&np.isfinite(response)


def residual_from_fixed_species_fit(frame:pd.DataFrame,component:dict,endpoint:str,
                                    selection:dict)->pd.Series:
    if component["status"]!="estimated":
        raise ValueError("Reference model is not fitted")
    predictors=selection["selected"]
    scales=selection["scales"]
    sample=prepare(frame,predictors)
    prediction=np.full(len(sample),component["reference_prediction"],float)
    for variable in predictors:
        for scale in ("__between","__within"):
            term=variable+scale
            expected=scales[variable]["mean"] if scale=="__between" else 0.
            prediction+=(
                (sample[term].to_numpy(float)-expected)/scales[variable]["sd"]
                *component["effects"][term]["raw_units_per_environment_sd"]
            )
    outcome=pd.to_numeric(sample[endpoint],errors="coerce")
    return pd.Series(outcome.to_numpy(float)-prediction,index=frame.index)


def unbalanced_observer_component(valid:pd.DataFrame, *, replicates:int=300)->dict:
    """MoM one-way observer variance inside species×cell, not RE/causal effect.

    For eligible strata containing >=2 observers and at least 2 within-observer
    residual degrees of freedom, MS_observer/MS_residual and the unbalanced
    effective n identify a within-cell observer-cluster component. Resampling
    5° cells, not individual repeated photographs, gives a sensitivity range.
    """
    required={"accepted_key","_cell","user_id","_residual"}
    if not required<=set(valid):
        raise ValueError("Missing original source observer/cell residual")
    if valid.user_id.isna().any() or not np.isfinite(valid._residual).all():
        raise ValueError("Observer variance input not finite and complete")
    keys=["accepted_key","_cell"]
    grouped=valid.groupby(keys+["user_id"],sort=True)["_residual"].agg(
        n="size",total="sum",sum_sq=lambda x:float((x.to_numpy(float)**2).sum())
    ).reset_index()
    grouped["group_mean"]=grouped.total/grouped.n
    blocks=grouped.groupby(keys,sort=True).agg(
        N=("n","sum"),G=("n","size"),overall=("total","sum"),
        sum_n2=("n",lambda x:float(np.square(x.to_numpy(float)).sum()))
    )
    grouped=grouped.join((blocks.overall/blocks.N).rename("cell_mean"),on=keys)
    grouped["ss_between"]=grouped.n*np.square(grouped.group_mean-grouped.cell_mean)
    grouped["ss_within"]=grouped.sum_sq-np.square(grouped.total)/grouped.n
    cellstats=grouped.groupby(keys,sort=True).agg(
        ss_b=("ss_between","sum"),
        ss_w=("ss_within","sum")
    ).join(blocks)
    cellstats["df_b"]=cellstats.G-1
    cellstats["df_w"]=cellstats.N-cellstats.G
    cellstats["effective_numer"]=cellstats.N-cellstats.sum_n2/cellstats.N
    retained=cellstats.loc[
        (cellstats.df_b>=1)&(cellstats.df_w>=2)
        &(cellstats.effective_numer>0)
    ].copy()
    if len(retained)<15:
        return {"status":"insufficient_multiauthor_source_cells",
                "taxon_cell_strata_with_repeat_and_two_authors":len(retained),
                "primary_reportable_random_camera_variance":False}
    arr=retained[["ss_b","ss_w","df_b","df_w","effective_numer"]].to_numpy(float)
    def estimate(block):
        s=block.sum(axis=0)
        if s[2]<=0 or s[3]<=0 or s[4]<=0:
            return None
        msb=s[0]/s[2]
        msw=max(0.,s[1]/s[3])
        n0=s[4]/s[2]
        obs=max(0.,(msb-msw)/n0)
        frac=obs/(obs+msw) if obs+msw>0 else None
        return {"observer_nested_variance":float(obs),
                "within_observer_variance":float(msw),
                "observer_nested_fraction":float(frac) if frac is not None else None,
                "boundary_zero":bool(msb<=msw),
                "effective_n":float(n0)}
    point=estimate(arr)
    rng=np.random.default_rng(BOOTSTRAP_SEED)
    fractions=[]
    for _ in range(replicates):
        draw=estimate(arr[rng.integers(0,len(arr),size=len(arr))])
        if draw is not None and draw["observer_nested_fraction"] is not None:
            fractions.append(draw["observer_nested_fraction"])
    if len(fractions)<.95*replicates:
        raise ValueError("Too many invalid variance re-samples")
    q=np.quantile(fractions,[.025,.5,.975]).tolist()
    return {
        "status":"estimated_nested_observer_taxon_cell_residual_variance_not_pure_camera_effect",
        "taxon_cell_strata_with_repeat_and_two_authors":int(len(retained)),
        "taxa_supported":int(retained.reset_index().accepted_key.nunique()),
        "observation_rows_in_valid_cells":int(retained.N.sum()),
        "within_observer_extra_residual_df":int(retained.df_w.sum()),
        "unbalanced_MoM":point,
        "block_resampling_replicates":replicates,
        "source_cell_bootstrap_fraction_quantiles":[float(v) for v in q],
        "scope":"Residual observer-within-species-cell component after frozen climate fixed effects. Observer may mean camera, individual, site and time; not pure camera measurement error or true crossed observer REML."
    }


def effect_line(component:dict,scales:dict)->dict:
    if component.get("status")!="estimated":
        raise ValueError("Orientation model did not fit")
    effect=component["effects"]["chelsa_bio12__between"]
    unit=100/scales["chelsa_bio12"]["sd"]
    within=component["effects"]["chelsa_bio12__within"]
    return {
        "observations":int(component["n"]),"accepted_taxa":int(component["n_species"]),
        "between_rain_deg_per_100mm":float(effect["raw_units_per_environment_sd"]*unit),
        "between_species_pointwise_95_interval_deg_per_100mm":[float(v*unit) for v in effect["raw_ci95"]],
        "within_rain_deg_per_100mm":float(within["raw_units_per_environment_sd"]*unit),
        "within_species_pointwise_95_interval_deg_per_100mm":[float(v*unit) for v in within["raw_ci95"]],
        "species_random_intercept_variance_ratio":float(component["variance_ratio"]),
        "caveat":"species random-intercept five-covariate observational fit; source/observer selection sensitivity, not causal rain effect"
    }


def observer_sensitivity(frame:pd.DataFrame,model:dict,mapping:list[dict],
                         *,replicates:int=300)->dict:
    if model.get("input_sha256")!=EXPECTED_SHA:
        raise ValueError("Wrong corrected native source model")
    selection=model["selection"]
    if selection["selected"]!=FIXED_PREDICTORS:
        raise ValueError("Five-climate model changed")
    scales=selection["scales"]
    prep=prepare(frame,FIXED_PREDICTORS)
    name=ENDPOINTS["orientation"]
    included=finite_endpoint(prep,name)
    if included.sum()<100000:
        raise ValueError("Native primary orientation source unexpectedly small")
    original=model["models"]["orientation"]["components"][name]
    if (int(original["n"]),int(original["n_species"]))!=(170066,244):
        raise ValueError("Already-published base estimate cohort changed")
    reference=effect_line(original,scales)
    if len(prep)!=313350 or prep.accepted_key.nunique()!=306:
        raise ValueError("Native frame changed")
    mask=~prep._mixed_accepted_key
    no_collision=prep.loc[mask].copy()
    if no_collision.accepted_key.nunique()!=299:
        raise ValueError("WCVP seven collision keys not excluded")
    ranked=frame.user_id.astype(str).value_counts()
    n_top=max(1,int(math.ceil(.01*len(ranked))))
    top=set(ranked.index[:n_top])
    outcome=prep.loc[included&mask].copy()
    if len(outcome)<20000:
        raise ValueError("Source concept-clean orientation sample too small")
    models={}
    candidates={
        "exclude_mixed_WCVP_keys":outcome,
        "one_observation_per_taxon_cell_original_observer":(
            outcome.sort_values("obs_id").drop_duplicates(
                ["accepted_key","_cell","user_id"],keep="first"
            )),
        "remove_top_one_percent_source_observers":outcome.loc[
            ~outcome.user_id.astype(str).isin(top)],
    }
    for scenario,sample in candidates.items():
        fit=fit_module(sample,[name],scales,FIXED_PREDICTORS,
                       required_environment=FIXED_PREDICTORS)
        if fit["status"]!="estimated":
            raise ValueError("Fitted observer sensitivity not identifiable: "+scenario)
        models[scenario]=effect_line(fit["components"][name],scales)
        models[scenario]["same_selected_five_predictors"]=True
        models[scenario]["observer_treatment"]=scenario
    # Directly assess residual repetition *after* species×cell conditioning,
    # then estimate the nested observer variance floor for 3 distinct endpoint
    # families. The style component remains confounded with plant/site.
    nested={}
    for title,endpoint in ENDPOINTS.items():
        e=finite_endpoint(prep,endpoint)&mask
        d=prep.loc[e].copy()
        component_name={
            "orientation":"orientation",
            "colour_chroma":"colour",
            "outline_aspect":"outline"
        }[title]
        comp=model["models"][component_name]["components"][endpoint]
        d["_residual"]=residual_from_fixed_species_fit(d,comp,endpoint,selection)
        nested[title]=unbalanced_observer_component(d,replicates=replicates)
    photo_counts=pd.to_numeric(outcome.source_photo_count,errors="coerce")
    multi=int(photo_counts.ge(2).sum())
    return {
        "status":"observer_balanced_five_covariate_sensitivity_estimated_no_same_head_validation",
        "source_sha256":EXPECTED_SHA,
        "reference":reference,
        "scenario_effects":models,
        "source_original_observers":int(prep.user_id.nunique()),
        "original_top_one_percent_observers_removed":n_top,
        "original_user_influence_weighted_by_full_pre_detector_native_roster":True,
        "source_photo_repeat_opportunity":{
            "orientation_rows_with_multiple_original_photo_links":multi,
            "total_orientation_rows_excluding_concept_collisions":len(outcome),
            "distinct_same_head_repeated_measurements_confirmed":0,
            "classification":"multiple photos per iNaturalist observation, NOT proof of repeated same head or independent geometry"
        },
        "nested_observer_within_species_cell_MoM":nested,
        "limitations":[
            "Balanced source/taxon-cell-observer one-per-stratum comparison is NOT an observer random-intercept likelihood fit.",
            "Nested repeated-observer residual variance includes different plants, individual sites and microclimate within 5-degree cells.",
            "Source photo count != paired images of exactly same capitulum; source database retains observation-level aggregated angles.",
            "Unmeasured value-dependent missingness and camera-pose bias remain possible.",
            "Equal-species/observer deletion and five-exposure conditional fits have different effective samples and conditional estimands.",
            "Among-species rainfall gradient does not identify rain shielding, evolutionary transitions or local adaptation.",
        ]
    }


def execute(table:Path,observers:Path,models:Path,mapping:Path,output:Path)->dict:
    sha=hashlib.sha256(table.read_bytes()).hexdigest()
    if sha!=EXPECTED_SHA:
        raise ValueError("Wrong original corrected native input bytes")
    imported=json.loads(models.read_text(encoding="utf-8"))
    f=pd.read_csv(table,dtype={"obs_id":str,"accepted_key":str},low_memory=False)
    o=pd.read_csv(observers,dtype={"obs_id":str,"user_id":str},low_memory=False)
    mapped=authority_mapping(mapping)
    f=original_observer_cohort(f,o,mapped)
    result=observer_sensitivity(f,imported,mapped)
    result["model_json_sha256"]=hashlib.sha256(models.read_bytes()).hexdigest()
    output.mkdir(parents=True,exist_ok=True)
    (output/"observer_rain_robustness.json").write_text(
        json.dumps(result,ensure_ascii=False,indent=2,allow_nan=False)+"\n")
    print(json.dumps({
        "status":result["status"],"reference":result["reference"]["between_rain_deg_per_100mm"],
        "by_original_observer_sensitivity":{
            k:v["between_rain_deg_per_100mm"] for k,v in result["scenario_effects"].items()},
        "same_observation_multiple_photos":result["source_photo_repeat_opportunity"]["orientation_rows_with_multiple_original_photo_links"],
    }))
    return result


def main()->None:
    p=argparse.ArgumentParser(description=__doc__)
    for name in ("table","observers","models","mapping","output"):
        p.add_argument("--"+name,required=True,type=Path)
    a=p.parse_args()
    execute(a.table,a.observers,a.models,a.mapping,a.output)


if __name__=="__main__":
    main()
