"""Quantify S1 photograph-trait climatic coverage bias when linking to EAzami294.

Uses every one of 244 original photographed WCVP keys (and all 237 clean keys).
The 113 exact Moreyra public-run matches versus remaining 131 keys are NOT
random samples. Equal-species OLS HC3 is an observation-level ecological
robustness diagnostic, not S1's original observation-level estimator,
nor phylogenetic adaptation or event history.

Primary target: annual precipitation BIO12 interaction after fixed five
original image-eligible CHELSA species means (BIO1, BIO4, BIO12, BIO15, RSDS).
Image axis is camera-relative, not gravity-referenced botanical nodding.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np

ORIGINAL_S1=244
CLEAN=237
EXACT_NUCLEAR_RUN_CANDIDATES=113
FIVE=("chelsa_bio01","chelsa_bio04","chelsa_bio12",
      "chelsa_bio15","chelsa_rsds_mean")
RESPONSE=("S1_mean_EXIF_angle_deg","S1_median_EXIF_angle_deg")


def read_source(path,required):
    with Path(path).open(newline="",encoding="utf-8-sig") as f:
        rows=list(csv.DictReader(f))
    if not rows or not required<=set(rows[0]):
        raise ValueError("Original S1/EAzami taxon source missing columns")
    return rows


def cohort(source,bridge):
    if len(source)!=ORIGINAL_S1 or len(bridge)!=ORIGINAL_S1:
        raise ValueError("S1 full 244 taxon cohort missing")
    lookup={x["accepted_key"]:x for x in bridge}
    if len(lookup)!=ORIGINAL_S1:
        raise ValueError("Duplicate WCVP source accepted key in bridge")
    rows=[]
    for x in source:
        key=x["accepted_key"]
        b=lookup.get(key)
        if not b or b["source_names"]!=x["source_names"]:
            raise ValueError("Original S1 taxon identity mismatch against EAzami bridge")
        n=int(x["eligible_photo_observations"])
        if n!=int(b["n_original_photo"]):
            raise ValueError("EAzami bridge does not preserve S1 images")
        collided=x["source_concept_collision"].lower()=="true"
        if collided != (b["WCVP_source_alias_collision"].lower()=="true"):
            raise ValueError("WCVP source collisions inconsistent")
        matched=int(b["n_exact_run_BioSample_linked_to_EAzami294"])>0
        source_name_candidate=x["Moreyra2025_SRA_taxon_name_candidate"].lower()=="true"
        if matched!=source_name_candidate:
            raise ValueError("Moreyra 113 source matched flagged differently in bridge")
        out={
            "accepted_key":key,"source_names":x["source_names"],
            "n_images":n,"collided":collided,"Moreyra_exact_run_candidate":matched,
        }
        for var in (*FIVE,*RESPONSE):
            field=(var+"_photo_eligible_mean" if var in FIVE else var)
            out[var]=float(x[field])
            if not np.isfinite(out[var]):
                raise ValueError("Nonfinite original taxon outcome or environment")
        if not 0<=out[RESPONSE[0]]<=180 or not 0<=out[RESPONSE[1]]<=180:
            raise ValueError("Invalid original EXIF orientation image angle")
        if out["chelsa_bio12"]<=0:
            raise ValueError("Nonpositive original CHELSA annual precipitation")
        rows.append(out)
    if sum(x["n_images"] for x in rows)!=170066:
        raise ValueError("Original measured images not 170066")
    if sum(x["collided"] for x in rows)!=7:
        raise ValueError("All S1 source collisions must remain 7")
    if sum(x["Moreyra_exact_run_candidate"] for x in rows)!=EXACT_NUCLEAR_RUN_CANDIDATES:
        raise ValueError("Exactly 113 Moreyra public-run candidate WCVP keys required")
    if any(x["collided"] and x["Moreyra_exact_run_candidate"] for x in rows):
        raise ValueError("Collided source accepted concepts cannot be promoted")
    return rows


def hc3(X,y):
    """Heteroskedastic robust equal-species OLS using an explicit design."""
    n,p=X.shape
    if n<=p+2 or np.linalg.matrix_rank(X)<p or np.linalg.cond(X)>1e8:
        raise ValueError("Unidentified five-climate and source-coverage interaction")
    inv=np.linalg.inv(X.T@X)
    b=inv@(X.T@y)
    u=y-X@b
    h=np.sum((X@inv)*X,axis=1)
    if np.any(h>=1-1e-9):
        raise ValueError("Leverage one: HC3 inference undefined")
    correction=(u/(1-h))**2
    V=inv@((X.T*correction)@X)@inv
    return b,V


def run_panel(rows,response):
    n=len(rows)
    y=np.array([x[response] for x in rows],float)
    orig=np.array([[x[e] for e in FIVE] for x in rows],float)
    mean=orig.mean(axis=0);sd=orig.std(axis=0)
    if np.any(sd<=1e-10):
        raise ValueError("Five original CHELSA axes not identifiable")
    z=(orig-mean)/sd
    match=np.array([int(x["Moreyra_exact_run_candidate"]) for x in rows])
    X=np.column_stack([np.ones(n),z,match,match*z[:,2]])
    b,V=hc3(X,y)
    def linear(label,weights):
        w=np.asarray(weights,float)
        effect=float(w@b)*(100./sd[2])
        se=float(np.sqrt(max(0.,w@V@w)))*(100./sd[2])
        return {"name":label,"image_deg_per_100mm_BIO12":effect,
                "hc3_95_interval":[effect-1.96*se,effect+1.96*se]}
    unmatched=np.zeros(len(b));unmatched[3]=1
    matched=unmatched.copy();matched[-1]=1
    interaction=np.zeros(len(b));interaction[-1]=1
    strata=[]
    for m in (0,1):
        selected=[i for i,x in enumerate(rows) if int(x["Moreyra_exact_run_candidate"])==m]
        ys=y[selected]
        strata.append({
            "source_Moreyra_public_run_candidate":bool(m),
            "n_taxa":len(selected),
            "n_source_photos":int(sum(rows[i]["n_images"] for i in selected)),
            "mean_angle_deg":float(ys.mean()),
            "mean_BIO12_mm":float(orig[selected,2].mean()),
            "mean_BIO15":float(orig[selected,3].mean()),
            "mean_BIO1":float(orig[selected,0].mean()),
        })
    return {
       "response":response,
       "n_taxa":n,"source_matched_taxa":int(match.sum()),
       "source_unmatched_taxa":int(n-match.sum()),
       "equal_species_weight":True,
       "fixed_original_five_climate_axes":list(FIVE),
       "measured_BIO12_species_SD_mm":float(sd[2]),
       "source_level_contrast":strata,
       "annual_precipitation_slopes":[
            linear("no_Moreyra_SRA_candidate",unmatched),
            linear("Moreyra_SRA_source_candidate",matched),
            linear("matched_minus_unmatched_slope",interaction),
       ],
       "not_geographically_or_phylogenetically_adjusted":True,
       "S1_original_photo_observation_level_model_not_refit":True,
       "adaptation_or_independent_trait_origins_inferred":False,
    }


def analyse(original,bridge):
    source=cohort(original,bridge)
    clean=[x for x in source if not x["collided"]]
    if len(clean)!=CLEAN:
        raise ValueError("Source clean 237-WCVP corpus drifted")
    result={
        "schema":"azami_S1_all244_to_EAzami294_113_SRA_coverage_gradient_diagnostic_v2",
        "status":"SOURCE_VERIFIED_ECOLOGICAL_COHORT_REPRESENTATIVENESS_NOT_EVOLUTIONARY_CAUSATION",
        "total_original_photo_taxa":ORIGINAL_S1,
        "total_original_photos":170066,
        "WCVP_alias_collisions_excluded_in_primary":7,
        "primary_clean_taxa":CLEAN,
        "matched_SRA_name_candidate_taxa":EXACT_NUCLEAR_RUN_CANDIDATES,
        "matched_sample_runs_verified_not_accepted_WCVP_concept":True,
        "panels":{
            scope:{response:run_panel(subset,response) for response in RESPONSE}
            for scope,subset in (("all244_sensitivity",source),("clean237_primary",clean))
        },
        "measurement_depth_sensitivity":{},
        "inference_limit":{
            "source_run_name_match_is_nonrandom_phylogenetic_coverage":True,
            "original_S1_was_170066_observation_level_not_equal_species":True,
            "EXIF_photo_angle_not_gravity_relative_botanical_posture":True,
            "original_S1_between_species_effect_not_replaced_by_these_models":True,
            "OU_BM_from_independent_294_tip_tree_not_executed":True,
            "environment_and_observer_spatial_confounding_unremoved":True,
        },
    }
    for threshold in (5,10):
        retained=[x for x in clean if x["n_images"]>=threshold]
        n_matched=sum(x["Moreyra_exact_run_candidate"] for x in retained)
        n_unmatched=len(retained)-n_matched
        if min(n_matched,n_unmatched)<10 or len(retained)<25:
            result["measurement_depth_sensitivity"][f"n_ge_{threshold}"]={
                "status":"blocked_not_enough_taxa_in_both_SRA_coverage_groups",
                "n_taxa":len(retained),"n_matched":n_matched,"n_unmatched":n_unmatched,
            }
        else:
            result["measurement_depth_sensitivity"][f"n_ge_{threshold}"]={
                "n_taxa":len(retained),"n_matched":n_matched,"n_unmatched":n_unmatched,
                "mean_photo_angle":run_panel(retained,RESPONSE[0]),
                "median_photo_angle":run_panel(retained,RESPONSE[1]),
            }
    return result


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--s1-census",type=Path,required=True)
    p.add_argument("--eazami-bridge",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    a=p.parse_args()
    census=read_source(a.s1_census,{"accepted_key","source_names",
        "eligible_photo_observations","source_concept_collision",
        "Moreyra2025_SRA_taxon_name_candidate",
        *RESPONSE,*{e+"_photo_eligible_mean" for e in FIVE}})
    bridge=read_source(a.eazami_bridge,{"accepted_key","source_names",
        "n_original_photo","WCVP_source_alias_collision",
        "n_exact_run_BioSample_linked_to_EAzami294"})
    result=analyse(census,bridge)
    result["source_sha256"]={
        "original_full244":hashlib.sha256(a.s1_census.read_bytes()).hexdigest(),
        "EAzami294_sample_bridge":hashlib.sha256(a.eazami_bridge.read_bytes()).hexdigest(),
    }
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(result,ensure_ascii=False,indent=2,
                     sort_keys=True,allow_nan=False)+"\n",encoding="utf-8")
    print("T2_S1_244_113_NUCLEAR_COVERAGE_ECOLOGICAL_DIAGNOSTIC="+json.dumps({
        "n_total":result["total_original_photo_taxa"],
        "n_clean":result["primary_clean_taxa"],
        "n_Moreyra":result["matched_SRA_name_candidate_taxa"],
        "clean237_mean_angle":result["panels"]["clean237_primary"][RESPONSE[0]],
        "n_ge_5_mean_angle":result["measurement_depth_sensitivity"]["n_ge_5"],
        "n_ge_10_mean_angle":result["measurement_depth_sensitivity"]["n_ge_10"],
    },ensure_ascii=False,sort_keys=True),flush=True)


if __name__=="__main__":
    main()
