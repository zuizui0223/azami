"""Test whether seven WCVP name mergers alter all-S1 BIO12 image-angle gradient.

The S1 base contains 244 WCVP keys and 170066 observations. Seven keys
each collapse two different recorded names. The source-verified 251
(accepted_key, source_name) observations give a strictly reversible aggregate:
no photographic observations are generated or discarded.

This contrasts original-key versus source-name operational units at the same
five fixed CHELSA variables, with equal unit weighting, WCVP-key cluster
robust standard errors and paired WCVP-key bootstrap for the full population.
It does NOT turn synonyms into accepted biological species, use a nuclear
tree or infer adaptive nodding.
"""
from __future__ import annotations
import argparse
from collections import defaultdict
import json
from pathlib import Path

import numpy as np
import pandas as pd

from analysis.v4.joint_s1_t2_phylogeny import EXPOSURES

N_PHOTOS=170066
N_KEYS=244
N_UNITS=251
SEED=20261022
BOOT=999


def read_units(path):
    df=pd.read_csv(path,dtype={"accepted_key":str})
    columns={"accepted_key","original_source_name","original_photo_count",
             "mean_image_angle_deg",*EXPOSURES}
    if not columns<=set(df) or len(df)!=N_UNITS:
        raise ValueError("Verified original S1 source-name census absent")
    if (df.accepted_key.nunique()!=N_KEYS or
        df.original_photo_count.sum()!=N_PHOTOS or
        df.duplicated(["accepted_key","original_source_name"]).any()):
        raise ValueError("Original 244 / 251 / 170066 identity contract drift")
    if df.groupby("accepted_key").size().gt(1).sum()!=7:
        raise ValueError("Original seven WCVP mixed-key pairs changed")
    if df[list(EXPOSURES)+["mean_image_angle_deg"]].isna().any().any():
        raise ValueError("Missing original trait/climate")
    return df


def make_accepted(df):
    rows=[]
    for key,g in df.groupby("accepted_key"):
        n=g.original_photo_count.to_numpy(float)
        def weighted(col):
            return float(np.average(g[col].to_numpy(float),weights=n))
        rows.append({"accepted_key":key,"original_photo_count":int(n.sum()),
             "mean_image_angle_deg":weighted("mean_image_angle_deg"),
             **{e:weighted(e) for e in EXPOSURES}})
    result=pd.DataFrame(rows)
    if len(result)!=N_KEYS or int(result.original_photo_count.sum())!=N_PHOTOS:
        raise ValueError("WCVP accepted-key means fail reversible photo accounting")
    return result


def slope_fit(df):
    if len(df)<30:
        raise ValueError("Insufficient original taxa for fixed five-climate analysis")
    raw=df[list(EXPOSURES)].to_numpy(float)
    mean=raw.mean(axis=0); scale=raw.std(axis=0)
    if (scale<=1e-9).any() or not np.isfinite(raw).all():
        raise ValueError("Unidentified original environmental design")
    X=np.column_stack((np.ones(len(raw)),(raw-mean)/scale))
    y=df.mean_image_angle_deg.to_numpy(float)
    if np.linalg.matrix_rank(X)!=X.shape[1]:
        raise ValueError("Nonidentifiable five-axis climate")
    beta=np.linalg.lstsq(X,y,rcond=None)[0]
    resid=y-X@beta
    inv=np.linalg.inv(X.T@X)
    groups=df.accepted_key.astype(str).to_numpy()
    uniq=np.unique(groups)
    if len(uniq)<15:
        raise ValueError("Insufficient independent accepted-key clusters")
    bread=np.zeros((len(beta),len(beta)))
    for g in uniq:
        v=X[groups==g].T@resid[groups==g]
        bread+=np.outer(v,v)
    n=len(df);p=len(beta);G=len(uniq)
    variance=inv@bread@inv*(G/(G-1))*((n-1)/(n-p))
    coef=100*beta[3]/scale[2]
    sd=100*np.sqrt(max(0.,variance[3,3]))/scale[2]
    return {"n_operational_rows":int(n),"independent_WCVP_key_clusters":int(G),
            "BIO12_image_angle_deg_per_100mm":float(coef),
            "WCVP_key_cluster_robust_95_interval":[float(coef-1.96*sd),
                                                    float(coef+1.96*sd)],
            "equal_operational_unit_weight":True,
            "exif_angle_not_gravity_posture":True,
            "spatial_and_true_nuclear_phylogeny_not_adjusted":True}


def point_slope(df):
    """Same equal-unit five-CHELSA point estimate, without expensive robust SE."""
    raw=df[list(EXPOSURES)].to_numpy(float)
    mean=raw.mean(axis=0); sd=raw.std(axis=0)
    if (sd<=1e-9).any():
        raise ValueError("Degenerate source CHELSA environmental design")
    x=np.column_stack((np.ones(len(df)),(raw-mean)/sd))
    if np.linalg.matrix_rank(x)!=len(EXPOSURES)+1:
        raise ValueError("No full rank five-environment phylo-coverage design")
    b=np.linalg.lstsq(x,df.mean_image_angle_deg.to_numpy(float),rcond=None)[0]
    return float(100*b[3]/sd[2])


def bootstrap_difference(oper,accepted,draws=BOOT,seed=SEED):
    if draws<100:
        raise ValueError("Require enough paired resamples")
    keys=sorted(set(oper.accepted_key))
    if len(keys)!=N_KEYS or set(accepted.accepted_key)!=set(keys):
        raise ValueError("Paired 244 accepted keys required for the bootstrap")
    oper_by={key:g for key,g in oper.groupby("accepted_key")}
    key_by={key:g for key,g in accepted.groupby("accepted_key")}
    rng=np.random.default_rng(seed)
    values=[]
    for _ in range(draws):
        sampled=rng.choice(keys,size=len(keys),replace=True)
        # Regression point estimators only; resampled accepted-key clusters
        # are paired. Duplicated labels are resampling draws, not new species.
        many=pd.concat([oper_by[k] for k in sampled],ignore_index=True)
        few=pd.concat([key_by[k] for k in sampled],ignore_index=True)
        try:
            diff=point_slope(many)-point_slope(few)
        except ValueError:
            continue
        values.append(diff)
    if len(values)<draws*.95:
        raise ValueError("Too many unidentifiable bootstrap draws")
    a=np.array(values)
    return {"resampling_unit":"one_original_WCVP_accepted_key",
            "paired_resamples":len(a),"seed":seed,
            "source_name_minus_accepted_key_BIO12_slope_95_percentile":
                [float(np.quantile(a,.025)),float(np.quantile(a,.975))],
            "number_of_simulated_new_species":0,
            "original_sample_inference_not_caused_by_taxonomy":True}


def compare(path):
    operational=read_units(path)
    accepted=make_accepted(operational)
    results={}
    for k in (1,5,10):
        unit=operational.loc[operational.original_photo_count>=k].copy()
        keys=accepted.loc[accepted.original_photo_count>=k].copy()
        results[f"min_images_{k}"]={
           "WCVP_original_accepted_key":slope_fit(keys),
           "operational_original_source_name":slope_fit(unit),
        }
        results[f"min_images_{k}"]["delta_operational_minus_key"]=(
            results[f"min_images_{k}"]["operational_original_source_name"]["BIO12_image_angle_deg_per_100mm"]-
            results[f"min_images_{k}"]["WCVP_original_accepted_key"]["BIO12_image_angle_deg_per_100mm"])
    boot=bootstrap_difference(operational,accepted)
    return {
       "schema":"azami_V4_S1_244_key_251_source_name_BIO12_sensitivity_v1",
       "status":"REAL_ORIGINAL_S1_7_NAME_COLLAPSES_EFFECT_ON_PHOTO_BIO12_SLOPE",
       "original_source_photos":N_PHOTOS,
       "original_WCVP_keys":N_KEYS,
       "source_name_nested_in_key_units":N_UNITS,
       "mixed_original_WCVP_keys":7,
       "fixed_original_5_env_variables":list(EXPOSURES),
       "output_measurement":"EXIF_relative_photo_head_angle_degrees",
       "nuclear_294_tip_tree_claimed":False,
       "biologically_distinct_251_species_claimed":False,
       "five_environment_regression":results,
       "paired_all_WCVP_key_bootstrap":boot,
    }


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--source-251",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    result=compare(a.source_251)
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(result,ensure_ascii=False,indent=2,
                     sort_keys=True,allow_nan=False)+"\n")
    print("T2_S1_251_VS_244_ORIGINAL_BIO12_SOURCE_COLLAPSE_SENSITIVITY="+json.dumps({
       k:{"accepted_BIO12":v["WCVP_original_accepted_key"]["BIO12_image_angle_deg_per_100mm"],
          "source_name_BIO12":v["operational_original_source_name"]["BIO12_image_angle_deg_per_100mm"],
          "difference":v["delta_operational_minus_key"],
          "n_accepted":v["WCVP_original_accepted_key"]["n_operational_rows"],
          "n_source":v["operational_original_source_name"]["n_operational_rows"]}
       for k,v in result["five_environment_regression"].items()},sort_keys=True),flush=True)
    print("T2_S1_SOURCE_COLLAPSE_PAIRED_BOOTSTRAP="+json.dumps(result["paired_all_WCVP_key_bootstrap"],sort_keys=True),flush=True)


if __name__=="__main__":
    main()
