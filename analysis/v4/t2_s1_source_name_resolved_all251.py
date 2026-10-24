"""Preserve all ORIGINAL Cirsium source names hidden inside 244 WCVP accepted keys.

The original S1 sample is 170066 authentic photo-angle observations across
244 WCVP accepted keys. Seven keys pool TWO different original scientific
names each. A 244-key model can be retained, but it is not a 244 biological
species phylogeny. We reconstruct 251 distinct (accepted_key, source_name)
operational units directly from the original native table, without merging
e.g. brevicaule and irumtiense or pretending every source synonym is a
distinct accepted species. The seven original accepted keys stay intact.

This is a phenotype/climate aggregation, not a taxonomic revision or a nuclear
tree; EXIF image axis is not true gravity-referenced botanical posture.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

from analysis.v4.native_integration_audit import eligible
from analysis.v4.joint_s1_t2_phylogeny import (
    ENDPOINT, EXPOSURES, MODEL_PATH, COVERAGE_PATH, MAPPING_PATH,
    actual_native_roster,
)
from analysis.v4.s1_continuous_phylogenetic_history import SHA

EXPECTED=(170066,244,251,7)
FOCAL={"Cirsium_brevicaule","Cirsium_irumtiense"}
FOCAL_KEY="206745753"


def source_rows(frame, mapping):
    required={"accepted_key","taxon_name",ENDPOINT,ENDPOINT+"__eligible",*EXPOSURES}
    if not required<=set(frame):
        raise ValueError("The original source taxon and five-exposure fields are not present")
    mask=eligible(frame,[ENDPOINT]) & np.isfinite(frame[list(EXPOSURES)]).all(axis=1)
    valid=frame.loc[mask,["accepted_key","taxon_name",ENDPOINT,*EXPOSURES]].copy()
    if valid[ENDPOINT].isna().any():
        raise ValueError("Missing original source image angle")
    valid["accepted_key"]=valid.accepted_key.astype(str)
    valid["source_name"]=valid.taxon_name.astype(str).str.strip().str.replace(" ","_",regex=False)
    if valid.source_name.eq("").any() or valid.source_name.eq("nan").any():
        raise ValueError("Original source-name missing and cannot be assigned to a species")
    if not valid[ENDPOINT].between(0,180).all():
        raise ValueError("Original source photo angle outside EXIF-relative range")
    if not valid.accepted_key.isin(mapping).all():
        raise ValueError("Unknown WCVP accepted concept in source")
    produced=[]
    for (key,name),group in valid.groupby(["accepted_key","source_name"],sort=True):
        record=mapping[key]
        raw_names={str(x).replace(" ","_") for x in record["source_names"]}
        if name not in raw_names:
            raise ValueError("Original source name absent from its recorded WCVP concept")
        measurements=group[ENDPOINT].to_numpy(float)
        count=len(measurements)
        tree_names={str(x) for x in record["tree_matches"]}
        produced.append({
            "accepted_key":key,
            "original_source_name":name,
            "original_photo_count":count,
            "mean_image_angle_deg":float(np.mean(measurements)),
            "median_image_angle_deg":float(np.median(measurements)),
            "exact_source_name_has_frozen_placement_tip":name in tree_names,
            "frozen_tree_tip":name if name in tree_names else "",
            "WCVP_key_collides_source_names":len(raw_names)>1,
            "biological_species_taxonomic_revision_verified":False,
            "botanical_gravity_posture_inferred_from_photo_angle":False,
            **{x:float(group[x].mean()) for x in EXPOSURES},
        })
    details=pd.DataFrame(produced)
    if len(details)!=details[["accepted_key","original_source_name"]].drop_duplicates().shape[0]:
        raise ValueError("Original source operational units are not unique")
    return valid,details


def reconstruct(table: Path):
    if hashlib.sha256(table.read_bytes()).hexdigest()!=SHA:
        raise ValueError("Original native decoded S1 SHA changed")
    models=json.loads(MODEL_PATH.read_text(encoding="utf-8"))
    coverage=json.loads(COVERAGE_PATH.read_text(encoding="utf-8"))[ENDPOINT]
    if models["input_sha256"]!=SHA or tuple(models["selection"]["selected"])!=EXPOSURES:
        raise ValueError("Original five CHELSA environmental axes differ from S1")
    mapping_rows=json.loads(MAPPING_PATH.read_text(encoding="utf-8"))
    mapping={str(x["accepted_key"]):x for x in mapping_rows}
    if len(mapping)!=len(mapping_rows):
        raise ValueError("Duplicate WCVP mapping key")
    frame=pd.read_csv(table,low_memory=False)
    complete_by_key=actual_native_roster(frame,models,coverage)
    original,detail=source_rows(frame,mapping)
    actual=(len(original),original.accepted_key.nunique(),len(detail),
            sum(detail.groupby("accepted_key").size()>1))
    if actual!=EXPECTED or len(complete_by_key)!=244:
        raise ValueError(f"Original 170066/244/251/7 scope changed: {actual}")
    original_collided={key for key,group in complete_by_key.items()
                       if len(group["source_names"])>1}
    observed_collided={str(key) for key,group in detail.groupby("accepted_key")
                      if len(group)>1}
    if observed_collided!=original_collided or len(original_collided)!=7:
        raise ValueError("Collided accepted source taxon identities changed")
    for key,group in detail.groupby("accepted_key"):
        expected=complete_by_key[key]
        if int(group.original_photo_count.sum())!=int(expected["n"]):
            raise ValueError(f"Photo observations lost after source-name split for {key}")
        weighted=float(np.average(group.mean_image_angle_deg,
                                   weights=group.original_photo_count))
        if abs(weighted-expected["mean_image_angle_deg"])>1e-9:
            raise ValueError(f"Grouped photo mean changed after source-name split for {key}")
        for env in EXPOSURES:
            weighted_env=float(np.average(group[env],weights=group.original_photo_count))
            if abs(weighted_env-expected["mean_frozen_five_exposures"][env])>1e-7:
                raise ValueError(f"Original five-exposure mean was altered: {env} {key}")
    pairs=[]
    for key,part in detail.loc[detail.WCVP_key_collides_source_names].groupby("accepted_key"):
        if len(part)!=2:
            raise ValueError("Original source collision is not a two-source-name pair")
        left,right=part.sort_values("original_source_name").to_dict("records")
        pairs.append({
            "accepted_key":key,
            "source_taxon_A":left["original_source_name"],
            "source_taxon_B":right["original_source_name"],
            "photos_A":int(left["original_photo_count"]),
            "photos_B":int(right["original_photo_count"]),
            "mean_photo_angle_B_minus_A_deg":
                float(right["mean_image_angle_deg"]-left["mean_image_angle_deg"]),
            "mean_BIO12_B_minus_A_mm":
                float(right["chelsa_bio12"]-left["chelsa_bio12"]),
            "any_placement_same_source_tip":
                bool(left["exact_source_name_has_frozen_placement_tip"] or
                     right["exact_source_name_has_frozen_placement_tip"]),
            "source_pair_is_not_taxonomic_species_split_proof":True,
        })
    focal=detail.loc[(detail.accepted_key==FOCAL_KEY) &
                     detail.original_source_name.isin(FOCAL)].copy()
    if set(focal.original_source_name)!=FOCAL or len(focal)!=2:
        raise ValueError("Original distinct C. brevicaule and C. irumtiense source names lost")
    if not focal.exact_source_name_has_frozen_placement_tip.all():
        raise ValueError("Distinct source names no longer correspond to original tree placements")
    summary={
       "schema":"azami_v4_S1_original_source_name_resolved_251_units_v1",
       "status":"REAL_ORIGINAL_170066_IMAGES_SPLIT_SEVEN_MERGED_WCVP_KEYS_BY_SOURCE_NAME",
       "original_S1_observation_n":170066,
       "original_WCVP_accepted_key_count":244,
       "original_single_source_WCVP_keys":237,
       "original_multi_source_WCVP_keys":7,
       "retained_source_name_within_accepted_key_units":251,
       "number_biological_accepted_species_proven":False,
       "source_name_pairs_not_assumed_independent_species":True,
       "same_original_five_climate_axes":list(EXPOSURES),
       "camera_relative_image_angle_not_botanical_nodding":True,
       "S1_accepted_key_244_estimand_preserved":True,
       "original_seven_source_mixture_pairs":pairs,
       "C_brevicaule_C_irumtiense_distinct_source_measurements":
          focal[["original_source_name","original_photo_count",
                 "mean_image_angle_deg","median_image_angle_deg",
                 *EXPOSURES]].sort_values("original_source_name").to_dict("records"),
       "source_name_matched_original_52_placement_tip_count":
          int(detail.exact_source_name_has_frozen_placement_tip.sum()),
       "nuclear_species_phylogeny_added":False,
    }
    return detail,summary


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--native-table",type=Path,required=True)
    p.add_argument("--out-dir",type=Path,required=True)
    a=p.parse_args()
    detail,summary=reconstruct(a.native_table)
    a.out_dir.mkdir(parents=True,exist_ok=True)
    detail.to_csv(a.out_dir/"S1_ALL_251_SOURCE_NAME_BY_ACCEPTED_KEY_MEASUREMENTS.csv",
                  index=False,float_format="%.12g")
    (a.out_dir/"S1_ALL244_ACCEPTED_KEYS_SEVEN_MIXED_SOURCE_SPLITS.json").write_text(
        json.dumps(summary,indent=2,ensure_ascii=False,sort_keys=True,
                   allow_nan=False)+"\n",encoding="utf-8")
    print("T2_S1_FULL_251_SOURCE_NAME_RESOLUTION="+json.dumps(summary,
       ensure_ascii=False,sort_keys=True),flush=True)


if __name__=="__main__":
    main()
