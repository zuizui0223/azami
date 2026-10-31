"""Original-S1 observer-disjoint split-half geographic photographic QC.

Global hash partitions of the archived, AES-GCM authenticated original user_id:
no observer can contribute to both halves anywhere in the original S1 sample.
Only taxon × cell pairs with >=2 source photos from EACH observer half may be
compared; ≥8 cells and ≥40 original source photos/taxon remain required.
This is observational signal reliability, not camera-error identification.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd

from analysis.v4.native_integration import FEATURES
from analysis.v4.native_observers import original_observer_roster, EXPECTED_NATIVE_OBSERVERS
from analysis.v4.rf_photographic_split_half_reliability import (
    SOURCE_N, REPLICATES, GRAINS, GROUP, MIN_CELLS_TAXON, MIN_PHOTOS_TAXON,
    HUE_MIN_VECTOR_NORM, _spatial_photo_cohort, _normalise,
    _digest_distribution, _taxon_half_stats, sha256,
)

OBSERVER_SALT = "azami-original-S1-global-observer-disjoint-half-20261031-v1"


def _observer_side(user_ids: pd.Series, replicate: int) -> np.ndarray:
    # The same original observer always falls on the same side globally,
    # across taxa, 1° / 5° cells and all photographic endpoints.
    unique=user_ids.drop_duplicates().astype(str)
    sides={k: hashlib.blake2b(f"{OBSERVER_SALT}|{replicate}|{k}".encode(),
           digest_size=8).digest()[0] & 1 for k in unique}
    return user_ids.astype(str).map(sides).to_numpy(dtype=np.int8)


def _observer_cell_pairs(photo: pd.DataFrame, original_users: pd.DataFrame,
                         trait: str, *, rep:int) -> tuple[pd.DataFrame,dict]:
    if trait not in FEATURES or rep < 0 or rep >= REPLICATES:
        raise ValueError("OBSERVER_HALF_UNREGISTERED_TRAIT_OR_SEED")
    if (len(original_users) == 0 or not original_users.obs_id.is_unique or
        original_users.obs_id.isna().any() or original_users.user_id.isna().any()):
        raise ValueError("OBSERVER_HALF_UNAUTHENTICATED_OR_DUPLICATE_ORIGINAL_USER")
    frame=photo.merge(original_users[["obs_id","user_id"]],on="obs_id",
                      validate="many_to_one",how="left",indicator=True)
    if len(frame)!=len(photo) or not frame._merge.eq("both").all() or frame.user_id.isna().any():
        raise ValueError("ORIGINAL_PHOTO_OBSERVER_JOIN_LOST_IDS")
    frame["half"]=_observer_side(frame.user_id,rep)
    if frame.groupby("user_id").half.nunique().max()>1:
        raise ValueError("SAME_OBSERVER_ON_BOTH_EVALUATION_HALVES")
    names=FEATURES[trait]
    counts=frame.groupby(GROUP+["half"],sort=False,observed=True).size().unstack("half",fill_value=0)
    if 0 not in counts.columns or 1 not in counts.columns:
        return pd.DataFrame(columns=GROUP),{"cell_pairs_with_two_photos_per_observer_half":0,"taxa_with_repeated_cell_support":0}
    qualified=counts.loc[(counts[0]>=2)&(counts[1]>=2)].copy()
    n_before=int(len(qualified))
    n_photo=qualified.sum(axis=1).rename("photos")
    summ=n_photo.reset_index().groupby("taxon").agg(cells=("photos","size"),photos=("photos","sum"))
    keep=summ.index[(summ.cells>=MIN_CELLS_TAXON)&(summ.photos>=MIN_PHOTOS_TAXON)]
    qualified=qualified.reset_index()[GROUP]
    qualified=qualified.loc[qualified.taxon.isin(keep)]
    src=frame.merge(qualified,on=GROUP,how="inner",validate="many_to_one")
    if len(src)==0:
        return pd.DataFrame(columns=GROUP),{
            "cell_pairs_with_two_photos_per_observer_half":n_before,
            "taxa_with_repeated_cell_support":0,"qualified_cells_after_taxon_gate":0}
    mid=src.groupby(GROUP+["half"],observed=True)[names].median().unstack("half")
    result=mid.reset_index()[GROUP].copy()
    for half in (0,1):
        for n in names:
            result[f"{n}__half{half}"]=mid[(n,half)].to_numpy(float)
    hue_unstable=0
    if trait=="hue":
        retain=np.ones(len(result),dtype=bool)
        for half in (0,1):
            cols=[f"{n}__half{half}" for n in names]
            uv=result[cols].to_numpy(float)
            retain &= np.linalg.norm(uv,axis=1)>=HUE_MIN_VECTOR_NORM
            result.loc[:,cols]=_normalise(uv)
        hue_unstable=int((~retain).sum())
        result=result.loc[retain].copy()
    cells_by_taxon=result.groupby("taxon").size()
    qual_taxa=cells_by_taxon.index[cells_by_taxon>=MIN_CELLS_TAXON]
    result=result.loc[result.taxon.isin(qual_taxa)].reset_index(drop=True)
    return result,{
        "cell_pairs_with_two_photos_per_observer_half":n_before,
        "qualified_cells_after_taxon_gate":int(len(result)),
        "taxa_with_repeated_cell_support":int(result.taxon.nunique()),
        "unstable_half_hue_cells_skipped":hue_unstable,
        "globally_no_observer_is_in_both_photo_halves":True,
        "observer_user_ids_not_exported":True,
    }


def audit_observer_half(frame:pd.DataFrame,observer_roster:pd.DataFrame,
                        trait:str,*,degrees:int,reps:int=REPLICATES)->dict:
    if trait not in FEATURES or degrees not in GRAINS or reps!=REPLICATES:
        raise ValueError("UNREGISTERED_OR_CHANGED_OBSERVER_HALF_QC_DESIGN")
    photos,base=_spatial_photo_cohort(frame,trait,degrees)
    report={"construct":trait,"degrees":degrees,
            "photo_QC_initial_support":base,"n_registered_observer_splits":reps,
            "status":"source_authenticated_observer_QC_only",
            "observer_global_not_cellwise_partition":True,
            "each_observer_exclusive_to_one_half_per_split":True}
    if base["n_qualified_taxa"]<5:
        report["status"]="insufficient_original_photo_cell_support"
        return report
    estimates=[]; diagnostics=[]
    for rep in range(reps):
        pairs,support=_observer_cell_pairs(photos,observer_roster,trait,rep=rep)
        diagnostics.append(support)
        if support["taxa_with_repeated_cell_support"]<5:
            continue
        stats,skill,null_skill,macro,corr=_taxon_half_stats(
           pairs,trait,shuffle_seed=20261101+rep)
        if not stats or not np.isfinite([skill,null_skill,macro,corr]).all():
            continue
        estimates.append({"n_taxa":len(stats),"n_cells":len(pairs),
                          "skill":skill,"correlation":corr,
                          "permutation_gap":skill-null_skill})
    report["n_successful_observer_disjoint_splits"]=len(estimates)
    report["n_taxa_with_two_observer_halves_by_split"]=_digest_distribution(
        [x.get("taxa_with_repeated_cell_support",0) for x in diagnostics])
    report["n_cells_with_two_observer_halves_by_split"]=_digest_distribution(
        [x.get("qualified_cells_after_taxon_gate",0) for x in diagnostics])
    if not estimates:
        report["status"]="insufficient_geographical_observer_independent_replication"
        report["no_skill_claim_from_insufficient_support"]=True
        return report
    report.update({
        "status":"estimated_observer_disjoint_spatial_photo_QC_not_RF",
        "n_taxa_contributing_to_skill_by_split":_digest_distribution([x["n_taxa"] for x in estimates]),
        "n_cells_contributing_to_skill_by_split":_digest_distribution([x["n_cells"] for x in estimates]),
        "cross_observer_half_skill_vs_train_observer_half_taxon_median":_digest_distribution([x["skill"] for x in estimates]),
        "cross_observer_half_taxon_centered_cell_correlation":_digest_distribution([x["correlation"] for x in estimates]),
        "gain_over_within_taxon_shuffled_cell_correspondence":_digest_distribution([x["permutation_gap"] for x in estimates]),
        "n_splits_positive_skill":int(sum(x["skill"]>0 for x in estimates)),
        "n_splits_positive_within_taxon_median_correlation":int(sum(x["correlation"]>0 for x in estimates)),
        "random_observer_splits_not_replicate_biological_individuals_or_confidence_intervals":True,
        "not_an_RF_or_camera_error_estimate_or_adaptive_causal_result":True,
    })
    return report


def run(native:Path,receipt_path:Path,sealed_metadata:Path,output:Path,
        *,expected_n:int=SOURCE_N,constructs:tuple[str,...]|None=None,
        grains:tuple[int,...]=GRAINS)->dict:
    receipt=json.loads(Path(receipt_path).read_text(encoding="utf-8"))
    if (sha256(native)!=receipt.get("prepared_table_sha256") or
        receipt.get("native_observations")!=expected_n or
        receipt.get("environment_decoding_status")!="physical_units_gsp_invalid_encoding_masked_v1"):
        raise ValueError("OBSERVER_HALF_ORIGINAL_S1_SOURCE_AND_DECODE_NOT_AUTHENTICATED")
    keyhex=os.environ.get("AZAMI_METADATA_KEY","")
    if len(keyhex)!=64:
        raise ValueError("ORIGINAL_OBSERVER_METADATA_AES_KEY_NOT_AVAILABLE")
    choices=tuple(FEATURES) if constructs is None else tuple(constructs)
    if not choices or len(choices)!=len(set(choices)) or not set(choices)<=set(FEATURES):
        raise ValueError("OBSERVER_HALF_UNKNOWN_OR_DUPLICATE_CONSTRUCTS")
    if not grains or len(grains)!=len(set(grains)) or not set(grains)<=set(GRAINS):
        raise ValueError("OBSERVER_HALF_UNKNOWN_OR_DUPLICATE_SPATIAL_GRAINS")
    cols=["obs_id","accepted_key","latitude","longitude"]
    for construct in choices:
        for name in FEATURES[construct]:cols.extend((name,name+"__eligible"))
    photos=pd.read_csv(native,usecols=list(dict.fromkeys(cols)),
                       dtype={"obs_id":str},low_memory=False)
    if len(photos)!=expected_n or photos.obs_id.isna().any() or not photos.obs_id.is_unique:
        raise ValueError("OBSERVER_HALF_ORIGINAL_SOURCE_GRAIN_DRIFT")
    roster=original_observer_roster(Path(sealed_metadata).read_bytes(),
              bytes.fromhex(keyhex),photos.obs_id)
    if (len(roster)!=expected_n or roster.user_id.isna().any() or
        roster.user_id.nunique()!=EXPECTED_NATIVE_OBSERVERS):
        raise ValueError("ORIGINAL_OBSERVER_SOURCE_POPULATION_CHANGED")
    results=[audit_observer_half(photos,roster,c,degrees=d)
             for c in choices for d in grains]
    report={
        "schema":"azami_original_S1_disjoint_observer_8photo_traits_1deg_5deg_v1",
        "status":"verified_original_S1_source_and_observer_metadata_QC_only",
        "n_original_s1_observations":expected_n,
        "original_s1_decoded_sha256":sha256(native),
        "n_original_distinct_observers":int(roster.user_id.nunique()),
        "n_observer_respondent_ids_exported":0,
        "design":{"constructs":list(choices),"degrees":list(grains),
                  "seeds":REPLICATES,"min_photos_per_half_per_taxon_cell":2,
                  "min_cells_per_taxon":MIN_CELLS_TAXON,"min_original_photos_per_taxon":MIN_PHOTOS_TAXON,
                  "observer_assignment_globally_not_cellwise":True},
        "reference":"original-photo-ID-disjoint source split QC in PR #168; this is a STRICTER observer-disjoint followup with a changed support population",
        "not_equivalent_to_no_camera_or_individual_or_date_reuse":True,
        "not_an_ecological_causal_result_or_Random_Forest":True,
        "original_frozen_GEB_S1_S2_H3c_and_pr167_RF_untouched":True,
        "private_photos_locations_users_and_taxon_ids_not_exported":True,
        "results":results,
    }
    output.parent.mkdir(parents=True,exist_ok=True)
    output.write_text(json.dumps(report,indent=2,ensure_ascii=False,allow_nan=False)+"\n",encoding="utf-8")
    return report


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for k in ("native","receipt","sealed-metadata","output"):
        p.add_argument("--"+k,type=Path,required=True)
    args=p.parse_args()
    rep=run(args.native,args.receipt,args.sealed_metadata,args.output)
    print(json.dumps({"n_s1":rep["n_original_s1_observations"],
                      "n_observers":rep["n_original_distinct_observers"],
                      "diagnostics":len(rep["results"]),
                      "n_estimated":sum(x["status"]=="estimated_observer_disjoint_spatial_photo_QC_not_RF" for x in rep["results"])}),flush=True)

if __name__=="__main__":main()