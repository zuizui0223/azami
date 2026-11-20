"""Describe archived image-angle within-photo dispersion and usable-head support.

This is NOT measurement error against independent human or gravity labels.
A photo may contain different heads. Its between-head MAD can quantify
internal variation / apparent disagreement but not angular accuracy.

Inputs are private original source rows. Output is aggregate-only JSON.
"""
from __future__ import annotations
import numpy as np
import pandas as pd

COUNTER_NAMES=("n_usable_heads_orientation","orientation_angle_degrees_n_usable_heads")
MAD_NAMES=("orientation_angle_degrees_mad","orientation_angle_degrees_MAD",
           "orientation_angle_mad")
CONF_NAMES=("orientation_confidence","orientation_confidence_median",
            "orientation_confidence_mean")
FLIP_NAMES=("orientation_flip_delta_degrees",
            "orientation_flip_delta_degrees_median")

def _choose(columns: list[str], names:tuple[str,...])->str|None:
    return next((c for c in names if c in columns),None)

def summarise(selected_ids, original_measurements:pd.DataFrame)->dict:
    if "obs_id" not in original_measurements:
        raise ValueError("ORIGINAL_MEASUREMENTS_REQUIRE_OBS_ID")
    if not original_measurements.obs_id.is_unique:
        raise ValueError("ORIGINAL_MEASUREMENTS_DUPLICATE_OBS_ID")
    ids=pd.Series(selected_ids,dtype=str).drop_duplicates()
    if ids.empty or ids.isna().any():
        raise ValueError("EMPTY_OR_INVALID_COHORT")
    if len(ids)!=len(selected_ids):
        raise ValueError("COHORT_OBS_ID_DUPLICATES")
    cols=original_measurements.columns.to_list()
    ncol=_choose(cols,COUNTER_NAMES)
    madcol=_choose(cols,MAD_NAMES)
    confcol=_choose(cols,CONF_NAMES)
    flipcol=_choose(cols,FLIP_NAMES)
    joined=pd.DataFrame({"obs_id":ids}).merge(
        original_measurements[["obs_id", *[x for x in
            (ncol,madcol,confcol,flipcol) if x is not None]]],
        how="left",on="obs_id",validate="one_to_one",indicator=True)
    present=int(joined["_merge"].eq("both").sum())
    report={
        "cohort_records":len(ids),
        "original_measurement_matched":present,
        "source_observation_angle_data_fields":{
            "usable_head_count":ncol,
            "within_photo_head_angle_mad":madcol,
            "orientation_confidence":confcol,
            "horizontal_flip_delta":flipcol,
        },
        "repeatability_is_not_accuracy":True,
        "angles_not_referenced_to_gravity":True
    }
    if present!=len(ids):
        raise ValueError("PHOTO_COHORT_MISSING_SOURCE_MEASUREMENT_ROW")
    def profile(series,mask=None):
        v=pd.to_numeric(series,errors="coerce").replace(
            [np.inf,-np.inf],np.nan)
        if mask is not None:
            v=v.loc[mask]
        t=v.dropna()
        return {
            "denominator":int(len(v)),
            "finite":int(t.size),
            "median":round(float(t.median()),4) if len(t) else None,
            "q05":round(float(t.quantile(.05)),4) if len(t) else None,
            "q95":round(float(t.quantile(.95)),4) if len(t) else None,
        }
    if ncol:
        counts=pd.to_numeric(joined[ncol],errors="coerce")
        if counts.isna().any() or counts.lt(1).any():
            raise ValueError("MAIN_COHORT_WITHOUT_USABLE_HEAD")
        report["usable_heads"]={
            "single_head_n":int(counts.eq(1).sum()),
            "single_head_fraction":round(float(counts.eq(1).mean()),6),
            "multiple_heads_n":int(counts.ge(2).sum()),
            "multiple_heads_fraction":round(float(counts.ge(2).mean()),6),
            "at_least_3_heads_n":int(counts.ge(3).sum()),
            **profile(counts),
        }
    else:
        report["usable_heads"]={"status":"head_count_unavailable"}
    if madcol:
        m=pd.to_numeric(joined[madcol],errors="coerce").replace(
             [np.inf,-np.inf],np.nan)
        mask=pd.to_numeric(joined[ncol],errors="coerce").ge(2) if ncol else m.notna()
        qualified=m.loc[mask]
        report["head_angle_MAD_degrees_when_at_least_two_heads"]={
            "warning":"within-photograph between-head heterogeneity, not measurement accuracy",
            **profile(qualified),
            "gt10_n":int((qualified>10).sum()),
            "gt20_n":int((qualified>20).sum()),
            "gt30_n":int((qualified>30).sum()),
            "gt20_fraction_of_multiple_head_photos":round(float(
                (qualified>20).sum()/max(int(mask.sum()),1)),6),
        }
    else:
        report["head_angle_MAD_degrees_when_at_least_two_heads"]={
            "status":"not retained in archived observation-level source"}
    if confcol:
        report["orientation_confidence"]=profile(joined[confcol])
    else:
        report["orientation_confidence"]={
            "status":"not retained at observation level; source quality gate did require usable head status"}
    if flipcol:
        delta=pd.to_numeric(joined[flipcol],errors="coerce")
        report["horizontal_flip_delta"]={
            **profile(delta),
            "gt20_n":int(delta.gt(20).sum()),
        }
    else:
        report["horizontal_flip_delta"]={
            "status":"not retained at observation level; archived usable-head gate includes mirror consistency"}
    return report
