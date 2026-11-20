"""Recover and audit S1 V5 angle-risk prevalence from authenticated July 2026 source.

Fallback for missing original prepared-S1 measurement receipt artifact: use the
source-sealed WCVP observation metadata, the source-sealed head measurements
and SHA-verified original CHELSA v2.1 five-axis rasters already used in V5.
Exact original annual V5 cohort size and grid supports are hard audit gates.

No original source IDs, taxon names, coordinates or individual measurements
are exported, only diagnostic aggregate counts.
"""
from __future__ import annotations
from pathlib import Path
import os
import json
import pandas as pd
import numpy as np

from analysis.v4.native_introduced_opportunity import (
    decrypt_checked, ORIGINAL_SEALED_SHA, ORIGINAL_SHA, ORIGINAL_AAD,
    source_photos, source_roster, authority_source
)
from analysis.v4.native_cached_climate import VARIABLES
from analysis.v5.native_pca_env_fullcoverage import apply_original_five_CHELSA
from analysis.v5.native_pca_env import background_mesh
from analysis.v4.s1_climate_trait_hypervolume import geo_aggregate, AXES
from analysis.v5.audit_angle_suspects import (
    ANGLE, ASPECT, ANNUAL, RAIN, _profile, _rain_quintiles, EXPECTED_NATIVE,
    EXPECTED_BY_GRAIN
)

PRIVATE=Path("work/submission/private")
OUT=Path("work/submission/public")

def main():
    blob=decrypt_checked(
        (PRIVATE/"original-metadata.aesgcm").read_bytes(),
        os.environ["AZAMI_METADATA_KEY"],
        ORIGINAL_SEALED_SHA,ORIGINAL_SHA,ORIGINAL_AAD
    )
    authority=authority_source((PRIVATE/"authority350.zip").read_bytes())
    roster=source_roster(source_photos(blob),authority)
    del blob
    native=roster.loc[roster.native_primary].copy()
    if len(native)!=313350:
        raise ValueError("NATIVE_WCVP_COHORT_CHANGED")
    measured=pd.read_csv(
        PRIVATE/"original-measured.csv.gz",
        dtype={"obs_id":str},low_memory=False
    )
    if not measured.obs_id.is_unique:
        raise ValueError("ORIGINAL_OBSERVATION_MEASUREMENT_DUPLICATES")
    traits={
        "orientation_angle_degrees_median":ANGLE,
        "corolla_lab_chroma_median":"corolla_lab_chroma",
        "shape_aspect_ratio_median":ASPECT,
    }
    qc=["n_usable_heads_orientation","n_usable_heads_colour",
        "n_usable_heads_shape"]
    required=["obs_id",*traits,*qc]
    if not set(required).issubset(measured.columns):
        raise ValueError("SOURCE_HEAD_MEASUREMENT_FIELDS_CHANGED: "+
                         repr(sorted(set(required)-set(measured))))
    native=native[["obs_id","accepted_key","latitude","longitude"]].merge(
        measured[required],on="obs_id",how="left",validate="one_to_one")
    for k in [*traits,*qc]:
        native[k]=pd.to_numeric(native[k],errors="coerce")
    pre=(np.isfinite(native[list(traits)].to_numpy(float)).all(axis=1)
         & native["orientation_angle_degrees_median"].between(0,180))
    for c in qc:
        pre &= native[c].gt(0)
    eligible=native.loc[pre,[
        "accepted_key","latitude","longitude",*traits
    ]].copy().rename(columns=traits)
    print(json.dumps({
       "native_source_records":313350,
       "three_trait_preclimate_eligible":len(eligible),
    }))
    if len(eligible)<EXPECTED_NATIVE:
        raise ValueError("PRECLIMATE_FALLS_BELOW_V5_ORIGINAL")
    land=background_mesh(authority)
    eligible,land,_=apply_original_five_CHELSA(
        eligible,land,
        source_meta=Path("analysis/v4/environment_source_metadata.json"),
        work=PRIVATE/"s1_recovery_CHELSA5"
    )
    selected=eligible.loc[np.isfinite(eligible[list(VARIABLES)].to_numpy(float)).all(axis=1)].copy()
    print(json.dumps({"original_annual_climate_and_three_traits":len(selected)}))
    if len(selected)!=EXPECTED_NATIVE:
        raise ValueError("S1_SOURCE_RECOVERY_COUNT_DIFFERS_FROM_LOCKED_V5_"+str(len(selected)))

    for name in AXES+ANNUAL:
        if name not in selected.columns:
            raise ValueError("V5_ANNUAL_AXES_NOT_VERIFIED: "+name)

    grid={}
    for deg,expected in EXPECTED_BY_GRAIN.items():
        cells,_=geo_aggregate(selected,deg,ANNUAL)
        frame=selected.copy()
        frame["_latcell"]=np.floor((frame.latitude+90)/deg).astype(int)
        frame["_loncell"]=np.floor((frame.longitude+180)/deg).astype(int)
        valid=cells[["accepted_key","_latcell","_loncell"]]
        frame=frame.merge(valid,on=["accepted_key","_latcell","_loncell"],
                          how="inner",validate="many_to_one")
        if (len(frame),cells.accepted_key.nunique(),len(cells))!=(
            expected["records"],expected["taxa"],expected["taxon_cells"]
        ):
            raise ValueError(f"S1_RECOVERED_GRID_COHORT_DIFFERS_{deg}: "+
                             repr((len(frame),cells.accepted_key.nunique(),len(cells))))
        grid[f"{deg}deg"]={
            "taxa":expected["taxa"],"taxon_cells":expected["taxon_cells"],
            **_profile(frame)
        }
    report={
       "schema":"azami_v5_S1_recovered_originals_angle_proxy_flags_v1",
       "denominator":"Source-confirmed 157822 V5 native photo-derived observation records",
       "source":"July-2026 encrypted WCVP+original measurements+SHA-authenticated CHELSA v2.1 1981-2010",
       "quality_interpretation":"Rotated outline aspect is not PCA eigenvalue anisotropy or proven measurement error",
       "all_s1":_profile(selected),
       "by_grain":grid,
       "rain_quintile":_rain_quintiles(selected)
    }
    OUT.mkdir(parents=True,exist_ok=True)
    path=OUT/"v5_angle_suspect_s1_recovered.json"
    path.write_text(json.dumps(report,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    print("V5_ANGLE_S1_RECOVERED_AUDIT_AGGREGATED_JSON_START")
    print(json.dumps(report,ensure_ascii=False,separators=(",",":")))
    print("V5_ANGLE_S1_RECOVERED_AUDIT_AGGREGATED_JSON_END")

if __name__=="__main__":main()
