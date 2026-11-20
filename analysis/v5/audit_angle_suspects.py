"""Source-locked prevalence audit for potentially unstable photographed head angle.

Measure *flags* on the same 157,822 native image-derived records used by the
annual five-axis V5 analysis, not verified mismeasurements. All public results
are aggregated, with NO source IDs, exact coordinates or taxon identifiers.

IMPORTANT: the rotated bounding-box aspect ratio is NOT the PCA eigenvalue
anisotropy used by the orientation estimator. Near-circular outline flags
are proxies and do not establish that an angle is incorrect.

Examples:
    python -m analysis.v5.audit_angle_suspects \
      --s1-source work/bio18_geo/private/decoded/observations.csv.gz \
      --s1-receipt work/bio18_geo/private/decoded/measurement_receipt.json \
      --output work/bio18_geo/public/v5_angle_suspect_audit.json
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from analysis.v4.s1_climate_trait_hypervolume import load_source, geo_aggregate

ANGLE = "orientation_image_vertical_angle"
ASPECT = "capitulum_outline_aspect_ratio"
RAIN = "chelsa_bio12"
ANNUAL = ("chelsa_bio04", "chelsa_bio15")
EXPECTED_NATIVE = 157822
EXPECTED_BY_GRAIN = {
    1: {"records": 153843, "taxa": 63, "taxon_cells": 7117},
    5: {"records": 132998, "taxa": 31, "taxon_cells": 814},
}

def _round(f: float) -> float:
    return round(float(f), 6)

def _profile(frame: pd.DataFrame) -> dict:
    angle = frame[ANGLE].to_numpy(float)
    aspect = frame[ASPECT].to_numpy(float)
    assert np.isfinite(angle).all() and np.isfinite(aspect).all()
    assert (angle >= 0).all() and (angle <= 180).all()
    assert (aspect > 0).all()
    n = len(frame)
    if not n:
        raise ValueError("EMPTY_ANGLE_QUALITY_COHORT")
    # Symmetric near-one ratio; e.g. 1/1.1 < A < 1.1
    near = {
        str(cutoff): np.abs(np.log(aspect)) <= np.log(cutoff)
        for cutoff in [1.05, 1.10, 1.20, 1.30, 1.50]
    }
    bins = {
        "image_angle_gt90": angle > 90,
        "image_angle_ge110": angle >= 110,
        "image_angle_ge150": angle >= 150,
        "image_angle_60_to_120": (angle >= 60) & (angle <= 120),
        "image_angle_0_to_30": angle <= 30,
    }
    bins.update({f"near_axis_ratio_{key}": mask for key, mask in near.items()})
    for cut in ("1.05", "1.10", "1.20", "1.30"):
        for ending in ("image_angle_gt90", "image_angle_ge110", "image_angle_60_to_120"):
            bins[f"near_{cut}_and_{ending}"] = near[cut] & bins[ending]
    counts = {
        key: {"n": int(mask.sum()), "fraction": _round(mask.mean())}
        for key, mask in bins.items()
    }
    # A diagnostic only; does not imply that a sign should flip.
    def distrib(arr):
        return {str(q): _round(np.quantile(arr, float(q))) for q in (0.05,.25,.5,.75,.95)}
    return {
        "n_source_measurement_records": int(n),
        "image_angle_degrees_quantiles": distrib(angle),
        "outline_aspect_quantiles": distrib(aspect),
        "flag_counts": counts,
        "angle90plus_rate_by_outline_class": [
            {
                "outline_class": name,
                "n": int(mask.sum()),
                "n_angle_gt90": int((mask & bins["image_angle_gt90"]).sum()),
                "fraction_angle_gt90": (_round((mask & bins["image_angle_gt90"]).sum()/mask.sum())
                                       if mask.any() else None),
            }
            for name, mask in (
                ("near_one_aspect_le_1.10",near["1.1"] if "1.1" in near else near["1.10"]),
                ("outside_near_one_aspect",~near["1.10"]),
            )
        ],
    }

def _rain_quintiles(frame: pd.DataFrame) -> list[dict]:
    # Equal-count source record quintiles, not biological causal evidence.
    ranked = frame[RAIN].rank(method="first")
    idx = pd.qcut(ranked, 5, labels=False, duplicates="raise")
    out=[]
    for j in range(5):
        part = frame.loc[idx.eq(j)]
        entry = _profile(part)
        out.append({
            "rain_quintile": j+1,
            "rain_median_mm": _round(part[RAIN].median()),
            "n": len(part),
            "near_aspect_1p1_count": entry["flag_counts"]["near_axis_ratio_1.1"]["n"]
                                    if "near_axis_ratio_1.1" in entry["flag_counts"]
                                    else entry["flag_counts"]["near_axis_ratio_1.10"]["n"],
            "near_aspect_1p1_rate": entry["flag_counts"]["near_axis_ratio_1.10"]["fraction"],
            "near_aspect_1p1_and_angle_gt90_count":
                entry["flag_counts"]["near_1.10_and_image_angle_gt90"]["n"],
            "near_aspect_1p1_and_angle_gt90_rate":
                entry["flag_counts"]["near_1.10_and_image_angle_gt90"]["fraction"],
        })
    return out

def audit_s1(source: Path, receipt: Path) -> dict:
    selected, receipt_info = load_source(source,receipt,ANNUAL)
    if len(selected)!=EXPECTED_NATIVE:
        raise ValueError(f"ORIGINAL_V5_ANNUAL_N_CHANGED: {len(selected)}")
    all_prof=_profile(selected)
    grain_results={}
    for degree, expected in EXPECTED_BY_GRAIN.items():
        cells,_=geo_aggregate(selected, degree, ANNUAL)
        # Source observations retained by taxon+cell membership, not merely
        # those taxa overall. This is *not* rarefied cell-pair counts.
        frame=selected.copy()
        frame["_latcell"]=np.floor((frame.latitude+90)/degree).astype(int)
        frame["_loncell"]=np.floor((frame.longitude+180)/degree).astype(int)
        allowed=cells[["accepted_key","_latcell","_loncell"]].copy()
        frame=frame.merge(allowed,on=["accepted_key","_latcell","_loncell"],
                          how="inner",validate="many_to_one")
        if (len(frame),cells.accepted_key.nunique(),len(cells)) != (
            expected["records"],expected["taxa"],expected["taxon_cells"]
        ):
            raise ValueError(f"ORIGINAL_V5_GRID_COHORT_CHANGED: {degree}")
        grain_results[f"{degree}deg"] = {
            "taxa":expected["taxa"],"taxon_cells":expected["taxon_cells"],
            **_profile(frame),
        }
    return {
        "schema":"azami_v5_s1_image_angle_proxy_flags_v1",
        "denominator":"photo-derived observation-level rows, not individual heads or unique JPEGs",
        "quality_interpretation":(
            "near-one rotated outline aspect is a PROXY, not eigenvalue "
            "anisotropy, head roll, human-matched orientation or a measured "
            "probability of error; >90 is coded angle, not necessarily a defect"
        ),
        "source_verified":receipt_info["source_observations"]==313350,
        "annual_v5_eligible":all_prof,
        "source_grains":grain_results,
        "annual_precipitation_source_record_quintiles":_rain_quintiles(selected),
    }

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--s1-source",type=Path,required=True)
    p.add_argument("--s1-receipt",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    args=p.parse_args()
    result=audit_s1(args.s1_source,args.s1_receipt)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,ensure_ascii=False,indent=2)+"\n",
                           encoding="utf-8")
    # Single aggregated JSON is safe for job logs; no IDs/coordinates.
    print("V5_ANGLE_AUDIT_AGGREGATED_JSON_START")
    print(json.dumps(result,ensure_ascii=False,separators=(",",":")))
    print("V5_ANGLE_AUDIT_AGGREGATED_JSON_END")

if __name__=="__main__":
    main()
