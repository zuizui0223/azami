"""Reconstruct complete S2 photo-source and continent curves from public JSON.

No original photographs, coordinates, photographer IDs or private token needed:
use ONLY the successfully source-authenticated S2_full_joint_rain_angle_landscape
artifact. Outputs all-curve distributions, not chosen width/tail effects.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

SPECIES=("Cirsium arvense","Cirsium vulgare")
INPUT_STATUS="ORIGINAL_S2_JOINT_RAIN_AND_IMAGE_ANGLE_ALL_INTRODUCED"


def _vec(values, *, unit_mass:bool=False):
    a=np.asarray(values,dtype=float)
    if a.ndim!=1 or not len(a) or not np.isfinite(a).all() or np.any(a<0):
        raise ValueError("Invalid public aggregate probability vector")
    if unit_mass and not np.isclose(a.sum(),1,atol=1e-7):
        raise ValueError("Full original-source probability mass not one")
    return a


def _matrix(values):
    a=np.asarray(values,dtype=float)
    if (a.ndim!=2 or 0 in a.shape or
        not np.isfinite(a).all() or (a<0).any() or
        not np.isclose(a.sum(),1,atol=1e-7)):
        raise ValueError("Source-verified two-dimensional distribution invalid")
    return a


def reconstruct(report:dict)->dict:
    """Exactly reconstruct whole-range climatic and image-angle ECDFs.

    Figure construction downstream may inspect every original plotting knot.
    Suppressed low-support regions remain in the *primary* all-introduced
    curves; the restricted curve is a separate diagnostic only.
    """
    if (report.get("status")!=INPUT_STATUS or
        not report.get("all_introduced_remains_primary") or
        not report.get("visualization_only_no_trait_width_or_tail_indices_selected") or
        set(report.get("species",{}))!=set(SPECIES) or
        any(report.get("privacy",{}).values())):
        raise ValueError("This analysis requires authentic privacy-safe S2 public joint aggregates")
    cases={}
    for species,case in report["species"].items():
        rain_edges=_vec(case["display_climate_BIO12_edges_mm_year"])
        angle_edges=_vec(case["display_angle_edges_image_degrees"])
        if (len(rain_edges)<10 or len(angle_edges)<10 or
            (np.diff(rain_edges)<=0).any() or (np.diff(angle_edges)<=0).any()):
            raise ValueError("Invalid unselected original joint rainfall/image display axes")
        curves={}
        for range_name in ("native","all_introduced"):
            panel=case[range_name]
            source=_vec(panel["original_wild_source_rain_bin_mass"],unit_mass=True)
            photographed=_vec(panel["original_photo_rain_bin_mass"],unit_mass=True)
            joint=_matrix(panel["original_photo_joint_mass"])
            if joint.shape!=(len(rain_edges)-1,len(angle_edges)-1):
                raise ValueError("Rain and image-angle plot axes mismatch")
            if not np.allclose(joint.sum(axis=1),photographed,atol=1e-7):
                raise ValueError("Photograph source climate and joint distributions disagree")
            curves[range_name]={
                "original_observations":int(panel["original_photo_observations"]),
                "geographical_cells":int(panel["original_geographic_cells"]),
                "original_wild_source_BIO12_ECDF":np.cumsum(source).tolist(),
                "photo_eligible_BIO12_ECDF":np.cumsum(photographed).tolist(),
                "same_photo_image_angle_ECDF":np.cumsum(joint.sum(axis=0)).tolist(),
            }
        global_angle=np.asarray(curves["all_introduced"]["same_photo_image_angle_ECDF"])
        groups=case["introduced_continents_supplementary"]["groups"]
        total=np.zeros(len(angle_edges)-1,float)
        retained=np.zeros_like(total)
        retained_weight=0.
        regional={}
        for name,part in groups.items():
            weight=float(part["fraction_global_cell_weight"])
            matrix=_matrix(part["joint_mass_within_region"])
            if matrix.shape!=(len(rain_edges)-1,len(angle_edges)-1):
                raise ValueError("Original regional joint angle dimensions mismatch")
            if weight<=0 or weight>1:
                raise ValueError("Invalid region fraction of global cell-balanced source")
            curve=matrix.sum(axis=0)
            total+=weight*curve
            if name!="other_low_support":
                retained+=weight*curve
                retained_weight+=weight
            regional[name]={
                "original_photo_observations":int(part["photo_n"]),
                "geographical_cells":int(part["photo_cells"]),
                "global_original_cell_weight_fraction":weight,
                "full_image_angle_ECDF":np.cumsum(curve).tolist(),
                "full_photo_BIO12_ECDF":np.cumsum(matrix.sum(axis=1)).tolist(),
            }
        if not np.allclose(np.cumsum(total),global_angle,atol=1e-8):
            raise ValueError("All introduced picture could not be reassembled")
        curves["regional_photo_angle_ECDF"]=regional
        curves["introduced_supported_regions_only_ECDF_sensitivity"]=(
            np.cumsum(retained/retained_weight).tolist()
            if retained_weight>0 else None)
        curves["rain_edges_mm"]=rain_edges.tolist()
        curves["angle_edges_image_degrees"]=angle_edges.tolist()
        curves["same_rain_photo_angle_conditional_difference_or_null"]=[
            [None if a is None or b is None else float(b-a)
             for a,b in zip(row_native,row_intro)]
            for row_native,row_intro in zip(
                case["native"]["photo_angle_conditional_on_rain_or_null"],
                case["all_introduced"]["photo_angle_conditional_on_rain_or_null"])
        ]
        curves["regional_reconstruction_exact"]=True
        cases[species]=curves
    return {
        "schema":"S2_all_introduced_photo_coverage_and_regional_ECDF_v1",
        "status":"DERIVED_FROM_PREVERIFIED_SOURCE_AGGREGATE_ONLY",
        "primary_comparison":"native versus ALL WCVP explicitly introduced regions",
        "no_selected_quantile_width_or_extreme_angle_indicator":True,
        "all_photo_identity_and_coordinates_absent":True,
        "note":"Entire PDF/CDF curves only; regional and photo-sampling display diagnoses are not ecological significance tests",
        "species":cases,
    }


def main()->None:
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--joint-json",type=Path,required=True)
    parser.add_argument("--out",type=Path,required=True)
    args=parser.parse_args()
    report=json.loads(args.joint_json.read_text(encoding="utf-8"))
    derived=reconstruct(report)
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(derived,ensure_ascii=False,indent=2,allow_nan=False)+"\n",encoding="utf-8")
    print("S2 joint full distribution diagnostic: original provenance preserved; no selected indices")


if __name__=="__main__":
    main()
