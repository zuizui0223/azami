"""V5 I+II direct sensitivity: BIO05 to BIO10 while keeping BIO18 and seasonal axes."""
from __future__ import annotations
import json
from pathlib import Path
import pandas as pd
import numpy as np

from analysis.v5.background_paired_overlap import (
    _source_original, _historical_photo_pixel_cache
)
from analysis.v5.annual_warm_sources import sample_warm, PRIVATE, OUT, EXTRA
from analysis.v5.bio10_s1 import run as run_s1
from analysis.v5.bio10_pair_diagnostics import run as run_pair_selection_diagnostics
from analysis.v5.bio10_s2 import run as run_s2
from analysis.v5.bio10_s2_angle_chroma import run as run_s2_angle_chroma
from analysis.v5.bio10_sources import sample_bio10, PINNED_BIO10_SHA256

ROOT=Path("work/submission")
NAME="v5_part1_part2_bio05_vs_bio10_same_seasonality.json"


def main():
    OUT.mkdir(parents=True,exist_ok=True)
    original=PRIVATE/"decoded/observations.csv.gz"
    receipt=PRIVATE/"decoded/measurement_receipt.json"
    if not original.exists() or not receipt.exists():
        raise ValueError("ORIGINAL_VERIFIED_S1_NATIVE_SOURCE_REQUIRED")
    native=pd.read_csv(original,dtype={"obs_id":str},low_memory=False)
    if len(native)!=313350 or not native.obs_id.is_unique:
        raise ValueError("S1_ORIGINAL_NATIVE_CENSUS_DRIFT")
    wild,land,authority,source_five=_source_original()
    historic,cache_receipt=_historical_photo_pixel_cache(wild)
    wild,land,native,warm_receipt=sample_warm(wild,land,native)
    wild,land,native,bio10_receipt=sample_bio10(wild,land,native)
    if len(historic)!=len(wild) or not historic.obs_id.equals(wild.obs_id):
        raise ValueError("ORIGINAL_HISTORICAL_V5_PHOTO_SOURCE_ORDER_DRIFT")
    for var in (*EXTRA, 'chelsa_bio10'):
        historic[var]=wild[var].to_numpy(float)
    for var,meta in warm_receipt.items():
        if (meta["wild_finite"]<len(wild)*.98 or
            meta["native_s1_finite"]<len(native)*.98 or
            meta["background_finite"]<len(land)*.95):
            raise ValueError("ORIGINAL_PINNED_WARM_PIXEL_COVERAGE_INADEQUATE_"+var)
    # Deliberately keep both existing annual climate and new warm climate
    # on the SAME original observations and photo QC. No missing fill.
    s1=run_s1(native,original,receipt)
    s1["cross_taxon_pair_selection_diagnostics"]=run_pair_selection_diagnostics(
        native,original,receipt)
    print({"part":"I_pair_selection",
           "pair_overlap":[{k:r[k] for k in (
               "degree","model_a","model_b","taxon_pair_jaccard_mean",
               "cell_pair_jaccard_within_shared_taxon_pairs_mean")}
               for r in s1["cross_taxon_pair_selection_diagnostics"]["selected_taxon_pair_overlap"]]},
          flush=True)
    if s1["cohort"]["matched_annual_and_warm_source_photo_n"]!=157822:
        raise ValueError("BIO10_S1_COHORT_NOT_MATCHED_TO_OLD_V5")
    prev={(1,"annual_plus_seasonality"):0.01866528656793442,
          (1,"warm_plus_seasonality"):-0.029660276762971418,
          (5,"annual_plus_seasonality"):0.05407560661642653,
          (5,"warm_plus_seasonality"):-0.03188218217448156}
    for x in s1["cross"]:
        k=(x["degrees"],x["model"])
        if k in prev and abs(x["excess_trait_similarity"]-prev[k])>1e-9:
            raise ValueError("BIO10_S1_PREVIOUS_ANNUAL_BIO05_MODEL_DRIFT")
    print({"part":"I","source_photos":s1["cohort"]["matched_annual_and_warm_source_photo_n"],
           "cross":[{k:r[k] for k in ("degrees","model","excess_trait_similarity","p")}
                    for r in s1["cross"]]},flush=True)
    measured=pd.read_csv(PRIVATE/"original-measured.csv.gz",
                         usecols=["obs_id","orientation_angle_degrees_median",
                                  "corolla_lab_chroma_median","shape_aspect_ratio_median"],
                         dtype={"obs_id":str},low_memory=False)
    if measured.obs_id.duplicated().any():
        raise ValueError("S2_DECRYPTED_MEASUREMENTS_DUPLICATED")
    s2=run_s2(wild,land,authority,historic,measured)
    if s2["cohort"]["photo_eligible_n"] !=83872:
        raise ValueError("BIO10_S2_FIXED_PHOTO_COHORT_DRIFT")
    prev_rho={("Cirsium arvense","annual_plus_seasonality"):.4840274796049806,
              ("Cirsium arvense","warm_plus_seasonality"):.44296073660607793,
              ("Cirsium vulgare","annual_plus_seasonality"):.12232504337767498,
              ("Cirsium vulgare","warm_plus_seasonality"):-.021820284978179706}
    for x in s2["regional_overlap"]:
        if x["head_dimensions"]==3 and x["bandwidth_pc_sd"]==.5:
            k=(x["species"],x["model"])
            if k in prev_rho and abs(x["background_corrected_rho"]-prev_rho[k])>1e-9:
                raise ValueError("BIO10_S2_PREVIOUS_ANNUAL_BIO05_REGIONAL_RHO_DRIFT")
    s2["angle_and_chroma"]=run_s2_angle_chroma(historic,measured)
    print({"part":"II_angle_and_chroma",
           "primary_within_region":[{
               "species":r["species"],"climate":r["climate_model"],
               "native_beta":r["rain_partial_slope"]["native"]["beta"],
               "introduced_beta":r["rain_partial_slope"]["introduced"]["beta"],
               "difference_p_holm":r["holm_6_interaction"]
           } for r in s2["angle_and_chroma"]["results"]
             if r["model"]=="within_20degree_region" and r["trait"]=="angle"]},flush=True)
    print({"part":"II",
           "niches":[{k:r[k] for k in ("species","model",
                       "full_range_climate_overlap_D",
                       "introduced_corrected_mass_outside_native95")}
                     for r in s2["full_niche"]],
           "regional":[{k:r[k] for k in ("species","model","head_dimensions",
                       "background_corrected_rho",
                       "background_corrected_permutation_p")}
                     for r in s2["regional_overlap"]
                     if r["head_dimensions"]==3 and r["bandwidth_pc_sd"]==.5]},flush=True)
    report={
        "schema":"V5_PART_I_AND_II_BIO05_TO_BIO10_5AXIS_20261110_v1",
        "status":"BIO10_SOURCE_PENDING_PIN_REPEAT" if PINNED_BIO10_SHA256 is None else "BIO10_SOURCE_PINNED_AND_S1_S2_FITS_COMPLETE",
        "frozen_original_core": {
            "native_s1":313350,"wild_s2":297622,
            "all_five_old_climate_pixel_sources":source_five,
            "historical_photo_climate_cache":cache_receipt["status"],
            "warm_original_pinned_pixel_sources":warm_receipt,
            "bio10_original_newly_audited_raster":bio10_receipt
        },
        "part_I":s1,
        "part_II":s2,
        "not_a_prospective_prespecification":True,
        "bio10_mean_warmest_quarter_matches_BIO18_warmest_quarter_window":True,
        "not_a_test_of_genetic_evolution":True,
        "no_individual_photo_taxon_observer_coordinate_record_export":True,
    }
    dest=OUT/NAME
    dest.write_text(json.dumps(report,indent=2,ensure_ascii=False,allow_nan=False)+"\n",encoding="utf-8")
    print({"status":"V5_PARTS_I_II_BIO10_MODEL_COMPARISON_WRITTEN",
           "file":str(dest),"S1_n":s1["cohort"]["matched_annual_and_warm_source_photo_n"],
           "S2_n":s2["cohort"]["photo_eligible_n"]},flush=True)


if __name__=="__main__":main()
