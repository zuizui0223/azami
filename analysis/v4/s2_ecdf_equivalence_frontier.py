"""S2: quantify uncertainty needed to claim whole-angle-distribution equivalence.

Uses ONLY previously source-authenticated full ECDF curves for WCVP-native
versus ALL introduced Cirsium arvense/C. vulgare:
  GitHub Actions run 39592586205, artifact 12341964253,
  file SHA256 765742c0a2fdaef094f26de7999c9c322e5397a3c74ac0912edec4c2d09d9e4f.

This is a post hoc *margin-identification frontier*, not a newly
prespecified equivalence test. No ecological minimally important difference
(SESOI) is supplied. A 95% exploratory simultaneous native/introduction
cell-bootstrap ECDF band yields a conservative criterion:
  declare whole-grid equivalence at margin delta ONLY IF
  max(|simultaneous lower|, |simultaneous upper|) < delta,
  AND a domain-grounded delta was chosen without reference to this result.
The original full-ECDF grid is 401 pooled-rank knots, not every observation
jump, and the 399 cell bootstrap is exploratory not a biological CI.

At no point does this code convert photographic direction into true
gravity-referenced nodding, claim no evolution, or replace the positive
S2 wetter-introduced rainfall occupancy outcome.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

ORIGINAL_SOURCE_SHA256="765742c0a2fdaef094f26de7999c9c322e5397a3c74ac0912edec4c2d09d9e4f"
ORIGINAL_SOURCE_RUN=39592586205
ORIGINAL_SOURCE_ARTIFACT=12341964253
ORIGINAL_SOURCE_ARTIFACT_DIGEST="sha256:f892cd15cbf8be0232e9ec5286d98b168ae47f813685bb1480ed6c53e638a910"
SPECIES={
    "Cirsium arvense":{"rain":(30523,38406),"photo":(21870,25657)},
    "Cirsium vulgare":{"rain":(26161,49088),"photo":(16246,30287)}
}


def input_check(path: Path) -> dict:
    source=path.read_bytes()
    if hashlib.sha256(source).hexdigest()!=ORIGINAL_SOURCE_SHA256:
        raise ValueError("Refusing unverified S2 original full-rain/full-angle ECDF source")
    r=json.loads(source)
    if (r.get("status")!="ORIGINAL_NATIVE_ALL_INTRODUCED_FULL_EMPIRICAL_CURVES"
            or r.get("no_selected_outcome_threshold_or_effect_score") is not True
            or r.get("complete_ECDF_geographic_uncertainty_without_selected_metric") is not True
            or set(r.get("species",{}))!=set(SPECIES)
            or any(r.get("privacy",{}).values())):
        raise ValueError("S2 source scope, eligible species or privacy drift")
    return r


def inspect_curve(curve,*,n_native,n_introduced,unit):
    if (curve.get("native_original_observations")!=n_native
            or curve.get("all_introduced_original_observations")!=n_introduced
            or curve.get("unit")!=unit
            or curve.get("native_and_introduced_geographic_cells_equal_weight") is not True):
        raise ValueError("S2 original taxon/range denoms changed")
    b=curve["whole_curve_spatial_uncertainty"]
    if (b.get("status")!="EXPLORATORY_FULL_ECDF_GEOGRAPHIC_CELL_UNCERTAINTY"
            or b.get("resamples_requested")!=399 or b.get("resamples_valid")<380
            or b.get("simultaneous_confidence_envelope_exploratory_conditional_on_source") is not True):
        raise ValueError("Insufficient calibrated source-cell ECDF resamples")
    x=np.asarray(curve["display_only_unselected_common_knots"],float)
    f=np.asarray(curve["native_empirical_cdf"],float)
    g=np.asarray(curve["all_introduced_empirical_cdf"],float)
    d=np.asarray(b["introduced_minus_native_ECDF"],float)
    lo=np.asarray(b["simultaneous_95pct_spatial_sensitivity_lower"],float)
    hi=np.asarray(b["simultaneous_95pct_spatial_sensitivity_upper"],float)
    if (x.ndim!=1 or min(map(len,(x,f,g,d,lo,hi)))<100
        or not all(z.shape==x.shape for z in (f,g,d,lo,hi))
        or not all(np.isfinite(z).all() for z in (x,f,g,d,lo,hi))
        or not np.all(np.diff(x)>0)
        or not np.all(np.diff(f)>=-1e-10)
        or not np.all(np.diff(g)>=-1e-10)
        or not np.allclose(d,g-f,rtol=1e-6,atol=1e-7)
        or np.any(lo>d+1e-8) or np.any(hi<d-1e-8)
        or np.any(lo< -1-1e-8) or np.any(hi>1+1e-8)):
        raise ValueError("Reported full-curve uncertainty grid is inconsistent")
    obs_sup=float(np.max(np.abs(d)))
    # Supremum of upper absolute displacement among ALL gridwise members of
    # reported simultaneous 95% region: a conservative critical margin.
    upper=float(max(np.max(np.abs(lo)),np.max(np.abs(hi))))
    # Conservative lower bound on |cdf difference| from the same region.
    gap=np.maximum(np.where(lo>0,lo,0.),np.where(hi<0,-hi,0.))
    lower=float(np.max(gap))
    arg=int(np.argmax(np.abs(d)))
    if not (0<=lower<=obs_sup+1e-8<=upper+1e-8<=1+1e-8):
        raise ValueError("Equivalence supremum uncertainty bounds invalid")
    right_shifts=float(np.mean(d< -1e-10))
    zero_crossing=bool(np.all(lo<=0)&np.all(hi>=0))
    return {
      "type":"UNSELECTED_FULL_CDF_GRID_DIFFERENCE_NATIVE_VS_ALL_INTRODUCED",
      "unit":unit,
      "n_native_original_cell_balanced_observations":n_native,
      "n_introduced_original_cell_balanced_observations":n_introduced,
      "n_native_spatial_5deg_cells":b["native_5degree_cells"],
      "n_introduced_spatial_5deg_cells":b["introduced_5degree_cells"],
      "display_knots":len(x),
      "observed_sup_absolute_CDF_difference":obs_sup,
      "max_observed_difference_at_value":float(x[arg]),
      "observed_signed_CDF_difference_at_max":float(d[arg]),
      "simultaneous_95pct_envelope_lower_bound_on_sup_CDF_difference":lower,
      "simultaneous_95pct_envelope_upper_bound_on_sup_CDF_difference":upper,
      "minimum_delta_for_95pct_simultaneous_band_containment":upper,
      "delta_units":"absolute cumulative probability = 100 * delta percentage points",
      "hypothetical_margin_decision_rule":"If delta>upper and delta independently selected in advance, source-grid band would fit within ±delta; without a justified delta, no equivalence claim is made",
      "all_95pct_CDF_difference_bands_include_zero_at_reported_knots":zero_crossing,
      "fraction_of_display_knots_where_introduced_CDF_below_native":right_shifts,
      "ECDF_bootstrap_draws":b["resamples_valid"],
      "NOT_a_prospectively_chosen_SESOI":True,
      "NOT_a_calibrated_formal_equivalence_TOST":True,
      "NOT_a_proof_of_unchanged_biological_head_posture":True,
      "reported_Plotting_grid_only_not_exact_all_empirical_CDF_jumps":True,
      "source_bootstrap_uncertainty_exploratory_and_spatial_cells_not_independent_genetic_replicates":True,
    }


def evaluate(source_path:Path,output_path:Path):
    source=input_check(source_path)
    result={}
    for name,sample in SPECIES.items():
        entry=source["species"][name]
        if set(entry)!={"rainfall_environment_full_distribution","image_angle_full_distribution"}:
            raise ValueError("Original two S2 positive/neutral source distributions changed")
        result[name]={
          "introduced_rain_occurrence_shift":inspect_curve(
              entry["rainfall_environment_full_distribution"],
              n_native=sample["rain"][0],n_introduced=sample["rain"][1],
              unit="mm_per_year"),
          "introduced_photo_angle_equivalence_frontier":inspect_curve(
              entry["image_angle_full_distribution"],
              n_native=sample["photo"][0],n_introduced=sample["photo"][1],
              unit="image_degrees")
        }
    receipt={
      "schema":"V4_S2_ECDF_IMAGE_ANGLE_WHOLE_CURVE_EQUIVALENCE_MARGIN_FRONTIER_20261022",
      "status":"ACTUAL_SOURCE_GRID_UNCERTAINTY_FRONTIER_COMPUTED_WITHOUT_SELECTING_AN_EQUIVALENCE_MARGIN",
      "source_provenance":{
        "workflow_run":ORIGINAL_SOURCE_RUN,
        "artifact_id":ORIGINAL_SOURCE_ARTIFACT,
        "artifact_digest":ORIGINAL_SOURCE_ARTIFACT_DIGEST,
        "aggregate_source_sha256":ORIGINAL_SOURCE_SHA256,
        "filename":source_path.name,
        "same_S1_and_S2_original_species_not_independent_studies":True,
      },
      "answer_now":{
        "S2_primary_positive":"Both species occur at wetter annual-precipitation distributions in full introduced vs native WCVP ranges.",
        "S2_photo_orientation":"No full-curve CDF difference strong enough to exclude zero on the original 399-cell bootstrap simultaneous bands; NOT a proof of equivalence.",
        "whether_equivalence_conclusion_is_possible_without_SESOI":"NO",
        "SESOI_must_be_fixed_from_independent_biological_meaning_and_angle_reliability":True,
        "minimum_margin_given_for_FEASIBILITY_ONLY":True,
        "not_a_replacement_for_original_S1_positive_among_taxon_climate_association":True
      },
      "species":result,
      "limits":{
        "CDF_supremum_is_observation_distribution_not_evolutionary_nodding":True,
        "preselected_biological_SESOI_not_available":True,
        "cannot_certify_true_same_head_camera_angle_reliability":True,
        "two_species_do_not_establish_genus_wide_evolutionary_invariance":True,
        "bootstraps_not_corrected_for_sampling_areas_between_range_states":True,
        "missing_climate_for_original_cached_S2_transfer_does_not_affect_full_CDF_observed_distribution":True,
        "no_threshold_defined_from_favourable_observed_p_value":True,
      }
    }
    output_path.parent.mkdir(parents=True,exist_ok=True)
    output_path.write_text(json.dumps(receipt,indent=2,sort_keys=True,
                                      ensure_ascii=False,allow_nan=False)+"\n",encoding="utf-8")
    print("S2_ECDF_EQUIVALENCE_FRONTIER_SUMMARY="+json.dumps({
      "status":receipt["status"],
      "species":{k:{
        "photo_angle_frontier":v["introduced_photo_angle_equivalence_frontier"],
        "rain_source_distribution":{z:v["introduced_rain_occurrence_shift"][z]
           for z in ["observed_sup_absolute_CDF_difference",
                    "all_95pct_CDF_difference_bands_include_zero_at_reported_knots"]}
      } for k,v in result.items()}
    },sort_keys=True,ensure_ascii=False),flush=True)
    return receipt


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source",required=True,type=Path)
    parser.add_argument("--out",required=True,type=Path)
    args=parser.parse_args()
    evaluate(args.source,args.out)


if __name__=="__main__":
    main()
