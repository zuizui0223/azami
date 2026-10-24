"""S1 -> S2 -> T2 climate-niche bridge, without imposing a favored climate sign.

T2 uses two real, non-overlapping nuclear trees and independently described
botanical U/D (NA44, Japan12), not photographic angle thresholds.
An omnibus five-climate present-niche test uses source-fixed S1 axes and an
Mk prevalence-conditioned null. It is NOT a transition-time climate test.
S1/S2 original aggregate receipts remain unchanged and are not independent
replications of the same native observations.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from analysis.v4.t2_h1_real_nuclear_phylo_same_s1 import source_panels
from analysis.v4.t2_phylogenetic_climate_null import (
    fitted_rate_given_root, simulate_mk_tip_matrix
)

AXES = (
    "chelsa_bio12", "chelsa_bio15", "chelsa_bio01",
    "chelsa_bio04", "chelsa_rsds_mean"
)
FOCAL = ("Cirsium arvense", "Cirsium vulgare")
ROOTS = ("erect", "stationary")
SEED = 20261024
RIDGE = 0.15
SOURCE_S1S2 = Path(
    "analysis/v4/results/20261022_S1_S2_BETWEEN_WITHIN_NATIVE_INTRO_COUNTS_VERIFIED.json"
)
SOURCE_S2 = Path(
    "analysis/v4/results/20261021_unified_abcde/"
    "S2_NATIVE_INTRODUCED_CLIMATE_TRANSFER_20261021.json"
)


def evidence_bridge(receipt: dict, transfer: dict, source: pd.DataFrame) -> dict:
    """Verify that S2 native photos are a subset of S1, never replicate them."""
    assert receipt["S1"]["n_observations"] == 170066
    assert receipt["S1"]["n_accepted_keys"] == 244
    assert len(source) == 251 and source.accepted_key.nunique() == 244
    assert int(source.original_photo_count.sum()) == 170066
    assert set(transfer["species"]) == set(FOCAL)
    assert set(receipt["S2"]["same_species_comparisons"]) == set(FOCAL)
    b = float(receipt["S1"]["native_244_between_species_BIO12_degrees_per_100mm"])
    w = float(receipt["S1"]["source_native_244_within_species_BIO12_degrees_per_100mm"])
    bci = receipt["S1"]["between_95pct_interval"]
    wci = receipt["S1"]["within_95pct_interval"]
    if not (bci[0] < b < bci[1] and wci[0] < w < wci[1]):
        raise ValueError("S1 receipt interval mismatch")
    original = source.assign(
        normalized_name=source.original_source_name.str.replace("_", " ", regex=False).str.strip()
    )
    focal = {}
    for name in FOCAL:
        rows = original.loc[original.normalized_name == name]
        if len(rows) != 1:
            raise ValueError("S2 focal original taxon has ambiguous S1 source name: " + name)
        row = rows.iloc[0]
        rr = receipt["S2"]["same_species_comparisons"][name]
        tr = transfer["species"][name]
        native_n = int(rr["native"]["original_orientation_eligible_observations"])
        if int(row.original_photo_count) != native_n:
            raise ValueError("S1 native/S2 native original photos do not match for " + name)
        if int(tr["native"]["original_orientation"]) != native_n:
            raise ValueError("S2 source missing original native images for " + name)
        shift = rr["differences"]
        delta = float(shift["introduced_minus_native_median_BIO12_mm"])
        angle_delta = float(shift["introduced_minus_native_photo_angle_median_degrees"])
        if delta <= 0:
            raise ValueError("Previously source-verified wetter introduced median changed")
        values = {}
        for range_state in ("native", "introduced"):
            z = tr[range_state]
            ci = z["cluster_ci95"]
            effect = float(z["bio12_deg_100mm"])
            if not ci[0] <= effect <= ci[1]:
                raise ValueError("S2 regional coefficient/CI mismatch")
            values[range_state] = {
                "original_angle_eligible": int(z["original_orientation"]),
                "model_climate_complete": int(z["environment_complete"]),
                "model_balanced_observer_cell": int(z["balanced_observer_cell_n"]),
                "BIO12_slope_image_deg_per100mm": effect,
                "cell_cluster_95CI": ci,
                "model_CI_includes_zero": bool(ci[0] <= 0 <= ci[1]),
            }
        focal[name] = {
            "S1_exact_native_source_name": str(row.original_source_name),
            "S1_exact_native_photo_count": native_n,
            "S1_botanical_posture_code": str(row.botanical_U_D_state),
            "same_original_native_source_as_S1": True,
            "introduced_minus_native_occurrence_BIO12_median_mm": delta,
            "introduced_minus_native_image_angle_median_deg": angle_delta,
            "native_and_introduced_are_not_the_same_individual": True,
            "regional_slope_evidence": values,
            "native_model_prediction_skill_in_introduced":
                tr["predictive_transfer"]["native_trained_skill_relative_to_native_constant"],
            "original_S2_native_intro_slope_difference":
                tr["introduced_minus_native"],
        }
    return {
        "S1": {
            "n_native_source_images": 170066, "n_accepted_taxon_keys": 244,
            "same_five_CHELSA_predictors": list(AXES),
            "BIO12_between_taxa_image_deg_per100mm": b,
            "BIO12_between_taxa_95CI": bci,
            "BIO12_within_taxa_image_deg_per100mm": w,
            "BIO12_within_taxa_95CI": wci,
            "descriptive_between_minus_within": b-w,
            "descriptive_difference_is_NOT_jointly_tested_without_covariance": True,
        },
        "S2": {
            "n_supported_species": 2,
            "native_images_are_already_in_S1_not_independent_confirmation": True,
            "species": focal,
            "source_biased_to_two_introduced_species": True,
            "not_a_test_of_genetic_fixation_or_phylogenetic_origin": True,
        }
    }


def panel_inputs(panel: dict) -> tuple[np.ndarray, np.ndarray, np.ndarray, float]:
    frame = panel["data"]
    missing = set(AXES) - set(frame)
    if missing:
        raise ValueError("Original five S1 climate source columns missing: " + str(missing))
    x = frame.loc[:, list(AXES)].to_numpy(dtype=float)
    if not np.isfinite(x).all() or not (x[:, 0] > 0).all():
        raise ValueError("Nonfinite or physically invalid five climate source")
    n = len(x)
    y = np.asarray(panel["y"], dtype=np.int8)
    if n != len(y) or set(np.unique(y)) != {0, 1}:
        raise ValueError("Species/posture cohort mismatch")
    sd = x.std(axis=0)
    if np.any(sd <= 1e-12):
        raise ValueError("Climate axis constant within nuclear panel")
    z = (x - x.mean(axis=0)) / sd
    k = int(y.sum())
    d = z[y == 1].mean(axis=0) - z[y == 0].mean(axis=0)
    weight = 2 * k * (n - k) / n
    return z, y, d, weight


def shrink_inverse(cov: np.ndarray, ridge: float = RIDGE) -> np.ndarray:
    if cov.ndim != 2 or cov.shape[0] != cov.shape[1]:
        raise ValueError("Covariance not square")
    if not np.isfinite(cov).all():
        raise ValueError("Nonfinite null covariance")
    diag = np.diag(np.diag(cov))
    if np.any(np.diag(diag) <= 1e-14):
        raise ValueError("Degenerate null environment axis")
    v = (1-ridge) * cov + ridge * diag
    return np.linalg.inv(v)


def score_vectors(deltas: np.ndarray, center: np.ndarray, inverse: np.ndarray) -> np.ndarray:
    d = np.asarray(deltas, dtype=float) - center
    if d.ndim == 1:
        d = d[None, :]
    return np.einsum("bi,ij,bj->b", d, inverse, d)


def mk_conditional_draws(panel: dict, z: np.ndarray, y: np.ndarray,
                         root: str, target: int, seed: int) -> tuple[np.ndarray, dict]:
    tree = panel["tree"]
    labels = {name: int(v) for name, v in zip(panel["tips"], y)}
    rate, loglike, near_bound = fitted_rate_given_root(tree, labels, root)
    rng = np.random.default_rng(seed)
    n = len(y)
    k = int(y.sum())
    out = []
    accepted = 0
    proposed = 0
    max_proposals = max(250000, target * 150)
    while accepted < target and proposed < max_proposals:
        size = min(4096, max_proposals - proposed)
        simulated = simulate_mk_tip_matrix(tree, rate, size, rng, root)
        if simulated.shape != (size, n):
            raise ValueError("Mk simulated original-tip order drift")
        hit = simulated[simulated.sum(axis=1) == k]
        if len(hit):
            # z is centered within panel, so mean(D)-mean(U) = n*sum(D*z)/(k*(n-k)).
            contrasts = hit @ z * (n/(k*(n-k)))
            out.append(contrasts)
            accepted += len(hit)
        proposed += size
    if not out:
        raise RuntimeError("No source-count conditioned Mk tips")
    data = np.concatenate(out, axis=0)[:target]
    if len(data) < target:
        raise RuntimeError("Insufficient Mk conditioned draws; do not present a P value")
    return data, {
        "n_real_tip": n, "n_U": n-k, "n_D": k,
        "root_mode": root, "Mk_rate": float(rate), "Mk_loglikelihood": float(loglike),
        "rate_near_optimizer_boundary": bool(near_bound),
        "proposals": proposed, "count_conditioned_acceptances": accepted,
        "conditioned_used": target,
    }


def multi_axis_nuclear(panel_a: dict, panel_b: dict, draws: int, seed: int) -> dict:
    panels = (panel_a, panel_b)
    vectors = [panel_inputs(p) for p in panels]
    weights = np.array([v[3] for v in vectors], dtype=float)
    if len(set(panel_a["data"].original_source_name) &
           set(panel_b["data"].original_source_name)):
        raise ValueError("Two real nuclear panels are not independent taxa")
    observed = np.average(np.vstack([v[2] for v in vectors]), axis=0, weights=weights)
    results = {}
    for j, root in enumerate(ROOTS):
        draws_by_region = []
        details = {}
        for i, (p, (z,y,d,w)) in enumerate(zip(panels, vectors)):
            outcomes, meta = mk_conditional_draws(
                p,z,y,root,target=draws,seed=seed+j*100+i
            )
            draws_by_region.append(outcomes)
            details[p["name"]] = {
                **meta, "observed_D_minus_U_in_within_region_SD_units": d.tolist(),
                "weight_2nDnU_over_n": float(w)
            }
        null = (weights[0]*draws_by_region[0]+weights[1]*draws_by_region[1])/sum(weights)
        center = null.mean(axis=0)
        cov = np.cov(null,rowvar=False,ddof=1)
        inverse = shrink_inverse(cov)
        primary = float(score_vectors(observed,center,inverse)[0])
        distribution = score_vectors(null,center,inverse)
        p_omnibus = float((1+np.count_nonzero(distribution >= primary-1e-12))/(draws+1))
        se = null.std(axis=0,ddof=1)
        scores = np.abs((observed-center)/se)
        null_scores = np.abs((null-center)/se)
        max_null = np.max(null_scores,axis=1)
        axis_tests = []
        for idx, name in enumerate(AXES):
            ax = null_scores[:,idx]
            raw = float((1+np.count_nonzero(ax>=scores[idx]-1e-12))/(draws+1))
            adj = float((1+np.count_nonzero(max_null>=scores[idx]-1e-12))/(draws+1))
            axis_tests.append({
                "axis":name,"observed_D_minus_U_in_pooled_within_region_SD":float(observed[idx]),
                "signed_null_center":float(center[idx]),"two_sided_null_rank_P":raw,
                "five_axis_maxT_adjusted_P":adj
            })
        results[root] = {
            "observed_omnibus_mahalanobis_T":primary,
            "omnibus_nondirectional_Mk_conditional_mc_P":p_omnibus,
            "candidate_axis_contrasts_are_secondary_not_preferred_signs":axis_tests,
            "region_details":details,
            "n_Mk_conditioned_draws":draws,
        }
    return {
        "label":"T2_5D_CURRENT_NICHE_BOTANICAL_STATE_OMNIBUS_NOT_TRANSITION_HISTORY",
        "axes_fixed_from_existing_S1":list(AXES),
        "axes_selection_not_independent_of_original_S1_workflow": True,
        "five_axis_signed_D_minus_U_predeclared":False,
        "omnibus_preferred_sign_predeclared":False,
        "pooled_region_standardized_5D_observed_D_minus_U":observed.tolist(),
        "panels":[{"name":p["name"],"n":len(v[1]),"D":int(v[1].sum()),
                   "U":int(len(v[1])-v[1].sum())} for p,v in zip(panels,vectors)],
        "different_clade_phylogenies_preserved_separately":True,
        "no_combined_fake_global_tree":True,
        "root_sensitivities":results,
        "no_ancestral_precipitation_or_independent_gains_estimated":True,
        "confounders_not_removed_by_Mk":"lineage geography, species sampling, sources, photo errors",
    }


def main() -> None:
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--original251",type=Path,required=True)
    parser.add_argument("--out",type=Path,required=True)
    parser.add_argument("--draws",type=int,default=10000)
    args=parser.parse_args()
    if args.draws<2000:
        raise ValueError("At least 2000 conditioned Mk draws required")
    source=pd.read_csv(args.original251)
    bridge=evidence_bridge(
        json.loads(SOURCE_S1S2.read_text(encoding="utf-8")),
        json.loads(SOURCE_S2.read_text(encoding="utf-8")),
        source)
    na,jp=source_panels(args.original251)
    for taxon,detail in bridge["S2"]["species"].items():
        name=taxon.replace(" ","_")
        detail["in_NA44_or_Japan12_true_nuclear_posture_cohort"] = any(
            name in {str(x).replace(" ","_") for x in p["data"].original_source_name}
            for p in (na,jp))
    multi=multi_axis_nuclear(na,jp,args.draws,SEED)
    result={
        "schema":"azami_v4_S1_S2_T2_five_environment_nondirectional_bridge_v1",
        "source_status":"COMPUTED_FROM_FROZEN_REAL_S1_S2_AND_INDEPENDENT_NUCLEAR_U_D",
        "S1_S2":bridge,
        "T2_56_authentic_nuclear_botanical_states":multi,
        "T3_cross_tribe": {
            "status":"NOT_INGESTED_NO_FALSE_CROSS_TRIBE_REPLICATION",
            "needs":"rooted nuclear phylogeny plus independently confirmed U and D during anthesis",
        },
        "interpretive_boundaries":{
            "S2_native_photo_records_subset_of_S1":True,
            "image_vertical_angle_not_botanical_gravity_orientation":True,
            "multiniche_test_uses_present_taxon_centroid_environments_not_origin_climate":True,
            "all_climate_effects_observational_not_causal":True,
            "five_environment_axes_were_selected_upstream_in_S1":True,
            "no_species_response_or_genetic_invariance_inferred_from_S2_null":True
        }
    }
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(result,ensure_ascii=False,indent=2,sort_keys=True)+"\n",
                        encoding="utf-8")
    a=result["T2_56_authentic_nuclear_botanical_states"]["root_sensitivities"]
    print("S1_S2_T2_MULTINICHE_RESULT="+json.dumps({
        "S1_between_BIO12":bridge["S1"]["BIO12_between_taxa_image_deg_per100mm"],
        "S1_within_BIO12":bridge["S1"]["BIO12_within_taxa_image_deg_per100mm"],
        "S2_introduced_rain_shift_mm": {
            k:v["introduced_minus_native_occurrence_BIO12_median_mm"]
            for k,v in bridge["S2"]["species"].items()},
        "T2_nuc56_multiaxis_omnibus_P_by_root":{
            k:v["omnibus_nondirectional_Mk_conditional_mc_P"] for k,v in a.items()
        },
        "T2_five_climate_signed_effect":multi["pooled_region_standardized_5D_observed_D_minus_U"],
        "T2_axis_order":list(AXES),
    },sort_keys=True))
if __name__=="__main__":
    main()
