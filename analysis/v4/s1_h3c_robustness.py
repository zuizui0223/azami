"""Post-positive H3c sensitivity: threshold, trait, sector and KDE bandwidth.

The original H3c positive result and its Holm-2 family remain frozen.
All combinations are exploratory diagnostics, not new independent discoveries.
Uses the same source-verified S1, its photographic QC and equal-cell sampling.
Export only pooled counts and effect estimates, never individual source rows.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from analysis.v4.s1_climate_trait_hypervolume import (
    SEED, load_source, geo_aggregate, standardise_cells
)
from analysis.v4.s1_h3_climate_analogue_alignment import (
    CLIMATE_ANALOG_Z, KDE_BANDWIDTH, RAREFIED_CELLS,
    eligible_taxa, _sample_taxa, matched_pairs,
    observed_and_permutations, kde_overlap_and_permutations
)

DRAWS = 12
PERMUTATIONS = 99
MIN_DRAWS = 6
MIN_PAIRS = 40
MIN_SECTORS = 1
EXPECTED_SOURCE_SHA = "342d78e9a8a960e8d9bdcf0bc2cdacd13bee48214a8b962419d5a71d5aef46ed"


def declared_scenarios(resolution: int):
    if resolution not in (1, 5):
        raise ValueError("UNREGISTERED_SPATIAL_RESOLUTION")
    scenarios = [
        {"label": "narrow_climate_0p7", "threshold": .7},
        {"label": "wide_climate_1p1", "threshold": 1.1},
        {"label": "without_angle", "omit_axis": 0},
        {"label": "without_chroma", "omit_axis": 1},
        {"label": "without_outline_aspect", "omit_axis": 2},
    ]
    if resolution == 1:
        scenarios += [
            {"label": "without_western_sector", "omit_sector": 0},
            {"label": "without_central_sector", "omit_sector": 1},
            {"label": "without_eastern_sector", "omit_sector": 2},
            {"label": "kde_h0p50", "bandwidth": .50},
            {"label": "kde_h1p00", "bandwidth": 1.00},
        ]
    return scenarios


def effect_against_null(observed, null):
    """One contrast across draws; repetitions are not independent taxa."""
    obs = np.asarray(observed, dtype=float)
    simulations = np.asarray(null, dtype=float)
    if obs.ndim != 1 or simulations.ndim != 2 or len(obs) != len(simulations):
        raise ValueError("INVALID_H3C_SENSITIVITY_NULL_SHAPE")
    average_obs = float(obs.mean())
    permuted = simulations.mean(axis=0)
    average_null = float(permuted.mean())
    advantage = average_obs - average_null
    count = int((np.abs(permuted - average_null) >= abs(advantage)).sum())
    return {
        "observed": average_obs,
        "null_mean": average_null,
        "observed_minus_null_mean": advantage,
        "unadjusted_permutation_p_diagnostic_only": (1 + count) / (1 + len(permuted)),
        "sd_across_repeated_cell_samplings": float(obs.std()),
        "original_H3c_Holm_2_not_recomputed": True,
    }


def examine(taxa, scenario, *, seed=SEED, draws=DRAWS, perms=PERMUTATIONS):
    subset = [x for x in taxa if x["sector"] != scenario.get("omit_sector", -1)]
    sectors = {int(x["sector"]) for x in subset}
    if len(subset) < 20 or len(sectors) < 2:
        return {"status": "insufficient_taxa_or_sectors", "n_taxa": len(subset)}
    threshold = float(scenario.get("threshold", CLIMATE_ANALOG_Z))
    bandwidth = float(scenario.get("bandwidth", KDE_BANDWIDTH))
    if threshold <= 0 or bandwidth <= 0:
        raise ValueError("INVALID_PREDECLARED_SENSITIVITY")
    required_sectors = MIN_SECTORS if "omit_sector" in scenario else 2
    draw_rng = np.random.default_rng(seed)
    perm_rng = np.random.default_rng(seed+314)
    distance, distance_null, overlap, overlap_null = [], [], [], []
    source_pairs, analog_pairs = [], []
    failures = {}
    for _ in range(draws):
        climate, phenotype, section = _sample_taxa(subset, draw_rng)
        if "omit_axis" in scenario:
            phenotype[:,:,int(scenario["omit_axis"])] = 0.0
        # Matching on scaled climate coordinates is exactly equivalent to
        # comparing the unscaled distance to the fixed scenario threshold.
        scaled = climate * (CLIMATE_ANALOG_Z / threshold)
        pair, meta = matched_pairs(scaled, phenotype, section)
        if (pair is None or meta["n_pairs"] < MIN_PAIRS
                or meta["n_sector_pair_types"] < required_sectors):
            code = meta["status"] if pair is None else "insufficient_supported_pairs"
            failures[code] = failures.get(code,0)+1
            continue
        obs, null = observed_and_permutations(pair, len(subset), perm_rng, reps=perms)
        density, density_null = kde_overlap_and_permutations(
            pair, phenotype, perm_rng, reps=perms, bandwidth=bandwidth)
        distance.append(obs)
        distance_null.append(null)
        overlap.append(density["conditional_minus_unconditional_mean"])
        overlap_null.append(density_null)
        source_pairs.append(meta["n_pairs"])
        analog_pairs.append(meta["n_analog_cell_pairs"])
    if len(distance) < min(MIN_DRAWS,draws):
        return {
            "status": "insufficient_repeated_analog_support",
            "n_taxa":len(subset),"successful_draws":len(distance),
            "attempted_draws":draws,"failures":failures,
            "required_supported_taxon_pairs":MIN_PAIRS,
            "source_biological_effect_not_estimated":True,
        }
    return {
        "status": "estimated_exploratory_robustness",
        "n_taxa":len(subset),
        "n_geographic_sectors":len(sectors),
        "successful_draws":len(distance),
        "requested_draws":draws,
        "mean_supported_cross_sector_taxon_pairs":float(np.mean(source_pairs)),
        "mean_supported_climate_analog_cell_pairs":float(np.mean(analog_pairs)),
        "distance_advantage":effect_against_null(distance,distance_null),
        "conditional_KDE_overlap":effect_against_null(overlap,overlap_null),
        "climate_distance_threshold":threshold,
        "kde_bandwidth":bandwidth,
        "omitted_trait_axis":scenario.get("omit_axis"),
        "omitted_geographic_sector":scenario.get("omit_sector"),
        "not_independent_confirmatory_hypothesis_testing":True,
    }


def execute(source: Path, receipt: Path, output: Path):
    frame,provenance = load_source(source,receipt)
    if (provenance["decoded_sha256"] != EXPECTED_SOURCE_SHA
        or provenance["source_observations"] != 313350
        or provenance["quality_eligible_matched_observations"] != 157822):
        raise ValueError("H3C_FROZEN_SOURCE_MISMATCH")
    result = {}
    for resolution in (1,5):
        cells, support = geo_aggregate(frame,resolution)
        if support["source_taxa_with_cell_and_photo_support"]<25:
            result[f"{resolution}degree"] = {
                "status":"insufficient_eligible_taxa","support":support}
            continue
        z,_ = standardise_cells(cells)
        taxa = eligible_taxa(z)
        conditions={}
        for index, scenario in enumerate(declared_scenarios(resolution)):
            conditions[scenario["label"]] = examine(
                taxa,scenario,seed=SEED+resolution*500+index*1000)
        result[f"{resolution}degree"] = {
            "status":"source_matched_sensitivity",
            "n_taxa":len(taxa),
            "conditions":conditions,
            "geographic_source_support": {
                "eligible_taxon_cell_pairs":support["eligible_taxon_cell_pairs"],
                "eligible_taxa":support["source_taxa_with_cell_and_photo_support"]
            }
        }
    public = {
        "schema":"azami_h3c_post_positive_robustness_20261029_v1",
        "status":"SOURCE_VERIFIED_EXPLORATORY_DOES_NOT_AMEND_ORIGINAL_H3C",
        "source_provenance":provenance,
        "fixed_earlier_H3c":{
            "original_positive_1degree_distance_advantage":0.040963661993538125,
            "original_positive_5degree_distance_advantage":0.0943310892554725,
            "original_24_rarefaction_199_permutation_Holm_2_unmodified":True
        },
        "diagnostic_design":{
            "n_cell_draws":DRAWS,"permutations_per_draw":PERMUTATIONS,
            "original_climate_threshold":CLIMATE_ANALOG_Z,
            "original_bandwidth":KDE_BANDWIDTH,
            "same_taxon_complete_3D_phenotype_vector_resampled_as_unit":True,
            "all_post_positive_conditions_displayed_no_result_selection":True
        },
        "results":result,
        "limits":{
            "post_positive_exploratory_not_independent_validation":True,
            "no_phylogenetic_convergence_test":True,
            "no_space_climate_causal_disentangling":True,
            "frozen_S1_80_S2_128_tests_untouched":True,
            "anonymous_aggregate_only":True
        },
    }
    out=Path(output);out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(public,indent=2,ensure_ascii=False,allow_nan=False)+"\n",
                   encoding="utf-8")
    print(json.dumps({
       "status":public["status"],"source":provenance,
       "results":{k:{
           "n_taxa":v.get("n_taxa",0),
           "conditions":{n:{
              "status":x["status"],
              "delta_distance_vs_null":x.get("distance_advantage",{}).get("observed_minus_null_mean"),
              "delta_overlap_vs_null":x.get("conditional_KDE_overlap",{}).get("observed_minus_null_mean"),
              "source_pairs":x.get("mean_supported_cross_sector_taxon_pairs")
           } for n,x in v.get("conditions",{}).items()}
       } for k,v in result.items()}
    },ensure_ascii=False,allow_nan=False),flush=True)
    return public


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--source",type=Path,required=True)
    ap.add_argument("--receipt",type=Path,required=True)
    ap.add_argument("--output",type=Path,required=True)
    a=ap.parse_args()
    execute(a.source,a.receipt,a.output)


if __name__=="__main__":
    main()
