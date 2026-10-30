"""H3c original-specification 6-seed rarefaction stability audit.

Reuses the unchanged H3c climate analogue and KDE estimator. Only the random
sequences of eight-cell rarefactions and within-taxon full-vector permutations
vary, with every seed fixed in the design before this additional source run.
No new independently significant results are claimed from the six runs.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from analysis.v4.s1_climate_trait_hypervolume import (
    load_source, geo_aggregate, standardise_cells
)
from analysis.v4.s1_h3_climate_analogue_alignment import analogue_analysis

SEEDS = {r: [20261029+10000*(i+1)+r for i in range(6)] for r in (1,5)}
DRAWS = 24
PERMS = 99
SOURCE_SHA = "342d78e9a8a960e8d9bdcf0bc2cdacd13bee48214a8b962419d5a71d5aef46ed"


def one_panel(cells, resolution, draws=DRAWS, permutations=PERMS, seeds=None):
    seeds=SEEDS[resolution] if seeds is None else list(seeds)
    rows=[]
    for seed in seeds:
        original=analogue_analysis(cells,draws=draws,
                                  permutations=permutations,seed=int(seed))
        if original["status"]!="source_matched_exploratory_analogue_alignment":
            rows.append({"status":original["status"],
                         "seed":int(seed),
                         "successful_rarefaction_draws":original.get("successful_rarefaction_draws",0)})
            continue
        kde=original["climate_analogue_conditional_KDE_L2_overlap"]
        rows.append({
            "status":"estimated",
            "seed":int(seed),
            "n_taxa":original["n_taxa"],
            "successful_rarefaction_draws":original["successful_rarefaction_draws"],
            "mean_cross_sector_taxon_pairs":original["mean_eligible_cross_sector_taxon_pairs"],
            "matched_distance_advantage_observed":original["primary_trait_distance_advantage_of_climate_matched_cells"],
            "matched_distance_advantage_null_mean":original["null_mean"],
            "matched_distance_advantage_observed_minus_null":float(
                original["primary_trait_distance_advantage_of_climate_matched_cells"]
                -original["null_mean"]),
            "matched_distance_permutation_p_unadjusted_diagnostic":
                original["two_sided_taxon_block_cell_permutation_p_exploratory"],
            "KDE_conditional_minus_unconditional_observed":
                kde["conditional_minus_all_cells_overlap"],
            "KDE_conditional_minus_unconditional_null_mean":kde["null_mean"],
            "KDE_conditional_overlap_observed_minus_null":float(
                kde["conditional_minus_all_cells_overlap"]-kde["null_mean"]),
            "KDE_overlap_permutation_p_unadjusted_diagnostic":
                kde["two_sided_taxon_block_permutation_p"],
        })
    valid=[x for x in rows if x["status"]=="estimated"]
    if len(valid)!=len(seeds):
        return {
            "status":"incomplete_not_all_seed_sequences_estimable",
            "n_seeds_registered":len(seeds),
            "n_valid_seeds":len(valid),
            "runs":rows,
            "not_a_negative_or_equivalence_test":True,
        }
    def summary(effect,pvalue):
        values=np.array([r[effect] for r in valid],float)
        ps=np.array([r[pvalue] for r in valid],float)
        return {
            "n_independent_rng_sequences_not_independent_source_samples":len(values),
            "n_positive_relative_to_permutation_null":int((values>0).sum()),
            "min":float(values.min()),
            "median":float(np.median(values)),
            "max":float(values.max()),
            "range_of_unadjusted_permutation_p_diagnostic":[float(ps.min()),float(ps.max())],
            "number_of_nominal_P_under_0p05_not_independent_confirmations":int((ps<.05).sum()),
        }
    return {
        "status":"estimated_exploratory_random_sequence_stability",
        "n_taxa":valid[0]["n_taxa"],
        "n_seeds_registered":len(seeds),
        "n_cell_draws_per_seed":draws,
        "n_permutations_per_draw":permutations,
        "summary_distance":summary(
            "matched_distance_advantage_observed_minus_null",
            "matched_distance_permutation_p_unadjusted_diagnostic"),
        "summary_KDE_overlap":summary(
            "KDE_conditional_overlap_observed_minus_null",
            "KDE_overlap_permutation_p_unadjusted_diagnostic"),
        "all_seed_results":rows,
        "random_sequences_are_monte_carlo_sensitivity_not_phylogenetic_replications":True,
    }


def execute(original_source, original_receipt, output):
    frame,provenance=load_source(original_source,original_receipt)
    if (provenance["decoded_sha256"]!=SOURCE_SHA
        or provenance["source_observations"]!=313350
        or provenance["quality_eligible_matched_observations"]!=157822):
        raise ValueError("H3C_SEED_AUDIT_MUST_USE_ORIGINAL_S1_SOURCE")
    output_panels={}
    for resolution in (1,5):
        cells,gate=geo_aggregate(frame,resolution)
        if gate["source_taxa_with_cell_and_photo_support"]<25:
            output_panels[f"{resolution}degree"]={"status":"inadequate_taxon_support"}
            continue
        cells,_=standardise_cells(cells)
        output_panels[f"{resolution}degree"]=one_panel(cells,resolution)
    output_obj={
        "schema":"azami_h3c_original_24_rarefaction_by_6_fixed_seeds_v1",
        "status":"SOURCE_VERIFIED_MONTE_CARLO_STABILITY_NOT_NEW_CONFIRMATION",
        "source_provenance":provenance,
        "original_H3c_frozen_reference":{
           "1degree_distance_advantage":0.040963661993538125,
           "5degree_distance_advantage":0.0943310892554725,
           "1degree_Holm2":.010,
           "5degree_Holm2":.010,
           "original_H3c_results_not_refitted_or_replaced":True,
        },
        "registered_sensitivity":{
           "seeds":{str(k):v for k,v in SEEDS.items()},
           "draws":DRAWS,"permutations":PERMS,
           "single_unchanged_climate_threshold":.9,
           "single_unchanged_KDE_bandwidth":.75,
           "no_outcome_selected_hypotheses_added":True,
        },
        "panels":output_panels,
        "interpretation":{
           "repeated_rarefactions_not_independent_taxa_or_studies":True,
           "nominal_P_across_seeds_not_combined_into_new_Holm_family":True,
           "original_S1_80_and_S2_128_unchanged":True,
           "no_convergence_adaptation_or_causal_climate_test":True,
           "source_identifiers_coordinates_photos_absent":True,
        }
    }
    dest=Path(output);dest.parent.mkdir(parents=True,exist_ok=True)
    dest.write_text(json.dumps(output_obj,indent=2,allow_nan=False,
                               ensure_ascii=False)+"\n",encoding="utf-8")
    print(json.dumps({
        "status":output_obj["status"],
        "source":provenance,
        "panels":{k:{"status":v["status"],
                      "n_taxa":v.get("n_taxa"),
                      "distance":v.get("summary_distance"),
                      "overlap":v.get("summary_KDE_overlap")}
                  for k,v in output_panels.items()}
    },allow_nan=False,ensure_ascii=False),flush=True)
    return output_obj


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--source",type=Path,required=True)
    parser.add_argument("--receipt",type=Path,required=True)
    parser.add_argument("--output",type=Path,required=True)
    args=parser.parse_args()
    execute(args.source,args.receipt,args.output)


if __name__=="__main__":
    main()
