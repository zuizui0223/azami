"""S1-continuous-only V4 primary bridge: S1 -> S2 -> T2 -> T3 readiness.

The response is ALWAYS the original 0–180 degree camera/image-vertical
head-presentation angle. Never turn it into binary floristic U/D, infer true
gravity tilt, or count independent nodding origins.

Recalculate source-matched equal-taxon FIVE-CHELSA contrasts on all251 and
an authentic 109-tip GenBank ITS genealogy; check full S1 observation-derived
native photos against S2 native focal species. T3 is blocked until separate
Asteraceae photo measurements of the SAME definition exist.

This reuses prior source-authenticated files; it does not claim to refit a
global nuclear species tree or rerun the original 244-concept observation model.
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np

from analysis.v4.t2_s1_109_vs251_phylogenetic_sampling_explanation import (
    FIVE, EXPECTED_251_SHA256, EXPECTED_109_ML_SHA256, evaluate, read_source,
)
from analysis.v4.t2_s1_s2_multiniche_phylo_bridge import SOURCE_S1S2, SOURCE_S2

S2_FOCALS = ("Cirsium arvense", "Cirsium vulgare")
ENDPOINT = "orientation_image_vertical_angle"
T3_MIN_N_SOURCE_OBSERVATIONS_PER_TAXON = 10
T3_MIN_INDEPENDENT_TRIBES = 2


def crosscheck_s1_s2(source, original_s1_s2, original_s2):
    if original_s1_s2["S1"]["n_observations"] != 170066:
        raise ValueError("Original S1 photo observation universe changed")
    if original_s1_s2["S1"]["n_accepted_keys"] != 244:
        raise ValueError("Original S1 accepted taxon keys changed")
    byname = {x["original_source_name"].replace("_", " "): x for x in source}
    if len(byname) != 251 or set(original_s2["species"]) != set(S2_FOCALS):
        raise ValueError("Original S1 source-name or S2 two-species cohort changed")
    records = {}
    for name in S2_FOCALS:
        x = byname[name]
        sr = original_s1_s2["S2"]["same_species_comparisons"][name]
        so = original_s2["species"][name]
        native_s1 = int(x["original_photo_count"])
        n_native = int(sr["native"]["original_orientation_eligible_observations"])
        if native_s1 != n_native or n_native != int(so["native"]["original_orientation"]):
            raise ValueError("S2 native photo angles do not match S1 original native observations")
        if sr["differences"]["introduced_minus_native_median_BIO12_mm"] <= 0:
            raise ValueError("S2 wetter introduced occurrence median changed")
        for status in ("native", "introduced"):
            e = so[status]
            if e["cluster_ci95"][0] > e["bio12_deg_100mm"] or e["bio12_deg_100mm"] > e["cluster_ci95"][1]:
                raise ValueError("S2 native/introduced five-climate angle gradient drift")
        records[name] = {
            "same_native_original_image_observations_as_S1":native_s1,
            "introduced_image_observations":int(so["introduced"]["original_orientation"]),
            "native_source_photo_mean_angle_deg":float(x["mean_image_angle_deg"]),
            "introduced_minus_native_BIO12_occurrence_median_mm":
                float(sr["differences"]["introduced_minus_native_median_BIO12_mm"]),
            "introduced_minus_native_image_angle_median_deg":
                float(sr["differences"]["introduced_minus_native_photo_angle_median_degrees"]),
            "native_BIO12_slope_image_deg_per100mm":
                float(so["native"]["bio12_deg_100mm"]),
            "introduced_BIO12_slope_image_deg_per100mm":
                float(so["introduced"]["bio12_deg_100mm"]),
            "native_and_introduced_slopes_cluster95_cover_zero":
                all(so[status]["cluster_ci95"][0] <= 0 <= so[status]["cluster_ci95"][1]
                    for status in ("native","introduced")),
            "native_model_predictive_transfer_R2":
                float(so["predictive_transfer"]["native_trained_skill_relative_to_native_constant"]),
        }
    return {
        "S1_original_source_focal_taxa_are_same_as_S2_native":True,
        "not_two_independent_replicates":True,
        "S2_not_evidence_for_evolutionary_fixation":True,
        "species":records,
    }


def validate_t3_continuous_source_manifest(rows):
    """A source-level admission gate, NEVER construction of another-genus photos.

    Synthetic unit tests of this guard are not biological evidence.
    """
    mandatory={"taxon_name","tribe","n_original_photo_observations",
               "measurement_endpoint","camera_vertical_EXIF_aligned",
               "mean_image_angle_deg","source_photo_digest",
               "source_climate_digest","original_phenotype_observed"}
    if not rows:
        return {"status":"BLOCKED_NO_MEASURED_EXTERNAL_TRIBE_PHOTOS",
                "reason":"No T3 observations using exact S1 image measurement"}
    if any(set(x)!=mandatory for x in rows):
        raise ValueError("T3 lacks a source-backed S1-equivalent measurement manifest")
    if len({x["taxon_name"] for x in rows})!=len(rows):
        raise ValueError("T3 duplicated a taxon concept")
    tribes=set()
    valid=0
    for x in rows:
        if (x["measurement_endpoint"]!=ENDPOINT
            or x["camera_vertical_EXIF_aligned"]!="True"
            or x["original_phenotype_observed"]!="True"
            or not x["source_photo_digest"]
            or not x["source_climate_digest"]):
            raise ValueError("T3 cannot use categorical flora states or missing photo provenance")
        value=float(x["mean_image_angle_deg"])
        n=int(x["n_original_photo_observations"])
        if not np.isfinite(value) or not 0<=value<=180 or n<1:
            raise ValueError("Invalid independently measured photographic angle or count")
        if n>=T3_MIN_N_SOURCE_OBSERVATIONS_PER_TAXON:
            valid+=1
            tribes.add(x["tribe"])
    return {
        "status":("ELIGIBLE_FOR_T3_TAXON_COVERAGE_REVIEW"
                  if len(tribes)>=T3_MIN_INDEPENDENT_TRIBES
                  else "BLOCKED_INSUFFICIENT_TRIBE_MEASUREMENT_COVERAGE"),
        "n_photographically_measured_taxa":len(rows),
        "n_taxa_at_min10_photographs":valid,
        "n_distinct_tribes_at_min10":len(tribes),
        "note":"This gate does not certify nuclear trees, gravity tilt or climate response",
    }


def run(original251,its109):
    source=read_source(original251)
    results=evaluate(original251,its109)
    receipt=json.loads(SOURCE_S1S2.read_text())
    transfer=json.loads(SOURCE_S2.read_text())
    s2=crosscheck_s1_s2(source,receipt,transfer)
    cohorts=results["cohorts"]
    full=cohorts["S1_full251"]
    geno=cohorts["GenBank_ITS109"]
    missing=cohorts["unsequenced_S1_142"]
    def digest_summary(row):
        fit=row["equal_taxon_weight_HC3_photo_angle_vs_BIO12"]["mean_image_angle_deg"]
        return {
            "n_taxon_names":row["n_original_source_names"],
            "n_source_photographs":row["n_original_owner_photo_angle_observations"],
            "PHOTO_mean_BIO12_slope_deg_per100mm":fit["BIO12_effect_original_photo_angle_degrees_per100mm"],
            "PHOTO_mean_BIO12_HC3_95CI":fit["HC3_95percent_normal_ci"],
            "model_is_equal_taxon_weight_not_S1_native_original_mixed":True,
        }
    w=results["ITS_nuclear_source_coverage_by_wet_quartile"]
    supports=results["ITS_ML_branch_support"]
    assert geno["n_original_owner_photo_angle_observations"]==162146
    if not results["claim_stop"]["S1_109_species_are_NOT_random_sample_of_ALL_251"]:
        raise ValueError("GenBank ITS selection limitations lost")
    native=receipt["S1"]
    if (native["native_244_between_species_BIO12_degrees_per_100mm"]<=0
        or not native["between_95pct_interval"][0]<
                 native["native_244_between_species_BIO12_degrees_per_100mm"]<
                 native["between_95pct_interval"][1]):
        raise ValueError("Original S1 observer-level positive between pattern changed")
    return {
        "schema":"azami_V4_primary_S1_continuous_photo_only_S1_S2_T2_T3_v1",
        "status":"EXECUTED_REAL_SOURCE_S1_PHOTO_NATIVE_ITS109_AND_S2_JOINT_AUDIT",
        "endpoint":ENDPOINT,
        "endpoint_units":"degrees_relative_to_EXIF_oriented_IMAGE_vertical_0_to_180",
        "never_botanical_U_D":True,
        "U_D_flora_codes_not_used_in_any_inference":True,
        "no_fake_ancestor_nodding_gains":True,
        "S1_frozen_original_model":{
            "n_original_native_images":170066,
            "n_original_WCVP_accepted_keys":244,
            "between_BIO12_image_deg_per100mm":
                float(native["native_244_between_species_BIO12_degrees_per_100mm"]),
            "between_95CI":native["between_95pct_interval"],
            "within_BIO12_image_deg_per100mm":
                float(native["source_native_244_within_species_BIO12_degrees_per_100mm"]),
            "within_95CI":native["within_95pct_interval"],
            "original_5_axes":list(FIVE),
            "flexible_spatial_adjustment_no_independent_rain_causal_impact":True,
        },
        "S2_same_continuous_endpoint":s2,
        "T2_real_ITS_gene_tree_continuous_source":{
            "source_cohorts":{
                "all251":digest_summary(full),
                "ITS109":digest_summary(geno),
                "without_ITS142":digest_summary(missing),
            },
            "n_gene_tree_taxa":109,
            "n_original_S1_photo_observations_covered":162146,
            "fraction_original_S1_taxa_covered":109/251,
            "fraction_original_S1_photos_covered":162146/170066,
            "source_BIO12_wettest_quartile_n":w["n_S1_top_quartile"],
            "ITS_taxa_of_wettest_quartile_n":w["n_true_ITS_taxa_among_wettest_quartile"],
            "ITS_internal_supported_nodes":supports["supported_internal_shalrt_and_UFBoot_count"],
            "ITS_internals_with_support_labels":supports["original_ML_support_labeled_internal_branches"],
            "ITS_is_multicopy_single_locus_gene_tree_not_species_tree":True,
            "original_genuine_multilocus_NA42_complementary_existing_run":
                "https://github.com/zuizui0223/azami/blob/codex/v4-native-integration-workflow/.github/workflows/v4-t2-na42-nuclear-bm-ou.yml",
            "Q20_110_reference_copy_multi_locus_original111_QC_failed_and_not_promoted":True,
            "photo_angle_phylogenetic_lambda_exploratory_existing_109":0.0,
            "logBIO12_lambda_nugget5pct_exploratory_existing_109":0.6807,
            "lambda_and_BM_OU_are_prior_real_analyses_not_recomputed_here":True,
            "source_coverage_is_nonrandom_and_same_species_weights_not_observation_model":True,
        },
        "T3":validate_t3_continuous_source_manifest([]),
        "statistical_limits":{
            "S1_S2_nonindependent_native_original_photos":True,
            "T2_ITS_coverage_biased_against_wet_S1_tail":True,
            "T3_current_131_ITS_ETS_has_no_S1_comparable_photo_angles":True,
            "no_reinterpretation_as_gravitational_nodding":True,
            "no_U_to_D_Mk_transition_model":True,
            "no_historical_climate_or_fitness_inferred":True,
        },
        "source_sha256":{"original_S1_all251":EXPECTED_251_SHA256,
                         "real_GenBank_ITS109_tree":EXPECTED_109_ML_SHA256},
    }


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--original251",type=Path,required=True)
    p.add_argument("--genbank-its109-tree",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    o=run(a.original251,a.genbank_its109_tree)
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(o,indent=2,ensure_ascii=False,sort_keys=True)+"\n")
    print("S1_ONLY_T2_T3_CONTINUOUS="+json.dumps({
        "full251_BIO12":o["T2_real_ITS_gene_tree_continuous_source"]["source_cohorts"]["all251"]["PHOTO_mean_BIO12_slope_deg_per100mm"],
        "ITS109_BIO12":o["T2_real_ITS_gene_tree_continuous_source"]["source_cohorts"]["ITS109"]["PHOTO_mean_BIO12_slope_deg_per100mm"],
        "photos_ITS109":o["T2_real_ITS_gene_tree_continuous_source"]["n_original_S1_photo_observations_covered"],
        "T3":o["T3"]["status"],
    },sort_keys=True))


if __name__=="__main__":
    main()
