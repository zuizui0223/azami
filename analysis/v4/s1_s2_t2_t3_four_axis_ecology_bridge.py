"""Verified S1 -> S2 -> T2 -> T3 ECOLOGICAL four-climate research continuum.

This is a CROSS-STAGE SYNTHESIS of four existing real source artifacts
(not a fifth observational dataset or an independent replication).
The one biological question is whether among-Cirsium species
BIO18–continuous-image-angle association is observed (S1), carries over
to native/introduced comparisons (S2), survives genuinely nuclear tip
sampling and phylogenetic correction (T2), and recurs in other tribes
using the SAME validated photo angle (T3).

Never imply a T3 angle from a photo URL or a nuclear ITS linked gene tree,
or infer gravitational downward U/D from EXIF vertical angle.
"""
from __future__ import annotations
import argparse,hashlib,json
from pathlib import Path

SHA={
 "S1_S2":"080c3c66228516834f62a32e8a7849199ff6580818030de19cc038b4c73a45bb",
 "T2_true_nuclear":"abab1bdf1e5ea74fac207b49a1dc3a277bce13a71ac5edc328720f27819b75f4",
 "T2_BIO18_gap":"894997bf2faa371b9664b55564c8be527aa699bfd7f379a95ceafe4dd451510b",
 "T3_extant_public_photos":"ffb90a009b031522480f3b5618db2c7ad70050b993e6eee70f0a428e05cd20d2",
}
AXES=("chelsa_bio18","chelsa_bio05","chelsa_rsds_mean","elevation_m")
SPECIES=("Cirsium arvense","Cirsium vulgare")
PANELS=("NA42_genuine_published_nuclear","Japan12_genuine_176_loci_nuclear")

def source_file(path,expected):
    raw=Path(path).read_bytes()
    digest=hashlib.sha256(raw).hexdigest()
    if digest!=SHA[expected]:
        raise ValueError("Original stage result SHA256 mismatch: "+expected)
    return json.loads(raw)

def synthesis(s,t,g,x):
    """Everything in this output is a deterministic transformation of sources."""
    if s["status"]!="COMPLETE_REAL_ORIGINAL_S1_AND_S2_SITE_EXPOSURES_FOUR_NEW_AXES" or tuple(s["variables"])!=AXES:
        raise ValueError("The current S1/S2 4-axis original result is not verified")
    if tuple(t["new_four_axes"])!=AXES or t["original_S1_native_image_universe"]!=170066:
        raise ValueError("T2 must use exactly same four original S1/S2 climate axes")
    if x["original_genbank_ITS_biological_taxa"]!=131 or tuple(x["four_axes_exactly_S1_S2"])!=AXES:
        raise ValueError("T3 source must use exact four S1/S2 environment exposures")
    if g["n_wet_global_species_MISSING_from_two_true_trees"]!=55:
        raise ValueError("Original true nuclear wet-tail sampling bias unverified")
    if t["no_plastid_or_flora_U_D_or_failed_111_Q20_tree"] is not True:
        raise ValueError("T2 chosen genuine nuclear species tree scope drifted")
    if x["n_actual_qualified_T3_S1_equivalent_continuous_head_angle_measurements"]!=0:
        raise ValueError("T3 photography phenotyping result changed; synthesize actual validated angles separately")

    y=s["S1"]["between_taxon"]["min_photo_1"]
    b=y["fit"]["coefficients"]["chelsa_bio18"]
    l,h=y["fit"]["HC3_95CI"]["chelsa_bio18"]
    sig=(b>0 and l>0)
    if y["n_source_names"]!=251 or not sig:
        raise ValueError("Original verified S1 BIO18 head-angle pattern lost")
    y10=s["S1"]["between_taxon"]["min_photo_10"]
    if y10["n_source_names"]!=149 or y10["fit"]["HC3_95CI"]["chelsa_bio18"][0]<=0:
        raise ValueError("Same-trait original S1 ten-photos support changed")
    intr=s["S1"]["within_taxon"]
    if intr["n_original_source_images_complete"]!=169020 or intr["n_source_taxon_clusters"]!=251:
        raise ValueError("Unmodified S1 within-species original source lost")
    within_ci=intr["cluster_taxon_95CI"]["chelsa_bio18"]
    if not within_ci[0]<=0<=within_ci[1]:
        raise ValueError("S1 within-taxonomic BIO18 uncertainty status changed")
    s2={}
    for species in SPECIES:
        area=s["S2"]["species"][species]
        if not all(area[k]["status"]=="FOUR_AXIS_S2_SOURCE_REAL_WCVP_RANGE" for k in ("native","introduced")):
            raise ValueError("S2 independent introduction comparison became unsupported")
        rows={}
        for k in ("native","introduced"):
            m=area[k]["original_S2_four_axis_cell_model"]
            if m["status"]!="S2_TRUE_IMAGE_ANGLE_FOUR_AXIS_CELL_CLUSTER_REGRESSION":
                raise ValueError("One of four S2 photo models did not fit")
            coef=m["coefficients_degrees_per_within_range_SD"]["chelsa_bio18"]
            ci=m["cell_cluster_95CI"]["chelsa_bio18"]
            rows[k]={"n_4climate_photos":area[k]["n_photo_four_axis_complete"],
                     "BIO18_photo_head_angle_coef_deg_per_within_range_SD":coef,
                     "BIO18_photo_head_angle_cluster_CI95":ci,
                     "supports_BIO18_positive_effect":bool(ci[0]>0)}
        shift=area["introduced_minus_native_full_wild_site_median"]["chelsa_bio18"]
        s2[species]={"introduced_minus_native_wild_BIO18_median_mm":shift,
                     "no_S2_photo_BIO18_positive_effect_confirmed":
                           not any(v["supports_BIO18_positive_effect"] for v in rows.values()),
                     "by_WCVP_range":rows}
    t2={}
    for name in PANELS:
        group=t["nuclear_panels"][name]["min_photos_1"]
        record=group["PGLS"]["estimated_Pagel_lambda"]
        ci=record["conditional_t_95CI_degrees_per_SD"][1]
        coef=record["coefficients_degrees_per_SD"][1]
        t2[name]={"original_S1_source_photographs":group["n_original_measured_photo_observations"],
                  "n_authentic_multilocus_nuclear_species":group["n_authentic_source_taxa"],
                  "BIO18_slope_deg_per_within_tree_SD":coef,"conditional_CI95":ci,
                  "n_true_nuclear_species_in_global_S1_BIO18_top_quartile":g["genuine_nuclear_species_tree_panels"][name]["n_global_S1_BIO18_wet_quartile_species"],
                  "BIO18_on_nuclear_panel_confirmed_positive":bool(ci[0]>0)}
    if any(v["BIO18_on_nuclear_panel_confirmed_positive"] for v in t2.values()):
        raise ValueError("Genuine S1 BIO18 T2 result now needs fresh interpretation")
    if x["n_source_species_with_photos"]!=48 or x["n_records_4_environment_complete"]!=720:
        raise ValueError("Real independently sourced T3 photography-climate census drifted")
    return {
      "schema":"azami_S1_S2_T2_T3_one_4axis_ecological_trajectory_20261024_v1",
      "status":"SOURCE_VERIFIED_S1_POSITIVE_S2_UNCONFIRMED_T2_SAMPLING_LIMIT_T3_PHOTO_QC_PENDING",
      "shared_exposures":list(AXES),
      "shared_trait":"S1_EXIF_camera_vertical_head_axis_angle_0_to_180_degrees",
      "central_test":"When an among-Cirsium species BIO18 photo-angle pattern is detectable, is it reproduced within populations/introductions, nuclear clades, and other Asteraceae tribes on the SAME quantitative trait?",
      "S1_between_species":{
         "n_photos_original_new4_complete":169020,
         "n_taxon_source_names":251,
         "BIO18_coeff_deg_per_all251_taxon_SD":b,
         "HC3_95CI": [l,h],
         "p":y["fit"]["HC3_two_sided_z_p"]["chelsa_bio18"],
         "n10_photo_species":149,
         "n10_coeff":y10["fit"]["coefficients"]["chelsa_bio18"],
         "n10_HC3_95CI":y10["fit"]["HC3_95CI"]["chelsa_bio18"],
         "positive_pattern_supported_by_S1_source":sig,
      },
      "S1_within_taxon":{
         "n_source_photo_records":169020,"n_original_taxon_clusters":251,
         "BIO18_coef":intr["coefficients_degrees_per_full_S1_photo_SD"]["chelsa_bio18"],
         "CI95":within_ci,"positive_within_taxon_effect_confirmed":False
      },
      "S2_two_species_introduction_gradient":{
         "native_photos_overlap_with_S1_not_independent_replication":True,
         "original_focal_all_wild_records":297622,
         "original_total_focal_photo_records":173737,
         "species":s2,
      },
      "T2_true_multilocus_nuclear_two_disjoint_geographic_panels":{
         "panels":t2,
         "n_global_S1_BIO18_top_quartile_source_taxa":63,
         "n_true_nuclear_top_quartile_present":8,
         "n_true_nuclear_top_quartile_missing":55,
         "North_America_tree_wet_source_taxa":0,
         "Japan_tree_wet_source_taxa":8,
         "Japan_top_wet_taxa_with_one_S1_photo":4,
         "conclusion":"No positive BIO18 nuclear PGLS established; critical nonrandom clade and wet-gradient sample support, not a full genus phylogenetic null.",
         "true_nuclear_species_trees_not_joined_into_one_global_tree":True,
      },
      "T3_independent_other_tribe_four_climate_photo_feasibility":{
          "n_original_GenBank_ITS_source_taxa":131,
          "n_exact_GBIF_taxon_matches":x["n_GBIF_exact_taxon_matches"],
          "n_photo_species":x["n_source_species_with_photos"],
          "n_geolocated_public_photo_records":x["n_georeferenced_photo_records_before_all4_climate"],
          "n_photo_records_with_all4_original_environment_pixels":720,
          "n_species_with_all4_climates":x["n_species_4_environment_complete"],
          "n_species_ready_for_MANUAL_pedicel_angle_QC":x["n_species_qualify_for_MANUAL_photo_angle_audit"],
          "n_genera_with_photography":x["n_genera_with_photo_records"],
          "Ligularia_species_qualified_for_manual_QC":x["by_genus"]["Ligularia"]["n_taxa_geographically_replicated_photo_pilot_candidate"],
          "measured_S1_equivalent_0_180_head_angles":0,
          "trait_environment_cross_tribe_replication_estimated":False,
          "true_multilocus_phylogeny_for_this_other_tribe":False,
      },
      "ecological_boundary":{
         "do_not_claim_independent_BIO18_clade_recurrence":True,
         "do_not_claim_T3_external_4axis_angletest_from_location_only":True,
         "do_not_infer_independent_genetic_adaptation_or_nodding":True,
         "true_gravity_angle_not_measured":True,
         "published_ITS_ETS_rDNA_linked_genealogy_not_validated_species_tree":True,
         "no_botanical_flora_U_D_substitution":True,
         "user_selected_posthoc_four_axes_instead_of_original_five":True,
      },
      "next_empirical_gate":"T2 extend ORIGINAL photographed wet BIO18 species across validated many-nuclear-locus species tree; T3 validate cross-genus capitulum+pedicel keypoint measurement on original licensed photos before angle~BIO18/5/RSDS/elevation and cross-tribe inference.",
      "sources_sha256":SHA,
      "no_new_unreproduced_numeric_estimates":True,
    }

def main():
    p=argparse.ArgumentParser()
    for k in ("s1-s2","t2","t2-gap","t3","out"):
        p.add_argument("--"+k,type=Path,required=True)
    a=p.parse_args()
    s=source_file(a.s1_s2,"S1_S2")
    t=source_file(a.t2,"T2_true_nuclear")
    g=source_file(a.t2_gap,"T2_BIO18_gap")
    x=source_file(a.t3,"T3_extant_public_photos")
    report=synthesis(s,t,g,x)
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(report,ensure_ascii=False,sort_keys=True,indent=2)+"\n")
    print("S1_S2_T2_T3_FOUR_AXIS_ONE_RESEARCH_TRAJECTORY="+json.dumps({
      "S1_BIO18_positive":report["S1_between_species"]["positive_pattern_supported_by_S1_source"],
      "S2_native_intro_BIO18_positive_confirmed":False,
      "T2_nuclear_wet_missing":55,
      "T3_n_other_tribe_4axis_photos":720,
      "T3_photo_angle_measured":0},sort_keys=True))
if __name__=="__main__":main()
