"""Revalidate and freeze V4 S1/S2 4var-BIO18 Holm & observable LEGO.

Use the exact original-code commit and GitHub Actions artifact; this script
requires only standard-library Python and published aggregate ZIP bytes.
It fails closed if the source bytes, variables, significant endpoints,
biological hypothesis family, or matched-cohort LEGO statistics drift.

A frozen scientific *interpretation* is not protected Git history, genetic
modularity, causal adaptation, or a permission to claim molecular G123.
"""
from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import math
from pathlib import Path
from statistics import median
import zipfile

SOURCE_REPO = "zuizui0223/azami"
SOURCE_COMMIT = "d2c45b00018a0c04e7c69a29857207dc1ac06d0c"
SOURCE_RUN = 40062126006
SOURCE_ARTIFACT = 12547876704
ARTIFACT_ZIP_SHA256 = "6e081dc081bb8ee2311a05ac5b7fcddd3f6309be2fa909e87d432789f4e0e205"
PINNED_MEMBERS = {
 "V4_S1_S2_four_vs_vif5_vif10_alltraits.json":
     "3738292e0a64da183dfe219c854c16cda199a7f9769f4928da8ab07657d60b70",
 "S1_S2_extreme_LEGO_results.json":
     "bb01d3b5da03148d84b2521c9c92cf14f7f6177435e145118f6ff6df8846733f",
 "V4_S1_S2_four_model_comparison.json":
     "d466f8931b4ed6ef55f7e7ce24cceedd00fcd002dc88cf83487094fe2afe9e7e",
}
FINAL_ENV = ("chelsa_bio18","chelsa_bio05",
             "chelsa_rsds_mean","topo_elevation")
FINAL_MODEL = "four_bio18"
SOURCE_CORE_N = 156683
SOURCE_CORE_TAXA = 174
MODULES = {
    "angle":0, "lightness":1,"chroma":1,"hue":1,
    "aspect":2,"circularity":2,"solidity":2,"width_profile":2,
}
EXPECTED_HOLM = {
    "S1": {
        ("angle","orientation_image_vertical_angle",
         "chelsa_bio18","between_equal_taxa"):
          (168907,0.04544684307686103,0.00011507938226118898,
           0.009206350580895119),
    },
    "S2":{
        ("Cirsium arvense/introduced/lightness",
         "corolla_lab_lightness","topo_elevation","within_species_range"):
          (26171,-0.17373187423223405,1.0252616543029317e-05,
           0.0013123349175077526),
        ("Cirsium vulgare/introduced/chroma",
         "corolla_lab_chroma","chelsa_rsds_mean","within_species_range"):
          (27200,0.1345037874025455,0.00011470971522737719,
           0.014568133833876904),
        ("Cirsium vulgare/introduced/lightness",
         "corolla_lab_lightness","chelsa_rsds_mean","within_species_range"):
          (27200,-0.14280073287481068,0.00015281061713213603,
           0.019254137758649138),
    }
}

def require(flag,reason):
    if not flag:
        raise ValueError(reason)

def close(a,b,reason,rel=1.e-9):
    require(math.isclose(float(a),float(b),rel_tol=rel,abs_tol=1.e-12),reason)

def holm(p):
    require(all(math.isfinite(x) and 0<=x<=1 for x in p),
            "INVALID_HOLM_INPUT")
    n=len(p); order=sorted(range(n),key=lambda i:p[i]); out=[None]*n
    max_so_far=0.
    for rank,i in enumerate(order):
        max_so_far=max(max_so_far,(n-rank)*p[i])
        out[i]=min(1.,max_so_far)
    return out

def exact_label_reference(rows):
    names=tuple(MODULES)
    require(len(rows)==28,"LEGO_28_EDGES_MISSING")
    require({tuple(sorted((r["left"],r["right"]))) for r in rows}
        =={tuple(sorted(c)) for c in itertools.combinations(names,2)},
        "LEGO_28_PAIR_REGISTER_DRIFT")
    def contrast(mapping,column):
        same=[r[column] for r in rows
              if mapping[r["left"]]==mapping[r["right"]]]
        cross=[r[column] for r in rows
               if mapping[r["left"]]!=mapping[r["right"]]]
        require((len(same),len(cross))==(9,19),
                "LEGO_BLOCK_SIZE_NOT_1_3_4")
        return median(same)-median(cross)
    labelings=[]
    for special in names:
        remainder=[v for v in names if v!=special]
        for colour in itertools.combinations(remainder,3):
            labelings.append({v:(
                0 if v==special else 1 if v in colour else 2) for v in names})
    require(len(labelings)==280,"LEGO_EXACT_REFERENCE_NOT_280")
    result={}
    for name in ("within_rv","among_rv"):
        o=contrast(MODULES,name)
        null=[contrast(m,name) for m in labelings]
        count=sum(v>=o-1.e-12 for v in null)
        result[name]={"observed_same_minus_cross_median":o,
                      "n_allocations":280,
                      "allocations_at_least_observed":count,
                      "conditional_label_reference_fraction":count/280,
                      "NOT_evolutionary_or_genetic_inferential_p_value":True}
    require(result["within_rv"]["allocations_at_least_observed"]==1,
            "WITHIN_LABEL_REFERENCE_DRIFT")
    require(result["among_rv"]["allocations_at_least_observed"]==3,
            "AMONG_LABEL_REFERENCE_DRIFT")
    return result

def check_families(d):
    result={}
    for study in ("S1","S2"):
        key=study+"/"+FINAL_MODEL
        family=d["multiplicity"][key]
        require(family["n_tests"]==({"S1":80,"S2":128}[study]),
                "HOLM_DISCOVERY_FAMILY_CHANGED "+study)
        rows=family["rows"]
        require(len(rows)==family["n_tests"],
                "HOLM_ROW_COUNT_CHANGED")
        require(all(r["coordinate"] not in ("corolla_hue_sin","corolla_hue_cos")
                    for r in rows),"HUE_DOUBLE_COUNTED")
        check=holm([r["p_original"] for r in rows])
        for i,r in enumerate(rows):
            close(check[i],r["p_Holm"],"INDEPENDENT_HOLM_DISAGREEMENT")
        passed={}
        for r in rows:
            if r["p_Holm"]>=.05:
                continue
            key=(r["group"],r["coordinate"],r["environment"],r["scale"])
            require(key not in passed,"DUPLICATE_SIGNIFICANT_ENDPOINT")
            expected=EXPECTED_HOLM[study].get(key)
            require(expected is not None,"UNREGISTERED_HOLM_POSITIVE "+str(key))
            n,beta,p,p_adj=expected
            require(r["n"]==n and r["trait_cohort_vif_pass"],
                    "SIGNIFICANT_SOURCE_OBSERVATIONS_CHANGED")
            close(r["beta_std"],beta,"SIGNIFICANT_BETA_CHANGED")
            close(r["p_original"],p,"SIGNIFICANT_RAW_P_CHANGED")
            close(r["p_Holm"],p_adj,"SIGNIFICANT_HOLM_CHANGED")
            require(r["ci95_wald"] is not None
                    and r["ci95_wald"][0]<beta<r["ci95_wald"][1],
                    "SIGNIFICANT_95CI_MISSING")
            passed[key]={"n":n,"beta_std":beta,
                         "ci95_wald":r["ci95_wald"],
                         "p_original":p,"p_BH":r["p_BH"],
                         "p_Holm":p_adj}
        require(set(passed)==set(EXPECTED_HOLM[study]),
                "EXPECTED_HOLM_FINDING_MISSING "+study)
        result[study]={"tested_biological_endpoints":len(rows),
                       "Holm_0_05":len(passed),
                       "significant":[{"group":k[0],"coordinate":k[1],
                                       "environment":k[2],"scale":k[3],
                                       **v} for k,v in passed.items()]}
    return result

def check_lego(d,old):
    s=d["S1"]["joint_8_constructs_28_pairs"]
    require((s["n"],s["taxa"],s["n_pairs"])==
            (SOURCE_CORE_N,SOURCE_CORE_TAXA,28),
            "FROZEN_SOURCE_COMMON_TRAIT_COHORT_CHANGED")
    base=s["models"][FINAL_MODEL]
    require(base["status"]=="ESTIMATED"
            and tuple(base["predictors"])==FINAL_ENV,
            "FROZEN_CORE_MODEL_NOT_BIO18_4VAR")
    require((base["n_among_stronger_raw"],
             base["n_among_stronger_adjusted"])==(21,24),
            "FROZEN_RV_BETWEEN_WITHIN_COUNTER_CHANGED")
    within=base["module_association"]["within"]
    among=base["module_association"]["among"]
    for k,v in (("within",within),("among",among)):
        require(v["median_same_module_RV"]>v["median_cross_module_RV"],
                "LEGO_CROSS_EXCEEDS_SAME "+k)
    labels=exact_label_reference(base["pairs"])
    scenario=old["S1"]["scenario_results"][
        "warmest_quarter_BIO18_hot_BIO05_elevation_reduced"]
    require(scenario["status"]=="ESTIMATED_POSTHOC_EXPLORATORY"
            and tuple(scenario["predictors"])==FINAL_ENV,
            "ORIGINAL_SCENARIO_NOT_IDENTICAL_FOUR_VAR")
    # Independently matched analysis code and model cohorts must agree.
    expected_pairs={(r["left"],r["right"]):r for r in scenario["pairs"]}
    for r in base["pairs"]:
        e=expected_pairs[(r["left"],r["right"])]
        for key in ("within_rv","among_rv"):
            close(r[key],e[key],"LEGO_4VAR_COVARIANCE_DIFFERED")
    configurations=scenario["three_way_phenotype_after_environment"]
    require(configurations["taxa_eligible"]==135
            and configurations["taxa_with_all_8_visible_configurations"]==123,
            "ORIGINAL_EIGHT_CONFIGURATION_CHECK_CHANGED")
    triad=scenario["within_and_among_three_domain_configurations"]
    a=triad["within_and_among"]["within"]
    b=triad["within_and_among"]["among"]
    require(a["n_nonempty_triple_configurations"]==8
            and b["n_nonempty_triple_configurations"]==8,
            "TRIAD_CONFIGURATIONS_MISSING")
    r=d["S1"]["bract_to_all_8_constructs"][
        "bridges_to_all_eight_original_constructs"]["chroma"]
    br=r["models"][FINAL_MODEL]
    require(r["n"]==1007 and r["taxa"]==38
            and not br["trait_cohort_vif_pass"]
            and br["max_vif_on_trait_cohort"]>5,
            "BRACT_LIMITED_COHORT_SENSITIVITY_CHANGED")
    s2={}
    for name,stratum in d["S2"]["strata"].items():
        item=stratum["joint_28_pairs"]
        require(item["status"]=="S2_WITHIN_RANGE_PHENOTYPE_COVARIATION"
                and item["models"][FINAL_MODEL]["status"]=="ESTIMATED"
                and len(item["models"][FINAL_MODEL]["pairs"])==28,
                "S2_28_PAIRS_INCOMPLETE "+name)
        edges=item["models"][FINAL_MODEL]["pairs"]
        same=[x["within_rv"] for x in edges
              if x["left_module"]==x["right_module"]]
        cross=[x["within_rv"] for x in edges
               if x["left_module"]!=x["right_module"]]
        s2[name]={"same_module_median_RV":median(same),
                  "cross_module_median_RV":median(cross),
                  "observed_n":item["n"],"cells":item["cells"],
                  "ONLY_within_range_photo_covariance":True}
    require(len(s2)==4,"S2_WCVP_FOUR_STRATA_MISSING")
    return {
        "S1":{"photos":s["n"],"taxa":s["taxa"],
             "n_visible_constructs":s["n_visible_constructs"],
             "n_pairs":s["n_pairs"],
             "raw_among_RV_stronger_pairs":21,
             "climate_adjusted_among_RV_stronger_pairs":24,
             "registered_group_sizes":{"presentation":1,"colour":3,"outline":4},
             "module_RV_medians":{"within_taxon":within,"among_taxa":among},
             "conditional_280_group_label_reference":labels,
             "n_taxa_eligible_for_above_median_triad":135,
             "n_taxa_with_all_8_above_median_triad_states":123,
             "mean_centered_phenotype_triad_within":a,
             "mean_centered_phenotype_triad_among":b,
             "above_median_triad_parity_excess":configurations["three_way_excess_parity"],
             "above_median_triad_parity_resample95":
                configurations["cluster_resampling_95pct_parity_interval"],
             "bract_chroma_original_photos":r["n"],
             "bract_chroma_original_taxa":r["taxa"],
             "bract_chroma_between_RV_four_BIO18":br["among_rv"],
             "bract_cohort_max_VIF":br["max_vif_on_trait_cohort"]},
        "S2":{"within_range_photo_covariation_by_focal_species_state":s2,
              "S2_original_bract_measurement_unavailable":True},
        "claim":"Visible module-like phenotype covariance & accessible configurations; genetic or functional LEGO untested.",
        "higher_order_conclusion":"configuration threshold changes third-order evidence: not robust proof of functional 3-way LEGO",
    }

def check_zip(path):
    data=Path(path).read_bytes()
    require(hashlib.sha256(data).hexdigest()==ARTIFACT_ZIP_SHA256,
            "FROZEN_SOURCE_ARTIFACT_ZIP_SHA256_MISMATCH")
    with zipfile.ZipFile(path) as z:
        parsed={}
        for name,expected in PINNED_MEMBERS.items():
            raw=z.read(name)
            require(hashlib.sha256(raw).hexdigest()==expected,
                    "PINNED_ORIGINAL_AGGREGATE_DRIFT "+name)
            parsed[name]=json.loads(raw)
    d=parsed["V4_S1_S2_four_vs_vif5_vif10_alltraits.json"]
    old=parsed["S1_S2_extreme_LEGO_results.json"]
    comparison=parsed["V4_S1_S2_four_model_comparison.json"]
    require(d["status"]=="SOURCE_VERIFIED_S1_S2_FOUR_MODELS_ALL_AVAILABLE_PHENOTYPES"
            and d["V4_NOT_V3"] is True,
            "NOT_V4_S1_S2_FINAL_SOURCE")
    require(d["S1"]["original_native_n"]==313350
            and d["S2"]["original_wild_source_n"]==297622,
            "ORIGINAL_SOURCE_OBSERVATION_NUMBERS_CHANGED")
    v=d["selection"]["fixed_model_set_diagnostics"][FINAL_MODEL]
    require(tuple(v["columns"])==FINAL_ENV
            and v["source_vif_pass"] and v["max_vif"]<5,
            "FINAL_BIO18_VIF_FIVE_FAILED")
    require(d["S2"]["optional_bract_source"]["status"]
            =="S2_BRACT_NOT_VERIFIABLY_QC_MEASURED",
            "S2_BRACT_EVIDENCE_CHANGED")
    require(comparison["status"]==
            "VERIFIED_FOUR_MODEL_CORRECTED_VS_UNCORRECTED",
            "SENSITIVITY_AUDIT_CHANGED")
    return {"state":"FINAL_V4_S1_S2_BIO18_HOLM_SCIENTIFIC_SNAPSHOT_VERIFIED",
       "repository":SOURCE_REPO,"scientific_analysis_commit":SOURCE_COMMIT,
       "authenticated_original_run":SOURCE_RUN,"artifact_id":SOURCE_ARTIFACT,
       "artifact_zip_sha256":ARTIFACT_ZIP_SHA256,
       "verified_JSON_members_sha256":PINNED_MEMBERS,
       "primary_model":{"name":FINAL_MODEL,
                        "climate_exposures":list(FINAL_ENV),
                        "original_full_context_max_VIF":v["max_vif"],
                        "not_outcome_driven_vifstep":True},
       "claim_gate":"Holm stepdown adjusted P<0.05 separately within 80 S1 / 128 S2 registered biological tests",
       "Holm_primary":check_families(d),
       "visible_phenotype_LEGO":check_lego(d,old),
       "other_models":"BIO12 four-variable, VIFstep5, VIFstep10 remain documented SI/sensitivity, NOT replacements for final primary claims.",
       "scientific_claim_limits":{"CIRSIN_not_heritable_genetic_LEGO":True,
          "no_physical_spine_or_mucilage_from_photo_proxy":True,
          "no_gravity_calibrated_nodding":True,
          "no_photograph_colour_calibration_or_pigment_biochemistry":True,
          "no_phylogenetic_or_geography_independent_adaptation_proven":True,
          "other_model_and_other_outcome_significance_does_not_upgrade_Holm_primary":True,
          "final_source_is_original_V4_not_V3":True},
       "no_private_observation_identifiers_emitted":True}

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--artifact-zip",type=Path,required=True)
    parser.add_argument("--output",type=Path,required=True)
    args=parser.parse_args()
    verified=check_zip(args.artifact_zip)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(verified,ensure_ascii=False,indent=2,
                                      allow_nan=False)+"\n",encoding="utf-8")
    print(json.dumps({"status":verified["state"],
      "holm_S1":verified["Holm_primary"]["S1"]["Holm_0_05"],
      "holm_S2":verified["Holm_primary"]["S2"]["Holm_0_05"],
      "registered_label_nulls":verified["visible_phenotype_LEGO"]["S1"][
          "conditional_280_group_label_reference"]},ensure_ascii=False))

if __name__=="__main__":
    main()
