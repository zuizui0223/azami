"""Source-pinned V4 S1/S2 BH-FDR primary reinterpretation.

No original/private observation rows needed. Reuses original-source-verified
aggregate from successful Actions run 40059949131, artifact 12546882413.
The photographic hue is one two-coordinate CIRCULAR trait: use the 2-df
cluster-sandwich Wald test rather than treating sin and cos as independent
biological traits in the discovery family.

Main and SI are separate posthoc, exploratory contrasts; S1 and S2 are
separate population/model families. Conditional observational associations
are NOT causal, adaptation, hereditary modularity or independent replication.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

ORIGINAL_SHA256 = "4089eda26a917971a043f7edcecf139e44dae302fefa7c4ea18ade38b7a3998f"
ORIGINAL_RUN = 40059949131
ORIGINAL_ARTIFACT = 12546882413
EXPECTED_OLD = {"S1/main_vif5":120, "S1/si_vif10":168,
                "S2/main_vif5":200, "S2/si_vif10":280}
EXPECTED_NEW = {"S1/main_vif5":100, "S1/si_vif10":140,
                "S2/main_vif5":160, "S2/si_vif10":224}
HUE_COORDINATE_DIAGNOSTICS = {"corolla_hue_sin", "corolla_hue_cos"}


def bh(pvalues):
    """Benjamini-Hochberg step-up q-values; zero/one and ties included."""
    p = list(pvalues)
    if not p:
        return []
    if any(not isinstance(v, (float, int)) or not math.isfinite(v)
           or not 0. <= v <= 1. for v in p):
        raise ValueError("INVALID_PRIMARY_DISCOVERY_P_VALUE")
    rank = sorted(range(len(p)), key=lambda i: (p[i], i))
    out = [0.]*len(p)
    smallest = 1.
    for rev in range(len(p)-1,-1,-1):
        i = rank[rev]
        smallest = min(smallest, p[i]*len(p)/(rev+1))
        out[i] = smallest
    return out


def paired_source_effect(table, family, row):
    """Use source-reported Wald intervals rather than fabricating CI from P."""
    coordinate = row["coordinate"]
    if coordinate == "hue_joint_sin_cos":
        return None
    model_name = family.split("/",1)[1]
    if family.startswith("S1/"):
        entry=table["S1"]["traits"][row["group"]]
    else:
        parts=row["group"].split("/")
        if len(parts)!=3:
            raise ValueError("SOURCE_S2_STRATUM_NOT_MAPPED")
        entry=table["S2"]["strata"][parts[0]+"/"+parts[1]][
            "source_measured_visible_traits"][parts[2]]
    model=entry["models"][model_name]
    if model["status"]!="ESTIMATED":
        raise ValueError("NON_ESTIMATED_SOURCE_FAMILY")
    for effect in model["component_effects"][coordinate]:
        if (effect["environment"]==row["environment"] and
            effect["scale"]==row["scale"]):
            if abs(effect["p_uncorrected"]-row["p_original"])>1.e-10:
                raise ValueError("SOURCE_COEFFICIENT_P_MISMATCH")
            if abs(effect["beta_std"]-row["beta_std"])>1.e-10:
                raise ValueError("SOURCE_COEFFICIENT_ESTIMATE_MISMATCH")
            return effect.get("ci95_wald")
    raise ValueError("SOURCE_COEFFICIENT_NOT_FOUND")


def handle(table):
    if (table.get("status") !=
        "SOURCE_VERIFIED_S1_S2_MAIN_SI_ALL_AVAILABLE_PHENOTYPES" or
        table.get("V4_NOT_V3") is not True or
        table.get("S1", {}).get("original_native_n")!=313350 or
        table.get("S2", {}).get("original_wild_source_n")!=297622):
        raise ValueError("UNVERIFIED_OR_NON_V4_AGGREGATE")
    outcomes={}
    for family,old_n in EXPECTED_OLD.items():
        old=table["multiplicity"][family]
        if old["n_tests"]!=old_n or len(old["rows"])!=old_n:
            raise ValueError("CHANGED_SOURCE_EXPLORATION_FAMILY")
        rows=[r for r in old["rows"]
              if r["coordinate"] not in HUE_COORDINATE_DIAGNOSTICS]
        if len(rows)!=EXPECTED_NEW[family]:
            raise ValueError("EXPECTED_JOINT_HUE_TRAIT_FAMILY_MISSING")
        names=[(r["group"],r["coordinate"],r["environment"],r["scale"])
               for r in rows]
        if len(names)!=len(set(names)):
            raise ValueError("DUPLICATED_TRAIT_ENVIRONMENT_FAMILY_TEST")
        if not any(r["coordinate"]=="hue_joint_sin_cos" for r in rows):
            raise ValueError("JOINT_CIRCULAR_HUE_TEST_NOT_PRESERVED")
        qs=bh([r["p_original"] for r in rows])
        updated=[]
        for r,q in zip(rows,qs):
            if r.get("test_kind") not in (
                "coordinate_exploratory","joint_circular_hue_2df"):
                raise ValueError("UNRECOGNISED_TEST_KIND")
            entry={k:v for k,v in r.items()
                   if k not in ("p_BH","p_Holm","Holm_lt_0_05",
                                "BH_lt_0_05")}
            if r["coordinate"]=="hue_joint_sin_cos":
                if r["beta_std"] is not None:
                    raise ValueError("JOINT_HUE_HAS_NO_SIGNED_SCALAR_BETA")
                entry["ci95_wald"]=None
            else:
                entry["ci95_wald"]=paired_source_effect(table,family,r)
            entry.update(q_BH_FDR_primary=float(q),
                         q_lt_0_05=bool(q<.05),
                         biological_endpoint=("circular_hue" if
                            r["coordinate"]=="hue_joint_sin_cos" else
                            r["group"].split("/")[-1]),
                         joint_hue_not_separate_sin_cos=
                           r["coordinate"]=="hue_joint_sin_cos")
            updated.append(entry)
        found=[r for r in updated if r["q_lt_0_05"]]
        outcomes[family]={"status":"BH_PRIMARY_SOURCE_AUTHENTICATED",
            "source_component_tests_including_redundant_hue":old_n,
            "n_primary_trait_environment_tests":len(rows),
            "n_q_lt_0_05":len(found),
            "n_q_lt_0_05_and_trait_cohort_VIF_pass":sum(
                r["trait_cohort_vif_pass"] for r in found),
            "findings":sorted(found,key=lambda x:x["q_BH_FDR_primary"]),
            "rows":updated}
    both={}
    for level in ("S1","S2"):
        main=outcomes[level+"/main_vif5"]["rows"]
        si=outcomes[level+"/si_vif10"]["rows"]
        def keyed(rows):
            return {(r["group"],r["coordinate"],r["environment"],r["scale"]):r
                    for r in rows if r["q_lt_0_05"]}
        m,s=keyed(main),keyed(si)
        shared=[]
        for key, row in m.items():
            if key not in s:
                continue
            other=s[key]
            consistent=(row["beta_std"] is None or other["beta_std"] is None or
                        row["beta_std"]*other["beta_std"]>0)
            if consistent:
                shared.append({"trait_environment":list(key),
                    "main_beta_std":row["beta_std"],
                    "si_beta_std":other["beta_std"],
                    "main_q":row["q_BH_FDR_primary"],
                    "si_q":other["q_BH_FDR_primary"],
                    "hue_joint_no_scalar_direction":row["beta_std"] is None})
        both[level]={"n_joint_FDR_positive":len(shared),
                     "comparisons":shared,
                     "not_independent_replication":True}
    return {"analysis_id":"azami_V4_S1_S2_alltraits_BH_FDR_primary_20261024",
        "status":"BH_PRIMARY_FROM_VERIFIED_ORIGINAL_V4_ARTIFACT",
        "source_github_actions_run":ORIGINAL_RUN,
        "source_artifact_id":ORIGINAL_ARTIFACT,
        "source_json_sha256":ORIGINAL_SHA256,
        "correction":"Benjamini-Hochberg FDR q<0.05; no Holm primary gate",
        "family_definition":"S1/S2 populations separate; main5/SI7 models separately; \
all available phenotype constructs incl bract in S1; circular hue once via joint 2df test",
        "families":outcomes,"both_main_and_SI":both,
        "limitations":{"exploratory_posthoc":True,
            "BH_FDR_interpretation_depends_on_test_dependency":True,
            "each_S1_S2_native_observation_not_independent":True,
            "SI_conditioned_on_added_RSDS_BIO13":True,
            "bract_is_image_projection_not_true_spines_or_mucilage":True,
            "image_angle_not_gravity_calibrated":True,
            "strong_geographic_phylogenetic_camera_confounding_unresolved":True,
            "no_adaptive_or_genetic_LEGO_causal_inference":True,
            "old_Holm_original_artifact_retained_only_as_historical_audit":True}
        }


def main():
    a=argparse.ArgumentParser()
    a.add_argument("--source-aggregate",type=Path,required=True)
    a.add_argument("--output",type=Path,required=True)
    args=a.parse_args()
    src=args.source_aggregate.read_bytes()
    if hashlib.sha256(src).hexdigest()!=ORIGINAL_SHA256:
        raise ValueError("SOURCE_AGGREGATE_SHA256_DRIFT")
    output=handle(json.loads(src))
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(output,indent=2,ensure_ascii=False,
                                       allow_nan=False)+"\n",encoding="utf-8")
    print(json.dumps({"status":output["status"],
        "counts":{k:[v["n_primary_trait_environment_tests"],v["n_q_lt_0_05"]]
                  for k,v in output["families"].items()},
        "same_association_in_main_and_SI":{
            k:v["n_joint_FDR_positive"] for k,v in
            output["both_main_and_SI"].items()}},ensure_ascii=False))

if __name__=="__main__":
    main()
