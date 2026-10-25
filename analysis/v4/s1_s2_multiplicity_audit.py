"""Posthoc multiplicity audit of executed V4 S1/S2 candidate climate slopes.

Not a preregistration and not a spatial/phylogenetic/causal validation. Every
originally tested S1 candidate/trait-coordinate/within-or-between term, and
each estimable S2 trait-by-range hydric coefficient, enters ONE exploratory
family per study layer. No cherry-picking a promising predictor, region,
species or visible trait before FDR/Holm correction. Hue sin/cos coordinates
are conservatively both counted; this does not establish two biological hues.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import numpy as np
from statsmodels.stats.multitest import multipletests


def validated_probability(x):
    v=float(x)
    if not np.isfinite(v) or v<0 or v>1:
        raise ValueError("UNTRUSTWORTHY_TRAIT_TEST_P_VALUE")
    return v


def collect_s1(d):
    rows=[]
    for scenario,model in d.get("S1",{}).get("scenario_results",{}).items():
        if model.get("status")!="ESTIMATED_POSTHOC_EXPLORATORY":
            continue
        for construct,entry in model["trait_environment_component_gradients"].items():
            for coordinate in entry["components"]:
                for row in coordinate["components"]:
                    rows.append({
                        "study":"S1",
                        "scenario":scenario,
                        "construct":construct,
                        "coordinate":coordinate["trait_coordinate"],
                        "exposure":row["environment"],
                        "scale":row["scale"],
                        "p_original":validated_probability(row["p_component"]),
                        "estimate_std":float(row["beta_std"]),
                        "n":model["sample_n"],
                        "taxa":model["taxa_n"],
                    })
    return rows


def collect_s2(d):
    rows=[]
    strata=d.get("S2_environmental_models",{}).get("strata",{})
    for stratum,entry in strata.items():
        for construct,trait in entry["original_visible_trait_environment_models"].items():
            for scenario,model in trait.get("range_models",{}).items():
                for coordinate,outcome in model["component_diagnostics_not_independent_endpoints"].items():
                    if outcome.get("status")!="EXPLORATORY_OBSERVED_S2_CLIMATE_GRADIENT":
                        continue
                    rows.append({
                        "study":"S2",
                        "stratum":stratum,
                        "scenario":scenario,
                        "construct":construct,
                        "coordinate":coordinate,
                        "exposure":outcome["hydric_predictor"],
                        "scale":"within_species_range",
                        "p_original":validated_probability(outcome["hydric_cluster_p_uncorrected"]),
                        "estimate_std":float(outcome["hydric_coefficient_per_sample_SD"]),
                        "n":outcome["n"],
                        "cells":outcome["cells"],
                        "photographers":outcome["observers"],
                    })
    return rows


def correct_family(rows):
    if not rows:
        return {"status":"NO_ESTIMABLE_ORIGINAL_TESTS","n_tests":0,"rows":[]}
    p=[x["p_original"] for x in rows]
    bh=multipletests(p,method="fdr_bh")[1]
    holm=multipletests(p,method="holm")[1]
    out=[]
    for r,b,h in zip(rows,bh,holm):
        out.append({**r,"p_BH_exploratory_family":float(b),
                    "p_Holm_exploratory_family":float(h),
                    "BH_0_05":bool(b<.05),
                    "Holm_0_05":bool(h<.05)})
    return {"status":"COMPLETED_EXPLORATORY_FAMILY",
            "n_tests":len(rows),
            "n_BH_005":int(sum(x["BH_0_05"] for x in out)),
            "n_Holm_005":int(sum(x["Holm_0_05"] for x in out)),
            "top30_by_raw_p":sorted(out,key=lambda x:x["p_original"])[:30],
            "rows":out}


def audit(d):
    require=d.get("status")=="SOURCE_VERIFIED_V4_EXTENSION_EXPLORATORY"
    if not require or d.get("not_V3_reanalysis") is not True:
        raise ValueError("CANNOT_REUSE_NON_V4_AGGREGATE")
    s1=correct_family(collect_s1(d))
    s2=correct_family(collect_s2(d))
    if s1["n_tests"]<100 or s2["n_tests"]<30:
        raise ValueError("S1_S2_NEW_FAMILY_INCOMPLETE")
    return {
      "analysis_id":"azami_V4_s1_s2_posthoc_multicandidate_multiplicity_20261024",
      "scope":"Separately corrected complete originally performed S1 and S2 exploratory candidate searches",
      "S1":s1,"S2":s2,
      "s1_coordinates_not_biologically_independent":True,
      "s2_same_original_photo_measurements_reused_across_scenarios":True,
      "BH_or_Holm_not_spatial_phylogenetic_or_adaptive_validation":True,
      "not_preregistered":True,
    }


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--aggregate",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    args=p.parse_args()
    data=audit(json.loads(args.aggregate.read_text()))
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(data,indent=2,allow_nan=False)+"\n")
    print(json.dumps({
        name:{k:data[name][k] for k in ("n_tests","n_BH_005","n_Holm_005")}
        for name in ("S1","S2")}))

if __name__=="__main__":
    main()
