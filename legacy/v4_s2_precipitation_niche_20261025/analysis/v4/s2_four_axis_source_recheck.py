"""Verify S2's four frozen exposure margins and all four-exposure angle tests.

Reads ONLY the pinned public aggregate ZIP from Actions 40062126006 and
exports only anonymous summary rows. An aggregate cannot identify a joint
four-dimensional niche shift, source-population genetics or equivalence.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
from pathlib import Path
import zipfile

ZIP_SHA="6e081dc081bb8ee2311a05ac5b7fcddd3f6309be2fa909e87d432789f4e0e205"
MEMBERS={
 "S1_S2_extreme_LEGO_results.json":"bb01d3b5da03148d84b2521c9c92cf14f7f6177435e145118f6ff6df8846733f",
 "V4_S1_S2_four_vs_vif5_vif10_alltraits.json":"3738292e0a64da183dfe219c854c16cda199a7f9769f4928da8ab07657d60b70"
}
ENV=("chelsa_bio18","chelsa_bio05","chelsa_rsds_mean","topo_elevation")
SPECIES=("Cirsium arvense","Cirsium vulgare")
SOURCE_N={
 ("Cirsium arvense","native"):78646,("Cirsium arvense","introduced"):70964,
 ("Cirsium vulgare","native"):62368,("Cirsium vulgare","introduced"):85644}
EXPECTED_BIO18={"Cirsium arvense":48.2,"Cirsium vulgare":49.7}

def audited(artifact:Path)->tuple[dict,dict]:
    if hashlib.sha256(artifact.read_bytes()).hexdigest()!=ZIP_SHA:
        raise ValueError("SOURCE_ZIP_SHA256_MISMATCH")
    with zipfile.ZipFile(artifact) as z:
        data={}
        for name,expected in MEMBERS.items():
            content=z.read(name)
            if hashlib.sha256(content).hexdigest()!=expected:
                raise ValueError("SOURCE_MEMBER_SHA256_MISMATCH: "+name)
            data[name]=json.loads(content)
    return tuple(data.values())

def verify(aggregate:dict, models:dict)->tuple[list[dict],list[dict]]:
    if aggregate["S2"]["source_wild_two_species_observations"]!=297622:
        raise ValueError("S2_SOURCE_POPULATION_CHANGED")
    if models["S2"]["original_wild_source_n"]!=297622:
        raise ValueError("S2_MODEL_POPULATION_CHANGED")
    rows=[]
    for sp in SPECIES:
        strata=aggregate["S2"]["taxa"][sp]
        for state in ("native","introduced"):
            if strata[state]["observations"]!=SOURCE_N[(sp,state)]:
                raise ValueError("WCVP_RANGE_SOURCE_DENOMINATOR_CHANGED")
        for v in ENV:
            native=strata["native"]["climate"][v]
            introduced=strata["introduced"]["climate"][v]
            rows.append({
                "species":sp,"variable":v,
                "native_observed_site_n":native["n"],
                "introduced_observed_site_n":introduced["n"],
                "native_median":native["median"],"introduced_median":introduced["median"],
                "median_introduced_minus_native":introduced["median"]-native["median"],
                "native_q05":native["q05"],"native_q95":native["q95"],
                "introduced_q05":introduced["q05"],"introduced_q95":introduced["q95"]})
        d=rows[-len(ENV)]["median_introduced_minus_native"]
        if abs(d-EXPECTED_BIO18[sp])>0.011:
            raise ValueError("FROZEN_BIO18_RANGE_MEDIAN_CHANGED")
    tests=[]
    family=models["multiplicity"]["S2/four_bio18"]
    if family["n_tests"]!=128:
        raise ValueError("HOLM_FAMILY_CHANGED")
    for row in family["rows"]:
        if row["coordinate"]=="orientation_image_vertical_angle":
            if row["environment"] not in ENV:
                raise ValueError("NOT_THE_FROZEN_4VAR_ANGLE_MODEL")
            tests.append({
                "stratum":row["group"],"exposure":row["environment"],
                "n_photo_eligible":row["n"],"standardized_beta":row["beta_std"],
                "raw_p":row["p_original"],"holm_p":row["p_Holm"]})
    if len(tests)!=16 or any(row["holm_p"]<.05 for row in tests):
        raise ValueError("ORIGINAL_S2_ANGLE_HOLM_RESULT_CHANGED")
    return rows,tests

def main()->None:
    p=argparse.ArgumentParser()
    p.add_argument("--artifact",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    a.out.mkdir(parents=True,exist_ok=True)
    climate,angle=verify(*audited(a.artifact))
    def write_csv(name,rows):
        with (a.out/name).open("w",encoding="utf-8",newline="") as f:
            wr=csv.DictWriter(f,fieldnames=list(rows[0]))
            wr.writeheader()
            wr.writerows(rows)
    write_csv("S2_four_axis_site_occupancy_margins.csv",climate)
    write_csv("S2_four_axis_angle_coefficients.csv",angle)
    (a.out/"S2_four_axis_source_recheck.json").write_text(
        json.dumps({
            "status":"ORIGINAL_FOUR_AXIS_MARGINAL_SITE_OCCUPANCY_AND_ANGLE_VERIFIED",
            "provenance":{"artifact_id":12547876704,"sha256":ZIP_SHA},
            "four_axis_marginal_occupancy":climate,
            "all_16_original_angle_coefficients":angle,
            "all_4D_joint_climate_niche_overlap_not_identifiable_from_aggregates":True,
            "range_main_effect_and_interaction_not_measured_by_stratum_slopes":True,
            "unmeasured_climatic_availability_and_introduction_history":True,
            "zero_Holm_angle_effects_does_not_prove_angle_equivalence":True,
            "data_privacy":"no individual obs/photographer or coordinates exported"
        },ensure_ascii=False,indent=2,allow_nan=False)+"\n",encoding="utf-8")
    print("OK verified four marginals in both species, 16 angle coefficients, 0 Holm positives")
if __name__=="__main__":main()
