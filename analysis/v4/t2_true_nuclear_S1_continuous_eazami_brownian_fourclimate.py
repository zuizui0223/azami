"""EAzami original Brownian ancestral-branch test for S1 CONTINUOUS angle on TRUE nuclear species trees.

The original 251 global native image angles and same S1/S2 CHELSA BIO18,
BIO05, RSDS and WorldClim elevation are rejoined to a genuine published
North American nuclear 63-tip species tree (NA42) and separate independent
Japan38 176-locus nuclear species tree (JP12), never grafts or ITS-only.
The EAzami Luu Brownian squared-change conditional-ancestral operator and
signed tip-to-node head angle delta x environment delta are used directly.
Permutation conditions on each genuine tree AND full climate niche vector.
This is NOT the EAzami discrete botanical U/D CTMC posterior; no camera
angle is dichotomised, no U/D transition or historical adaptation inferred.
"""
from __future__ import annotations
import argparse,copy,json
from pathlib import Path
from analysis.v4 import t2_s1_continuous_genuine_nuclear_two_tree_clustering as dna
from analysis.v4 import t2_world_s1_fouraxis_eazami_brownian_branch_regime as world

FOUR=world.FOUR
N_PERM=9999
PANELS={
 "NorthAmerica_true_nuclear42":42,
 "Japan_true_176loci_nuclear12":12,
}

def independent_tree(tree,tipmap):
    t=copy.deepcopy(tree)
    original={x.name for x in t.get_terminals()}
    if not set(tipmap)<=original:raise ValueError("No real molecular source-tip matching")
    for tip in t.get_terminals():
        if tip.name in tipmap:tip.name=tipmap[tip.name]
    if len({x.name for x in t.get_terminals()})!=len(t.get_terminals()):
        raise ValueError("Original nuclear DNA tip aliases collide")
    return t

def evaluate(original251,newfour,*,nperm=N_PERM):
    by,sets,_=world.source_panels(original251,newfour)
    raw=[dna.upstream(dna.NA_REPO,dna.NA_REF,dna.NA_TREE),
         dna.upstream(dna.J_REPO,dna.J_REF,dna.J_TREE),
         dna.upstream(dna.J_REPO,dna.J_REF,dna.J_MAP),
         dna.upstream(dna.J_REPO,dna.J_REF,dna.J_MANIFEST)]
    full,real,provenance=dna.read_sources(original251,*raw)
    if len(full)!=251 or len(by)!=251:raise ValueError("Original world S1 source population changed")
    res={}
    for i,(_, (tree,mapping)) in enumerate(real.items()):
        label=list(PANELS)[i]
        names=sorted(set(mapping.values()))
        if len(names)!=PANELS[label]:raise ValueError("Original authentic tree exact source tip coverage changed")
        for name in names:
            if name not in by or full[name]["WCVP_key_collides_source_names"]!="False":
                raise ValueError("S1 original source missing, collided or 4 site climate absent")
        original=independent_tree(tree,mapping)
        result={}
        for scope in ("all_edges","internal_edges"):
            r=world.evaluate_panel(by,names,[(label,original)],nperm,
                   world.stable_seed(label),edge_scope=scope)
            if r["n_original_tree_placement_scenarios"]!=1:
                raise ValueError("A true multilocus nuclear tree was treated as correlated scenarios")
            result[scope]=r
        sample10=[n for n in names if int(by[n]["n_photo_complete_four_axis"])>=10]
        if len(sample10)>=8:
            result["original_S1_photo_min10"]={
                "status":"COMPUTED_SAME_TRUE_NUCLEAR_TREE",
                **world.evaluate_panel(by,sample10,[(label,original)],
                    nperm,world.stable_seed(label+"_10"))}
        else:
            result["original_S1_photo_min10"]={
                "status":"STOP_TOO_FEW_GENUINE_MULTILOCUS_TIPS",
                "n_photo_min10_nuclear_tips":len(sample10),"minimum_for_branch_test":8}
        res[label]={"n_authentic_original_nuclear_S1_taxa":len(names),
             "original_source_fourclimate_photos":sum(int(by[n]["n_photo_complete_four_axis"]) for n in names),
             "n_photo10_tips":len(sample10),
             "source_type":"independently_published_authentic_NA63_species_tree" if i==0
                    else "EAzami_Japan38_176_original_nuclear_loci_species_tree",
             "comparison":result}
    if res["NorthAmerica_true_nuclear42"]["n_photo10_tips"]!=40:
        raise ValueError("Original true NA40 photos positive source gate changed")
    if res["Japan_true_176loci_nuclear12"]["n_photo10_tips"]!=6:
        raise ValueError("Japan12 S1 camera angle precision gate drift")
    return {
       "schema":"azami_S1_continuous_true_NA42_Japan12_EAzami_Brownian_branch_fourclimate_v1",
       "status":"TRUE_NUCLEAR_BROWNIAN_S1_PHOTO_BRANCH_FOUR_CLIMATE_PERMUTATION_EXECUTED",
       "S1_251_original_source_names":251,
       "same_original_S1_site_four_environment_complete_head_photos":169020,
       "axes":list(FOUR),
       "true_independent_multilocus_nuclear_species_tree_panels":2,
       "n_true_multi_locus_species_total":54,
       "not_one_joined_na_japan_species_tree":True,
       "S1_angular_outcome_is_CONTINUOUS_and_NOT_EAzami_discrete_U_D":True,
       "operator":"same_EAzami_Brownian_conditional_ancestor_squared_change_signed_branch_alignment",
       "conditional_tip_permutations":nperm,
       "BIO18_positive_sign_fixed_from_global_S1_S2_original_species_pattern":True,
       "no_plant_flora_botanical_U_D_or_chloroplast_or_ITS_primary_tree":True,
       "no_adaptive_origin_or_evolutionary_transition_counts_claimed":True,
       "source_provenance":provenance,
       "panels":res,
    }

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--original251",required=True,type=Path)
    p.add_argument("--newfour",required=True,type=Path)
    p.add_argument("--permutations",type=int,default=N_PERM)
    p.add_argument("--out",required=True,type=Path)
    x=p.parse_args()
    if x.permutations<999:raise ValueError("True nuclear S1 branch permutation needs >=999 replicates")
    v=evaluate(x.original251,x.newfour,nperm=x.permutations)
    x.out.parent.mkdir(parents=True,exist_ok=True)
    x.out.write_text(json.dumps(v,indent=2,sort_keys=True,ensure_ascii=False)+"\n")
    for lab,t in v["panels"].items():
        for scope,r in t["comparison"].items():
            if r.get("status","")=="STOP_TOO_FEW_GENUINE_MULTILOCUS_TIPS":continue
            print("TRUE_NUCLEAR_CONTINUOUS_S1_EAzami_BROWNIAN="+json.dumps({
              "panel":lab,"scope":scope,"taxa":r["n_original_source_names"],
              "nuc_tree":t["source_type"],
              "S1_BIO18_branch_score":r["scores_standardised_climatic_shift_per_signed_degree"]["chelsa_bio18"],
              "BIO18_positive_one_sided_fraction":r["BIO18_positive_direction_one_sided_fraction"],
              "BIO18_fourclimate_maxT_fraction":r["four_axis_maxT_fraction"]["chelsa_bio18"]
            },sort_keys=True),flush=True)
if __name__=="__main__":main()
