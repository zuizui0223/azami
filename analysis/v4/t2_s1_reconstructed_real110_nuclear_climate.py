"""S1 original251 nuclear-SRA DNA -> actual 110 tip technical ML -> BIO12.

A separately labelled POSTHOC one-source quality sensitivity after predeclared
111 source gate failed on one sparse C. tioganum. The 108 original selected
nuclear candidates are genuine 120k public Moreyra paired-read majority Q20
sites at shared 89/111 source reference positions. No reference-imputed bases
or pooling per-species samples. NOT validated single copy paralog-safe species
tree or dated phylogeny, and not an independent adaptive-origin analysis.

Botanical U/D is independently source-coded; current BIO12 from original S1
photo-taxon CHELSA centroids is environmental niche context, NOT a trait
orientation inferred from camera vertical nor historic selection.
"""
from __future__ import annotations
import argparse
import csv
from copy import deepcopy
import hashlib
import json
from pathlib import Path

import numpy as np
from Bio import Phylo

from analysis.v4.t2_na63_actual_tree_rain_gls import (
    brownian_shared_ancestry, fit_models, equal_variance_independent_OLS,
)
from analysis.v4.t2_real110_posthoc_one_sparse_source_sensitivity import (
    EXCLUDED_TIP, EXCLUDED_NAME, ORIGINAL_SHARED_SITES,
    ORIGINAL_LOCI, RECOVERED_LOCI, BEFORE_SAMPLES, AFTER_SAMPLES
)

def read_csv(p):
    with Path(p).open(newline="",encoding="utf-8-sig") as fh:
        return list(csv.DictReader(fh))


def analyze(s1_csv,source110_csv,source110_receipt,tree_path):
    s1=read_csv(s1_csv)
    source=read_csv(source110_csv)
    a=json.loads(Path(source110_receipt).read_text(encoding="utf-8"))
    if (len(s1)!=251 or len({r["accepted_key"] for r in s1})!=244
        or sum(int(r["original_photo_count"]) for r in s1)!=170066):
        raise ValueError("S1 original251 source/photos not conserved")
    if (a["original_all111_Q20_successful_real_SRA_BioSamples"]!=BEFORE_SAMPLES or
        a["posthoc_sensitivity_taxa"]!=AFTER_SAMPLES or
        a["observed_Q20_reference_positions_reused_unchanged"]!=ORIGINAL_SHARED_SITES or
        a["original_111_source_108_retained_loci"]!=RECOVERED_LOCI or
        a["original111_technical_gate_passed"] or
        a["original_111_reference_position_selection_recomputed_on110"] or
        a["true_Mk_nodding_transition_or_rainfall_adaptation_supported"]):
        raise ValueError("Original source source111-to110 post hoc selection was changed")
    if len(source)!=AFTER_SAMPLES or len({r["source_tip"] for r in source})!=AFTER_SAMPLES:
        raise ValueError("Do not infer false 110 independence of original SRA BioSamples")
    if len({r["source_BioSample"] for r in source})!=AFTER_SAMPLES:
        raise ValueError("One original biological sample used twice")
    if EXCLUDED_TIP in {r["source_tip"] for r in source} or EXCLUDED_NAME in {r["taxon_name"] for r in source}:
        raise ValueError("Failed quality specimen quietly restored without true Q20 source")
    matched={r["original_source_name"]:r for r in s1}
    seqnames={r["taxon_name"] for r in source}
    if not seqnames<=set(matched):
        raise ValueError("Molecular sources outside independently authenticated S1 original species roster")
    t=Phylo.read(str(tree_path),"newick")
    tips=[x.name for x in t.get_terminals()]
    byid={r["source_tip"]:r for r in source}
    if len(tips)!=AFTER_SAMPLES or len(set(tips))!=AFTER_SAMPLES or set(tips)!=set(byid):
        raise ValueError("Original 110 authenticated real DNA taxon IDs and IQTREE mismatch")
    if any(node.branch_length is None or node.branch_length<0
           for node in t.find_clades() if node is not t.root):
        raise ValueError("Tree must be built with actual real DNA branch lengths")
    traits={}
    for tip in tips:
        molecular=byid[tip]
        biological=matched[molecular["taxon_name"]]
        if biological["accepted_key"] not in tip or (
            biological["original_source_name"]!=molecular["taxon_name"]):
            raise ValueError("Exact S1 WCVP/source original binomial sample identity changed")
        if biological["botanical_U_D_state"] in ("U","D"):
            traits[tip]=(molecular,biological)
    d=sum(b["botanical_U_D_state"]=="D" for a,b in traits.values())
    u=sum(b["botanical_U_D_state"]=="U" for a,b in traits.values())
    result={
        "schema":"azami_originalS1_110_public_nuclear_candidate_loci_sensitivity_vs_climate_v1",
        "status":"OWN_SRA_110_TAXON_Q20_REAL_NUCLEAR_REFERENCE_ML_POSTHOC_TECHNICAL_TREE",
        "full_original_S1_species_units_retained":251,
        "original_S1_accepted_WCVP_keys":244,
        "original_public_SRA_BioSamples_available":BEFORE_SAMPLES,
        "Q20_posthoc_quality_omission":EXCLUDED_NAME,
        "n_original_S1_source_tip_Q20_sequences_in_tree":AFTER_SAMPLES,
        "original_common_Q20_reference_loci":RECOVERED_LOCI,
        "selected_reference_Q20_sites":ORIGINAL_SHARED_SITES,
        "n_authority_scored_binary_tree_tips":d+u,
        "n_botanical_D_tree_tips":d,
        "n_botanical_U_tree_tips":u,
        "n_botanical_unknown_tree_tips":AFTER_SAMPLES-d-u,
        "source_original_111_predeclared_technical_gate_succeeded":False,
        "same_original_111_Q20_alignment_mask_selected_pre_110_exclusion":True,
        "IQTREE_unrooted_branch_length_original_Q20_tree_SHA256":hashlib.sha256(Path(tree_path).read_bytes()).hexdigest(),
        "nuclear_loci_are_reference_copy_candidates_not_ortholog_curated":True,
        "posthoc_quality_exclusion_not_a_preregistered_110_species_model":True,
        "DNA_tree_reconstruction_done_independently_from_any_published_NA63_Japan38_phylogeny":True,
        "true_110_biological_species_tree_completed":False,
        "full_global_244_species_tree_completed":False,
        "historic_U_to_D_gain_count_inferred":False,
        "photos_used_to_choose_D_vs_U":False,
        "annual_rain_present_S1_distribution_not_origin_time_rain":True,
        "original_source_sra_receipt_SHA256":hashlib.sha256(Path(source110_receipt).read_bytes()).hexdigest(),
        "covariance_model":"root at ITS/SRA ML midpoint for sensitivity, NOT historical biological root",
        "comparison_estimated":False,
    }
    if d>=3 and u>=5:
        try:
            t.root_at_midpoint()
            part=deepcopy(t)
            for node in list(part.get_terminals()):
                if node.name not in traits:part.prune(node)
            if part.root.branch_length is not None:
                part.root.branch_length=0.
            ids=[x.name for x in part.get_terminals()]
            C=brownian_shared_ancestry(part,ids)
            heads=np.array([traits[x][1]["botanical_U_D_state"]=="D" for x in ids],float)
            rain=np.array([float(traits[x][1]["chelsa_bio12"]) for x in ids])
            tbi=np.array([float(traits[x][1]["chelsa_bio01"]) for x in ids])
            X=np.column_stack([np.ones(len(ids)),heads])
            Xt=np.column_stack([X,(tbi-tbi.mean())/tbi.std()])
            result["comparison_estimated"]=True
            result["S1_BIO12_mean_D_minus_U_mm"]=float(rain[heads==1].mean()-rain[heads==0].mean())
            result["same_original_S1_taxon_real_Q20_midpoint_PGLS"]={
                "BIO12_only":fit_models(rain,X,C),
                "BIO1_adjusted_BIO12":fit_models(rain,Xt,C)}
            result["same_original_S1_taxon_equal_variance_IID_OLS"]={
                "BIO12_only":equal_variance_independent_OLS(rain,X),
                "BIO1_adjusted_BIO12":equal_variance_independent_OLS(rain,Xt)}
            result["historical_direction_root_verified"]=False
        except (ValueError,np.linalg.LinAlgError) as e:
            result["comparison_block_reason"]=str(e)[:280]
    else:result["comparison_block_reason"]="Too few strict independently botanically scored U/D S1 tips"
    return result


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--s1",type=Path,required=True)
    p.add_argument("--source110",type=Path,required=True)
    p.add_argument("--receipt",type=Path,required=True)
    p.add_argument("--tree",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    out=analyze(a.s1,a.source110,a.receipt,a.tree)
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(out,sort_keys=True,ensure_ascii=False,indent=2,allow_nan=False)+"\n",encoding="utf-8")
    print("OWN_S1_REAL110_SRA_NUCLEAR_CANDIDATE_TREE_RAIN="+json.dumps({
        "actual_seq":out["n_original_S1_source_tip_Q20_sequences_in_tree"],
        "D":out["n_botanical_D_tree_tips"],"U":out["n_botanical_U_tree_tips"],
        "comparison":out["comparison_estimated"],
        "effect":out.get("S1_BIO12_mean_D_minus_U_mm"),
        "blocked":out.get("comparison_block_reason"),
    },ensure_ascii=False,sort_keys=True),flush=True)


if __name__=="__main__":
    main()
