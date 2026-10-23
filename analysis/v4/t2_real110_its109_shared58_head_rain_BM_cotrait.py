"""T2 source-real S1 original photo head presentation vs present rainfall contrasts.

The SAME source 58 Cirsium species are present in TWO newly reconstructed
original DNA genealogies: original SRA Q20 108-candidate-nuclear-locus
posthoc110 (original111 technical gate FAILED) and GenBank ITS109 (multi-copy
rDNA); no EAzami historic tree or generic grafted placements are used.

Fit Brownian ancestral shared-branch covariance to log(BIO12) and
S1 OWN MEASURED image-vertical head-angle mean on each original gene tree,
whiten residual trait vectors for phylogeny and an intercept, and compute
a SAME-TAXON cross-trait correlation with conditional t approximation.
Compare 0 and declared 5% tip photo/phylogeny residual nugget. These
nuggets are NOT independently calibrated. Neither regression is a
historic precipitation-at-gain model or a biological gravity-verified
downward-presentation angle. No adaptive radiation or within-day motion.
"""
from __future__ import annotations
import argparse
import json
import math
from pathlib import Path
import numpy as np
from scipy import linalg, stats

from analysis.v4.t2_110_multilocus_vs_109_its_photo_niche import (
    verify_sources, original_true_covariance, phenotype,
    N_MULTILOCUS, N_ITS, N_EXACT_OVERLAP,
    EXPECTED_S1_SHA, EXPECTED_ML108_SHA, EXPECTED_ITS_SHA
)

SOURCE_SCOPES=(("nuclear108_posthoc110","Q20","full"),
               ("GenBankITS109","ITS","full"),
               ("nuclear108_same58","Q20","shared"),
               ("GenBankITS_same58","ITS","shared"))
NUISANCE_TIP_FRACTIONS=(0.,.05)


def phylogenetic_pair_correlation(names,K,byname,nugget):
    if not (12<=len(names)<=110 and 0<=nugget<=.20):
        raise ValueError("Original source species or noise assumptions outside contract")
    if len(set(names))!=len(names) or K.shape!=(len(names),len(names)):
        raise ValueError("One-to-one original source S1 photo and gene tree taxon tips required")
    n=len(names)
    C=K+float(nugget)*np.median(np.diag(K))*np.eye(n)
    L=linalg.cholesky(C,lower=True,check_finite=True)
    q=linalg.solve_triangular(L,np.ones(n),lower=True)
    def independent_residual(original):
        value=linalg.solve_triangular(L,original,lower=True)
        return value - q*(q @ value)/(q @ q)
    rain=phenotype(names,byname,"log_BIO12")
    angle=phenotype(names,byname,"S1_photo_mean_image_head_angle")
    if np.std(rain)<=1e-12 or np.std(angle)<=1e-12:
        raise ValueError("Original two measured phenotypes have no variation")
    a=independent_residual(rain)
    b=independent_residual(angle)
    r=float(a@b/np.sqrt((a@a)*(b@b)))
    if abs(r)>=1:
        raise ValueError("Source species independent contrasts perfectly collinear")
    # Fixed-covariance approximate t after jointly eliminating the intercept.
    # Neither tree uncertainty nor non-phylogenetic spatial and measurement
    # uncertainty enter this t or CI.
    t=r*np.sqrt((n-2)/(1-r*r))
    p=float(2*stats.t.sf(abs(t),n-2))
    z=math.atanh(r)
    zcut=stats.norm.ppf(.975)/math.sqrt(n-3)
    ci=[float(math.tanh(z-zcut)),float(math.tanh(z+zcut))]
    return {
        "S1_exact_species_count":n,
        "S1_original_own_photo_vs_logBIO12_BM_whitened_correlation":r,
        "phylogeny_conditioned_approx_two_sided_t_p":p,
        "conditional_Fisher_z_95CI":ci,
        "tip_independent_noise_fraction_of_median_real_gene_tree_depth":float(nugget),
        "correlated_tip_anthesis_not_directly_measured":True,
        "head_angle_is_original_S1_camera_relative_not_gravity":True,
        "present_CHELSA_annual_rain_is_not_historical_gain_rain":True,
        "tree_root_midpoint_only_and_gene_copy_uncertainty_unpropagated":True,
        "spatial_and_biogeographic_confounders_not_adjusted":True,
        "not_BM_joint_evolutionary_rate_or_independent_gain_count":True,
    }


def analyze(s1,source110,receipt,real110,its109):
    byname,nuclear,ITS,nuclear_map,shared,source=verify_sources(
        s1,source110,receipt,real110,its109)
    its_map={x.name:x.name for x in ITS.get_terminals()}
    results={}
    for label,which,level in SOURCE_SCOPES:
        tree,mapping=(nuclear,nuclear_map) if which=="Q20" else (ITS,its_map)
        keep=(shared if level=="shared" else None)
        names,K=original_true_covariance(tree,mapping,keep=keep)
        results[label]={
            "source_true_gene_tree":"108_nuclear_candidate_loci_Q20_110_posthoc" if which=="Q20"
                                    else "ITS_multicopy_rDNA_109",
            "same_original_S1_photographic_taxon_names":len(names),
            "original_S1_photo_measurements_represented":
                sum(int(byname[x]["original_photo_count"]) for x in names),
            "actual_same_taxon_names_exact_with_other_tree":(level=="shared"),
            "sensitivity":{
                ("nugget5pct" if e else "nugget0pct"):
                phylogenetic_pair_correlation(names,K,byname,e)
                for e in NUISANCE_TIP_FRACTIONS
            }
        }
    lengths=[results[x]["same_original_S1_photographic_taxon_names"]
             for x,_,_ in SOURCE_SCOPES]
    if lengths!=[N_MULTILOCUS,N_ITS,N_EXACT_OVERLAP,N_EXACT_OVERLAP]:
        raise ValueError("Two genuine nuclear gene trees and exact 58-tip shared S1 cohort changed")
    return {
        "schema":"azami_T2_original_photo_real110_multilocus_ITS109_and_identical58_BM_whitened_cotrait_v1",
        "status":"EMPIRICAL_TWO_PUBLIC_NUCLEAR_GENEALOGIES_S1_OWN_HEAD_VS_PRESENT_RAIN_NO_COEVOLUTION_CLAIM",
        "original_frozen_S1_251_name_concepts":251,
        "original_S1_170066_photograph_measurements":170066,
        "original_111_Q20_primary_technical_gate_success":False,
        "Q20_original_nuclear_candidate_loci":108,
        "n_public_real_nuclear_ITS_taxa":109,
        "n_public_real_posthoc110_nuclear_candidates":110,
        "shared_58_original_S1_same_biological_taxon_names":sorted(shared),
        "source_biological_traits":["log_original_species_annual_BIO12",
                                    "original_S1_own_continuous_photo_image_head_angle"],
        "source_cross_trait_result_has_no_verified_erect_vs_nodding_gravity_codes":True,
        "source_tip_photo_error_five_pct_is_assumption_not_camera_calibration":True,
        "nuclear_copies_not_verified_independent_orthologs":True,
        "genealogy_ITS_multi_copy_rDNA_not_species_tree":True,
        "phylogenetic_ancestral_nodding_origins_inferred":False,
        "rain_shield_fitness_causality_estimated":False,
        "source_original_digests":{
            "S1_original_photo":EXPECTED_S1_SHA,
            "Q20_original110":EXPECTED_ML108_SHA,
            "ITS_original109":EXPECTED_ITS_SHA,
        },
        "results":results,
    }


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ("s1","source110","receipt","Q20-tree","ITS-tree","out"):
        p.add_argument("--"+name,type=Path,required=True)
    a=p.parse_args()
    result=analyze(a.s1,a.source110,a.receipt,a.Q20_tree,a.ITS_tree)
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(result,ensure_ascii=False,sort_keys=True,
                                indent=2,allow_nan=False)+"\n",encoding="utf-8")
    print("T2_ACTUAL_Q20_VS_ITS_SAME58_BM_HEAD_RAIN="+json.dumps({
        name:{
            "n":data["same_original_S1_photographic_taxon_names"],
            "r_nugget0":data["sensitivity"]["nugget0pct"][
                "S1_original_own_photo_vs_logBIO12_BM_whitened_correlation"],
            "r_nugget5":data["sensitivity"]["nugget5pct"][
                "S1_original_own_photo_vs_logBIO12_BM_whitened_correlation"],
            "approx_p_nugget5":data["sensitivity"]["nugget5pct"][
                "phylogeny_conditioned_approx_two_sided_t_p"],
        } for name,data in result["results"].items()
    },ensure_ascii=False,sort_keys=True),flush=True)


if __name__=="__main__":
    main()
