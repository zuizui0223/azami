"""Explain why original S1 interspecific BIO12 slope and newly inferred ITS109
phylogenetic comparison differ, without equating either to adaptive radiation.

Compare SAME original 5 CHELSA axes + original camera-frame head angle with
taxon equal-weight, HC3 heteroscedasticity-robust uncertainty for all S1
251 (accepted key x name), the 109 exact GenBank ITS tree tips, and 142
absent from that gene tree. The original 244-key, image-level mixed model
is an entirely different estimand and is never overwritten.

Also inspect the nonrandom environmental coverage of nuclear ITS accessions,
and corroborate the original new ML tree's SH-aLRT / UFBoot support. Does
not infer biological nodding, within-day plasticity, convergence, radiation,
diversification rate, fitness, or a dated species-tree trait event.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
from pathlib import Path

import numpy as np
from Bio import Phylo
from scipy import stats

EXPECTED_251_SHA256="47b32972b20564ade163907ef6fdf3eb2a69f8b42c8d87b8cf4ea75ca1682dc1"
EXPECTED_109_ML_SHA256="8fd10d299b4d520ce08f2552993b02d3cc1b10350f0fb68ea24874dbb9ccfe92"
EXPECTED_IMAGES=170066
NUCLEAR_TIPS=109
SOURCE_NAMES=251
ACCEPTED_KEYS=244
FIVE=("chelsa_bio01","chelsa_bio04","chelsa_bio12","chelsa_bio15","chelsa_rsds_mean")
RESPONSES=("mean_image_angle_deg","median_image_angle_deg")


def read_source(path):
    data=Path(path).read_bytes()
    if hashlib.sha256(data).hexdigest()!=EXPECTED_251_SHA256:
        raise ValueError("Frozen original S1 owner-photo and CHELSA source SHA mismatch")
    rows=list(csv.DictReader(io.StringIO(data.decode("utf-8-sig"))))
    if (len(rows)!=SOURCE_NAMES or len({r["accepted_key"] for r in rows})!=ACCEPTED_KEYS
        or len({r["original_source_name"] for r in rows})!=SOURCE_NAMES
        or sum(int(r["original_photo_count"]) for r in rows)!=EXPECTED_IMAGES
        or any(r["botanical_gravity_posture_inferred_from_photo_angle"]!="False" for r in rows)):
        raise ValueError("Original 251 source names, 244 accepted keys or 170066 measured photo angles changed")
    return rows


def read_ITS_tree(path):
    data=Path(path).read_bytes()
    if hashlib.sha256(data).hexdigest()!=EXPECTED_109_ML_SHA256:
        raise ValueError("Original independently reconstructed GenBank ITS tree SHA mismatch")
    tree=Phylo.read(io.StringIO(data.decode("utf-8-sig")),"newick")
    tips=[t.name for t in tree.get_terminals()]
    if len(tips)!=NUCLEAR_TIPS or len(set(tips))!=NUCLEAR_TIPS:
        raise ValueError("New ITS source ML tree not original 109 exact taxa")
    info=[]
    for node in tree.get_nonterminals():
        if not node.name or "/" not in node.name:continue
        left,right=node.name.split("/")
        info.append((float(left),float(right)))
    if len(info)!=106:
        raise ValueError("New ITS source tree ML support labels unexpectedly missing")
    return set(tips), info


def fit_equal_species_5axis_HC3(rows, endpoint):
    if len(rows)<12 or endpoint not in RESPONSES:
        raise ValueError("Insufficient original photo species for 5-axis model")
    env=np.array([[float(r[v]) for v in FIVE] for r in rows])
    y=np.array([float(r[endpoint]) for r in rows])
    if np.any(~np.isfinite(env)) or np.any(~np.isfinite(y)) or np.any(env[:,2]<=0):
        raise ValueError("Original native CHELSA/photo angles missing")
    sd=env.std(axis=0,ddof=0)
    if np.any(sd<=0):
        raise ValueError("Original S1 five climate variables are invariant")
    cov=(env-env.mean(axis=0))/sd
    cov[:,2]=(env[:,2]-env[:,2].mean())/100
    X=np.column_stack([np.ones(len(rows)),cov])
    if np.linalg.matrix_rank(X)!=6:
        raise ValueError("Original S1 five climate model rank invalid")
    inv=np.linalg.inv(X.T@X)
    beta=inv@(X.T@y)
    residual=y-X@beta
    hat=np.sum((X@inv)*X,axis=1)
    if np.any(1-hat<=1e-10):
        raise ValueError("An original S1 taxon completely determines the fitted model")
    stress=(X* (residual/(1-hat))[:,None])
    robust=inv@(stress.T@stress)@inv
    sderr=np.sqrt(np.diag(robust))
    z=stats.norm.ppf(.975)
    j=3
    if not np.isfinite(sderr[j]) or sderr[j]<=0:
        raise ValueError("Robust S1 rain/photo association is not estimable")
    return {
        "n_equal_weighted_original_S1_source_name_units":len(rows),
        "continuous_S1_photo_response":endpoint,
        "BIO12_effect_original_photo_angle_degrees_per100mm":float(beta[j]),
        "HC3_se":float(sderr[j]),
        "HC3_95percent_normal_ci":[float(beta[j]-z*sderr[j]),float(beta[j]+z*sderr[j])],
        "HC3_two_sided_z_p":float(2*stats.norm.sf(abs(beta[j]/sderr[j]))),
        "source_covariates":[*FIVE],
        "not_original_244_accepted_key_observation_mixed_model":True,
        "photo_angle_is_EXIF_image_vertical_not_botanical_nodding":True,
    }


def evaluate(s1,tree):
    rows=read_source(s1)
    tips,support=read_ITS_tree(tree)
    master={r["original_source_name"]:r for r in rows}
    if not tips<=set(master):
        raise ValueError("ITS gene tree contains species absent from original S1 photo measures")
    matched=[r for r in rows if r["original_source_name"] in tips]
    missing=[r for r in rows if r["original_source_name"] not in tips]
    if len(matched)!=109 or len(missing)!=142:
        raise ValueError("Actual 109 S1 sequence taxon matched out of original 251 is incomplete")
    n_photos=sum(int(r["original_photo_count"]) for r in matched)
    if n_photos!=162146:
        raise ValueError("Original S1 photo measurement represented by own ITS gene-tree changed")
    threshold=float(np.quantile([float(r["chelsa_bio12"]) for r in rows],.75))
    wet=lambda rows:sum(float(r["chelsa_bio12"])>=threshold for r in rows)
    if wet(rows)!=63 or wet(matched)!=16:
        raise ValueError("Frozen all-S1 highly wet tail / ITS source representation drift")
    # Exact 2x2 Fisher test of documented nonrandom selection, exploratory only.
    table=[[wet(matched),len(matched)-wet(matched)],
           [wet(missing),len(missing)-wet(missing)]]
    odds,p=stats.fisher_exact(table,alternative="two-sided")
    models={}
    for cohort,items in (("S1_full251",rows),("GenBank_ITS109",matched),
                         ("unsequenced_S1_142",missing)):
        data={
            "n_original_source_names":len(items),
            "n_original_owner_photo_angle_observations":
                sum(int(r["original_photo_count"]) for r in items),
            "mean_current_S1_BIO12_mm":float(np.mean([float(r["chelsa_bio12"]) for r in items])),
            "mean_original_S1_photo_angle_image_deg":float(
                np.mean([float(r["mean_image_angle_deg"]) for r in items])),
            "equal_taxon_weight_HC3_photo_angle_vs_BIO12":{
                res:fit_equal_species_5axis_HC3(items,res) for res in RESPONSES
            },
        }
        models[cohort]=data
    high=sum(sh>=80 and uf>=95 for sh,uf in support)
    return {
        "schema":"azami_T2_S1_full251_vs_own109_ITS_phylo_why_mismatch_source_qc_20261023_v1",
        "status":"ORIGINAL_S1_PHOTO_SPECIES_SELECTION_AND_GENE_TREE_RESOLUTION_AUDIT",
        "same_original_five_CHELSA_axes":list(FIVE),
        "master_S1_original_source_taxa":SOURCE_NAMES,
        "master_WCVP_accepted_taxon_keys":ACCEPTED_KEYS,
        "all_S1_170066_owner_measured_photo_angles_kept":EXPECTED_IMAGES,
        "actual_GenBank_ITS_ML_original_taxa":NUCLEAR_TIPS,
        "matched_photo_observation_fraction":n_photos/EXPECTED_IMAGES,
        "ITS_nuclear_source_coverage_by_wet_quartile":{
            "BIO12_wettest_quartile_threshold_mm":threshold,
            "n_S1_top_quartile":wet(rows),
            "n_true_ITS_taxa_among_wettest_quartile":wet(matched),
            "n_S1_top_quartile_without_ITS_tree":wet(missing),
            "exploratory_Fisher_source_ITS_coverage_odds_ratio":float(odds),
            "exploratory_Fisher_p":float(p),
            "nonrandom_ITS_coverage_detected_not_evolutionary_selection":True,
        },
        "ITS_ML_branch_support":{
            "source_original_ITS_alignment_sites":645,
            "source_parsimony_informative_sites":103,
            "supported_internal_shalrt_and_UFBoot_count":high,
            "original_ML_support_labeled_internal_branches":len(support),
            "supported_branch_definition":"SH-aLRT>=80 AND UFBoot>=95",
            "not_confirmed_Cirsium_multilocus_species_phylogeny":True,
        },
        "cohorts":models,
        "claim_stop":{
            "original_S1_244_species_observation_level_mixed_effect_not_recomputed":True,
            "S1_109_species_are_NOT_random_sample_of_ALL_251":True,
            "photo_head_angle_not_gravity_or_head_nodding_trait":True,
            "low_phylogenetic_signal_does_not_test_short_term_phenotypic_plasticity":True,
            "independent_adaptation_convergence_or_adaptive_radiation_proven":False,
            "genetic_stability_of_botanical_cirsium_nodding_verified":False,
            "broad_geographic_biome_effect_not_adjusted_in_this_repeat_source_frame":True,
        },
        "source":{
            "S1_original_all251_source_SHA256":EXPECTED_251_SHA256,
            "independent_S1_GenBank109_IQTREE_original_Newick_SHA256":EXPECTED_109_ML_SHA256,
            "S1_original_photo_source_Action":"https://github.com/zuizui0223/azami/actions/runs/39769949628",
            "GenBank_original_ITS_ML_Action":"https://github.com/zuizui0223/azami/actions/runs/39914024548",
        }
    }


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--original251",required=True,type=Path)
    parser.add_argument("--new_its109_tree",required=True,type=Path)
    parser.add_argument("--out",required=True,type=Path)
    options=parser.parse_args()
    answer=evaluate(options.original251,options.new_its109_tree)
    options.out.parent.mkdir(parents=True,exist_ok=True)
    options.out.write_text(json.dumps(answer,indent=2,ensure_ascii=False,
                                      sort_keys=True,allow_nan=False)+"\n",encoding="utf-8")
    print("S1_109_VS_251_PHYLO_TRAIT_DISCORDANCE="+json.dumps({
        "n_photo_genbank_ITS":answer["cohorts"]["GenBank_ITS109"]["n_original_owner_photo_angle_observations"],
        "n_taxa_matched":answer["cohorts"]["GenBank_ITS109"]["n_original_source_names"],
        "n_wet_quartile_matched":answer["ITS_nuclear_source_coverage_by_wet_quartile"]["n_true_ITS_taxa_among_wettest_quartile"],
        "high_supported_ITS_internal_nodes":answer["ITS_ML_branch_support"]["supported_internal_shalrt_and_UFBoot_count"],
        "original251_BIO12_beta":answer["cohorts"]["S1_full251"]["equal_taxon_weight_HC3_photo_angle_vs_BIO12"]["mean_image_angle_deg"]["BIO12_effect_original_photo_angle_degrees_per100mm"],
        "new_ITS109_BIO12_beta":answer["cohorts"]["GenBank_ITS109"]["equal_taxon_weight_HC3_photo_angle_vs_BIO12"]["mean_image_angle_deg"]["BIO12_effect_original_photo_angle_degrees_per100mm"],
    },sort_keys=True),flush=True)


if __name__=="__main__":
    main()
