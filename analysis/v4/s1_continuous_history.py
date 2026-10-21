"""EAzami-style continuous evolution on *all* eligible azami S1 orientations.

Primary phenotype is the already-measured EXIF-relative angle, never inferred
botanical erect/nodding labels. Every comparison uses the frozen S1 five
CHELSA predictors (BIO1, BIO4, BIO12, BIO15, RSDS) and the *same exact*
taxon-cohort from the 52 archived GBOTB/V.PhyloMaker2 placement scenarios.

Analyses:
  - 244 native S1 accepted taxa and all observation-level source-count audit.
  - All individually matched, collision-free tree species (up to 194).
  - Brownian root/ancestral continuous trait reconstruction from PIC contrasts.
  - Independent-contrast 5-predictor regressions, across 52 TREE SCENARIOS;
    mean angle primary, median angle observational sensitivity.
  - Tip coverage >= 5, >= 10 observations sensitivity; no threshold-coded
    biological downward events or invented biological ancestry.

All tree-specific p values/intervals condition on GRAFTED placement, branch
lengths, Brownian evolution, known taxon means, and no camera error; they are
exploratory, not evolutionary origin or rain-adaptation evidence.
"""
from __future__ import annotations

import argparse
from copy import deepcopy
import hashlib
import json
from pathlib import Path

from Bio import Phylo
import numpy as np
import pandas as pd
from scipy.stats import t as student_t

from analysis.v4.joint_s1_t2_phylogeny import (
    ENDPOINT, EXPOSURES, MODEL_PATH, MAPPING_PATH, COVERAGE_PATH,
    actual_native_roster, source_concept_index,
)
from analysis.v4.native_phylogeny import restore_trees
from analysis.v4.t2_pgls_placement_screen import check_placement_tree

TRAITS=("image_angle_mean_deg","chelsa_bio01","chelsa_bio04",
        "chelsa_bio12","chelsa_bio15","chelsa_rsds_mean")
EXTRA_TRAIT="image_angle_median_deg"
MIN_EFFECTIVE_BRANCH=1e-8
SOURCE_SHA="342d78e9a8a960e8d9bdcf0bc2cdacd13bee48214a8b962419d5a71d5aef46ed"


def taxon_roster(native_table, *, expected_source_sha=SOURCE_SHA):
    data=Path(native_table)
    if hashlib.sha256(data.read_bytes()).hexdigest()!=expected_source_sha:
        raise ValueError("Native 244-taxon input SHA does not match S1's original decoded source")
    models=json.loads(MODEL_PATH.read_text())
    if (models["input_sha256"]!=expected_source_sha or
        tuple(models["selection"]["selected"])!=EXPOSURES):
        raise ValueError("Frozen S1 source/co-adjusted five climate axes drift")
    coverage=json.loads(COVERAGE_PATH.read_text())[ENDPOINT]
    sample=pd.read_csv(data,low_memory=False)
    source=actual_native_roster(sample,models,coverage)
    if len(source)!=244 or sum(r["n"] for r in source.values())!=170066:
        raise ValueError("S1 native-eligible sample does not reproduce 244/170066")
    mapped=source_concept_index(json.loads(MAPPING_PATH.read_text()))
    if sum(r["status"]=="matched" for r in mapped.values())!=194:
        raise ValueError("194 original tree-matched keys changed")
    rows=[]
    for key,record in sorted(source.items()):
        names=record["source_names"]
        exact_tip=None
        status="not_on_52_placement_trees"
        if len(names)!=1:
            status="excluded_collided_or_multiple_source_names"
        elif key in mapped and mapped[key]["status"]=="matched":
            match=mapped[key]["tree_matches"]
            candidate=names[0].replace(" ","_")
            if match=={candidate}:
                exact_tip=candidate
                status="exact_independent_source_name_tip_match"
            else:
                status="excluded_taxon_key_tip_mismatch"
        climate=record["mean_frozen_five_exposures"]
        if set(climate)!=set(EXPOSURES) or not all(np.isfinite(float(x)) for x in climate.values()):
            raise ValueError("Missing frozen climate exposures in native species centroid")
        rows.append({
            "accepted_key":key,
            "taxon_source_names":names,
            "tree_tip":exact_tip,
            "join_status":status,
            "S1_observations":record["n"],
            "image_angle_mean_deg":record["mean_image_angle_deg"],
            "image_angle_median_deg":record["median_image_angle_deg"],
            **climate,
        })
    if len(rows)!=244 or len(set(x["accepted_key"] for x in rows))!=244:
        raise ValueError("Duplicate or missing S1 native taxon keys")
    good=[r for r in rows if r["tree_tip"] is not None]
    if len(set(r["tree_tip"] for r in good))!=len(good):
        raise ValueError("Distinct S1 accepted keys collapse onto the same tree tip")
    return rows, {"full_244_species":len(rows),"full_S1_observations":sum(r["S1_observations"] for r in rows),
                  "original_S1_52_tree_keys":194,
                  "exact_noncollided_phylogeny_tips":len(good),
                  "missing_due_to_tree_or_concept":244-len(good),
                  "exclusion_counts":dict(sorted((key,sum(r["join_status"]==key for r in rows))
                                            for key in {r["join_status"] for r in rows})),
                  "source_names_and_climate_from_authenticated_native_S1":True}


def _branch_length(child):
    b=child.branch_length
    if b is None or not np.isfinite(b) or b<0:
        raise ValueError("Missing/negative/nonfinite originally supplied branch length")
    return max(float(b),MIN_EFFECTIVE_BRANCH)


def ancestral_PIC(tree, by_tip):
    """Brownian Gaussian pruning and independent contrasts for k-ary nodes.

    Every node is estimated from its DESCENDANT tips (partial likelihood).
    Root is the GLS ancestral estimate with diffuse root prior. Contrasts are
    orthogonal under a conditional Brownian model; multifurcations are
    deterministically paired by child tip name (orthogonal reparameterization).
    """
    tipset={x.name for x in tree.get_terminals()}
    if tipset!=set(by_tip):
        raise ValueError("PIC tree tips and actual S1 phenotype cohort differ")
    k=len(next(iter(by_tip.values())))
    if not 2<=k<=8:
        raise ValueError("Expected continuous response with selected climates")
    means={}
    variances={}
    contrasts=[]
    adjustments=0
    original_zero=0
    for node in tree.find_clades(order="postorder"):
        if node.is_terminal():
            value=np.asarray(by_tip[node.name],dtype=float)
            if value.shape!=(k,) or not np.isfinite(value).all():
                raise ValueError("Nonfinite or missing true S1 tip value")
            means[node]=value
            variances[node]=0.0
            continue
        children=sorted(node.clades,key=lambda c:min(t.name for t in c.get_terminals()))
        if len(children)<2:
            raise ValueError("Unary phylogenetic node: prune and clean tree first")
        first=children[0]
        center=means[first].copy()
        v=variances[first]+_branch_length(first)
        original_zero+=int(first.branch_length==0)
        for child in children[1:]:
            raw_length=child.branch_length
            original_zero+=int(raw_length==0)
            child_v=variances[child]+_branch_length(child)
            delta=(center-means[child])/np.sqrt(v+child_v)
            contrasts.append(delta)
            merged_v=v*child_v/(v+child_v)
            center=(center*child_v+means[child]*v)/(v+child_v)
            v=merged_v
            adjustments+=1
        means[node]=center
        variances[node]=v
    pic=np.array(contrasts,dtype=float)
    if pic.shape!=(len(tipset)-1,k) or not np.isfinite(pic).all():
        raise ValueError("Felsenstein PIC dimension mismatch; no fake phylogenetic replicates")
    return {"PIC":pic,"root":means[tree.root].tolist(),
            "root_variance_branch_units":float(variances[tree.root]),
            "n_tip":len(tipset),"n_independent_contrasts_within_scenario":len(pic),
            "original_zero_length_branches":int(original_zero),
            "zero_branch_numeric_regularization":MIN_EFFECTIVE_BRANCH}


def _holm(p):
    p=np.asarray(p,float)
    order=np.argsort(p)
    sorted_adj=np.maximum.accumulate(np.minimum(1.0,(len(p)-np.arange(len(p)))*p[order]))
    out=np.empty(len(p))
    out[order]=sorted_adj
    return out.tolist()


def fit_PIC_multi(pic):
    """No-intercept contrasts regression with all fixed S1 predictors present."""
    a=np.asarray(pic,dtype=float)
    n,k=a.shape
    if k!=6 or n<25:
        raise ValueError("Too few taxon contrasts or wrong five-S1 climate design")
    y=a[:,0]
    X=a[:,1:]
    if np.linalg.matrix_rank(X)<5 or np.linalg.cond(X)>1e8:
        raise ValueError("Five-exposure independent-contrast design not identifiable")
    beta=np.linalg.lstsq(X,y,rcond=None)[0]
    resid=y-X@beta
    df=n-5
    variance=float(np.dot(resid,resid)/df)
    if not np.isfinite(variance) or variance<0:
        raise ValueError("Invalid Brownian conditional residual uncertainty")
    gram=np.linalg.inv(X.T@X)
    se=np.sqrt(np.diag(gram)*variance)
    if np.any(se<=0):
        raise ValueError("Degenerate contrast standard errors")
    tstat=beta/se
    p=2*student_t.sf(np.abs(tstat),df)
    crit=float(student_t.ppf(.975,df))
    return {
        "coefficients_per_1_cohort_species_SD_climate":{
            exp:float(x) for exp,x in zip(EXPOSURES,beta)},
        "pointwise_SE":{exp:float(x) for exp,x in zip(EXPOSURES,se)},
        "pointwise_model_P_conditional_on_grafted_tree":{
            exp:float(x) for exp,x in zip(EXPOSURES,p)},
        "within_scenario_five_variable_Holm_P":{
            exp:float(x) for exp,x in zip(EXPOSURES,_holm(p))},
        "pointwise_95pct_interval":{exp:[float(beta[i]-crit*se[i]),
                                       float(beta[i]+crit*se[i])]
                                   for i,exp in enumerate(EXPOSURES)},
        "df":int(df),
        "contrast_n":int(n),
        "Brownian_PIC_residual_variance":variance,
        "not_an_evolutionary_transition_rate_or_origin_count":True,
    }


def fit_full_species_OLS(roster):
    """All 244 S1 source taxon means, one observation per authority concept.

    This is a new descriptive taxon-level estimand, not a refit or replacement
    for the original observation-level mixed/cell-weighted S1 result.
    """
    if len(roster)!=244:
        raise ValueError("Full S1 species-level context requires all 244 keys")
    climates=np.array([[float(r[v]) for v in EXPOSURES] for r in roster])
    sd=climates.std(axis=0)
    if np.any(sd<=0):
        raise ValueError("Degenerate full source climate")
    X=np.column_stack([np.ones(len(roster)),(climates-climates.mean(axis=0))/sd])
    y=np.array([r["image_angle_mean_deg"] for r in roster])
    beta=np.linalg.lstsq(X,y,rcond=None)[0]
    return {
        "n_taxa":len(roster),
        "coefficient_deg_per_1_species_SD":{
            v:float(x) for v,x in zip(EXPOSURES,beta[1:])},
        "source":"same_S1_244_taxon_centroids_exploratory_equal_species_not_original_S1_observation_model",
    }


def analyze_scenario(name,original,eligible,*,min_obs):
    included=[r for r in eligible if r["S1_observations"]>=min_obs]
    by_name={r["tree_tip"]:r for r in included}
    if len(included)<35:
        return {"status":"not_enough_high_observation_taxa",
                "n_taxa":len(included)}
    tree=deepcopy(original)
    tree_tips=check_placement_tree(tree)
    if not set(by_name)<=tree_tips:
        raise ValueError("Original S1 source-name crosswalk not present in frozen tree")
    for leaf in list(tree.get_terminals()):
        if leaf.name not in by_name:
            tree.prune(leaf)
    if len(tree.root.clades)<2:
        raise ValueError("Pruning yielded unresolved unary root")
    tiporder=sorted(by_name)
    climate=np.array([[by_name[t][x] for x in EXPOSURES] for t in tiporder],float)
    mean=climate.mean(axis=0)
    sd=climate.std(axis=0)
    if np.any(sd<=1e-12):
        raise ValueError("No climate variance in 5-exposure eligible PIC cohort")
    z=(climate-mean)/sd
    mean_values={
        t:np.array([float(by_name[t]["image_angle_mean_deg"]),*z[i]],float)
        for i,t in enumerate(tiporder)}
    med_values={
        t:np.array([float(by_name[t]["image_angle_median_deg"]),*z[i]],float)
        for i,t in enumerate(tiporder)}
    mean_PIC=ancestral_PIC(tree,mean_values)
    median_PIC=ancestral_PIC(tree,med_values)
    main=fit_PIC_multi(mean_PIC["PIC"])
    median=fit_PIC_multi(median_PIC["PIC"])
    root=mean_PIC["root"]
    root_sd=np.sqrt(float(np.mean(mean_PIC["PIC"][:,0]**2)*
                            mean_PIC["root_variance_branch_units"]))
    return {
        "status":"conditional_Brownian_continuous_S1_photo_presentation",
        "scenario":name,
        "min_original_S1_photos_per_taxon":min_obs,
        "tree_tips_n":len(included),
        "PIC_contrasts_n":mean_PIC["n_independent_contrasts_within_scenario"],
        "root_photo_angle_deg_estimate":float(root[0]),
        "root_photo_angle_deg_Brownian_conditional_sd":float(root_sd),
        "root_five_climate_reconstruction_in_cohort_standardized_units":{
            e:float(x) for e,x in zip(EXPOSURES,root[1:])},
        "n_original_zero_length_branch_regularizations":mean_PIC["original_zero_length_branches"],
        "mean_image_angle":main,
        "median_image_angle_sensitivity":median,
        "no_binary_downward_origin_count":True,
        "root_not_gravity_referenced_or_dated":True,
    }


def execute(native_table, output_dir):
    """Reuse the original authenticated private data, archive ONLY species summaries."""
    roster,coverage=taxon_roster(native_table)
    eligible=[r for r in roster if r["tree_tip"] is not None]
    if coverage["exact_noncollided_phylogeny_tips"]<100:
        raise ValueError("Unexplained collapse from 194 original tree-eligible source taxa")
    trees=restore_trees()
    if len(trees)!=52:
        raise ValueError("Not the archived 52 PGLS sensitivity trees")
    full=fit_full_species_OLS(roster)
    all_results=[]
    for min_n in (1,5,10):
        for number,(name,phy) in enumerate(trees):
            entry=analyze_scenario(name,phy,eligible,min_obs=min_n)
            entry["tree_index"]=number
            all_results.append(entry)
    summary={
        "analysis":"EAzami_style_CONTINUOUS_S1_tip_phylogenetic_PIC_and_BM_ancestral_photo_angle",
        "source":"same_20261020_native_S1_244_170066_sha256_locked",
        "source_SHA256":SOURCE_SHA,
        "n_native_taxa":244,"n_native_photos":170066,
        "n_original_PGLS_matched_taxa":194,
        "primary_eligible_canonical_taxa":len(eligible),
        "n_placement_scenarios":52,
        "taxon_eligibility":coverage,
        "all244_equal_species_context":full,
        "covariates":list(EXPOSURES),
        "S1_original_analysis_not_refitted":True,
        "phenotype_original":"angle_relative_to_EXIF_image_vertical_NOT_true_gravitational_head_position",
        "continuous_trait_used_without_erect_nodding_threshold":True,
        "EAzami_like_phases_executed":[
            "continuous_Brownian_root_reconstruction",
            "phylogenetically_independent_contrasts",
            "joint_all_five_S1_climates",
            "same_tip_mean_vs_median_sensitivity",
            "52_placement_tree_sensitivity",
            "observation_count_eligibility_sensitivity"],
        "S1_244_not_all_tree_matched":True,
        "no_true_origin_event_count_or_absolute_dating":True,
        "no_evolutionary_rain_dependent_transition_rate_estimated":True,
        "not_prospective_EAzami_replication":True,
    }
    for min_n in (1,5,10):
        arr=[r for r in all_results if r["min_original_S1_photos_per_taxon"]==min_n
             and r["status"].startswith("conditional_")]
        summary[f"min_{min_n}_n_scenarios"]=len(arr)
        summary[f"min_{min_n}_n_taxa_range"]=(
            [min(r["tree_tips_n"] for r in arr),max(r["tree_tips_n"] for r in arr)]
            if arr else None)
        if not arr:
            continue
        summary[f"min_{min_n}_root_photo_angle_deg_range"]=[
            min(r["root_photo_angle_deg_estimate"] for r in arr),
            max(r["root_photo_angle_deg_estimate"] for r in arr)]
        for var in EXPOSURES:
            slopes=[r["mean_image_angle"]["coefficients_per_1_cohort_species_SD_climate"][var]
                    for r in arr]
            adjusted=[r["mean_image_angle"]["within_scenario_five_variable_Holm_P"][var]
                      for r in arr]
            ints=[r["mean_image_angle"]["pointwise_95pct_interval"][var]
                  for r in arr]
            median=[r["median_image_angle_sensitivity"]["coefficients_per_1_cohort_species_SD_climate"][var]
                    for r in arr]
            summary[f"min_{min_n}_{var}"]={
                "mean_angle_PIC_slope_deg_per_cohort_climate_SD_scenario_range":[
                    float(min(slopes)),float(max(slopes))],
                "median_angle_PIC_slope_deg_per_cohort_climate_SD_scenario_range":[
                    float(min(median)),float(max(median))],
                "positive_95pct_conditional_tree_intervals":sum(x[0]>0 for x in ints),
                "negative_95pct_conditional_tree_intervals":sum(x[1]<0 for x in ints),
                "scenario_Holm_P_below_0p05":sum(x<.05 for x in adjusted),
            }
    root=Path(output_dir);root.mkdir(parents=True,exist_ok=True)
    (root/"S1_244_taxon_continuous_source_roster.json").write_text(
        json.dumps({"contract":"nonpersonal_single_taxon_S1_mean_median_and_five_CHELSA",
                    "coverage":coverage,"taxa":roster},indent=2,ensure_ascii=False,
                   allow_nan=False)+"\n")
    (root/"S1_194_continuous_EAzami_style_history.json").write_text(
        json.dumps({"summary":summary,"tree_results":all_results},indent=2,
                   ensure_ascii=False,allow_nan=False)+"\n")
    print("V4_S1_CONTINUOUS_EAZAMI_STYLE_SUMMARY="+json.dumps(summary,sort_keys=True,
          ensure_ascii=False,allow_nan=False),flush=True)
    return summary


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--native-table",required=True)
    p.add_argument("--output-dir",required=True)
    args=p.parse_args()
    execute(args.native_table,args.output_dir)


if __name__=="__main__":
    main()
