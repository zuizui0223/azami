"""T2 primary *continuous* S1 image-axis phylogenetic history on all 194 tree taxa.

No categorical inference from image vertical. Reuse the same authenticated
170,066-observation / 244-accepted-key S1 sample and original five exposures.
On 52 archived GBOTB/V.PhyloMaker2 placement scenarios, fit a Gaussian
Brownian ancestry model to all 194 image-measured taxa. Report directional
degree-shift sensitivity, NOT numbers of true evolutionary nodding gains.
"""
from __future__ import annotations

import argparse
from copy import deepcopy
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import Phylo
from scipy.linalg import cho_factor, cho_solve
from scipy.stats import norm, t as student_t

from analysis.v4.joint_s1_t2_phylogeny import (
    ENDPOINT, EXPOSURES, MODEL_PATH, MAPPING_PATH, COVERAGE_PATH,
    actual_native_roster, source_concept_index,
)
from analysis.v4.native_phylogeny import restore_trees

SHA="342d78e9a8a960e8d9bdcf0bc2cdacd13bee48214a8b962419d5a71d5aef46ed"
INC=np.array([1., 100., 100., 10., 1.])
SHIFT_DEGREES=(5.,10.,20.)

def verified_species(table:Path):
    models=json.loads(MODEL_PATH.read_text())
    coverage=json.loads(COVERAGE_PATH.read_text())[ENDPOINT]
    if hashlib.sha256(table.read_bytes()).hexdigest()!=SHA or models["input_sha256"]!=SHA:
        raise ValueError("S1 original decoded source digest changed")
    if models["selection"]["selected"]!=list(EXPOSURES):
        raise ValueError("S1 five covariates cannot be reselected post-hoc")
    frame=pd.read_csv(table,low_memory=False)
    roster=actual_native_roster(frame,models,coverage)
    mapping=source_concept_index(json.loads(MAPPING_PATH.read_text()))
    if len(roster)!=244:
        raise ValueError("S1 has no longer 244 orientation-eligible accepted taxon keys")
    if sum(z["status"]=="matched" for z in mapping.values())!=194:
        raise ValueError("Historical 194-tree mapping drift")
    items=[]
    for key,group in sorted(roster.items()):
        hist=mapping.get(key)
        if not hist or hist["status"]!="matched":
            continue
        tip=list(hist["tree_matches"])
        if len(tip)!=1:
            raise ValueError("Tree tip identity not single-matched")
        angle=float(group["mean_image_angle_deg"])
        median=float(group["median_image_angle_deg"])
        if not np.isfinite(angle) or not np.isfinite(median):
            raise ValueError("Missing original image-angle phenotype")
        if not (0<=angle<=180 and 0<=median<=180):
            raise ValueError("Original image-angle outcome fell outside signed axis [0,180]")
        climate=[float(group["mean_frozen_five_exposures"][e]) for e in EXPOSURES]
        if not np.isfinite(climate).all():
            raise ValueError("Nonfinite source-cohort CHELSA mean")
        items.append({
         "tip":tip[0], "accepted_key":key,
         "taxon_source_names":sorted(group["source_names"]),
         "collision":len(group["source_names"])!=1,
         "n_eligible_S1_images":int(group["n"]),
         "mean_S1_image_axis_deg":angle,
         "median_S1_image_axis_deg":median,
         "S1_five_climate_means":climate,
        })
    if len(items)!=194 or len(set(x["tip"] for x in items))!=194:
        raise ValueError("S1 194 orientation measured PGLS tips changed")
    return items


def audited_pruned_tree(original,names):
    """Preserve inherited stem lengths, regularise zero edges numerically only."""
    tree=deepcopy(original)
    original_tips=[z.name for z in tree.get_terminals()]
    if len(original_tips)!=len(set(original_tips)) or not set(names)<=set(original_tips):
        raise ValueError("S1 species to historical tree mismatch or duplicate")
    for tip in list(tree.get_terminals()):
        if tip.name not in names:
            tree.prune(tip)
    lengths=[node.branch_length for node in tree.find_clades() if node is not tree.root]
    if any(l is None or not np.isfinite(l) or l<0 for l in lengths):
        raise ValueError("Missing/negative/nonfinite original phylogeny edge")
    pos=[l for l in lengths if l>0]
    if not pos:
        raise ValueError("Cannot infer ancestry on a branchless tree")
    epsilon=min(pos)*1e-7
    zeros=0
    for node in tree.find_clades():
        if node is not tree.root and node.branch_length<=0:
            node.branch_length=epsilon
            zeros+=1
    if len(tree.get_terminals())!=len(names):
        raise ValueError("Pruning eliminated a measured S1 tip")
    return tree,{"zero_length_edges_regularised":zeros,
                 "numerical_epsilon_not_biological_length":float(epsilon)}


def shared_path_matrix(tree,tips):
    """All nodes x tip shared-root-path length without repeated MRCA searches."""
    ordered=list(tree.find_clades(order="preorder"))
    index={name:i for i,name in enumerate(tips)}
    if len(index)!=len(tips):
        raise ValueError("Invalid duplicate S1 tip")
    descendants={}
    for node in tree.find_clades(order="postorder"):
        if node.is_terminal():
            descendants[node]=[index[node.name]]
        else:
            descendants[node]=[i for child in node.clades for i in descendants[child]]
    path={}
    depth={tree.root:0.0}
    path[tree.root]=np.zeros(len(tips),dtype=float)
    parent={}
    for node in ordered:
        for child in node.clades:
            parent[child]=node
            edge=float(child.branch_length)
            depth[child]=depth[node]+edge
            row=path[node].copy()
            row[descendants[child]]+=edge
            path[child]=row
    terminal_by_name={n.name:n for n in tree.get_terminals()}
    if set(terminal_by_name)!=set(tips):
        raise ValueError("Source tree taxa differ from authenticated S1 cohort")
    C=np.vstack([path[terminal_by_name[t]] for t in tips])
    if not np.allclose(C,C.T,atol=1e-8):
        raise ValueError("BM shared-root covariance non-symmetric")
    if min(np.linalg.eigvalsh(C))<=0:
        # With positive terminal lengths, a proper BM covariance must be PD.
        raise ValueError("Historical branch covariance not positive definite")
    return C, ordered, path, parent


def continuous_BM(tree,tips,values,*,nugget_fraction=.05,thresholds=SHIFT_DEGREES,cached_paths=None):
    """Fixed-tree continuous Brownian reconstruction with a tip-level nugget.

    Likelihood-derived root and rate. Per-edge conditional marginal probabilities
    are correlated; summed probabilities are expected edge *magnitude exceedances*
    under this particular diffusion model, NOT independently originated nodding.
    """
    if nugget_fraction<0 or nugget_fraction>1:
        raise ValueError("Unsupported tip-level residual nugget")
    y=np.asarray(values,dtype=float)
    if len(y)!=len(tips) or len(y)<3 or not np.isfinite(y).all():
        raise ValueError("Invalid complete continuous S1 source trait")
    C,nodes,path,parent=(cached_paths if cached_paths is not None
                         else shared_path_matrix(tree,tips))
    residual_nugget=nugget_fraction*float(np.mean(np.diag(C)))
    V=C+np.eye(len(tips))*residual_nugget
    ch=cho_factor(V)
    ones=np.ones(len(tips))
    inv1=cho_solve(ch,ones)
    mu=float((inv1@y)/(ones@inv1))
    residual=y-mu
    invres=cho_solve(ch,residual)
    rate=float(residual@invres / (len(tips)-1))
    if rate<=0 or not np.isfinite(rate):
        raise ValueError("Degenerate Brownian diffusion variance")
    edges=[]
    for node in nodes[1:]:
        ancestral=parent[node]
        contrast=path[node]-path[ancestral]
        increment=float(contrast@invres)
        remaining=float(node.branch_length - contrast@cho_solve(ch,contrast))
        if remaining< -1e-6:
            raise ValueError("Invalid Brownian ancestral conditional variance")
        uncertainty=float(np.sqrt(rate*max(0,remaining)))
        probs={}
        for degrees in thresholds:
            threshold=float(degrees)
            p=(float(norm.sf((threshold-increment)/uncertainty))
               if uncertainty>1e-10 else float(increment>threshold))
            probs[str(int(threshold))]=p
        edges.append(dict(
           terminal_edge=bool(node.is_terminal()),
           delta_image_vertical_deg_conditional_mean=increment,
           uncertainty_sd_deg=uncertainty,
           positive_degree_shift_probability=probs,
        ))
    internal=[z for z in edges if not z["terminal_edge"]]
    def aggregate(rows):
        return {str(int(t)):{
            "conditional_expected_edges_increase_by_at_least_deg":
               float(sum(z["positive_degree_shift_probability"][str(int(t))] for z in rows)),
            "posterior_marginal_p_ge_0_95_count":
               int(sum(z["positive_degree_shift_probability"][str(int(t))]>=.95 for z in rows)),
        } for t in thresholds}
    return {
      "n_tips":len(tips),"n_edges":len(edges),
      "n_internal_edges":len(internal),
      "fitted_root_image_axis_deg":mu,
      "Brownian_rate_image_degree_sq_per_branch_unit":rate,
      "tip_observation_nugget_fraction_of_mean_root_to_tip":nugget_fraction,
      "all_edges":aggregate(edges),"internal_edges_only":aggregate(internal),
      "interpretation":"Conditional expected counts of branch magnitude exceedances, not independent discrete origins, not actual gravity-referenced nodding transitions and not trait selection",
    }


def equal_species_five_climate_GLS(tree,tips,trait,climate,*,nugget_fraction=.05,cached_paths=None):
    """Distinct species-mean sensitivity; never substitute for S1 observation fit."""
    C,_,_,_=(cached_paths if cached_paths is not None
             else shared_path_matrix(tree,tips))
    cov=np.asarray(climate,dtype=float)
    y=np.asarray(trait,dtype=float)
    if cov.shape!=(len(tips),5) or len(y)!=len(tips):
        raise ValueError("Five fixed S1 climates not species-aligned")
    X=np.column_stack([np.ones(len(tips)),(cov-cov.mean(axis=0))/INC])
    V=C+np.eye(len(tips))*nugget_fraction*np.mean(np.diag(C))
    ch=cho_factor(V)
    ViX=cho_solve(ch,X)
    covb=cho_solve(cho_factor(X.T@ViX),np.eye(X.shape[1]))
    beta=covb@ViX.T@y
    residual=y-X@beta
    sigma2=float(residual@cho_solve(ch,residual)/(len(y)-X.shape[1]))
    se=np.sqrt(np.diag(covb)*sigma2)
    crit=float(student_t.ppf(.975,len(y)-X.shape[1]))
    return {
      "n_tips":len(tips),"df":len(y)-6,
      "status":"same_source_equal_species_Brownian_GLS_five_climate_DIAGNOSTIC_not_primary_S1",
      "per_exposure_increment":dict(zip(EXPOSURES,INC.tolist())),
      "coefficients_image_axis_deg_per_increment":{
          name:{"estimate":float(beta[i+1]),"lower95":float(beta[i+1]-crit*se[i+1]),
                "upper95":float(beta[i+1]+crit*se[i+1])}
          for i,name in enumerate(EXPOSURES)},
      "not_observation_weighted_original_S1":True,
      "not_a_causes_or_origins_test":True,
    }


def analyze(items,trees,*,thresholds=SHIFT_DEGREES,expected_taxa=194):
    names=sorted(x["tip"] for x in items)
    byname={x["tip"]:x for x in items}
    if len(names)!=expected_taxa or len(trees)!=52 or len(set(names))!=len(names):
        raise ValueError("Wrong source-cohort size, duplicate tips, or 52 tree scenarios")
    climate=np.array([byname[n]["S1_five_climate_means"] for n in names])
    mean=np.array([byname[n]["mean_S1_image_axis_deg"] for n in names])
    median=np.array([byname[n]["median_S1_image_axis_deg"] for n in names])
    source_collisions=[n for n in names if byname[n]["collision"]]
    results=[]
    for idx,(label,original) in enumerate(trees):
        pruned,diag=audited_pruned_tree(original,names)
        # Same exact 194 taxa and topology; compute and validate covariance
        # only once for the continuous mean, median and 5-environment GLS.
        cached=shared_path_matrix(pruned,names)
        fit=continuous_BM(pruned,names,mean,thresholds=thresholds,cached_paths=cached)
        medfit=continuous_BM(pruned,names,median,thresholds=thresholds,cached_paths=cached)
        spec=equal_species_five_climate_GLS(pruned,names,mean,climate,cached_paths=cached)
        results.append({"scenario":label,"scope":f"{expected_taxa}_measured_S1_species_image_axis",
              "branch_length":diag,"species_mean":fit,
              "species_median":medfit,"five_S1_climate_GLS":spec})
    def bounds(vals):
        return [float(min(vals)),float(max(vals))]
    summary={
      "status":f"EMPIRICAL_{expected_taxa}_S1_CONTINUOUS_IMAGE_ORIENTATION_PHYLOGENY",
      "n_original_S1_eligible_taxon_keys":244,
      "n_tree_taxon_keys_with_original_continuous_orientation":expected_taxa,
      "n_original_images_in_S1":170066,
      "n_PGLS_historical_placement_scenarios":len(trees),
      "continuous_phenotype":"EXIF_oriented_image_vertical_angle_degrees_NOT_gravity",
      "measurement_is_continuous_not_binary":True,
      "S1_original_five_exposures":list(EXPOSURES),
      "source_name_collisions_within_modeled_cohort":len(source_collisions),
      "source_name_collisions_within_194":len(source_collisions) if expected_taxa==194 else 0,
      "source_name_collision_tips":source_collisions,
      "source_collision_biological_interpretation":(
          "3 original alias collisions retained ONLY to reproduce legacy S1 placement sensitivity"
          if source_collisions else "no source-name collision tips retained in sensitivity"),
      "n_BM_nonindependent_branch_shifts_cannot_equal_discrete_origins":True,
      "evolutionary_nodding_gain_count_claimed":False,
      "study_scope":"original 52 GBOTB grafted placements, not independently sequenced global Cirsium species trees",
      "reused_data_no_frozen_S1_refit":True,
    }
    for state in ("species_mean","species_median"):
        for class_ in ("all_edges","internal_edges_only"):
            for threshold in thresholds:
                key=str(int(threshold))
                summary[f"{state}_{class_}_increase_ge_{key}deg_conditional_expected_edge_scenario_range"]=bounds([
                    r[state][class_][key]["conditional_expected_edges_increase_by_at_least_deg"]
                    for r in results])
                summary[f"{state}_{class_}_increase_ge_{key}deg_high_probability_branch_count_range"]=bounds([
                    r[state][class_][key]["posterior_marginal_p_ge_0_95_count"]
                    for r in results])
    for env in ("chelsa_bio12","chelsa_bio15","chelsa_bio01"):
        for field in ("estimate","lower95","upper95"):
            summary[f"diagnostic_same{expected_taxa}_GLS_{env}_{field}_range"]=bounds([
              r["five_S1_climate_GLS"]["coefficients_image_axis_deg_per_increment"][env][field]
              for r in results])
    return {"summary":summary, "52_scenario_details":results,
            "source_cohort_194_species":[{"tip":x["tip"],"accepted_key":x["accepted_key"],
                    "n_eligible_images":x["n_eligible_S1_images"],
                    "image_mean_deg":x["mean_S1_image_axis_deg"],
                    "image_median_deg":x["median_S1_image_axis_deg"],
                    "S1_five_climate_means":dict(zip(EXPOSURES,[float(z) for z in x["S1_five_climate_means"]])),
                    "source_concept_collision":x["collision"]}
                for x in sorted(items,key=lambda z:z["tip"])]}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--table",type=Path,required=True)
    parser.add_argument("--out",type=Path,required=True)
    args=parser.parse_args()
    items=verified_species(args.table)
    trees=restore_trees()
    results=analyze(items,trees)
    # Collision-free biological-taxon sensitivity: reuse ALL 52 original
    # placements, but do not treat original WCVP-merged names as single species.
    clean=[record for record in items if not record["collision"]]
    if len(clean)!=191:
        raise ValueError("Expected 194 mapped source keys minus exactly three unresolved aliases")
    sensitivity=analyze(clean,trees,expected_taxa=191)
    if sensitivity["summary"]["source_name_collisions_within_modeled_cohort"]!=0:
        raise ValueError("Collision-free ancestral sensitivity still includes merged names")
    results["collision_free_191_taxon_sensitivity"]=sensitivity["summary"]
    # The original S1 cohort is unevenly sampled. Species means supported by
    # only 1-4 photographs must never be treated as equally precise as those
    # with hundreds; this exclusion ladder is a *sensitivity*, not a new S1.
    sample_sensitivity={}
    for minimum, expected_n in ((5,171),(10,143)):
        selected=[record for record in clean
                  if record["n_eligible_S1_images"]>=minimum]
        if len(selected)!=expected_n:
            raise ValueError(f"Original S1 image-count >= {minimum} cohort changed")
        subset=analyze(selected,trees,expected_taxa=expected_n)
        sample_sensitivity[str(minimum)]={
          "n_original_S1_continuous_tips":expected_n,
          "n_source_name_collisions":0,
          "n_original_S1_images":int(sum(r["n_eligible_S1_images"] for r in selected)),
          "conditioning":"same pinned S1 source and same 52 grafted scenario trees",
          "not_independent_replication":True,
          "summary":subset["summary"],
        }
    results["minimum_S1_image_count_sensitivity"]=sample_sensitivity
    results["source_audit"]={
       "source_sha256":SHA,"native_observation_source_authenticated":True,
       "source_selection_sha256":hashlib.sha256(MODEL_PATH.read_bytes()).hexdigest(),
       "pgls_tree_bundle":"SHA-verified Zenodo 22703215 (existing native_phylogeny.restore_trees)",
       "tip_photometric_gravity_calibration":False,
       "nonindependent_52_tree_ensemble":True,
       "all_194_biological_nodding_state_verified":False,
    }
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(results,indent=2,ensure_ascii=False,allow_nan=False)+"\n",
                        encoding="utf-8")
    public_summary={**results["summary"],
       "collision_free_191_taxon_sensitivity":sensitivity["summary"]}
    print("T2_ALL_S1_CONTINUOUS_194_SUMMARY="+json.dumps(public_summary,
          sort_keys=True,ensure_ascii=False),flush=True)


if __name__=="__main__":
    main()
