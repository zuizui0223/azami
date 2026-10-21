"""One S1/T2 source-lane-controlled five-climate present-state comparison.

Source is the immutable, source-authenticated S1 244→194→46 taxon ledger.
Binary posture comes ONLY from Siniscalchi 2023 or EAzami authorities; angle
is always EXIF image-relative. Grafted PGLS trees are sensitivity scenarios,
not independent evidence of biological origins. No S1 numerical refit.
"""
from __future__ import annotations

import argparse
from copy import deepcopy
import hashlib
import io
import itertools
import json
from pathlib import Path
import urllib.request

from Bio import Phylo
import numpy as np
from scipy.linalg import cho_factor, cho_solve

from analysis.v4.joint_s1_t2_phylogeny import EXPOSURES
from analysis.v4.native_phylogeny import tree_covariance
from analysis.v4.t2_pgls_placement_screen import (
    SOURCE_REPO, SOURCE_COMMIT, AUTHOR_NUCLEAR_TREE_PATH,
    AUTHOR_NUCLEAR_TREE_GIT_BLOB, git_blob_sha, check_placement_tree,
)

RECEIPT=Path("analysis/v4/results/20261021_unified_abcde/S1_T2_JOINT_244_194_VERIFIED_20261021.json")
SOURCE_ARTIFACT_SHA="sha256:aa5a9ac298b8e36c9e96da6f1eb96319733136e2e3a105305092f7da2ca97bff"
NA="Siniscalchi2023_author_scored"
EA="EAzami_FDT4_existing_authority"
AXES=("chelsa_bio12","chelsa_bio15")


def authenticate_input(raw, receipt):
    inv=raw.get("T2_independent_state_inventory",{})
    scope=raw.get("model_scope",{})
    provenance=raw.get("provenance",{})
    if (raw.get("status")!="empirical_full_244_native_join_source_sha_verified"
        or raw.get("fixed_covariates")!=list(EXPOSURES)
        or scope.get("orientation_eligible_accepted_keys")!=244
        or scope.get("orientation_eligible_observations")!=170066
        or scope.get("placement_tree_covered_keys")!=194
        or scope.get("placement_tree_covered_observations")!=169912
        or scope.get("grafted_placement_scenarios")!=52
        or receipt.get("workflow_run_id")!=39584744194
        or receipt.get("workflow_artifact_digest")!=SOURCE_ARTIFACT_SHA
        or provenance.get("S1_model_source_sha256")!=receipt["source_provenance"]["native_sha256"]
        or provenance.get("EAzami_source_git_blob")!=receipt["source_provenance"]["EAzami_fdt4_git_blob"]
        or provenance.get("NA_source_git_blob")!=receipt["source_provenance"]["NA_head_position_git_blob"]):
        raise ValueError("Authenticated exact S1 244/194 observation-and-source contract changed")
    rows=[r for r in inv["taxon_source_audit"] if r.get("join_status")=="exact_S1_and_tree"]
    if len(rows)!=46:
        raise ValueError("Expected 46 external botany states matched to S1 and PGLS")
    keys=[str(r.get("accepted_key")) for r in rows]
    names=[r["source_taxon"].replace(" ","_") for r in rows]
    if len(set(keys))!=46 or len(set(names))!=46:
        raise ValueError("Multiple source names or keys used as one independent tip")
    bysource={(src,state):sum(r["source"]==src and r["state"]==state for r in rows)
              for src in (NA,EA) for state in ("nodding","erect")}
    if bysource!={(NA,"nodding"):3,(NA,"erect"):39,
                  (EA,"nodding"):3,(EA,"erect"):1}:
        raise ValueError("NA/EA fixed phenotype source composition changed")
    for row in rows:
        v=row.get("s1_five_climate_centroids",{})
        if (not row.get("s1_phylogeny_194_key_match") or
           set(v)!=set(EXPOSURES) or
           not all(np.isfinite(float(v[k])) for k in EXPOSURES)):
            raise ValueError("Missing exact published five-climate S1 taxon source")
    return rows


def cohort(rows):
    """Source lane is a fixed effect; NA-only uses no all-zero extra column."""
    rows=sorted(rows,key=lambda row:row["source_taxon"])
    names=[r["source_taxon"].replace(" ","_") for r in rows]
    y=np.array([float(r["state"]=="nodding") for r in rows],dtype=float)
    lane=np.array([int(r["source"]==EA) for r in rows],dtype=int)
    exposure=np.array([[float(r["s1_five_climate_centroids"][axis])
                        for axis in EXPOSURES] for r in rows])
    sd=exposure.std(axis=0)
    if np.any(sd<1e-12):
        raise ValueError("Degenerate original five-climate S1 predictors")
    z=(exposure-exposure.mean(axis=0))/sd
    has_region=bool(lane.any() and (1-lane).any())
    X=np.column_stack(([np.ones(len(rows)),lane,z] if has_region
                       else [np.ones(len(rows)),z]))
    if np.linalg.matrix_rank(X)!=X.shape[1] or np.linalg.cond(X)>1e6:
        raise ValueError("Unidentifiable five-environment source-adjusted model")
    return {"names":names,"y":y,"lane":lane,"X":X,
            "predictor_offset":2 if has_region else 1,
            "n":len(rows),"n_nodding":int(y.sum()),
            "NA":{"n":int((lane==0).sum()),"D":int(y[lane==0].sum())},
            "EA":{"n":int((lane==1).sum()),"D":int(y[lane==1].sum())}}


def exact_source_maps(data):
    """All exchangeable 0/1 assignments preserving source lane D counts."""
    y=data["y"]; group=data["lane"]
    strata=[np.where(group==s)[0].tolist() for s in sorted(set(group))]
    bins=[list(itertools.combinations(index,int(sum(y[index]))))
          for index in strata]
    total=int(np.prod([len(c) for c in bins]))
    if total>100000:
        raise ValueError("Exact region/source conditioned permutation universe too large")
    null=np.zeros((total,len(y)),dtype=float)
    for i,combos in enumerate(itertools.product(*bins)):
        for select in combos:
            null[i,list(select)]=1.
    if np.sum(np.all(null==y,axis=1))!=1:
        raise ValueError("Observed botanical composition must be in null exactly once")
    return null, total


def score(data,labels,covariance=None):
    """GLS coefficients are descriptive binary LPM; permutation is sharp null.

    Important: condition on NA/EA D prevalences; not a *conditional* test
    for one climate axis given a nonzero effect of the other climate axes.
    """
    X=data["X"]; n=len(X)
    if covariance is None:
        VI_X=X.copy();model="equal_species_fixed_source"
    else:
        A=np.asarray(covariance,dtype=float)
        if A.shape!=(n,n) or not np.isfinite(A).all():
            raise ValueError("Invalid phylogeny covariance for matched taxa")
        A=(A+A.T)/2+np.eye(n)*0.05
        A_ch=cho_factor(A)
        VI_X=cho_solve(A_ch,X)
        model="Brownian_covariance_fixed_5pct_nugget"
    W=cho_solve(cho_factor(X.T@VI_X),VI_X.T)
    columns=[data["predictor_offset"]+list(EXPOSURES).index(k) for k in AXES]
    weights=W[columns]
    observed=weights@data["y"]
    null=labels@weights.T
    centered=null.mean(axis=0);std=null.std(axis=0)
    if np.any(std<1e-12):
        raise ValueError("Degenerate exact source-constrained null")
    zobs=np.abs((observed-centered)/std)
    znull=np.abs((null-centered)/std)
    raw=(znull>=zobs).mean(axis=0)
    corrected=(znull.max(axis=1)[:,None]>=zobs[None,:]).mean(axis=0)
    return {"model":model,"n_taxa":data["n"],
            "n_downward":data["n_nodding"],
            "source_prevalence":{"NA":data["NA"],"EA":data["EA"]},
            "two_registered_precipitation_axes":list(AXES),
            "all_five_covariates_in_model":list(EXPOSURES),
            "estimated_linear_probability_coefficient_per_1SD":observed.tolist(),
            "source_stratified_exact_sharp_null_two_sided_p":raw.tolist(),
            "two_axis_maxT_source_stratified_exact_p":corrected.tolist(),
            "null_maps":len(null),
            "interpretation":"Conditional on source-lane counts under global sharp null of state/environment exchangeability; NOT axis-specific partial p if other exposure effects exist; no adaptation inference"}


def tree_covariance_for(tree,names):
    phy=deepcopy(tree)
    if not set(names)<=check_placement_tree(phy):
        raise ValueError("An independently-scored S1 taxon absent from tree")
    for node in list(phy.get_terminals()):
        if node.name not in names:
            phy.prune(node.name)
    return tree_covariance(phy,names)


def author_nuclear_tree():
    url=f"https://raw.githubusercontent.com/{SOURCE_REPO}/{SOURCE_COMMIT}/{AUTHOR_NUCLEAR_TREE_PATH}"
    with urllib.request.urlopen(url,timeout=90) as response:
        raw=response.read()
    if git_blob_sha(raw)!=AUTHOR_NUCLEAR_TREE_GIT_BLOB:
        raise ValueError("Author nuclear species tree Git blob drift")
    tree=Phylo.read(io.StringIO(raw.decode("utf-8-sig")),"newick")
    if len(tree.get_terminals())!=63:
        raise ValueError("Expected published 63-taxon nuclear reference")
    return tree


def run(joint,output,*,include_scenarios=True):
    from analysis.v4.native_phylogeny import restore_trees
    data=json.loads(Path(joint).read_text(encoding="utf-8"))
    receipt=json.loads(RECEIPT.read_text(encoding="utf-8"))
    rows=authenticate_input(data,receipt)
    whole=cohort(rows)
    na=cohort([r for r in rows if r["source"]==NA])
    if (whole["n"],whole["n_nodding"],na["n"],na["n_nodding"])!=(46,6,42,3):
        raise ValueError("Source-joined T2 cohort not as frozen")
    labels,nmaps=exact_source_maps(whole)
    nlabs,n_na_maps=exact_source_maps(na)
    if nmaps!=45920 or n_na_maps!=11480:
        raise ValueError("Unexpected source-count conditioned exact label universe")
    nuclear=author_nuclear_tree()
    nuclear42=score(na,nlabs,tree_covariance_for(nuclear,na["names"]))
    base46=score(whole,labels)
    base42=score(na,nlabs)
    results=[]
    if include_scenarios:
        trees=restore_trees()
        if len(trees)!=52:
            raise ValueError("Expected exact frozen 52 S1 historical placements")
        first=None
        for name,tree in trees:
            tipset=check_placement_tree(tree)
            if first is None:first=tipset
            elif tipset!=first:raise ValueError("Historical placement tip lists drift")
            results.append({"scenario":name,
                  "joined46":score(whole,labels,tree_covariance_for(tree,whole["names"])),
                  "sameNA42":score(na,nlabs,tree_covariance_for(tree,na["names"]))})
    summary={
       "status":"REAL_SOURCE_BLOCKED_FIVE_S1_COVARIATES_EMPIRICAL_ASSOCIATION",
       "S1_original_image_angle_native_taxa":244,
       "S1_placement_tree_original_cohort":194,
       "independent_binary_states_matched_to_S1_and_tree":46,
       "matched_NA":42,"matched_EA":4,"nodding_NA":3,"nodding_EA":3,
       "S1_five_fixed_exposures":list(EXPOSURES),
       "primary_precipitation":"chelsa_bio12",
       "secondary_precipitation":"chelsa_bio15",
       "NA_EA_source_blocked_null_maps":nmaps,
       "NA_only_null_maps":n_na_maps,
       "whole46_equal_species":base46,
       "NA42_equal_species":base42,
       "NA42_independent_published_nuclear":nuclear42,
       "PGLS_scenario_count":len(results),
       "source_confounded_region_response":"East Asia 3 D among 4 taxa, NA 3 D among 42; no independent regional replication or interaction power",
       "S1_frozen_model_refitted":False,
       "independent_evolutionary_origins_or_dates_estimated":False,
       "climate_dependent_transition_rates_estimated":False,
       "null_type":"sharp no-climate-state-association conditioned on original source counts, NOT single-axis partial null",
       "not_prospective_independent_EAzami_validation":True}
    if results:
        for key in ("joined46","sameNA42"):
            for outcome in ("estimated_linear_probability_coefficient_per_1SD",
                            "two_axis_maxT_source_stratified_exact_p"):
                for axis,i in (("BIO12",0),("BIO15",1)):
                    values=[r[key][outcome][i] for r in results]
                    summary[f"{key}_{axis}_{outcome}_scenario_range"]=[
                         float(min(values)),float(max(values))]
    doc={"summary":summary,"scenario_details":results,
         "authenticated_joint_source_run":39584744194,
         "joint_artifact_SHA":SOURCE_ARTIFACT_SHA,
         "NA_nuclear_source_DOI":"10.1086/724310",
         "S1_native_observation_input_SHA":receipt["source_provenance"]["native_sha256"]}
    path=Path(output);path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(doc,indent=2,sort_keys=True,allow_nan=False)+"\n",encoding="utf-8")
    return summary


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--joint",required=True);p.add_argument("--out",required=True)
    p.add_argument("--no-placement-scenarios",action="store_true")
    args=p.parse_args()
    summary=run(args.joint,args.out,include_scenarios=not args.no_placement_scenarios)
    print("V4_T2_SOURCE_BLOCKED_FIVE_CLIMATE_SUMMARY="+
          json.dumps(summary,sort_keys=True,ensure_ascii=False),flush=True)


if __name__=="__main__":
    main()
