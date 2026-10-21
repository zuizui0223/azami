"""T2 regional 5-climate contrast on the authenticated S1 joint cohort.

Independent botanical U/D states are kept separate from photo-vertical angles.
This is a *present-state* diagnostic, NOT a climate-dependent Mk transition
model. The 52 GBOTB placements are alternative grafts, not 52 replications.
Source-blocked exact assignment ranks are DESCRIPTIVE, not calibrated
phylogenetic p values (shuffling labels does not preserve tree covariance).
"""
from __future__ import annotations

import argparse
from copy import deepcopy
import itertools
import json
from math import comb
from pathlib import Path
import urllib.request

from Bio import Phylo
import numpy as np
from scipy.linalg import cho_factor, cho_solve

from analysis.v4.joint_s1_t2_phylogeny import EXPOSURES
from analysis.v4.native_phylogeny import tree_covariance
from analysis.v4.t2_pgls_placement_screen import (
    AUTHOR_NUCLEAR_TREE_GIT_BLOB, AUTHOR_NUCLEAR_TREE_PATH,
    SOURCE_COMMIT, SOURCE_REPO, check_placement_tree, git_blob_sha,
)

RECEIPT_PATH = Path(
    "analysis/v4/results/20261021_unified_abcde/"
    "S1_T2_JOINT_244_194_VERIFIED_20261021.json")
RECEIPT_DIGEST = "sha256:aa5a9ac298b8e36c9e96da6f1eb96319733136e2e3a105305092f7da2ca97bff"
EA_SOURCE = "EAzami_FDT4_existing_authority"
NA_SOURCE = "Siniscalchi2023_author_scored"
TEST_AXES = ("chelsa_bio12", "chelsa_bio15")
NUCLEAR_DOI = "10.1086/724310"
NUGGET = 0.05


def validate_joint(source: dict, receipt: dict) -> list[dict]:
    scope = source.get("model_scope", {})
    if (source.get("status") != "empirical_full_244_native_join_source_sha_verified"
        or tuple(source.get("fixed_covariates", ())) != EXPOSURES
        or (scope.get("orientation_eligible_accepted_keys"),
            scope.get("orientation_eligible_observations"),
            scope.get("placement_tree_covered_keys"),
            scope.get("placement_tree_covered_observations"),
            scope.get("grafted_placement_scenarios"))
           != (244, 170066, 194, 169912, 52)):
        raise ValueError("Frozen S1 five-climate 244/194 source cohort changed")
    evidence = receipt.get("results", {}).get("botanical_state_inventory", {})
    if (receipt.get("workflow_run_id") != 39584744194
        or receipt.get("workflow_artifact_digest") != RECEIPT_DIGEST
        or evidence.get("NA_exact_S1_194_tree_matches") != 42
        or evidence.get("EastAsia_exact_S1_194_tree_matches") != 4):
        raise ValueError("Authenticated archived joint-source receipt mismatch")
    provenance = source.get("provenance", {})
    locked = receipt["source_provenance"]
    if (provenance.get("S1_model_source_sha256") != locked["native_sha256"]
        or provenance.get("EAzami_source_git_blob") != locked["EAzami_fdt4_git_blob"]
        or provenance.get("NA_source_git_blob") != locked["NA_head_position_git_blob"]
        or provenance.get("native_decoded_source_authenticated") is not True):
        raise ValueError("Original independent source identities changed")
    rows = [r for r in source["T2_independent_state_inventory"]["taxon_source_audit"]
            if r.get("s1_phylogeny_194_key_match") is True]
    tip_names = [r["source_taxon"].replace(" ", "_") for r in rows]
    concept_ids = [str(r.get("accepted_key")) for r in rows]
    if (len(rows) != 46 or len(set(tip_names)) != 46
        or len(set(concept_ids)) != 46
        or sum(r["source"] == NA_SOURCE for r in rows) != 42
        or sum(r["source"] == EA_SOURCE for r in rows) != 4
        or sum(r["state"] == "nodding" for r in rows) != 6
        or sum(r["state"] == "nodding" and r["source"] == EA_SOURCE for r in rows) != 3):
        raise ValueError("Biological score coverage or taxonomy drift")
    for r in rows:
        if (r["join_status"] != "exact_S1_and_tree"
            or r["state"] not in {"erect", "nodding"}
            or int(r.get("s1_observations", 0)) < 1):
            raise ValueError("Unreviewed morphology/taxon pair")
        climates = r.get("s1_five_climate_centroids")
        if (not isinstance(climates, dict)
            or set(climates) != set(EXPOSURES)
            or not all(np.isfinite(float(climates[e])) for e in EXPOSURES)):
            raise ValueError("Incomplete fixed source five-climate centroid")
    return rows


def cohort(rows):
    """One unit per exact taxon; covariates z-standardized in THIS cohort."""
    rows = sorted(rows, key=lambda x: x["source_taxon"])
    tips = [r["source_taxon"].replace(" ", "_") for r in rows]
    states = np.array([int(r["state"] == "nodding") for r in rows], dtype=np.int8)
    groups = np.array([int(r["source"] == EA_SOURCE) for r in rows], dtype=np.int8)
    values = np.array([[float(r["s1_five_climate_centroids"][v]) for v in EXPOSURES]
                       for r in rows])
    if len(set(tips)) != len(tips) or not np.all(np.isfinite(values)):
        raise ValueError("Repeated tip or invalid climate centroid")
    scales = values.std(axis=0)
    if np.any(scales <= 1e-10):
        raise ValueError("No variation in an original climatic predictor")
    climate = (values-values.mean(axis=0)) / scales
    has_two_sources = len(set(groups)) == 2
    columns = [np.ones(len(rows))]
    if has_two_sources:
        columns.append(groups.astype(float))
    climate_start = len(columns)
    columns.extend([climate[:,j] for j in range(len(EXPOSURES))])
    X = np.column_stack(columns)
    if len(rows) <= X.shape[1] + 2 or np.linalg.matrix_rank(X) != X.shape[1]:
        raise ValueError("Nonidentifiable five-climate plus source model")
    if np.linalg.cond(X) > 1e8:
        raise ValueError("Ill-conditioned five-climate source model")
    src = {}
    for src_name,code in ((NA_SOURCE,0),(EA_SOURCE,1)):
        if np.any(groups == code):
            membership = groups == code
            src[src_name] = {"taxa":int(membership.sum()),
                             "nodding":int(states[membership].sum()),
                             "erect":int((1-states[membership]).sum()),
                             "mean_climate_by_state":{
                                  st:{v:float(values[membership & (states==val),j].mean())
                                      if np.any(membership & (states==val)) else None
                                      for j,v in enumerate(EXPOSURES)}
                                  for st,val in (("erect",0),("nodding",1))}}
    return {"tips":tips,"y":states.astype(float),"groups":groups,"X":X,
            "climate_start":climate_start,"n":len(tips),
            "n_nodding":int(states.sum()),"source_counts":src,
            "cohort_means":values.mean(axis=0).tolist(),
            "cohort_sds":scales.tolist()}


def source_blocked_state_maps(c, max_maps=100000):
    """Exact cross-product, retaining original nodding count in EACH source."""
    y = c["y"].astype(np.int8)
    groups = c["groups"]
    positions = [np.flatnonzero(groups == v) for v in sorted(set(groups))]
    k = [int(y[p].sum()) for p in positions]
    total = int(np.prod([comb(len(p), n) for p,n in zip(positions,k)]))
    if total < 2 or total > max_maps:
        raise ValueError("Degenerate or too-large source-stratified assignment universe")
    labels = np.zeros((total,len(y)),dtype=float)
    choices = [itertools.combinations(p, n) for p,n in zip(positions,k)]
    # combinations iterators must be materialized for Cartesian product reuse
    choices = [list(cset) for cset in choices]
    for i,parts in enumerate(itertools.product(*choices)):
        for part in parts:
            labels[i,list(part)] = 1
    observed = int(np.sum(np.all(labels == y,axis=1)))
    if observed != 1:
        raise ValueError("Observed label assignment not unique in exact universe")
    return labels, {"maps":total,"n_groups":len(positions),
                     "group_sizes":[int(len(p)) for p in positions],
                     "group_nodding":[int(n) for n in k],
                     "observed_assignment_in_universe":observed}


def weights(c, covariance=None, nugget=NUGGET):
    n = c["n"]
    X = c["X"]
    if covariance is None:
        VI = np.eye(n)
        method = "identity_covariance"
    else:
        a = np.asarray(covariance,dtype=float)
        if (a.shape != (n,n) or not np.all(np.isfinite(a))
            or nugget <= 0):
            raise ValueError("Invalid frozen PGLS/author nuclear covariance")
        a = 0.5*(a+a.T)
        # For all placements, branch lengths are only relative; the nugget is
        # explicitly a regularization sensitivity, not biological measurement.
        inv = cho_factor(a+nugget*np.eye(n))
        VI = cho_solve(inv,np.eye(n))
        method = "Brownian_shared_paths_fixed_relative_nugget_0.05"
    information = X.T@VI@X
    ch = cho_factor(information)
    W = cho_solve(ch, X.T@VI)
    idx = [c["climate_start"]+list(EXPOSURES).index(v) for v in TEST_AXES]
    return W[idx,:], method


def exact_scores(c, labels, covariance=None):
    """Combinatorial label-rank DIAGNOSTIC, not valid phylogenetic p-values."""
    w,method = weights(c,covariance)
    obs = w@c["y"]
    sims = (w@labels.T).T
    center = sims.mean(axis=0)
    sd = sims.std(axis=0)
    if np.any(sd <= 1e-10):
        raise ValueError("Degenerate conditional reference contrasts")
    target = np.abs((obs-center)/sd)
    reference = np.abs((sims-center)/sd)
    ranks = (reference >= target).mean(axis=0)
    two_axis = (reference.max(axis=1)[:,None] >= target).mean(axis=0)
    return {"status":"exploratory_source_stratified_combinatorial_ranks",
            "method":method,"n_taxa":c["n"],"n_nodding":c["n_nodding"],
            "same_five_predictors":list(EXPOSURES),
            "test_axes":list(TEST_AXES),
            "source_fixed_effect":len(set(c["groups"]))==2,
            "coefficient_is_linear_probability_change_per_cohort_SD":obs.tolist(),
            "reference_null_mean":center.tolist(),
            "reference_null_sd":sd.tolist(),
            "combinatorial_reference_rank":ranks.tolist(),
            "two_axis_maxT_reference_rank":two_axis.tolist(),
            "exact_map_count":int(len(labels)),
            "p_value_status":"NOT_CALIBRATED_FOR_PHYLOGENETIC_DEPENDENCE",
            "interpretation":"Source-blocked fixed-prevalence assignment diagnostic; NOT inferential P for climate-specific, tree-dependent or causal effects"}


def pruned_covariance(tree,names):
    phy = deepcopy(tree)
    if not set(names).issubset(check_placement_tree(phy)):
        raise ValueError("Source tip missing from phylogeny")
    for tip in list(phy.get_terminals()):
        if tip.name not in names:
            phy.prune(tip)
    return tree_covariance(phy,names)


def published_nuclear_tree():
    url=(f"https://raw.githubusercontent.com/{SOURCE_REPO}/{SOURCE_COMMIT}/"
         f"{AUTHOR_NUCLEAR_TREE_PATH}")
    with urllib.request.urlopen(url,timeout=90) as response:
        raw = response.read()
    if git_blob_sha(raw) != AUTHOR_NUCLEAR_TREE_GIT_BLOB:
        raise ValueError("Independently published nuclear tree Git blob drift")
    import io
    tree = Phylo.read(io.StringIO(raw.decode("utf-8-sig")),"newick")
    if len(tree.get_terminals()) != 63:
        raise ValueError("Independently published 63-tip nuclear tree changed")
    return tree, __import__("hashlib").sha256(raw).hexdigest()


def analyze(joint_path,output_path,*,placement=True,author_tree=True):
    evidence = json.loads(Path(joint_path).read_text(encoding="utf-8"))
    receipt = json.loads(RECEIPT_PATH.read_text(encoding="utf-8"))
    rows = validate_joint(evidence,receipt)
    full = cohort(rows)
    north = cohort([r for r in rows if r["source"]==NA_SOURCE])
    assert (full["n"],full["n_nodding"],north["n"],north["n_nodding"]) == (46,6,42,3)
    full_maps,full_counts = source_blocked_state_maps(full)
    north_maps,north_counts = source_blocked_state_maps(north)
    full_simple = exact_scores(full,full_maps)
    north_simple = exact_scores(north,north_maps)
    nuclear = None
    nuclear_source_sha = None
    if author_tree:
        phy, nuclear_source_sha = published_nuclear_tree()
        nuclear = exact_scores(north,north_maps,pruned_covariance(phy,north["tips"]))
    placements = []
    if placement:
        from analysis.v4.native_phylogeny import restore_trees
        trees = restore_trees()
        if len(trees) != 52:
            raise ValueError("Exactly 52 frozen phylogenetic placements required")
        for scenario,tree in trees:
            placements.append({
                "scenario":scenario,
                "all_46":exact_scores(full,full_maps,pruned_covariance(tree,full["tips"])),
                "same_NA_42":exact_scores(north,north_maps,pruned_covariance(tree,north["tips"])),
            })
    def scenario_range(group, field, col):
        vals = [row[group][field][col] for row in placements]
        return [float(min(vals)),float(max(vals))]
    summary = {
        "status":"ACTUAL_S1_FIVE_EXPOSURES_AND_BOTANICAL_STATE_REGIONAL_SENSITIVITY",
        "source_verified_S1_run_id":39584744194,
        "joint_reference_digest":RECEIPT_DIGEST,
        "total_exact_joint_taxa":46,
        "total_botanical_nodding":6,
        "north_america_taxa":42,
        "north_america_nodding":3,
        "east_asia_taxa":4,
        "east_asia_nodding":3,
        "fixed_exposures":list(EXPOSURES),
        "primary_axis":"chelsa_bio12",
        "secondary_axis":"chelsa_bio15",
        "source_stratified_maps_46":full_counts["maps"],
        "source_stratified_maps_NA42":north_counts["maps"],
        "scenarios_completed":len(placements),
        "source_imbalance_warning":"East Asian 3/4 nodding vs North American 3/42; regional climate interactions not identified",
        "rank_interpretation":"Descriptive exact assignments conditioned on source lane, NOT phylogeny-calibrated p-values",
        "nugget_fixed_relative_to_tree_covariance":NUGGET,
        "S1_original_models_unchanged":True,
        "image_angle_as_true_gravitational_state":False,
        "real_independent_nodding_origins_inferred":False,
        "environment_dependent_evolutionary_transition_rates_fitted":False,
        "unstandardized_coefficients_directly_comparable_to_S1":False,
        "full46_nonphylogenetic":full_simple,
        "NA42_nonphylogenetic":north_simple,
        "NA42_independent_nuclear_covariance":nuclear,
        "NA63_nuclear_tree_sha256":nuclear_source_sha,
    }
    if placements:
        for group in ("all_46","same_NA_42"):
            for col,axis in enumerate(("BIO12","BIO15")):
                summary[f"{group}_{axis}_GLS_coefficient_scenario_range"] = scenario_range(
                    group,"coefficient_is_linear_probability_change_per_cohort_SD",col)
                summary[f"{group}_{axis}_two_axis_reference_rank_scenario_range"] = scenario_range(
                    group,"two_axis_maxT_reference_rank",col)
    result = {"summary":summary,
              "source_blocking_46":full_counts,
              "source_blocking_NA42":north_counts,
              "source_descriptive_climate_by_state":full["source_counts"],
              "tree_placement_scenarios":placements,
              "original_joint_workflow_artifact":receipt["workflow_artifact_id"],
              "original_joint_workflow_digest":RECEIPT_DIGEST}
    path=Path(output_path);path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(result,indent=2,sort_keys=True,allow_nan=False)+"\n",
                    encoding="utf-8")
    return result


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--joint",required=True)
    ap.add_argument("--output",required=True)
    ap.add_argument("--no-52",action="store_true")
    args=ap.parse_args()
    report=analyze(args.joint,args.output,placement=not args.no_52)
    print("V4_T2_REGIONAL_CLIMATE_SUMMARY="+json.dumps(
        report["summary"],ensure_ascii=False,sort_keys=True),flush=True)


if __name__=="__main__":
    main()
