"""T2 exploratory re-use of the exact 52 GBOTB/V.PhyloMaker2 PGLS trees.

NOT a dated Cirsium molecular phylogeny, and NOT T2's independent primary
historical test. This screen keeps published author-scored head posture separate
from image-vertical angles and from the prospective double-blind scoring layer.
All 52 branch-length trees are checksum-restored from the existing V3 archive.
"""
from __future__ import annotations

import argparse
from collections import Counter
from copy import deepcopy
import csv
import hashlib
import io
import json
from pathlib import Path
import urllib.request

import numpy as np
from Bio import Phylo

from analysis.v4.rain_history import rate_fit, stochastic_maps

SOURCE_COMMIT = "c925d4af86160064448f794a9715f70b917a2786"
SOURCE_REPO = "carol-siniscalchi/cirsium"
SOURCE_PATH = "traits/traits.csv"
SOURCE_GIT_BLOB = "fa127f0587b676e6e004ea121424d26069586bcb"
AUTHOR_NUCLEAR_TREE_PATH = "traits/cirsium.phyluce95p.ultrametric.tre"
AUTHOR_NUCLEAR_TREE_GIT_BLOB = "6417a38793427bc7c8562dc0e2cab59b74d8c316"
SOURCE_DOI = "10.1086/724310"
EXPECTED_AUTHOR_TAXA = 63
EXPECTED_NODDING = 5

def git_blob_sha(raw: bytes) -> str:
    return hashlib.sha1(b"blob " + str(len(raw)).encode() + b"\0" + raw).hexdigest()


def read_author_states(raw: bytes) -> dict[str, int]:
    if git_blob_sha(raw) != SOURCE_GIT_BLOB:
        raise ValueError("Author's Git blob does not match pinned publication evidence")
    rows = list(csv.DictReader(io.StringIO(raw.decode("utf-8-sig"))))
    if len(rows) != EXPECTED_AUTHOR_TAXA:
        raise ValueError("Expected 63 original author-scored taxa")
    out = {}
    for row in rows:
        name, state = row.get("Species", ""), row.get("HeadPosition", "")
        if not name.startswith("Cirsium_") or name in out or state not in ("0", "1"):
            raise ValueError("Duplicate taxon, changed label or invalid published head state")
        out[name] = int(state)
    if sum(out.values()) != EXPECTED_NODDING:
        raise ValueError("Unexpected published nodding-state frequency")
    return out


def download_author_states() -> tuple[dict[str, int], str]:
    url = f"https://raw.githubusercontent.com/{SOURCE_REPO}/{SOURCE_COMMIT}/{SOURCE_PATH}"
    with urllib.request.urlopen(url, timeout=90) as handle:
        raw = handle.read()
    return read_author_states(raw), hashlib.sha256(raw).hexdigest()


def read_author_nuclear_tree(raw: bytes, expected_states: dict[str,int]):
    if git_blob_sha(raw) != AUTHOR_NUCLEAR_TREE_GIT_BLOB:
        raise ValueError("Published independent nuclear tree's pinned Git blob mismatch")
    tree=Phylo.read(io.StringIO(raw.decode("utf-8-sig")), "newick")
    tips=check_placement_tree(tree)
    if tips != set(expected_states):
        raise ValueError("Published nuclear tree tip set differs from original scored taxa")
    heights=list(tree.depths()[t] for t in tree.get_terminals())
    if max(heights)-min(heights)>0.001:
        raise ValueError("Published reference tree is not effectively ultrametric")
    return tree


def nuclear_reference_summary(tree, states, common, draws=150, seed=20261021):
    """Sequence-tree comparison for identical 44 (or N) taxa, not a global tree."""
    matched=set(common)
    if not matched <= set(states):
        raise ValueError("Nuclear comparison contains unknown morphological states")
    if matched != set(t.name for t in tree.get_terminals()):
        phy=deepcopy(tree)
        for tip in list(phy.get_terminals()):
            if tip.name not in matched:
                phy.prune(tip)
    else:
        phy=deepcopy(tree)
    matched_states={t:states[t] for t in matched}
    n_down=sum(matched_states.values())
    result={"n_taxa":len(matched),"n_nodding":n_down,"n_erect":len(matched)-n_down,
            "published_nuclear_tree_based":True,"cross_tree_same_tip_cohort":True,
            "not_independent_morphology_reannotation":True,
            "chronogram_absolute_dates_inferred":False}
    result.update(parsimony(phy,matched_states))
    if n_down>=2 and (len(matched)-n_down)>=2:
        q,lnlike,boundary=rate_fit(phy,matched_states)
        result.update(fitted_symmetric_Mk_rate=q, Mk_log_likelihood=lnlike,
                      optimizer_near_boundary=boundary,
                      conditional_mapping=stochastic_maps(phy,matched_states,q,
                                                          draws=draws,seed=seed))
    return result


def parsimony(tree, observed) -> dict[str, int]:
    """Minimal total changes with unknown / fixed root, for each rooted scenario."""
    costs = {}
    for node in tree.find_clades(order="postorder"):
        if node.is_terminal():
            state = observed[node.name]
            costs[node] = np.array([0, 10**8], dtype=np.int64) if state == 0 else np.array([10**8, 0], dtype=np.int64)
        else:
            children = [costs[c] for c in node.clades]
            costs[node] = np.array([sum(min(child[0] + (parent != 0), child[1] + (parent != 1))
                                        for child in children) for parent in (0, 1)], dtype=np.int64)
    root = costs[tree.root]
    return {"root_erect_min_steps": int(root[0]),
            "root_nodding_min_steps": int(root[1]),
            "unknown_root_min_steps": int(min(root))}


def check_placement_tree(tree):
    tips = [tip.name for tip in tree.get_terminals()]
    if any(not t for t in tips) or len(tips) != len(set(tips)):
        raise ValueError("PGLS tree tip identities not unique")
    if len(tree.root.clades) < 2:
        raise ValueError("PGLS tree has no auditable branching root")
    for node in tree.find_clades():
        if node is tree.root:
            continue
        length = node.branch_length
        if length is None or not np.isfinite(length) or length < 0:
            raise ValueError("PGLS tree missing / non-finite / negative branch lengths")
    return set(tips)


def evaluate_placements(trees, states, source_sha, output, *, draws=150, min_scored=12,
                        min_nodding=2, seed=20261021,
                        published_nuclear_tree=None, nuclear_tree_sha=None):
    if len(trees) != 52:
        raise ValueError("Expected the exact 52 archived PGLS placement scenarios")
    if len(states) != 63 or sum(states.values()) != 5:
        raise ValueError("Published author scores not source validated")
    results = []
    first_tips = None
    for number, (scenario, original) in enumerate(trees):
        tips = check_placement_tree(original)
        if first_tips is None:
            first_tips = tips
        elif first_tips != tips:
            raise ValueError("Historical scenario tip sets changed")
        matches = tips & states.keys()
        subset = {tip: states[tip] for tip in matches}
        counts = Counter(subset.values())
        record = {
            "scenario": scenario, "pgls_tree_tips": len(tips),
            "author_taxa_matching_pgls_exact": len(matches),
            "matched_erect": counts.get(0, 0), "matched_nodding": counts.get(1, 0),
            "excluded_published_author_taxa": sorted(states.keys() - matches),
            "source_status": "independent_published_taxon_morphology_NOT_dual_blind_scoring",
            "tree_status": "GBOTB_dated_backbone_with_VPhyloMaker2_genus_grafts_not_independent_nuclear_Cirsium_tree",
        }
        if len(matches) < 2:
            record.update(status="blocked_less_than_two_taxa")
            results.append(record)
            continue
        tree = deepcopy(original)
        for tip in list(tree.get_terminals()):
            if tip.name not in matches:
                tree.prune(tip)
        record.update(parsimony(tree, subset))
        if (len(matches) < min_scored or counts.get(1, 0) < min_nodding
                or counts.get(0, 0) < 2):
            record.update(status="descriptive_parsimony_only_insufficient_binary_coverage")
        else:
            q, lnlike, boundary = rate_fit(tree, subset)
            record.update(
                status="conditional_placement_Mk_sensitivity_NOT_independent_history",
                symmetric_Mk_rate_per_backbone_branch_unit=q,
                log_likelihood=lnlike,
                optimizer_near_boundary=boundary,
                mapping=stochastic_maps(tree, subset, q, draws=draws, seed=seed + number),
            )
        results.append(record)
    statuses = Counter(r["status"] for r in results)
    summary = {
        "stage":"T2_PGLS_52_TREE_EXPLORATORY",
        "scientific_status":"placement_sensitivity_only_NOT_independent_origins_or_dates_or_rain_causation",
        "author_evidence":{"doi":SOURCE_DOI, "repo":SOURCE_REPO, "commit":SOURCE_COMMIT,
            "path":SOURCE_PATH, "git_blob":SOURCE_GIT_BLOB,
            "source_sha256":source_sha, "n_taxa":63, "n_nodding":5,
            "annotation_basis":"original_publication_taxon_scored_head_position_not_camera_vertical"},
        "pgls_tree_evidence":{"doi":"GBOTB.extended.LCVP/V.PhyloMaker2 archived in azami Zenodo record 22703215",
            "backbone_used_for_placement":True, "n_scenarios":52,
            "not_independent_phylogenetic_replicates":True,
            "root_and_branch_scaling_not_revalidated_for_Cirsium_events":True},
        "n_scenario_results":len(results),
        "status_counts":dict(sorted(statuses.items())),
        "matched_taxa_n_range":[min(r["author_taxa_matching_pgls_exact"] for r in results),
                                max(r["author_taxa_matching_pgls_exact"] for r in results)],
        "matched_nodding_n_range":[min(r["matched_nodding"] for r in results),
                                   max(r["matched_nodding"] for r in results)],
        "not_independent_repeated_selection":True,
        "rain_temperature_evolutionary_test_executed":False,
        "empirical_transition_dates_estimated":False,
        "ancestral_state_independent_nuclear_confirmation":False,
        "source_photo_orientation_as_true_biological_posture":False,
    }
    eligible=[r for r in results if r["status"].startswith("conditional_placement")]
    if eligible:
        for key in ("unknown_root_min_steps", "root_erect_min_steps", "root_nodding_min_steps"):
            values=[r[key] for r in eligible]
            summary[key+"_scenario_range"]=[min(values),max(values)]
        for key in ("erect_to_nodding_q025_median_q975", "nodding_to_erect_q025_median_q975"):
            values=[r["mapping"][key][1] for r in eligible]
            summary[key+"_scenario_median_range"]=[float(min(values)),float(max(values))]
    elif results[0].get("unknown_root_min_steps") is not None:
        summary["descriptive_only_unknown_root_min_steps_scenario_range"]=[
            min(r["unknown_root_min_steps"] for r in results if "unknown_root_min_steps" in r),
            max(r["unknown_root_min_steps"] for r in results if "unknown_root_min_steps" in r)]
    if published_nuclear_tree is not None:
        if first_tips is None:
            raise ValueError("No PGLS taxa for nuclear comparison")
        if not nuclear_tree_sha:
            raise ValueError("Nuclear tree's original SHA256 must be stored")
        shared=set(states) & first_tips
        # Same original 63-tip publication tree is separately analysed with all
        # published tips, but ONLY shared tips compare against PGLS 52 scenarios.
        summary["independent_nuclear_reference"]={
            "source_repo":SOURCE_REPO,"commit":SOURCE_COMMIT,
            "path":AUTHOR_NUCLEAR_TREE_PATH,
            "git_blob":AUTHOR_NUCLEAR_TREE_GIT_BLOB,
            "source_sha256":nuclear_tree_sha,
            "source_doi":SOURCE_DOI,
            "common_tip_comparison":nuclear_reference_summary(
                published_nuclear_tree,states,shared,draws=draws,seed=seed+52),
            "complete_published_63_tip_context":nuclear_reference_summary(
                published_nuclear_tree,states,set(states),draws=draws,seed=seed+53),
            "not_global_cirsium_tree":True,
            "not_a_confirmed_dated_origin_analysis":True,
        }
    report = {"summary":summary, "placement_scenarios":results}
    target=Path(output)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(report,ensure_ascii=False,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", default="work/t2-pgls-placement-screen.json")
    parser.add_argument("--draws", type=int, default=150)
    parser.add_argument("--min-scored", type=int, default=12)
    parser.add_argument("--min-nodding", type=int, default=2)
    opts = parser.parse_args()
    if opts.draws < 1 or opts.draws > 10000:
        raise ValueError("Draw count out of range")
    from analysis.v4.native_phylogeny import restore_trees
    states, source_sha=download_author_states()
    trees=restore_trees()
    tree_url=(f"https://raw.githubusercontent.com/{SOURCE_REPO}/{SOURCE_COMMIT}/"
              f"{AUTHOR_NUCLEAR_TREE_PATH}")
    with urllib.request.urlopen(tree_url,timeout=90) as handle:
        tree_raw=handle.read()
    published_tree=read_author_nuclear_tree(tree_raw,states)
    results=evaluate_placements(trees,states,source_sha,opts.out,
                                published_nuclear_tree=published_tree,
                                nuclear_tree_sha=hashlib.sha256(tree_raw).hexdigest(),
                                draws=opts.draws,min_scored=opts.min_scored,
                                min_nodding=opts.min_nodding)
    print("T2_PGLS_PLACEMENT_SCREEN_AGGREGATE="+json.dumps(results["summary"],
        ensure_ascii=False,sort_keys=True),flush=True)


if __name__ == "__main__":
    main()
