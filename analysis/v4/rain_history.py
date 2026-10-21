"""T2/T3 rain-orientation history: audited evidence intake and conditional Mk pilot.

This module does NOT convert photo-vertical angles into gravity-relative posture.
No branch timing, rain-driven selection, or climate-causal inference is implemented.
Use only independently reviewed botanical states and rooted nuclear species trees.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
from copy import deepcopy
import csv
import hashlib
import io
import json
from pathlib import Path

from Bio import Phylo
import numpy as np
from scipy.optimize import minimize_scalar

BINARY = {"erect": 0, "nodding": 1}
STATES = set(BINARY) | {"variable", "unknown"}
SCORE_COLS = {
    "record_id", "accepted_taxon", "concept_id", "source_taxon",
    "source_id", "source_locator", "evidence_quote", "head_stage",
    "annotator_id", "annotation_role", "blind_climate", "blind_phylogeny",
    "state",
}
MAP_COLS = {"tip", "accepted_taxon", "concept_id", "decision", "tribe",
            "phenotype_reference", "head_stage", "homology_reviewed"}
CLIMATE_COLS = {"concept_id", "chelsa_bio12", "chelsa_bio01",
                "elevation_m", "region", "sampling_frame"}
FORBIDDEN_SCORE_FRAGMENTS = ("bio01", "bio12", "bio15", "climate", "rainfall", "phylo_score")


class EvidenceGateError(ValueError):
    """Missing or contradictory independent evidence: do not infer history."""


def load_csv(path, required):
    with Path(path).open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None or not required <= set(reader.fieldnames):
            raise EvidenceGateError(f"{path}: missing required columns {sorted(required - set(reader.fieldnames or []))}")
        rows = list(reader)
    return rows


def score_botanical_records(rows):
    """Resolve two climate/tree-blinded primary scores plus optional third review.

    One record = one botanical source passage. Conflicting *sources* are never
    silently majority-voted into a single species state.
    """
    if not rows:
        raise EvidenceGateError("No botanical evidence records; templates are not data")
    if any(any(s in k.lower() for s in FORBIDDEN_SCORE_FRAGMENTS)
           for k in (set(rows[0]) - SCORE_COLS)):
        raise EvidenceGateError("Scoring sheets must not contain climate or phylogeny estimates")
    records = defaultdict(list)
    for row in rows:
        if row["state"] not in STATES:
            raise EvidenceGateError(f"Bad state {row['state']!r}")
        if row["annotation_role"] not in {"primary", "adjudicator"}:
            raise EvidenceGateError("Unknown annotation role")
        if row["blind_climate"] != "yes" or row["blind_phylogeny"] != "yes":
            raise EvidenceGateError("All assessors must be blind to climate/tree inference")
        for field in ("record_id", "accepted_taxon", "concept_id", "source_taxon",
                      "source_id", "source_locator", "evidence_quote", "annotator_id"):
            if not row[field].strip():
                raise EvidenceGateError(f"Blank {field} in record {row.get('record_id')}")
        if row["head_stage"] not in {"anthesis", "bud", "postanthesis", "uncertain"}:
            raise EvidenceGateError("Developmental head stage must be explicit")
        records[row["record_id"]].append(row)
    consensus = defaultdict(list)
    record_audit = []
    for rid, entries in sorted(records.items()):
        shared = ("accepted_taxon", "concept_id", "source_taxon", "source_id",
                  "source_locator", "evidence_quote", "head_stage")
        if any(len({r[f] for r in entries}) != 1 for f in shared):
            raise EvidenceGateError(f"Source identity mismatch across scorers: {rid}")
        primary = [r for r in entries if r["annotation_role"] == "primary"]
        third = [r for r in entries if r["annotation_role"] == "adjudicator"]
        reviewers = [r["annotator_id"] for r in entries]
        if len(primary) != 2 or len(third) > 1 or len(set(reviewers)) != len(reviewers):
            raise EvidenceGateError(f"Need two different primary assessors and at most one independent adjudicator: {rid}")
        if primary[0]["state"] == primary[1]["state"]:
            resolved = primary[0]["state"]
            method = "double_agreement"
        elif third:
            resolved = third[0]["state"]
            method = "third_adjudication"
        else:
            resolved = None
            method = "awaiting_adjudication"
        if entries[0]["head_stage"] != "anthesis":
            resolved = None
            method = "excluded_nonanthesis"
        key = (entries[0]["concept_id"], entries[0]["accepted_taxon"])
        consensus[key].append(resolved)
        record_audit.append({
            "record_id": rid, "concept_id": key[0], "accepted_taxon": key[1],
            "status": method, "state": resolved,
        })
    decisions = {}
    for (concept, name), votes in consensus.items():
        if concept in decisions:
            raise EvidenceGateError(f"Concept ID linked to multiple taxon names: {concept}")
        informative = [v for v in votes if v not in (None, "unknown")]
        if len(set(informative)) > 1:
            state, status = None, "conflicting_botanical_sources"
        elif not informative:
            state, status = None, "unresolved_or_unknown"
        else:
            state = informative[0]
            status = "accepted_binary" if state in BINARY else "variable_nonbinary"
        decisions[concept] = {"accepted_taxon": name, "state": state, "status": status,
                              "source_records": len(votes)}
    return decisions, record_audit


def validate_tree(tree_path, manifest_path):
    """Require declared independently validated phylogeny; a grafted sensitivity is not one."""
    raw = Path(tree_path).read_bytes()
    meta = json.loads(Path(manifest_path).read_text(encoding="utf-8"))
    required = {"tree_sha256", "source_doi", "tree_type", "root_method",
                "branch_unit", "review_status"}
    if not required <= meta.keys():
        raise EvidenceGateError(f"Tree manifest missing: {sorted(required - meta.keys())}")
    if hashlib.sha256(raw).hexdigest() != meta["tree_sha256"]:
        raise EvidenceGateError("Tree bytes fail SHA256 verification")
    if (meta["tree_type"] != "independent_nuclear_species_tree"
            or meta["review_status"] != "independently_reviewed"
            or not str(meta["source_doi"]).startswith("10.")
            or len(str(meta["root_method"]).strip()) < 8):
        raise EvidenceGateError("Only independently reviewed rooted nuclear species trees qualify")
    if meta["branch_unit"] not in {"substitutions_per_site", "million_years"}:
        raise EvidenceGateError("Unknown branch units: placement/dummy unit not accepted")
    if meta["branch_unit"] == "million_years" and not all(
        meta.get(f) for f in ("calibration_doi", "calibration_method", "calibration_uncertainty")
    ):
        raise EvidenceGateError("Dated tree needs source calibration, method and uncertainty")
    try:
        tree = Phylo.read(io.StringIO(raw.decode("utf-8")), "newick")
    except Exception as exc:
        raise EvidenceGateError("Invalid rooted Newick input") from exc
    if len(tree.root.clades) != 2:
        raise EvidenceGateError("Tree root must have two explicitly reviewed descendant clades")
    tips = tree.get_terminals()
    names = [t.name for t in tips]
    if any(not name for name in names) or len(set(names)) != len(names):
        raise EvidenceGateError("Tip names must be present and unique")
    for node in tree.find_clades():
        if node is tree.root:
            continue
        t = node.branch_length
        if t is None or not np.isfinite(t) or t <= 0:
            raise EvidenceGateError("Every branch must have an independently supported positive length")
    return tree, meta


def audited_tips(tree, crosswalk, states, stage):
    """One concept per molecular tip, no grafted key collisions, only verified source states."""
    if stage not in ("T2", "T3"):
        raise EvidenceGateError("stage must be T2 or T3")
    mapped = {}
    concepts = set()
    for row in crosswalk:
        tip, concept = row["tip"], row["concept_id"]
        if not tip or not concept or tip in mapped or concept in concepts:
            raise EvidenceGateError("Tree-tip/concept crosswalk is not one-to-one")
        mapped[tip] = row
        concepts.add(concept)
    tree_names = {t.name for t in tree.get_terminals()}
    if not tree_names <= mapped.keys():
        raise EvidenceGateError("Every tip needs an explicit taxonomic review row")
    observed, taxa, blocked = {}, [], defaultdict(int)
    for tip in sorted(tree_names):
        row = mapped[tip]
        if row["decision"] != "verified_exact":
            blocked["unverified_taxon_mapping"] += 1
            continue
        state_row = states.get(row["concept_id"])
        if state_row is None or state_row["accepted_taxon"] != row["accepted_taxon"]:
            blocked["no_matching_blinded_source_state"] += 1
            continue
        if row["head_stage"] != "anthesis":
            blocked["nonanthesis_crosswalk"] += 1
            continue
        if stage == "T3" and (row["homology_reviewed"] != "yes"
                              or row["phenotype_reference"] != "botanical_head_posture"
                              or not row["tribe"]):
            blocked["unreviewed_cross_tribe_homology"] += 1
            continue
        if state_row["state"] not in BINARY:
            blocked["unknown_variable_or_conflict"] += 1
            continue
        observed[tip] = BINARY[state_row["state"]]
        taxa.append({"tip": tip, "concept_id": row["concept_id"],
                     "tribe": row["tribe"], "state": state_row["state"]})
    return observed, taxa, dict(blocked)


def transition_matrix(rate, length):
    x = np.exp(-2.0 * rate * length)
    return np.array([[(1 + x) / 2, (1 - x) / 2],
                     [(1 - x) / 2, (1 + x) / 2]])


def pruning_likelihood(tree, observed, rate):
    """Scaled two-state equal-rates Mk pruning, avoiding product underflow."""
    partials = {}
    scales = {}
    for node in tree.find_clades(order="postorder"):
        if node.is_terminal():
            state = observed.get(node.name)
            partials[node] = np.ones(2) if state is None else np.eye(2)[state]
            scales[node] = 0.0
            continue
        value = np.ones(2)
        scale = 0.0
        for child in node.clades:
            value *= transition_matrix(rate, child.branch_length) @ partials[child]
            scale += scales[child]
        factor = float(value.max())
        if factor <= 0:
            return -np.inf, {}
        partials[node] = value / factor
        scales[node] = scale + np.log(factor)
    root = tree.root
    log_likelihood = float(np.log(0.5 * partials[root].sum()) + scales[root])
    return log_likelihood, partials


def rate_fit(tree, observed):
    lengths = [c.branch_length for c in tree.find_clades() if c is not tree.root]
    scale = float(np.median(lengths))
    obj = lambda lr: -pruning_likelihood(tree, observed, np.exp(lr) / scale)[0]
    fit = minimize_scalar(obj, bounds=(-12, 8), method="bounded", options={"xatol": 1e-5})
    if not fit.success:
        raise EvidenceGateError("Mk rate optimization failed")
    return float(np.exp(fit.x) / scale), float(-fit.fun), bool(abs(fit.x + 12) < 0.005 or abs(fit.x - 8) < 0.005)


def conditional_poisson_events(mean, odd, rng):
    """Poisson flips on a branch conditional on endpoints' parity."""
    if mean < 1e-9:
        return 1 if odd else 0
    if mean > 25:
        while True:
            n = int(rng.poisson(mean))
            if (n % 2 == 1) == odd:
                return n
    # Stable parity-normalised inversion. No no-change assumption on long branches.
    p = np.exp(-mean) * (mean if odd else 1.0)
    n = 1 if odd else 0
    mass = -np.expm1(-2 * mean) / 2 if odd else (1 + np.exp(-2 * mean)) / 2
    target = rng.random() * mass
    subtotal = 0.0
    while n < 1000:
        subtotal += p
        if subtotal >= target:
            return n
        p *= mean * mean / ((n + 1) * (n + 2))
        n += 2
    raise EvidenceGateError("Failed to sample conditional branch transitions")


def stochastic_maps(tree, observed, rate, draws=1000, seed=144, tribe_by_tip=None):
    """Sample ancestral states and *all* symmetric Mk flips, not just node changes."""
    if draws < 1 or draws > 100000:
        raise EvidenceGateError("Posterior draw count out of supported bounds")
    _, partials = pruning_likelihood(tree, observed, rate)
    rng = np.random.default_rng(seed)
    counts = np.zeros((draws, 2), dtype=int)
    root_nodding = 0
    tribe_clades = {}
    tribe_counts = {}
    if tribe_by_tip is not None:
        for node in tree.find_clades(order="postorder"):
            if node.is_terminal():
                tribe_clades[node] = {tribe_by_tip[node.name]}
            else:
                tribe_clades[node] = set().union(*(tribe_clades[c] for c in node.clades))
        tribe_counts = {tribe: np.zeros(draws, dtype=int)
                        for tribe in set(tribe_by_tip.values())}
    for i in range(draws):
        root_weights = partials[tree.root] / partials[tree.root].sum()
        root_state = int(rng.choice(2, p=root_weights))
        root_nodding += root_state
        stack = [(tree.root, root_state)]
        while stack:
            parent, old = stack.pop()
            for child in parent.clades:
                probs = transition_matrix(rate, child.branch_length)[old] * partials[child]
                probs /= probs.sum()
                new = int(rng.choice(2, p=probs))
                n = conditional_poisson_events(rate * child.branch_length,
                                               old != new, rng)
                gains = (n + (old == 0)) // 2
                counts[i, 0] += gains  # erect -> nodding
                counts[i, 1] += (n + (old == 1)) // 2  # nodding -> erect
                if tribe_counts and len(tribe_clades[child]) == 1:
                    tribe = next(iter(tribe_clades[child]))
                    tribe_counts[tribe][i] += gains
                stack.append((child, new))
    result = {
        "draws": draws, "seed": seed,
        "root_nodding_p": root_nodding / draws,
        "erect_to_nodding_q025_median_q975": np.quantile(counts[:, 0], [0.025, 0.5, 0.975]).tolist(),
        "nodding_to_erect_q025_median_q975": np.quantile(counts[:, 1], [0.025, 0.5, 0.975]).tolist(),
        "p_at_least_two_erect_to_nodding": float(np.mean(counts[:, 0] >= 2)),
        "warning": "Conditional on one reviewed tree, binary-only tip scores, symmetric Mk and rate point estimate; NOT climate-dependent transitions, independent origin times, or fitness selection."
    }
    if tribe_counts:
        result["within_tribe_gain_distribution"] = {
            k: np.quantile(v, [0.025, 0.5, 0.975]).tolist()
            for k, v in sorted(tribe_counts.items())
        }
        result["p_gains_in_two_or_more_tribes_conditional_Mk"] = float(np.mean(
            np.sum(np.stack([v > 0 for v in tribe_counts.values()]), axis=0) >= 2
        ))
        result["cross_tribe_warning"] = (
            "Changes on mixed-tribe stem branches are NOT attributed; two within-tribe "
            "event histories under a symmetric Mk do NOT establish rain-driven convergence."
        )
    return result


def climate_coverage(rows, taxa):
    by_concept = defaultdict(list)
    for row in rows:
        by_concept[row["concept_id"]].append(row)
    complete = []
    for taxon in taxa:
        group = by_concept.get(taxon["concept_id"], [])
        if len(group) != 1:
            continue
        try:
            p = float(group[0]["chelsa_bio12"])
            t = float(group[0]["chelsa_bio01"])
            elev = float(group[0]["elevation_m"])
        except (ValueError, TypeError):
            continue
        if (np.isfinite(p) and np.isfinite(t) and np.isfinite(elev)
                and group[0]["sampling_frame"] and group[0]["region"]):
            complete.append(taxon["concept_id"])
    return {"complete_climate_and_elevation": len(complete),
            "state_scored_tips": len(taxa),
            "complete_concept_ids": sorted(complete)}


def execute(scores_path, crosswalk_path, tree_path, manifest_path, output_path,
            stage="T2", climate_path=None, draws=1000, min_scored=20):
    scores = load_csv(scores_path, SCORE_COLS)
    decisions, audit = score_botanical_records(scores)
    tree, manifest = validate_tree(tree_path, manifest_path)
    crosswalk = load_csv(crosswalk_path, MAP_COLS)
    observed, taxa, blocked = audited_tips(tree, crosswalk, decisions, stage)
    tribes = sorted({row["tribe"] for row in taxa if row["tribe"]})
    result = {
        "stage": stage, "status": "blocked_no_historical_inference",
        "source_sha256": {
            k: hashlib.sha256(Path(v).read_bytes()).hexdigest()
            for k, v in {"scores": scores_path, "crosswalk": crosswalk_path,
                         "tree": tree_path, "manifest": manifest_path}.items()},
        "source_tree_doi": manifest["source_doi"], "branch_unit": manifest["branch_unit"],
        "tree_taxa": len(tree.get_terminals()), "scored_taxa": len(observed),
        "n_erect": sum(s == 0 for s in observed.values()),
        "n_nodding": sum(s == 1 for s in observed.values()),
        "resolved_source_records": sum(x["status"] in ("double_agreement", "third_adjudication") for x in audit),
        "unresolved_source_records": sum(x["status"] not in ("double_agreement", "third_adjudication") for x in audit),
        "mapping_exclusions": blocked, "tribes": tribes,
        "gravitational_photographic_angle_validated": False,
        "rain_vs_temperature_evolutionary_association_fitted": False,
        "transition_ages_estimated": False,
        "convergent_rain_adaptation_established": False,
    }
    if climate_path:
        climate = load_csv(climate_path, CLIMATE_COLS)
        result["source_sha256"]["climate"] = hashlib.sha256(Path(climate_path).read_bytes()).hexdigest()
        result["climate_coverage"] = climate_coverage(climate, taxa)
    result["scored_tip_fraction"] = len(observed) / len(tree.get_terminals())
    result["unscored_tips_excluded_from_Mk_pilot"] = len(tree.get_terminals()) - len(observed)
    if stage == "T3" and len(tribes) < 2:
        result["stop_reason"] = "T3 requires reviewed anthesis posture in at least two tribes"
    elif len(observed) < min_scored or min_scored < 4:
        result["stop_reason"] = f"Requires at least {max(min_scored, 4)} reviewed binary tips"
    elif min(result["n_erect"], result["n_nodding"]) < 2:
        result["stop_reason"] = "Both binary states require at least two reviewed tips"
    else:
        # Remove unscored branches: otherwise unconstrained hidden flips on them
        # inflate the reported counts of independent state changes.
        inference_tree = deepcopy(tree)
        for leaf in list(inference_tree.get_terminals()):
            if leaf.name not in observed:
                inference_tree.prune(leaf.name)
        rate, lnL, boundary = rate_fit(inference_tree, observed)
        result["status"] = "conditional_trait_only_Mk_pilot"
        result["fitted_symmetric_Mk_rate_per_branch_unit"] = rate
        result["Mk_log_likelihood"] = lnL
        result["rate_optimizer_boundary"] = boundary
        result["transitions"] = stochastic_maps(
            inference_tree, observed, rate, draws=draws,
            tribe_by_tip={t["tip"]: t["tribe"] for t in taxa} if stage == "T3" else None
        )
        result["stop_reason"] = "No joint rain-temperature historical model, no camera-gravity or biological fitness calibration; Mk pilot is not a rain-adaptation test"
    destination = Path(output_path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(json.dumps(result, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
                           encoding="utf-8")
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for option in ("scores", "crosswalk", "tree", "manifest", "output"):
        parser.add_argument("--" + option, required=True)
    parser.add_argument("--stage", choices=["T2", "T3"], default="T2")
    parser.add_argument("--climate")
    parser.add_argument("--draws", type=int, default=1000)
    parser.add_argument("--min-scored", type=int, default=20)
    opts = parser.parse_args()
    result = execute(opts.scores, opts.crosswalk, opts.tree, opts.manifest,
                     opts.output, stage=opts.stage, climate_path=opts.climate,
                     draws=opts.draws, min_scored=opts.min_scored)
    print(json.dumps({"status": result["status"], "stop_reason": result["stop_reason"],
                      "scored_taxa": result["scored_taxa"]}, ensure_ascii=False))
    return 0 if result["status"] == "conditional_trait_only_Mk_pilot" else 2


if __name__ == "__main__":
    raise SystemExit(main())
