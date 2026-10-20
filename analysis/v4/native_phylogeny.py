"""Observation-level phylogenetic REML using exact grouped crossproducts."""
import argparse
import hashlib
import io
import json
import os
from pathlib import Path
import urllib.request
import zipfile

import numpy as np
import pandas as pd
from Bio import Phylo
from scipy.linalg import cho_factor, cho_solve
from scipy.optimize import minimize
from scipy.stats import norm

from analysis.v4.native_integration_audit import MODULES, eligible
from analysis.v4.native_trait_models import prepare
from reproducibility.run_current_analysis import INPUTS, ZENODO_URL, ZENODO_SHA, verified


def statistics(x, y, groups):
    labels, index = np.unique(groups, return_inverse=True)
    sizes = np.bincount(index).astype(float)
    sx = np.zeros((len(labels), x.shape[1]))
    sy = np.zeros(len(labels))
    np.add.at(sx, index, x)
    np.add.at(sy, index, y)
    return dict(labels=labels, sizes=sizes, sx=sx, sy=sy,
                xx=x.T @ x, xy=x.T @ y, yy=float(y @ y),
                n=len(y), p=x.shape[1])


def evaluate(s, a, ratios):
    """Profile REML; covariance is I + Z G Z' in residual-variance units."""
    g = ratios[0] * a + ratios[1] * np.eye(len(a))
    root = np.sqrt(s["sizes"])
    b = np.eye(len(a)) + root[:, None] * g * root[None, :]
    ch = cho_factor(b)
    # Woodbury in a form valid even when G is singular or zero.
    h = (np.eye(len(a)) - cho_solve(ch, np.eye(len(a)))) / root[:, None] / root[None, :]
    info = s["xx"] - s["sx"].T @ h @ s["sx"]
    rhs = s["xy"] - s["sx"].T @ h @ s["sy"]
    ic = cho_factor(info)
    beta = cho_solve(ic, rhs)
    sse = s["yy"] - s["sy"] @ h @ s["sy"] - beta @ rhs
    if sse <= 0:
        raise ValueError("Nonpositive residual sum of squares")
    variance = sse / (s["n"] - s["p"])
    objective = ((s["n"] - s["p"]) * np.log(variance)
                 + 2 * np.log(np.diag(ch[0])).sum()
                 + 2 * np.log(np.diag(ic[0])).sum())
    return float(objective), beta, variance * cho_solve(ic, np.eye(s["p"])), variance


def fit(s, a, phylogenetic=True):
    def objective(v):
        try:
            ratios = np.exp(v) if phylogenetic else [0., np.exp(v[0])]
            return evaluate(s, a, ratios)[0]
        except (ValueError, np.linalg.LinAlgError):
            return 1e100
    starts = [[-2., -2.], [1., -6.], [-6., 1.]] if phylogenetic else [[-2.], [1.]]
    options = [minimize(objective, start, method="L-BFGS-B",
                        bounds=[(-18., 14.)] * len(start),
                        options={"ftol": 1e-11, "maxiter": 300})
               for start in starts]
    good = [o for o in options if o.success and np.isfinite(o.fun) and o.fun < 1e90]
    if not good:
        raise ValueError("No converged variance fit")
    best = min(good, key=lambda o: o.fun)
    ratios = np.exp(best.x) if phylogenetic else np.array([0., np.exp(best.x[0])])
    # Explicit zero-variance boundaries, including independent-species model.
    candidates = [ratios, np.zeros(2)]
    if phylogenetic:
        for k in (0, 1):
            def one(v):
                r = np.zeros(2)
                r[k] = np.exp(v[0])
                return evaluate(s, a, r)[0]
            opt = minimize(one, [-2.], method="L-BFGS-B", bounds=[(-18., 14.)])
            if opt.success:
                r = np.zeros(2)
                r[k] = np.exp(opt.x[0])
                candidates.append(r)
    ratios = min(candidates, key=lambda r: evaluate(s, a, r)[0])
    obj, beta, cov, variance = evaluate(s, a, ratios)
    return beta, cov, {"objective": obj, "phylogenetic_variance": float(ratios[0]*variance),
                       "independent_species_variance": float(ratios[1]*variance),
                       "residual_variance": variance,
                       "boundary": bool(np.any(ratios < 2e-8) or np.any(ratios > np.exp(13.9)))}


def restore_trees():
    with urllib.request.urlopen(ZENODO_URL, timeout=180) as response:
        bundle = verified(response.read(), ZENODO_SHA, "Zenodo bundle")
    artifact, digest, members = INPUTS["historical"]
    with zipfile.ZipFile(io.BytesIO(bundle)) as z:
        raw = verified(z.read(f"inputs/artifact-{artifact}-historical.zip"), digest, "trees")
    trees = []
    with zipfile.ZipFile(io.BytesIO(raw)) as z:
        for member, (_, sha) in members.items():
            content = verified(z.read(member), sha, member).decode()
            for i, tree in enumerate(Phylo.parse(io.StringIO(content), "newick")):
                trees.append((Path(member).stem + f"_{i}", tree))
    if len(trees) != 52:
        raise ValueError("Expected 52 frozen placements")
    return trees


def tree_covariance(tree, names):
    tips = {t.name: t for t in tree.get_terminals()}
    if len(tips) != len(tree.get_terminals()):
        raise ValueError("Duplicate tree tips")
    for clade in tree.find_clades():
        if clade is not tree.root and (clade.branch_length is None or
                not np.isfinite(clade.branch_length) or clade.branch_length < 0):
            raise ValueError("Missing, negative or nonfinite branch length")
    depth = tree.depths()
    a = np.empty((len(names), len(names)))
    for i, left in enumerate(names):
        for j in range(i+1):
            a[i, j] = a[j, i] = depth[tree.common_ancestor(tips[left], tips[names[j]])]
    if np.mean(np.diag(a)) <= 0:
        raise ValueError("No positive tree depth")
    return a / np.mean(np.diag(a))


def execute(table, output):
    frozen = Path("analysis/v4/results/20261020_native_integration")
    models = json.loads((frozen/"native_trait_models.json").read_text())
    if hashlib.sha256(table.read_bytes()).hexdigest() != models["input_sha256"]:
        raise ValueError("Native table SHA mismatch")
    selection = models["selection"]
    predictors, scales = selection["selected"], selection["scales"]
    frame = prepare(pd.read_csv(table, low_memory=False), predictors)
    trees = restore_trees()
    tipsets = [set(t.name for t in tree.get_terminals()) for _, tree in trees]
    if any(tips != tipsets[0] for tips in tipsets):
        raise ValueError("Placement tip sets differ")
    # accepted_key is an authority identifier, NOT a scientific name.
    # Use only exact recorded scientific names; do not silently invent placements.
    names = frame.taxon_name.astype(str).str.strip().str.replace(" ", "_", regex=False)
    mapping, audit = {}, []
    for key, group in frame.assign(_name=names).groupby("accepted_key"):
        candidates = sorted(set(group._name) & tipsets[0])
        mapping[key] = candidates[0] if len(candidates) == 1 else None
        audit.append(dict(accepted_key=str(key), source_names=sorted(set(group._name)),
                          tree_matches=candidates, status="matched" if len(candidates)==1
                          else "ambiguous" if candidates else "absent"))
    assigned = [value for value in mapping.values() if value is not None]
    if len(assigned) != len(set(assigned)):
        raise ValueError("Multiple authority taxa map to one tree tip")
    frame["tree_key"] = frame.accepted_key.map(mapping)
    terms = [v+s for v in predictors for s in ("__between", "__within")]
    rows, coverage, failures = [], {}, []
    output.mkdir(parents=True, exist_ok=True)
    (output/"taxon_mapping.json").write_text(json.dumps(audit, indent=2))
    for module in ("orientation", "colour", "outline"):
        for endpoint in MODULES[module]:
            valid = eligible(frame, [endpoint]) & np.isfinite(frame[terms]).all(axis=1)
            all_sample = frame.loc[valid]
            sample = all_sample.loc[all_sample.tree_key.isin(tipsets[0])]
            coverage[endpoint] = dict(full_n=len(all_sample), covered_n=len(sample),
                full_species=int(all_sample.accepted_key.nunique()),
                covered_species=int(sample.accepted_key.nunique()),
                missing_taxa=sorted(set(all_sample.loc[all_sample.tree_key.isna(), "taxon_name"])))
            (output/"coverage.json").write_text(json.dumps(coverage, indent=2))
            if sample.accepted_key.nunique() < 30:
                failures.append({"endpoint": endpoint, "reason": "fewer than 30 tree-covered taxa"})
                continue
            x = np.column_stack([np.ones(len(sample))] +
                [sample[t].to_numpy()/scales[t.rsplit("__", 1)[0]]["sd"] for t in terms])
            y = sample[endpoint].to_numpy(float)
            center, scale = y.mean(), y.std()
            s = statistics(x, (y-center)/scale, sample.tree_key.to_numpy())
            for number, (name, tree) in enumerate(trees):
                print(endpoint, name, flush=True)
                a = tree_covariance(tree, s["labels"])
                variants = [("phylogenetic", True)]
                if number == 0:
                    variants.insert(0, ("same_cohort_independent", False))
                for variant, phy in variants:
                    try:
                        beta, cov, diagnostic = fit(s, a, phy)
                        for k, term in enumerate(terms, 1):
                            factor = scale/scales[term.rsplit("__", 1)[0]]["sd"]
                            se = np.sqrt(cov[k, k])
                            rows.append(dict(endpoint=endpoint, tree=name, model=variant,
                                term=term, n=len(sample), species=len(s["labels"]),
                                coefficient=beta[k]*factor, lower95=(beta[k]-1.96*se)*factor,
                                upper95=(beta[k]+1.96*se)*factor,
                                p=2*norm.sf(abs(beta[k]/se)), **diagnostic))
                    except (ValueError, np.linalg.LinAlgError) as error:
                        failures.append(dict(endpoint=endpoint, tree=name, model=variant, reason=str(error)))
                pd.DataFrame(rows).to_csv(output/"coefficients.csv", index=False)
    (output/"receipt.json").write_text(json.dumps(dict(
        run_id=os.environ.get("GITHUB_RUN_ID"), git_sha=os.environ.get("GITHUB_SHA"),
        input_sha256=models["input_sha256"], tree_bundle_sha256=ZENODO_SHA,
        failures=failures, intervals="model-based normal conditional on fitted variances",
        scope="Historical-placement sensitivity; no molecular-tree or joint-colour claim"),
        indent=2))
    if failures:
        raise SystemExit("Some fits unavailable; inspect partial results")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--table", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    execute(**vars(parser.parse_args()))
