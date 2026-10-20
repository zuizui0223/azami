"""Joint covariance-component estimation by PSD-constrained quadratic moments."""
import argparse
import hashlib
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import wishart

from analysis.v4.native_integration import FEATURES, PAIRS, association
from analysis.v4.native_integration_audit import eligible
from analysis.v4.native_trait_models import prepare
from analysis.v4.native_phylogeny import restore_trees, tree_covariance


def psd(matrix):
    w, v = np.linalg.eigh((matrix+matrix.T)/2)
    return (v*np.maximum(w, 0))@v.T


def square_root(matrix):
    w, v = np.linalg.eigh((matrix+matrix.T)/2)
    return v*np.sqrt(np.maximum(w, 0))


def kernels(a, counts):
    g = len(a)
    h = np.eye(g)-np.ones((g, g))/g
    c = np.array([h@a@h, h])
    gram = np.einsum("aij,bij->ab", c, c)
    condition = float(np.linalg.cond(gram))
    if condition > 1e8:
        raise ValueError("Phylogenetic and independent kernels not separable")
    noise = (h/np.asarray(counts)[None, :])@h
    return c, gram, np.einsum("aij,ji->a", c, noise), condition


def solve_components(s, gram, tolerance=1e-9, maxiter=10000):
    """Convex block coordinate descent; both component matrices remain PSD."""
    unconstrained = np.linalg.solve(gram, s.reshape(2, -1)).reshape(s.shape)
    p, u = psd(unconstrained[0]), psd(unconstrained[1])
    for iteration in range(maxiter):
        previous_p, previous_u = p, u
        p = psd((s[0]-gram[0, 1]*u)/gram[0, 0])
        u = psd((s[1]-gram[0, 1]*p)/gram[1, 1])
        change = max(np.linalg.norm(p-previous_p), np.linalg.norm(u-previous_u))
        if change < tolerance*(1+np.linalg.norm(p)+np.linalg.norm(u)):
            return np.array([p, u]), unconstrained, iteration+1
    raise ValueError("PSD covariance solver did not converge")


def estimate(means, r, kernel):
    c, gram, correction, _ = kernel
    s = np.einsum("gi,kgh,hj->kij", means, c, means, optimize=True)
    s -= correction[:, None, None]*r
    components, unbounded, iterations = solve_components(s, gram)
    return np.array([components[0], components[1], r]), unbounded, iterations


def observed_summaries(values, groups):
    labels, inv = np.unique(groups, return_inverse=True)
    counts = np.bincount(inv)
    means = np.zeros((len(labels), values.shape[1]))
    np.add.at(means, inv, values)
    means /= counts[:, None]
    residual = values-means[inv]
    df = len(values)-len(labels)
    return labels, counts, means, residual.T@residual/df, df


def metrics(matrices, indices):
    result = []
    total = np.diag(matrices.sum(axis=0))
    for matrix in matrices:
        fractions = np.diag(matrix)/total
        values = []
        for left, right in PAIRS:
            ix, iy = indices[left], indices[right]
            # Correlations near vanishing component variances are uninformative.
            if min(fractions[ix].sum(), fractions[iy].sum()) < 1e-6:
                values.append([np.nan, np.nan])
            else:
                values.append(association(matrix, ix, iy))
        result.append(values)
    return np.asarray(result)


def run_tree(means, r, counts, df, a, indices, rng, replicates):
    kernel = kernels(a, counts)
    matrices, unbounded, iterations = estimate(means, r, kernel)
    point = metrics(matrices, indices)
    boot = np.empty((replicates, *point.shape))
    boot_fractions = np.empty((replicates, 3, len(r)))
    ar = square_root(a)
    pr, ur, rr = [square_root(m) for m in matrices]
    for b in range(replicates):
        synthetic = (ar@rng.normal(size=means.shape)@pr.T
                     + rng.normal(size=means.shape)@ur.T
                     + (rng.normal(size=means.shape)@rr.T)/np.sqrt(counts)[:, None])
        rb = wishart.rvs(df=df, scale=r, random_state=rng)/df
        mb, _, _ = estimate(synthetic, rb, kernel)
        boot[b] = metrics(mb, indices)
        boot_fractions[b] = np.array([np.diag(m) for m in mb])/np.diag(mb.sum(axis=0))
    intervals = []
    for level in range(3):
        for pair, (left, right) in enumerate(PAIRS):
            row = dict(level=("phylogenetic", "independent_species", "within")[level],
                       left=left, right=right)
            for k, key in enumerate(("rv", "correlation")):
                values = boot[:, level, pair, k]
                finite = values[np.isfinite(values)]
                row[key] = float(point[level, pair, k]) if np.isfinite(point[level, pair, k]) else None
                row[key+"_finite"] = len(finite)
                bounds = np.quantile(finite, [.025, .975]) if len(finite)>=.95*replicates else [None, None]
                row[key+"_lower"], row[key+"_upper"] = bounds
            intervals.append(row)
    return intervals, dict(
        matrices=matrices.tolist(), unconstrained=unbounded.tolist(),
        eigenvalues=[np.linalg.eigvalsh(m).tolist() for m in matrices],
        unconstrained_eigenvalues=[np.linalg.eigvalsh(m).tolist() for m in unbounded],
        gram_condition=kernel[3], iterations=iterations,
        variance_fractions=(np.array([np.diag(m) for m in matrices])/
                            np.diag(matrices.sum(axis=0))).tolist(),
        variance_fraction_intervals=np.quantile(boot_fractions, [.025, .975], axis=0).tolist())


def execute(table, output, replicates=200):
    frozen = Path("analysis/v4/results/20261020_native_integration")
    model = json.loads((frozen/"native_trait_models.json").read_text())
    integration = json.loads((frozen/"native_integration.json").read_text())
    mapping_path = Path("analysis/v4/results/20261020_native_phylogeny/taxon_mapping.json")
    mapping = {str(r["accepted_key"]): r["tree_matches"][0]
               for r in json.loads(mapping_path.read_text()) if r["status"]=="matched"}
    if hashlib.sha256(table.read_bytes()).hexdigest() != model["input_sha256"]:
        raise ValueError("Native input hash mismatch")
    selection = model["selection"]
    predictors, scales = selection["selected"], selection["scales"]
    frame = prepare(pd.read_csv(table, low_memory=False), predictors)
    frame["tree_key"] = frame.accepted_key.astype(str).map(mapping)
    columns = sum(FEATURES.values(), [])
    mask = eligible(frame, columns) & np.isfinite(frame[predictors]).all(axis=1)
    counts = frame.loc[mask].groupby("accepted_key").size()
    mask &= frame.accepted_key.isin(counts[counts>=5].index)
    original_n = int(mask.sum())
    original_species = int(frame.loc[mask, "accepted_key"].nunique())
    if (original_n, original_species) != (integration["common_cohort"]["n"],
                                       integration["common_cohort"]["species"]):
        raise ValueError("Common cohort changed")
    sample = frame.loc[mask & frame.tree_key.notna()].copy()
    if sample.accepted_key.nunique()<30 or not sample.obs_id.is_unique:
        raise ValueError("Invalid tree-covered cohort")
    denom = np.array([integration["coordinate_scales"][v] for v in columns])
    raw = sample[columns].to_numpy(float)
    adjusted = raw.copy()
    components = {k: v for m in model["models"].values() for k, v in m["components"].items()}
    for k, endpoint in enumerate(columns):
        for variable in predictors:
            for suffix in ("__between", "__within"):
                term = variable+suffix
                value = sample[term].to_numpy()-(scales[variable]["mean"] if suffix=="__between" else 0)
                adjusted[:, k] -= value/scales[variable]["sd"] * components[endpoint]["effects"][term]["raw_units_per_environment_sd"]
    versions = {name: observed_summaries(v/denom, sample.tree_key.to_numpy())
                for name, v in (("raw", raw), ("adjusted", adjusted))}
    indices = {k: [columns.index(v) for v in values] for k, values in FEATURES.items()}
    trees = restore_trees()
    output.mkdir(parents=True, exist_ok=True)
    rows, diagnostics, failures = [], {}, []
    receipt = dict(run_id=os.environ.get("GITHUB_RUN_ID"), git_sha=os.environ.get("GITHUB_SHA"),
                   input_sha256=model["input_sha256"], n=len(sample),
                   species=int(sample.accepted_key.nunique()),
                   original_common_n=original_n, original_common_species=original_species,
                   columns=columns, replicates=replicates, seed=20261021,
                   estimator="PSD-constrained quadratic method of moments, not REML",
                   interval_scope="Conditional Gaussian parametric bootstrap, pointwise",
                   contract_sha256=hashlib.sha256(Path("analysis/v4/COVARIANCE_CONTRACT.md").read_bytes()).hexdigest())
    (output/"receipt.json").write_text(json.dumps(receipt, indent=2))
    for t, (name, tree) in enumerate(trees):
        a = tree_covariance(tree, versions["raw"][0])
        for version, (_, counts, means, r, df) in versions.items():
            print(name, version, flush=True)
            try:
                estimates, diagnostic = run_tree(means, r, counts, df, a, indices,
                    np.random.default_rng(20261021+t), replicates)
                rows.extend(dict(tree=name, version=version, **row) for row in estimates)
                diagnostics[name+"/"+version] = diagnostic
            except (ValueError, np.linalg.LinAlgError) as error:
                failures.append(dict(tree=name, version=version, reason=str(error)))
            pd.DataFrame(rows).to_csv(output/"covariance_pairs.csv", index=False)
            (output/"matrices.json").write_text(json.dumps(diagnostics, indent=2, allow_nan=False))
            (output/"failures.json").write_text(json.dumps(failures, indent=2))
    if failures:
        raise SystemExit("Some covariance decompositions failed; inspect saved evidence")


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--table", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--replicates", type=int, default=200)
    execute(**vars(p.parse_args()))
