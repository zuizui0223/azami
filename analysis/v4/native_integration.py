"""Descriptive paired phenotype integration with explicit conditional intervals."""
import argparse
import hashlib
import itertools
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd

from analysis.v4.native_integration_audit import eligible
from analysis.v4.native_trait_models import prepare

FEATURES = {
    "angle": ["orientation_image_vertical_angle"],
    "lightness": ["corolla_lab_lightness"],
    "chroma": ["corolla_lab_chroma"],
    "hue": ["corolla_hue_sin", "corolla_hue_cos"],
    "aspect": ["capitulum_outline_aspect_ratio"],
    "circularity": ["capitulum_outline_circularity"],
    "solidity": ["capitulum_outline_solidity"],
    "width_profile": ["capitulum_width_profile_cv"],
}
PAIRS = list(itertools.combinations(FEATURES, 2))


def summaries(values, groups):
    labels, inverse = np.unique(groups, return_inverse=True)
    means, covariances, indices = [], [], []
    for group in range(len(labels)):
        ix = np.flatnonzero(inverse == group)
        block = values[ix]
        mean = block.mean(axis=0)
        centred = block-mean
        means.append(mean)
        covariances.append(centred.T@centred/len(ix))
        indices.append(ix)
    return np.array(means), np.array(covariances), indices


def covariance_levels(means, covariances, weights):
    weights = np.asarray(weights, float)
    weights = weights/weights.sum()
    mean = weights@means
    centred = means-mean
    return (np.einsum("g,gij->ij", weights, covariances),
            (centred*weights[:, None]).T@centred)


def association(covariance, left, right):
    xx = covariance[np.ix_(left, left)]
    yy = covariance[np.ix_(right, right)]
    xy = covariance[np.ix_(left, right)]
    denominator = np.sqrt(np.sum(xx*xx)*np.sum(yy*yy))
    if denominator <= 1e-15:
        return np.array([np.nan, np.nan])
    rv = np.sum(xy*xy)/denominator
    signed = (xy[0, 0]/np.sqrt(xx[0, 0]*yy[0, 0])
              if len(left) == len(right) == 1 else np.nan)
    return np.array([rv, signed])


def interval(values):
    finite = values[np.isfinite(values)]
    if len(finite) < .95*len(values):
        return {"status": "insufficient_finite_replicates", "finite": len(finite)}
    return {"status": "estimated", "finite": len(finite),
            "lower": float(np.quantile(finite, .025)),
            "upper": float(np.quantile(finite, .975))}


def analyse(raw, adjusted, groups, columns, pairs, rng, replicates=1000, equal_n=False):
    summaries_by_version = [summaries(values, groups) for values in (raw, adjusted)]
    n_species = len(summaries_by_version[0][0])
    if n_species < 20:
        return {"status": "insufficient_species", "n": len(raw), "species": n_species}
    indices = {name: [columns.index(v) for v in FEATURES[name]] for pair in pairs for name in pair}
    point = np.empty((2, 2, len(pairs), 2))
    boot = np.empty((replicates, 2, 2, len(pairs), 2))
    for version, (means, cov, _) in enumerate(summaries_by_version):
        for scale, matrix in enumerate(covariance_levels(means, cov, np.ones(n_species))):
            for p, (left, right) in enumerate(pairs):
                point[version, scale, p] = association(matrix, indices[left], indices[right])
    for replicate in range(replicates):
        weights = np.bincount(rng.integers(0, n_species, n_species), minlength=n_species)
        for version, (means, cov, _) in enumerate(summaries_by_version):
            for scale, matrix in enumerate(covariance_levels(means, cov, weights)):
                for p, (left, right) in enumerate(pairs):
                    boot[replicate, version, scale, p] = association(matrix, indices[left], indices[right])
    rows = []
    for p, (left, right) in enumerate(pairs):
        row = {"left": left, "right": right, "estimates": {}}
        for version, version_name in enumerate(("raw", "adjusted")):
            for scale, scale_name in enumerate(("within", "between")):
                entry = {"rv": float(point[version, scale, p, 0]),
                         "rv_interval": interval(boot[:, version, scale, p, 0])}
                if np.isfinite(point[version, scale, p, 1]):
                    entry.update(correlation=float(point[version, scale, p, 1]),
                                 correlation_interval=interval(boot[:, version, scale, p, 1]))
                row["estimates"][version_name+"/"+scale_name] = entry
            row["estimates"][version_name+"/between_minus_within_rv"] = {
                "estimate": float(point[version, 1, p, 0]-point[version, 0, p, 0]),
                "interval": interval(boot[:, version, 1, p, 0]-boot[:, version, 0, p, 0])}
        for scale, scale_name in enumerate(("within", "between")):
            row["estimates"]["raw_minus_adjusted/"+scale_name] = {
                "estimate": float(point[0, scale, p, 0]-point[1, scale, p, 0]),
                "interval": interval(boot[:, 0, scale, p, 0]-boot[:, 1, scale, p, 0])}
        rows.append(row)
    result = {"status": "estimated", "n": len(raw), "species": n_species, "pairs": rows}
    if equal_n:
        diagnostics = np.empty((replicates, 2, len(pairs)))
        for rep in range(replicates):
            chosen = np.array([rng.choice(ix) for ix in summaries_by_version[0][2]])
            for version, values in enumerate((raw, adjusted)):
                means = summaries_by_version[version][0]
                centred = values[chosen]-means
                cov = np.cov(centred, rowvar=False, bias=True)
                for p, (left, right) in enumerate(pairs):
                    diagnostics[rep, version, p] = point[version, 1, p, 0]-association(
                        cov, indices[left], indices[right])[0]
        result["equal_n_diagnostic"] = {
            name: {"median_pair_delta_q025_q50_q975":
                   np.quantile(np.median(diagnostics[:, version], axis=1), [.025, .5, .975]).tolist(),
                   "scope": "Resampling range, not a confidence interval; between means fixed."}
            for version, name in enumerate(("raw", "adjusted"))}
    return result


def execute(table, models_path, output):
    model = json.loads(models_path.read_text())
    if hashlib.sha256(table.read_bytes()).hexdigest() != model["input_sha256"]:
        raise ValueError("Integration input differs from gradient model")
    selection = model["selection"]
    predictors, scales = selection["selected"], selection["scales"]
    frame = prepare(pd.read_csv(table, low_memory=False), predictors)
    if not frame.obs_id.is_unique:
        raise ValueError("Duplicate observation identity")
    columns = sum(FEATURES.values(), [])
    components = {endpoint: component for module in model["models"].values()
                  for endpoint, component in module["components"].items()}
    coordinate_scales = {v: float(frame.loc[eligible(frame, [v]), v].std(ddof=0)) for v in columns}
    hue_scale = np.sqrt(sum(coordinate_scales[v]**2 for v in FEATURES["hue"])/2)
    for v in FEATURES["hue"]:
        coordinate_scales[v] = hue_scale
    if any(not np.isfinite(s) or s <= 0 for s in coordinate_scales.values()):
        raise ValueError("Invalid phenotype coordinate scaling")
    adjusted = frame[columns].copy()
    for endpoint in columns:
        component = components[endpoint]
        if component["status"] != "estimated":
            raise ValueError("Missing environmental component")
        for variable in predictors:
            for suffix in ("__between", "__within"):
                term = variable+suffix
                value = frame[term]-(scales[variable]["mean"] if suffix == "__between" else 0)
                adjusted[endpoint] -= value/scales[variable]["sd"] * \
                    component["effects"][term]["raw_units_per_environment_sd"]
    env_complete = np.isfinite(frame[predictors]).all(axis=1)
    rng = np.random.default_rng(20261020)
    def run(endpoints, pairs, equal_n):
        mask = env_complete & eligible(frame, endpoints)
        counts = frame.loc[mask].groupby("accepted_key").size()
        mask &= frame.accepted_key.isin(counts[counts >= 5].index)
        denom = np.array([coordinate_scales[v] for v in endpoints])
        return analyse(frame.loc[mask, endpoints].to_numpy()/denom,
                       adjusted.loc[mask, endpoints].to_numpy()/denom,
                       frame.loc[mask, "accepted_key"].to_numpy(), endpoints, pairs, rng,
                       equal_n=equal_n)
    common = run(columns, PAIRS, True)
    pair_specific = {}
    for left, right in PAIRS:
        print("Integration:", left, right, flush=True)
        pair_specific[left+"/"+right] = run(FEATURES[left]+FEATURES[right], [(left, right)], False)
    result = {
        "run_id": os.environ.get("GITHUB_RUN_ID"), "git_sha": os.environ.get("GITHUB_SHA"),
        "input_sha256": model["input_sha256"],
        "model_sha256": hashlib.sha256(models_path.read_bytes()).hexdigest(),
        "contract_sha256": hashlib.sha256(Path("analysis/v4/INTEGRATION_CONTRACT.md").read_bytes()).hexdigest(),
        "coordinate_scales": coordinate_scales, "common_cohort": common,
        "pair_specific": pair_specific, "replicates": 1000, "seed": 20261020,
        "interval_scope": "Species-cluster bootstrap; adjusted intervals conditional on fixed "
                          "environment coefficients. No spatial or coefficient-estimation uncertainty.",
        "claim_boundary": "Descriptive phenotype association, not latent random-effect/genetic "
                          "covariance, adaptation, or a causal environment-explained fraction."}
    output.mkdir(parents=True, exist_ok=True)
    (output/"native_integration.json").write_text(json.dumps(result, indent=2, allow_nan=False)+"\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    for key in ("table", "models", "output"):
        parser.add_argument("--"+key, type=Path, required=True)
    args = parser.parse_args()
    execute(args.table, args.models, args.output)
