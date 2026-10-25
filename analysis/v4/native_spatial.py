"""Fixed spatial concentration and conditional-residual diagnostics."""
import argparse
import hashlib
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

from analysis.v4.native_integration_audit import MODULES, eligible
from analysis.v4.native_trait_models import prepare
from analysis.v4.observation_models import fit_module


def cells(latitude, longitude):
    lat, lon = np.asarray(latitude, float), np.asarray(longitude, float)
    if not np.isfinite(lat).all() or not np.isfinite(lon).all() or \
            np.any(abs(lat) > 90) or np.any(abs(lon) > 180):
        raise ValueError("Invalid coordinates")
    x = np.floor(((lon+180) % 360)/5).astype(int)
    y = np.clip(np.floor((np.sin(np.deg2rad(lat))+1)/.05).astype(int), 0, 39)
    return y*72+x


def conditional_residuals(sample, component, predictors, scales):
    prediction = np.full(len(sample), component["reference_prediction"])
    for variable in predictors:
        for suffix in ("__between", "__within"):
            term = variable+suffix
            value = sample[term].to_numpy()-(scales[variable]["mean"] if suffix == "__between" else 0)
            prediction += value/scales[variable]["sd"] * \
                component["effects"][term]["raw_units_per_environment_sd"]
    endpoint = sample.attrs["endpoint"]
    residual = pd.Series(sample[endpoint].to_numpy()-prediction, index=sample.index)
    group = residual.groupby(sample.accepted_key)
    ratio = component["variance_ratio"]
    random_intercept = ratio*group.transform("sum")/(1+ratio*group.transform("size"))
    return (residual-random_intercept).to_numpy()


def moran(cell_ids, residuals, permutations=999):
    means = pd.Series(residuals).groupby(np.asarray(cell_ids)).mean()
    if len(means) < 20 or means.std() <= 1e-12:
        return {"status": "insufficient_spatial_support", "cells": len(means)}
    ids = means.index.to_numpy()
    lat = np.arcsin(-1+(ids//72+.5)*.05)
    lon = np.deg2rad(-180+(ids % 72+.5)*5)
    xyz = np.column_stack([np.cos(lat)*np.cos(lon), np.cos(lat)*np.sin(lon), np.sin(lat)])
    neighbors = cKDTree(xyz).query(xyz, k=9)[1][:, 1:]
    edges = set()
    for i, row in enumerate(neighbors):
        for j in row:
            edges.add((min(i, int(j)), max(i, int(j))))
    left, right = np.array(sorted(edges)).T
    values = means.to_numpy()-means.mean()
    def statistic(v):
        return float(len(v)*np.sum(v[left]*v[right])/(len(left)*np.sum(v*v)))
    observed = statistic(values)
    expected = -1/(len(values)-1)
    rng = np.random.default_rng(20261020)
    exceed = sum(abs(statistic(rng.permutation(values))-expected) >= abs(observed-expected)
                 for _ in range(permutations))
    return {"status": "estimated", "cells": len(means), "I": observed,
            "permutation_p": (exceed+1)/(permutations+1), "permutations": permutations,
            "cell_mean_residual_sd": float(means.std()),
            "scope": "Cell-mean residual diagnostic; not proof of independent errors."}


def balanced_indices(sample, seed=20261020):
    ordered = sample.sort_values("obs_id")
    return ordered.groupby(["accepted_key", "_cell"], sort=True).sample(
        n=1, random_state=seed).index


def execute(table, models_path, output):
    models = json.loads(models_path.read_text())
    if hashlib.sha256(table.read_bytes()).hexdigest() != models["input_sha256"]:
        raise ValueError("Spatial input differs from main model")
    selection = models["selection"]
    predictors, scales = selection["selected"], selection["scales"]
    frame = prepare(pd.read_csv(table, low_memory=False), predictors)
    if not frame.obs_id.is_unique:
        raise ValueError("Duplicate observation identities")
    frame["_cell"] = cells(frame.latitude, frame.longitude)
    complete = np.isfinite(frame[predictors]).all(axis=1)
    results = {}
    for module, outcomes in MODULES.items():
        for endpoint in outcomes:
            print("Space:", endpoint, flush=True)
            sample = frame.loc[complete & eligible(frame, [endpoint])].copy()
            reference = models["models"][module]["components"][endpoint]
            if len(sample) != reference["n"] or reference["status"] != "estimated":
                raise ValueError("Endpoint sample/model mismatch")
            sample.attrs["endpoint"] = endpoint
            residuals = conditional_residuals(sample, reference, predictors, scales)
            subset = sample.loc[balanced_indices(sample)]
            fit = fit_module(subset, [endpoint], scales, predictors, required_environment=predictors)
            results[endpoint] = {
                "main_n": len(sample), "balanced_n": len(subset),
                "main_species": int(sample.accepted_key.nunique()),
                "balanced_species": int(subset.accepted_key.nunique()),
                "main_effects": reference["effects"],
                "residual_spatial_diagnostic": moran(sample["_cell"].to_numpy(), residuals),
                "balanced_fit": fit}
    result = {
        "run_id": os.environ.get("GITHUB_RUN_ID"), "git_sha": os.environ.get("GITHUB_SHA"),
        "input_sha256": models["input_sha256"], "endpoints": results,
        "model_sha256": hashlib.sha256(models_path.read_bytes()).hexdigest(),
        "contract_sha256": hashlib.sha256(Path("analysis/v4/SPACE_CONTRACT.md").read_bytes()).hexdigest(),
        "scope": "Fixed spatial concentration sensitivity; main model not replaced. "
                 "Species-cluster intervals remain conditional on that error treatment. "
                 "Neither native filtering nor thinning certifies spatial independence."}
    output.mkdir(parents=True, exist_ok=True)
    (output/"native_spatial.json").write_text(json.dumps(result, indent=2, allow_nan=False)+"\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    for key in ("table", "models", "output"):
        parser.add_argument("--"+key, type=Path, required=True)
    args = parser.parse_args()
    execute(args.table, args.models, args.output)
