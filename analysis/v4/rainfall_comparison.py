"""Same-cohort rainfall-only diagnostic; not a replacement primary analysis."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

from analysis.v4.observation_models import fit_module

RAIN = "chelsa_bio12"
TRAIT = "orientation_image_vertical_angle"


def comparison_frame(frame, predictors):
    """Lock the full-model cohort while retaining population-defined exposures."""
    frame = frame.copy()
    for variable in predictors:
        mean = frame.groupby("accepted_key")[variable].transform("mean")
        frame[variable+"__between"] = mean
        frame[variable+"__within"] = frame[variable]-mean
    eligible = (frame[TRAIT+"__eligible"].eq(True) & np.isfinite(frame[TRAIT])
                & np.isfinite(frame[predictors]).all(axis=1))
    return frame.loc[eligible].copy()


def execute(table, baseline_path, output):
    baseline = json.loads(baseline_path.read_text())
    if hashlib.sha256(table.read_bytes()).hexdigest() != baseline["input_sha256"]:
        raise ValueError("Comparison input differs from main fit")
    if baseline["input_receipt"].get("environment_decoding_status") != \
            "physical_units_gsp_invalid_encoding_masked_v1":
        raise ValueError("Corrected baseline required")
    predictors = baseline["predictors"]
    if RAIN not in predictors:
        raise ValueError("Main model has no annual rainfall term")
    frame = pd.read_csv(table, low_memory=False)
    if not frame.obs_id.is_unique:
        raise ValueError("Duplicate observation identities")
    sample = comparison_frame(frame, predictors)
    reference = baseline["models"]["orientation"]["components"][TRAIT]
    if len(sample) != reference["n"]:
        raise ValueError("Comparison cohort differs from main model")
    reduced = fit_module(sample, [TRAIT], baseline["scales"], [RAIN],
                         required_environment=predictors)
    if reduced["status"] != "estimated":
        raise ValueError("Rainfall-only comparison could not be estimated")
    component = reduced["components"][TRAIT]
    if (component["n"], component["n_species"]) != \
            (reference["n"], reference["n_species"]):
        raise ValueError("Comparison changed observation or species counts")
    rows = []
    for scenario, fit in (("rainfall_only", component), ("full_environment", reference)):
        for level in ("between", "within"):
            effect = fit["effects"][RAIN+"__"+level]
            factor = 100 / baseline["scales"][RAIN]["sd"]
            rows.append(dict(
                scenario=scenario, level=level, n=fit["n"],
                n_species=fit["n_species"],
                degrees_per_100mm=effect["raw_units_per_environment_sd"]*factor,
                lower95=effect["raw_ci95"][0]*factor,
                upper95=effect["raw_ci95"][1]*factor))
    result = dict(
        baseline_run=baseline["run_id"], input_sha256=baseline["input_sha256"],
        baseline_sha256=hashlib.sha256(baseline_path.read_bytes()).hexdigest(),
        full_predictors=predictors, rainfall_only_fit=reduced, comparison=rows,
        scope="Post-result same-cohort specification diagnostic, not a V3 replay, "
              "independent confirmation, causal effect, or replacement primary test. "
              "Response observations are not aggregated. Species means and scales "
              "retain the full native-population definitions. Pointwise intervals "
              "do not change the existing primary multiplicity family. This "
              "comparison cannot identify which covariate explains a difference.")
    output.mkdir(parents=True, exist_ok=True)
    (output/"rainfall_comparison.json").write_text(
        json.dumps(result, indent=2, allow_nan=False)+"\n")
    pd.DataFrame(rows).to_csv(output/"rainfall_comparison.csv", index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    for key in ("table", "baseline", "output"):
        parser.add_argument("--"+key, type=Path, required=True)
    args = parser.parse_args()
    execute(args.table, args.baseline, args.output)
