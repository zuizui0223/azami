"""Shared observation-level environment stage for the integration workflow."""
import argparse
import hashlib
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests

from analysis.v4.native_integration_audit import MODULES, eligible
from analysis.v4.observation_models import fit_module
from analysis.v4.environment_gradients import exposure_support, gradient_points
from analysis.v4.select_native_environment import diagnose


def prepare(frame, predictors):
    frame = frame.copy()
    for variable in predictors:
        frame[variable] = frame[variable].where(np.isfinite(frame[variable]))
        frame[variable+"__between"] = frame.groupby("accepted_key")[variable].transform("mean")
        frame[variable+"__within"] = frame[variable]-frame[variable+"__between"]
    return frame


def execute(table, receipt_path, selection_path, output):
    receipt = json.loads(receipt_path.read_text())
    selection = json.loads(selection_path.read_text())
    digest = hashlib.sha256(table.read_bytes()).hexdigest()
    if digest != receipt["prepared_table_sha256"] or digest != selection["input_sha256"]:
        raise ValueError("Input and environment selection do not match")
    if selection["status"] != "selected" or receipt.get("environment_decoding_status") != \
            "physical_units_gsp_invalid_encoding_masked_v1":
        raise ValueError("Passed environment and decoding gates required")
    if selection["contract_sha256"] != hashlib.sha256(
            Path("analysis/v4/ENVIRONMENT_CONTRACT.md").read_bytes()).hexdigest():
        raise ValueError("Environment contract changed since selection")
    predictors, scales = selection["selected"], selection["scales"]
    frame = pd.read_csv(table, low_memory=False)
    if not frame.obs_id.is_unique or len(frame) != receipt["native_observations"]:
        raise ValueError("Native input population changed")
    frame = prepare(frame, predictors)
    complete = np.isfinite(frame[predictors]).all(axis=1)
    models, diagnostics = {}, {}
    rows = []
    for module, outcomes in MODULES.items():
        print("Environmental component:", module, flush=True)
        model = fit_module(frame, outcomes, scales, predictors, required_environment=predictors)
        models[module] = model
        for endpoint in outcomes:
            sample = frame.loc[complete & eligible(frame, [endpoint])]
            support = exposure_support(sample, predictors)
            columns = [v+s for v in predictors for s in ("__between", "__within")]
            try:
                design = diagnose({"endpoint_design": sample}, columns)
            except ValueError as error:
                design = {"status": "not_identifiable", "reason": str(error)}
            diagnostics[endpoint] = {"n": len(sample), "design": design, "support": support}
            component = model.get("components", {}).get(endpoint, {})
            if component.get("status") != "estimated":
                continue
            if component["n"] != len(sample):
                raise ValueError("Endpoint-specific cohort mismatch")
            component["gradients"] = {}
            for term, effect in component["effects"].items():
                variable = term.rsplit("__", 1)[0]
                component["gradients"][term] = gradient_points(
                    effect, support[term], scales[variable]["sd"])
                increment = {"chelsa_bio01": 1., "chelsa_bio04": 100.,
                             "chelsa_bio12": 100., "chelsa_bio15": 10.,
                             "chelsa_rsds_mean": 1.}[variable]
                factor = increment/scales[variable]["sd"]
                rows.append(dict(module=module, endpoint=endpoint, term=term,
                                 n=component["n"], species=component["n_species"],
                                 environmental_increment=increment,
                                 trait_change=effect["raw_units_per_environment_sd"]*factor,
                                 lower95=effect["raw_ci95"][0]*factor,
                                 upper95=effect["raw_ci95"][1]*factor))
    families = {}
    for name, modules in (("main", ["orientation", "colour", "outline"]),
                          ("si", ["bract_si"])):
        tests = [models[m].get("tests", {}).get(v+s, {})
                 for m in modules for v in predictors for s in ("__between", "__within")]
        adjusted = multipletests([t.get("p", 1.) for t in tests], method="holm")[1]
        for test, value in zip(tests, adjusted):
            if test.get("status") == "estimated":
                test["holm_p"] = float(value)
                test["family"] = name
        families[name] = {"size": len(tests), "method": "Holm",
                          "unestimated_tests_retained_as_p1": True}
    result = {
        "run_id": os.environ.get("GITHUB_RUN_ID"), "git_sha": os.environ.get("GITHUB_SHA"),
        "input_sha256": digest, "selection": selection, "models": models,
        "diagnostics": diagnostics, "families": families,
        "scope": "Conditional observational gradients, not causal or adaptation tests. "
                 "Response-specific species random intercepts; species-cluster intervals. "
                 "No spatial/phylogenetic correction yet. This is the environment "
                 "component, not completion of the integration analysis.",
        "units": {"chelsa_bio01": "degC", "chelsa_bio04": "BIO4 index",
                  "chelsa_bio12": "mm/year", "chelsa_bio15": "CV percentage points",
                  "chelsa_rsds_mean": "MJ/m2/day"}}
    output.mkdir(parents=True, exist_ok=True)
    (output/"native_trait_models.json").write_text(json.dumps(result, indent=2, allow_nan=False)+"\n")
    pd.DataFrame(rows).to_csv(output/"native_trait_effects.csv", index=False)
    if any(m["status"] != "estimated" for m in models.values()):
        raise SystemExit("Some modules not identifiable; inspect saved partial results.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    for key in ("table", "receipt", "selection", "output"):
        parser.add_argument("--"+key, type=Path, required=True)
    args = parser.parse_args()
    execute(args.table, args.receipt, args.selection, args.output)
