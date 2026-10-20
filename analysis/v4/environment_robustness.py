"""Fixed species-abundance and broad-space checks of the corrected main fit."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

from analysis.v4.environment_response import MODULES
from analysis.v4.observation_models import fit_module, spherical_basis


def abundance_ranking(frame, outcomes, predictors):
    eligible = np.isfinite(frame[predictors]).all(axis=1)
    any_trait = np.zeros(len(frame), dtype=bool)
    for endpoint in outcomes:
        any_trait |= frame[endpoint+"__eligible"].eq(True) & np.isfinite(frame[endpoint])
    counts = frame.loc[eligible & any_trait].groupby("accepted_key").size()
    # Stable identity-based order resolves abundance ties, without looking at effects.
    return counts.sort_index().sort_values(ascending=False, kind="stable")


def execute(table, baseline_path, output):
    baseline = json.loads(baseline_path.read_text())
    if hashlib.sha256(table.read_bytes()).hexdigest() != baseline["input_sha256"]:
        raise ValueError("Sensitivity input differs from main fit")
    if baseline["input_receipt"].get("environment_decoding_status") != \
            "physical_units_gsp_invalid_encoding_masked_v1":
        raise ValueError("Corrected baseline required")
    frame = pd.read_csv(table, low_memory=False)
    predictors, scales = baseline["predictors"], baseline["scales"]
    if not frame.obs_id.is_unique:
        raise ValueError("Duplicate observation identities")
    for variable in predictors:
        frame[variable+"__between"] = frame.groupby("accepted_key")[variable].transform("mean")
        frame[variable+"__within"] = frame[variable]-frame[variable+"__between"]
    spatial = []
    spatial_status = "not_evaluable_missing_coordinates"
    if {"latitude", "longitude"} <= set(frame) and \
            np.isfinite(frame[["latitude","longitude"]]).all(axis=None):
        basis = spherical_basis(frame.latitude.to_numpy(), frame.longitude.to_numpy())
        sd = basis.std(axis=0)
        if np.all(sd > 1e-12):
            spatial = ["broad_space_"+str(i) for i in range(basis.shape[1])]
            frame[spatial] = (basis-basis.mean(axis=0))/sd
            spatial_status = "evaluated_eight_spherical_terms"
    results, rows = {}, []
    for module, outcomes in MODULES.items():
        ranking = abundance_ranking(frame, outcomes, predictors)
        scenarios = [(f"exclude_top_{n}", frame.loc[~frame.accepted_key.isin(ranking.index[:n])],
                      [], ranking.index[:n].tolist()) for n in (2,5)]
        if spatial:
            scenarios.append(("broad_space", frame, spatial, []))
        for scenario, subset, adjustments, excluded in scenarios:
            print(module, scenario, flush=True)
            fit = fit_module(subset, outcomes, scales, predictors,
                             adjustments=adjustments, required_environment=predictors)
            results[module+"/"+scenario] = {
                "excluded_species": [str(x) for x in excluded],
                "excluded_observations_with_any_module_trait": int(ranking.loc[excluded].sum()),
                "fit": fit}
            for endpoint, component in fit.get("components", {}).items():
                reference = baseline["models"][module]["components"][endpoint]
                if scenario == "broad_space" and component["n"] != reference["n"]:
                    raise ValueError("Spatial comparison changed endpoint observations")
                for term, effect in component.get("effects", {}).items():
                    if term not in reference.get("effects", {}):
                        continue
                    main = reference["effects"][term]
                    rows.append(dict(module=module, endpoint=endpoint, scenario=scenario,
                                     term=term, n=component["n"], species_n=component["n_species"],
                                     main_change_per_sd=main["raw_units_per_environment_sd"],
                                     diagnostic_change_per_sd=effect["raw_units_per_environment_sd"],
                                     diagnostic_lower=effect["raw_ci95"][0],
                                     diagnostic_upper=effect["raw_ci95"][1]))
    result = {
        "baseline_run": baseline["run_id"], "input_sha256": baseline["input_sha256"],
        "baseline_sha256": hashlib.sha256(baseline_path.read_bytes()).hexdigest(),
        "predictors": predictors, "scales": scales, "spatial_status": spatial_status,
        "scenarios": results,
        "scope": "Fixed diagnostic fits, not alternative preferred main models. "
                 "Abundance exclusions do not rule out influence of every individual species. "
                 "Eight spherical terms describe broad geography, not a complete spatial covariance "
                 "model or phylogenetic correction. Primary Holm family remains unchanged."
    }
    output.mkdir(parents=True, exist_ok=True)
    (output/"environment_robustness.json").write_text(json.dumps(result,indent=2,allow_nan=False)+"\n")
    pd.DataFrame(rows).to_csv(output/"environment_robustness_effects.csv",index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    for key in ("table","baseline","output"):
        parser.add_argument("--"+key,type=Path,required=True)
    args = parser.parse_args()
    execute(args.table,args.baseline,args.output)
