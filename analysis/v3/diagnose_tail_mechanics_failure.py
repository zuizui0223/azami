"""Reproduce the single retained 999-draw mechanics failure without changing acceptance.

Synthetic only. This diagnostic is pinned to the observed mechanics failure
heterogeneous_slopes / outer 7 / 2-degree / bootstrap replicate 562 and records
which module fails plus the original SciPy termination messages. It does not
redraw, alter thresholds, authorize ecology or contribute a replacement draw.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from .dependence_resampling import cohort_partition, crossed_draw, source_partition
from .joint_partial_pooling import PoolingNotEstimable
from .module_ecology import fit_matched_module
from .protected_artifacts import new_json, require
from .simulate_full_family_tail_pilot import _specs, generated_modules
from .workflow import canonical_digest

SCENARIO = "heterogeneous_slopes"
OUTER = 7
GRID = 2
REPLICATE = 562
EXPECTED_SEED = 2026110868


def run(out: Path):
    spec, _ = _specs()
    require(spec["scenarios"].index(SCENARIO) == 2, "Scenario ordering changed")
    seed = int(spec["seed"]) + 10000 * spec["scenarios"].index(SCENARIO) + OUTER
    require(seed == EXPECTED_SEED, "Pinned failure seed changed")
    responses, _, x, g, nuisance, components, lat, lon = generated_modules(SCENARIO, seed)
    source = source_partition(np.arange(len(x)), g, components, lat, lon, grid_degrees=GRID)
    positions, _, _ = cohort_partition(source, np.arange(len(x)))
    indices, copy_taxa = crossed_draw(source, positions, seed=seed, replicate=REPLICATE)

    modules = []
    for module, response in responses.items():
        row = {"module": module, "response_coordinates": int(response.shape[1])}
        try:
            fitted = fit_matched_module(response[indices], x[indices], copy_taxa, nuisance[indices])
            row.update(status="estimated", coefficients_finite=bool(np.isfinite(fitted["coefficients"]).all()),
                       optimizer_attempts=fitted["diagnostics"][0]["within_optimizer_attempts"])
        except PoolingNotEstimable as error:
            row.update(status="not_estimable", error_type=type(error).__name__, reason=str(error),
                       optimizer_attempts=error.optimizer_attempts)
        except (ValueError, RuntimeError, np.linalg.LinAlgError) as error:
            row.update(status="not_estimable", error_type=type(error).__name__, reason=str(error))
        modules.append(row)

    report = {
        "schema_version": 1,
        "status": "PINNED_TAIL_FAILURE_REPRODUCED_FOR_DIAGNOSIS_ONLY",
        "source_mechanics_run_id": 34239896034,
        "mechanics_contract_canonical_sha256": canonical_digest(spec),
        "scenario": SCENARIO,
        "outer_replicate": OUTER,
        "grid_degrees": GRID,
        "bootstrap_replicate": REPLICATE,
        "seed": seed,
        "sampled_observations": int(len(indices)),
        "sampled_taxon_copies": int(len(set(copy_taxa))),
        "modules": modules,
        "replacement_draw_generated": False,
        "acceptance_rule_changed": False,
        "empirical_trait_environment_values_read": 0,
        "ecological_models_executed": 0,
        "ecological_fitting_authorized": False,
    }
    require(sum(row["status"] == "not_estimable" for row in modules) >= 1,
            "Pinned failure no longer reproduces under unchanged acceptance")
    new_json(out, report)
    print(json.dumps(report, allow_nan=False), flush=True)
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    run(args.out)
