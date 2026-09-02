#!/usr/bin/env python3
"""Run the trait–environment association atlas for full and precise-coordinate cohorts.

The two cohorts are intentionally produced side by side:
- all_coordinate_eligible: all rows with complete predeclared environmental predictors.
- positional_accuracy_le_10km: rows with iNaturalist positional accuracy <=10 km.

The script invokes the report-safe hardened atlas unchanged for each cohort,
preventing the precision sensitivity comparison from using different modelling
rules while retaining a valid report/provenance output.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parent
HARDENED = ROOT / "46_run_trait_environment_association_atlas_reportfix.py"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run full and <=10 km coordinate-precision association atlases.")
    parser.add_argument("--environment-long", required=True)
    parser.add_argument("--environment-provenance", required=True)
    parser.add_argument("--plan", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--n-folds", type=int, default=5)
    parser.add_argument("--seed", type=int, default=20260706)
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def run_atlas(environment_long: Path, provenance: Path, plan: Path, out_dir: Path, n_folds: int, seed: int) -> None:
    command = [
        sys.executable,
        str(HARDENED),
        "--environment-long", str(environment_long),
        "--environment-provenance", str(provenance),
        "--plan", str(plan),
        "--out-dir", str(out_dir),
        "--n-folds", str(n_folds),
        "--seed", str(seed),
    ]
    subprocess.run(command, check=True)


def read_summary(path: Path) -> dict:
    return json.loads((path / "trait_environment_association_provenance.json").read_text(encoding="utf-8"))


def main() -> None:
    args = parse_args()
    if args.n_folds < 2:
        raise ValueError("--n-folds must be at least 2")
    input_path = Path(args.environment_long)
    provenance_path = Path(args.environment_provenance)
    plan_path = Path(args.plan)
    for path in [input_path, provenance_path, plan_path, HARDENED]:
        if not path.is_file():
            raise FileNotFoundError(path)

    frame = pd.read_csv(input_path, dtype=str, keep_default_na=False)
    required = {"obs_id", "trait_id", "coordinate_precision_tier"}
    missing = sorted(required.difference(frame.columns))
    if missing:
        raise ValueError(f"Environment-long table is missing {missing}")
    if frame.duplicated(["obs_id", "trait_id"]).any():
        raise ValueError("Environment-long table has duplicate obs_id/trait_id rows")

    output = Path(args.out_dir)
    output.mkdir(parents=True, exist_ok=True)
    cohorts = {
        "all_coordinate_eligible": frame,
        "positional_accuracy_le_10km": frame.loc[
            frame["coordinate_precision_tier"].isin(["high_le_1km", "moderate_1_to_10km"])
        ].copy(),
    }
    summaries: dict[str, dict] = {}
    for cohort_id, cohort in cohorts.items():
        if cohort.empty:
            raise ValueError(f"Cohort is empty: {cohort_id}")
        cohort_input = output / f"{cohort_id}_environment_long.csv"
        cohort.to_csv(cohort_input, index=False, encoding="utf-8-sig")
        cohort_output = output / cohort_id
        run_atlas(cohort_input, provenance_path, plan_path, cohort_output, args.n_folds, args.seed)
        summary = read_summary(cohort_output)
        summary["n_trait_rows_input_to_cohort"] = int(len(cohort))
        summary["n_observations_input_to_cohort"] = int(cohort["obs_id"].nunique())
        summaries[cohort_id] = summary

    comparison = pd.DataFrame([
        {
            "cohort": cohort_id,
            "n_observations": details["n_observations_input_to_cohort"],
            "n_trait_rows": details["n_trait_rows_input_to_cohort"],
            "n_state_outcomes_screened": details["n_trait_state_outcomes_screened"],
            "n_state_outcomes_eligible": details["n_trait_state_outcomes_eligible"],
            "n_coefficient_rows": details["n_coefficient_rows"],
            "n_validation_rows": details["n_validation_rows"],
        }
        for cohort_id, details in summaries.items()
    ])
    comparison.to_csv(output / "coordinate_precision_sensitivity_cohort_summary.csv", index=False, encoding="utf-8-sig")
    provenance = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "semantic_status": "Coordinate-precision sensitivity wrapper for the fully automated trait–environment association atlas. It changes only the observation cohort, not the trait model, predictors, outcome thresholds, or grouped-validation design.",
        "cohorts": {
            "all_coordinate_eligible": "All rows in the supplied environment-long table.",
            "positional_accuracy_le_10km": "Rows whose coordinate_precision_tier is high_le_1km or moderate_1_to_10km.",
        },
        "input_sha256": {
            "environment_long": sha256(input_path),
            "environment_provenance": sha256(provenance_path),
            "analysis_plan": sha256(plan_path),
        },
        "n_folds": int(args.n_folds),
        "seed": int(args.seed),
        "cohort_summaries": summaries,
    }
    (output / "coordinate_precision_sensitivity_provenance.json").write_text(
        json.dumps(provenance, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(json.dumps({"cohorts": summaries}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
