#!/usr/bin/env python3
"""Run the unchanged frozen-v2 atlas after restoring all 27 registered endpoints.

This lane changes the trait input only: five fields that were already measured at
historical head level but omitted from the old observation aggregation are restored.
It keeps the full frozen 46,276-observation cohort, all nine frozen environmental
predictors, the original min-5/min-2 rules, 9,999-label-permutation inference and
the original global BH implementation by invoking the existing v2 runner unchanged.

A fresh CHELSA reconstruction is accepted only if the common historical endpoints
reproduce their frozen among-taxon standardized coefficients to <=1e-10 before the
new full-27 run starts.  This is a value-identity gate, not a byte-identity claim.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from analysis.v3 import run_full_vifstep_sensitivity as audit

ROOT = Path(__file__).resolve().parents[2]
TRAIT_SHA = "d775794f2bce2dfd0c1f63c5c8e01778c518f6eeb327bf0d9944045143a02344"
CANONICAL_ENV_SHA = "e242aa7ce69d12b11937c1335e84b9638799c50b42ef36b95725e77190df98e7"
TRAIT_COLUMNS = [
    "obs_id", "taxon_name", "endpoint_id", "module", "analysis_tier",
    "measurement_available", "value",
]


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def as_bool(series: pd.Series) -> pd.Series:
    return series.astype(str).str.strip().str.lower().isin({"true", "1", "yes"})


def build_trait_input(args: argparse.Namespace) -> tuple[Path, pd.DataFrame, pd.DataFrame, dict]:
    if sha256(args.traits_long) != TRAIT_SHA:
        raise SystemExit("Frozen trait-universe SHA mismatch")
    contract = pd.read_csv(args.contract, dtype=str, keep_default_na=False)
    if len(contract) != 27 or contract.endpoint_id.duplicated().any():
        raise SystemExit("Expected exactly 27 unique registered endpoints")
    traits = pd.read_csv(args.traits_long, usecols=TRAIT_COLUMNS, low_memory=False)
    environment = pd.read_csv(args.environment, low_memory=False)
    traits["obs_id"] = traits.obs_id.astype(str)
    traits["taxon_name"] = traits.taxon_name.astype(str)
    environment["obs_id"] = environment.obs_id.astype(str)
    environment["taxon_name"] = environment.taxon_name.astype(str)
    if len(environment) != 46276 or environment.obs_id.nunique() != 46276 or environment.taxon_name.nunique() != 259:
        raise SystemExit("Full v2 environment cohort identity failed")
    for predictor in audit.PREDICTORS:
        if predictor not in environment.columns:
            raise SystemExit(f"Missing predictor: {predictor}")
        environment[predictor] = pd.to_numeric(environment[predictor], errors="coerce")
        if environment[predictor].notna().mean() < .98:
            raise SystemExit(f"Coverage below frozen threshold: {predictor}")
    traits["value"] = pd.to_numeric(traits.value, errors="coerce")
    available = as_bool(traits.measurement_available)
    measured = traits.loc[available & traits.value.notna()].copy()
    measured, recovery = audit.restore_historical_fields(measured, environment, contract, args.recovered_display)
    if measured.duplicated(["obs_id", "endpoint_id"]).any():
        raise SystemExit("Recovered trait input is not unique by obs_id/endpoint_id")
    recovered_counts = {
        endpoint: int(pd.to_numeric(measured.loc[measured.endpoint_id.eq(endpoint), "value"], errors="coerce").notna().sum())
        for endpoint in audit.RECOVERED
    }
    if any(count <= 0 for count in recovered_counts.values()):
        raise SystemExit(f"One or more recovered endpoints remain empty: {recovered_counts}")
    args.work_dir.mkdir(parents=True, exist_ok=True)
    trait_out = args.work_dir / "full27_recovered_traits_long.csv"
    measured[TRAIT_COLUMNS].to_csv(trait_out, index=False)
    return trait_out, contract, environment, {"recovery": recovery, "recovered_counts": recovered_counts}


def value_identity_gate(contract: pd.DataFrame, traits_path: Path, environment: pd.DataFrame,
                        frozen_among: Path, minimum: int) -> dict:
    traits = pd.read_csv(traits_path, low_memory=False)
    traits["value"] = pd.to_numeric(traits.value, errors="coerce")
    units = audit.inferential_units(contract)
    env_taxon = audit.taxon_environment(environment)
    result = audit.verify_frozen_univariate(units, traits, env_taxon, frozen_among, minimum)
    if int(result.get("n_comparable_rows", 0)) < 100:
        raise SystemExit(f"Too few frozen rows reproduced: {result}")
    if float(result.get("maximum_absolute_coefficient_error", np.inf)) > 1e-10:
        raise SystemExit(f"Frozen coefficient reproduction failed: {result}")
    return result


def compare_common(new_path: Path, old_path: Path, scope: str) -> dict:
    new = pd.read_csv(new_path, low_memory=False)
    old = pd.read_csv(old_path, low_memory=False)
    new = new.loc[new.scope.eq(scope)].copy()
    old = old.loc[old.scope.eq(scope)].copy()
    keys = ["unit_id", "predictor", "inferential_unit"]
    both = old.merge(new, on=keys, suffixes=("_old", "_new"))
    both = both.loc[both.status_old.eq("ok") & both.status_new.eq("ok")].copy()
    errors = []
    p_errors = []
    for row in both.to_dict("records"):
        if row["inferential_unit"] == "linear_endpoint":
            candidates = ["beta_std", "beta_std_among"]
            col = next((c for c in candidates if c + "_old" in both.columns), None)
            if col is not None and pd.notna(row.get(col + "_old")) and pd.notna(row.get(col + "_new")):
                errors.append(abs(float(row[col + "_old"]) - float(row[col + "_new"])))
        else:
            for col in ("beta_sine_std", "beta_cosine_std", "beta_sine_std_among", "beta_cosine_std_among"):
                if col + "_old" in both.columns and pd.notna(row.get(col + "_old")) and pd.notna(row.get(col + "_new")):
                    errors.append(abs(float(row[col + "_old"]) - float(row[col + "_new"])))
        if pd.notna(row.get("p_value_old")) and pd.notna(row.get("p_value_new")):
            p_errors.append(abs(float(row["p_value_old"]) - float(row["p_value_new"])))
    return {
        "scope": scope,
        "n_common_ok_rows": int(len(both)),
        "max_effect_error": float(max(errors) if errors else 0.0),
        "max_p_value_error": float(max(p_errors) if p_errors else 0.0),
        "note": "q values are not required to match because restoring five endpoints expands the BH family",
    }


def signal_rows(path: Path, scope: str) -> list[dict]:
    frame = pd.read_csv(path, low_memory=False)
    q = pd.to_numeric(frame.get("q_fdr_bh_global_family"), errors="coerce")
    hit = frame.loc[frame.scope.eq(scope) & frame.status.eq("ok") & q.lt(.05)].copy()
    columns = [c for c in ["scope", "unit_id", "module", "inferential_unit", "predictor", "n_taxa",
                              "beta_std", "beta_std_among", "p_value", "beta_sine_std_among",
                              "beta_cosine_std_among", "effect_magnitude", "effect_direction_degrees",
                              "q_fdr_bh_global_family"] if c in hit.columns]
    return hit[columns].sort_values(["q_fdr_bh_global_family", "unit_id", "predictor"]).to_dict("records")


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--contract", type=Path, required=True)
    p.add_argument("--analysis-contract", type=Path, required=True)
    p.add_argument("--traits-long", type=Path, required=True)
    p.add_argument("--environment", type=Path, required=True)
    p.add_argument("--recovered-display", type=Path, required=True)
    p.add_argument("--frozen-among", type=Path, required=True)
    p.add_argument("--frozen-within", type=Path, required=True)
    p.add_argument("--work-dir", type=Path, required=True)
    p.add_argument("--out-dir", type=Path, required=True)
    args = p.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    trait_out, contract, environment, input_report = build_trait_input(args)
    gate = value_identity_gate(contract, trait_out, environment, args.frozen_among, 5)
    actual_env_sha = sha256(args.environment)

    subprocess.run([
        sys.executable, str(ROOT / "analysis/run_geb_v2_full27_environment_atlas.py"),
        "--contract", str(args.contract),
        "--analysis-contract", str(args.analysis_contract),
        "--traits-long", str(trait_out),
        "--environment", str(args.environment),
        "--out-dir", str(args.out_dir),
    ], cwd=ROOT, check=True)

    new_among = args.out_dir / "v2_full27_environment_among.csv"
    new_within = args.out_dir / "v2_full27_environment_within.csv"
    common = [
        compare_common(new_among, args.frozen_among, "among_taxon_min5"),
        compare_common(new_among, args.frozen_among, "among_taxon_min2"),
    ]
    # Within uses the same keys but may have a single scope; compare all shared scope names.
    old_within = pd.read_csv(args.frozen_within, low_memory=False)
    new_within_frame = pd.read_csv(new_within, low_memory=False)
    for scope in sorted(set(old_within.scope).intersection(set(new_within_frame.scope))):
        common.append(compare_common(new_within, args.frozen_within, scope))
    if any(row["max_effect_error"] > 1e-10 or row["max_p_value_error"] > 1e-12 for row in common):
        raise SystemExit(f"Unchanged-v2 common rows failed reproduction: {common}")

    report = {
        "status": "FULL27_EXACT_V2_FLOW_COMPLETE",
        "meaning": "same v2 analysis flow; only restores five historically measured registered endpoints",
        "cohort": {"observations": 46276, "taxa": 259},
        "registered_endpoints": 27,
        "environment_predictors": audit.PREDICTORS,
        "environment_identity": {
            "canonical_sha256": CANONICAL_ENV_SHA,
            "actual_sha256": actual_env_sha,
            "canonical_byte_identity": actual_env_sha == CANONICAL_ENV_SHA,
            "value_identity_gate": gate,
        },
        "input_extension": input_report,
        "common_v2_reproduction": common,
        "fdr_signals": {
            "among_taxon_min5": signal_rows(new_among, "among_taxon_min5"),
            "among_taxon_min2": signal_rows(new_among, "among_taxon_min2"),
        },
        "claim_boundary": [
            "This is not a new model family: it invokes the existing v2 atlas runner unchanged.",
            "Restoring five endpoints expands the BH family, so q values for old endpoints may change even when their beta and raw p values are identical.",
            "The four colour-composition fractions are a closed composition and are not four independent biological discoveries.",
            "Environmental coefficients remain marginal associations along correlated observational gradients.",
        ],
    }
    (args.out_dir / "full27_exact_v2_flow_report.json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    print("FULL27_EXACT_V2_FLOW=PASS")
    print("FULL27_EXACT_V2_COMPACT_JSON=" + json.dumps({
        "environment_value_gate": gate,
        "common_reproduction": common,
        "among_min5_fdr_n": len(report["fdr_signals"]["among_taxon_min5"]),
        "among_min5_fdr": report["fdr_signals"]["among_taxon_min5"],
    }, separators=(",", ":"), allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
