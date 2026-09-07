"""Build and diagnose a phenotype-blind v3 abiotic exposure matrix.

This module must not read capitulum trait files. It samples CHELSA 2.1 long-term
monthly climatologies at observation coordinates using the observation month,
then reports coverage and environmental redundancy before any trait join.
"""
from __future__ import annotations

import argparse
from dataclasses import dataclass
import hashlib
import json
import math
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
DEFAULT_CONTRACT = ROOT / "analysis" / "v3" / "environment_exposure_contract.json"


@dataclass(frozen=True)
class RasterCandidate:
    output_id: str
    variable: str
    construct: str
    unit: str


def load_contract(path: Path = DEFAULT_CONTRACT) -> dict:
    data = json.loads(path.read_text(encoding="utf-8"))
    if data.get("status") != "environment_only_design_before_trait_join":
        raise ValueError("Environment exposure contract is not in environment-only design state")
    if not data.get("selection_is_phenotype_blind"):
        raise ValueError("Environment contract must prohibit phenotype-guided selection")
    return data


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def month_url(contract: dict, variable: str, month: int) -> str:
    if not 1 <= int(month) <= 12:
        raise ValueError(f"Invalid observation month: {month}")
    return contract["monthly_url_template"].format(variable=variable, month=int(month))


def validate_observation_frame(frame: pd.DataFrame, require_native: bool) -> pd.DataFrame:
    required = {"obs_id", "latitude", "longitude", "observation_month"}
    missing = sorted(required - set(frame.columns))
    if missing:
        raise ValueError("Missing observation columns: " + ", ".join(missing))
    out = frame.copy()
    if out["obs_id"].isna().any() or out["obs_id"].astype(str).duplicated().any():
        raise ValueError("obs_id must be present and unique")
    out["latitude"] = pd.to_numeric(out["latitude"], errors="coerce")
    out["longitude"] = pd.to_numeric(out["longitude"], errors="coerce")
    out["observation_month"] = pd.to_numeric(out["observation_month"], errors="coerce")
    invalid = (
        out["latitude"].isna()
        | out["longitude"].isna()
        | ~out["latitude"].between(-90, 90)
        | ~out["longitude"].between(-180, 180)
        | ~out["observation_month"].between(1, 12)
    )
    if invalid.any():
        raise ValueError(f"Invalid coordinate/month rows: {int(invalid.sum())}")
    out["observation_month"] = out["observation_month"].astype(int)
    if require_native:
        if "native_range_status" not in out.columns:
            raise ValueError("native_range_status is required for the primary environment cohort")
        out = out.loc[out["native_range_status"].astype(str).eq("native")].copy()
        if out.empty:
            raise ValueError("No native-range observations remain")
    return out.reset_index(drop=True)


def _apply_raster_scale(value: float, scale: float, offset: float) -> float:
    return float(value) * float(scale) + float(offset)


def sample_monthly_candidate(
    frame: pd.DataFrame,
    candidate: RasterCandidate,
    contract: dict,
    url_builder=month_url,
) -> pd.Series:
    """Sample one CHELSA variable, opening only months represented in frame."""
    try:
        import rasterio
        from rasterio.warp import transform
    except ImportError as exc:  # pragma: no cover - exercised in full environment only
        raise RuntimeError("rasterio is required for CHELSA sampling; install project [full]") from exc

    result = pd.Series(np.nan, index=frame.index, dtype="float64", name=candidate.output_id)
    env_options = {
        "GDAL_DISABLE_READDIR_ON_OPEN": "EMPTY_DIR",
        "CPL_VSIL_CURL_ALLOWED_EXTENSIONS": ".tif,.TIF",
        "GDAL_HTTP_MULTIRANGE": "YES",
        "VSI_CACHE": "TRUE",
        "VSI_CACHE_SIZE": "50000000",
    }
    with rasterio.Env(**env_options):
        for month, group in frame.groupby("observation_month", sort=True):
            source = url_builder(contract, candidate.variable, int(month))
            with rasterio.open(source) as dataset:
                lons = group["longitude"].astype(float).tolist()
                lats = group["latitude"].astype(float).tolist()
                if dataset.crs is None:
                    raise ValueError(f"Raster has no CRS: {source}")
                if str(dataset.crs).upper() not in {"EPSG:4326", "OGC:CRS84"}:
                    xs, ys = transform("EPSG:4326", dataset.crs, lons, lats)
                else:
                    xs, ys = lons, lats
                scale = dataset.scales[0] if dataset.scales else 1.0
                offset = dataset.offsets[0] if dataset.offsets else 0.0
                sampled = dataset.sample(zip(xs, ys), indexes=1, masked=True)
                values = []
                for item in sampled:
                    scalar = item[0]
                    if np.ma.is_masked(scalar):
                        values.append(np.nan)
                    else:
                        value = float(scalar)
                        values.append(_apply_raster_scale(value, scale, offset) if math.isfinite(value) else np.nan)
                result.loc[group.index] = values
    return result


def candidate_specs(contract: dict) -> list[RasterCandidate]:
    specs = []
    for row in contract["monthly_candidates"]:
        specs.append(
            RasterCandidate(
                output_id=str(row["id"]),
                variable=str(row["chelsa_variable"]),
                construct=str(row["construct"]),
                unit=str(row["unit"]),
            )
        )
    return specs


def coverage_table(environment: pd.DataFrame, variables: Iterable[str]) -> pd.DataFrame:
    n = len(environment)
    rows = []
    for variable in variables:
        values = pd.to_numeric(environment[variable], errors="coerce")
        finite = np.isfinite(values.to_numpy(dtype=float))
        rows.append(
            {
                "variable": variable,
                "n_total": n,
                "n_finite": int(finite.sum()),
                "coverage": float(finite.mean()) if n else np.nan,
                "mean": float(values[finite].mean()) if finite.any() else np.nan,
                "sd": float(values[finite].std(ddof=1)) if finite.sum() > 1 else np.nan,
            }
        )
    return pd.DataFrame(rows)


def correlation_long(environment: pd.DataFrame, variables: list[str], method: str) -> pd.DataFrame:
    corr = environment[variables].corr(method=method, min_periods=3)
    rows = []
    for i, left in enumerate(variables):
        for right in variables[i + 1 :]:
            rows.append({"method": method, "variable_a": left, "variable_b": right, "correlation": corr.loc[left, right]})
    return pd.DataFrame(rows)


def complete_standardized(environment: pd.DataFrame, variables: list[str]) -> tuple[pd.DataFrame, list[str]]:
    data = environment[variables].apply(pd.to_numeric, errors="coerce").replace([np.inf, -np.inf], np.nan).dropna()
    retained = [column for column in variables if data[column].nunique(dropna=True) > 1]
    data = data[retained]
    if data.empty or not retained:
        return data, retained
    sd = data.std(ddof=0)
    retained = [column for column in retained if float(sd[column]) > 0]
    data = data[retained]
    if not retained:
        return data, retained
    return (data - data.mean()) / data.std(ddof=0), retained


def vif_table(environment: pd.DataFrame, variables: list[str]) -> tuple[pd.DataFrame, dict]:
    standardized, retained = complete_standardized(environment, variables)
    report = {
        "n_complete_rows": int(len(standardized)),
        "n_candidate_variables": len(variables),
        "n_nonconstant_complete_variables": len(retained),
        "matrix_rank": 0,
        "condition_number": None,
    }
    if len(retained) == 0 or len(standardized) < 3:
        return pd.DataFrame(columns=["variable", "vif"]), report
    matrix = standardized.to_numpy(dtype=float)
    report["matrix_rank"] = int(np.linalg.matrix_rank(matrix))
    singular = np.linalg.svd(matrix, compute_uv=False)
    report["condition_number"] = float(singular[0] / singular[-1]) if singular[-1] > np.finfo(float).eps else float("inf")
    rows = []
    for idx, variable in enumerate(retained):
        y = matrix[:, idx]
        others = np.delete(matrix, idx, axis=1)
        if others.shape[1] == 0:
            vif = 1.0
        else:
            design = np.column_stack([np.ones(len(others)), others])
            coef, *_ = np.linalg.lstsq(design, y, rcond=None)
            fitted = design @ coef
            ss_res = float(np.sum((y - fitted) ** 2))
            ss_tot = float(np.sum((y - y.mean()) ** 2))
            r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else 1.0
            vif = float("inf") if 1.0 - r2 <= 1e-12 else 1.0 / (1.0 - r2)
        rows.append({"variable": variable, "vif": vif})
    return pd.DataFrame(rows), report


def redundancy_components(correlations: pd.DataFrame, variables: list[str], threshold: float) -> list[list[str]]:
    parent = {v: v for v in variables}

    def find(x: str) -> str:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a: str, b: str) -> None:
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[rb] = ra

    pearson = correlations.loc[correlations["method"].eq("pearson")]
    for row in pearson.itertuples(index=False):
        if pd.notna(row.correlation) and abs(float(row.correlation)) >= threshold:
            union(str(row.variable_a), str(row.variable_b))
    groups: dict[str, list[str]] = {}
    for variable in variables:
        groups.setdefault(find(variable), []).append(variable)
    return sorted((sorted(values) for values in groups.values()), key=lambda x: (len(x), x), reverse=True)


def diagnose_environment(environment: pd.DataFrame, variables: list[str], threshold: float) -> tuple[dict, dict[str, pd.DataFrame]]:
    coverage = coverage_table(environment, variables)
    pearson = correlation_long(environment, variables, "pearson")
    spearman = correlation_long(environment, variables, "spearman")
    correlations = pd.concat([pearson, spearman], ignore_index=True)
    vif, matrix = vif_table(environment, variables)
    clusters = redundancy_components(correlations, variables, threshold)
    report = {
        "status": "ENVIRONMENT_ONLY_DIAGNOSTICS_COMPLETE",
        "n_rows": int(len(environment)),
        "variables": variables,
        "matrix": matrix,
        "absolute_correlation_redundancy_threshold": float(threshold),
        "redundancy_components": clusters,
        "trait_columns_read": 0,
        "representation_frozen": False,
        "note": "Diagnostics flag redundancy only. Final representation must be chosen from biological meaning and environment-only diagnostics before trait fitting.",
    }
    return report, {"coverage": coverage, "correlations": correlations, "vif": vif}


def build_matrix(observations: pd.DataFrame, contract: dict, require_native: bool, variables: list[str] | None = None) -> tuple[pd.DataFrame, list[str]]:
    frame = validate_observation_frame(observations, require_native=require_native)
    specs = candidate_specs(contract)
    if variables:
        wanted = set(variables)
        specs = [spec for spec in specs if spec.output_id in wanted]
        unknown = wanted - {spec.output_id for spec in specs}
        if unknown:
            raise ValueError("Unknown environment variables: " + ", ".join(sorted(unknown)))
    output = frame.copy()
    for spec in specs:
        output[spec.output_id] = sample_monthly_candidate(output, spec, contract)
    return output, [spec.output_id for spec in specs]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--observations", type=Path, required=True, help="CSV with obs_id, latitude, longitude, observation_month and optionally native_range_status")
    parser.add_argument("--contract", type=Path, default=DEFAULT_CONTRACT)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--allow-nonnative-pilot", action="store_true", help="Engineering pilot only; cannot freeze final environment representation")
    parser.add_argument("--variables", nargs="*", default=None)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    contract = load_contract(args.contract)
    source_sha = sha256_file(args.observations)
    observations = pd.read_csv(args.observations, low_memory=False)
    matrix, variables = build_matrix(observations, contract, require_native=not args.allow_nonnative_pilot, variables=args.variables)
    threshold = float(contract["redundancy_rule"]["default_absolute_correlation_flag"])
    report, tables = diagnose_environment(matrix, variables, threshold)
    report.update(
        {
            "contract_id": contract["contract_id"],
            "source_observations_sha256": source_sha,
            "native_only": not args.allow_nonnative_pilot,
            "engineering_pilot_only": bool(args.allow_nonnative_pilot),
            "representation_may_be_frozen_from_this_run": not args.allow_nonnative_pilot,
        }
    )
    args.out_dir.mkdir(parents=True, exist_ok=False)
    matrix.to_csv(args.out_dir / "environment_candidate_matrix.csv", index=False)
    tables["coverage"].to_csv(args.out_dir / "environment_coverage.csv", index=False)
    tables["correlations"].to_csv(args.out_dir / "environment_correlations.csv", index=False)
    tables["vif"].to_csv(args.out_dir / "environment_vif.csv", index=False)
    (args.out_dir / "environment_diagnostics.json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({"status": report["status"], "n_rows": report["n_rows"], "variables": variables, "pilot": report["engineering_pilot_only"]}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
