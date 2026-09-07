"""Build and diagnose a phenotype-blind v3 abiotic exposure matrix.

This module must not read capitulum trait files. It samples CHELSA 2.1 long-term
monthly climatologies at observation coordinates using the observation month and
can sample broader annual/seasonal representations at the same coordinates. It
then reports coverage and environmental redundancy before any trait join. If the
source cohort carries ``equal_taxon_weight``, the primary correlation/VIF matrix
uses those weights so photo-rich taxa cannot dominate environmental selection.
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
    temporal_mode: str = "monthly"
    url: str | None = None


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
    if "equal_taxon_weight" in out.columns:
        out["equal_taxon_weight"] = pd.to_numeric(out["equal_taxon_weight"], errors="coerce")
        if out["equal_taxon_weight"].isna().any() or (out["equal_taxon_weight"] <= 0).any():
            raise ValueError("equal_taxon_weight must be finite and positive")
    if require_native:
        if "native_range_status" not in out.columns:
            raise ValueError("native_range_status is required for the primary environment cohort")
        out = out.loc[out["native_range_status"].astype(str).eq("native")].copy()
        if out.empty:
            raise ValueError("No native-range observations remain")
    return out.reset_index(drop=True)


def _apply_raster_scale(value: float, scale: float, offset: float) -> float:
    return float(value) * float(scale) + float(offset)


def _rasterio_modules():
    try:
        import rasterio
        from rasterio.warp import transform
    except ImportError as exc:  # pragma: no cover - full environment only
        raise RuntimeError("rasterio is required for CHELSA sampling; install project [full]") from exc
    return rasterio, transform


def _sample_dataset(dataset, frame: pd.DataFrame, transform) -> list[float]:
    lons = frame["longitude"].astype(float).tolist()
    lats = frame["latitude"].astype(float).tolist()
    if dataset.crs is None:
        raise ValueError("Raster has no CRS")
    if str(dataset.crs).upper() not in {"EPSG:4326", "OGC:CRS84"}:
        xs, ys = transform("EPSG:4326", dataset.crs, lons, lats)
    else:
        xs, ys = lons, lats
    scale = dataset.scales[0] if dataset.scales else 1.0
    offset = dataset.offsets[0] if dataset.offsets else 0.0
    values: list[float] = []
    for item in dataset.sample(zip(xs, ys), indexes=1, masked=True):
        scalar = item[0]
        if np.ma.is_masked(scalar):
            values.append(np.nan)
            continue
        value = float(scalar)
        values.append(_apply_raster_scale(value, scale, offset) if math.isfinite(value) else np.nan)
    return values


def _raster_env_options() -> dict[str, str]:
    return {
        "GDAL_DISABLE_READDIR_ON_OPEN": "EMPTY_DIR",
        "CPL_VSIL_CURL_ALLOWED_EXTENSIONS": ".tif,.TIF",
        "GDAL_HTTP_MULTIRANGE": "YES",
        "VSI_CACHE": "TRUE",
        "VSI_CACHE_SIZE": "50000000",
    }


def sample_monthly_candidate(
    frame: pd.DataFrame,
    candidate: RasterCandidate,
    contract: dict,
    url_builder=month_url,
) -> pd.Series:
    """Sample one CHELSA monthly candidate, opening only represented months."""
    rasterio, transform = _rasterio_modules()
    result = pd.Series(np.nan, index=frame.index, dtype="float64", name=candidate.output_id)
    with rasterio.Env(**_raster_env_options()):
        for month, group in frame.groupby("observation_month", sort=True):
            source = url_builder(contract, candidate.variable, int(month))
            with rasterio.open(source) as dataset:
                result.loc[group.index] = _sample_dataset(dataset, group, transform)
    return result


def sample_static_candidate(frame: pd.DataFrame, candidate: RasterCandidate) -> pd.Series:
    """Sample one annual/seasonal CHELSA candidate at all observation points."""
    if not candidate.url:
        raise ValueError(f"Static candidate {candidate.output_id} has no URL")
    rasterio, transform = _rasterio_modules()
    with rasterio.Env(**_raster_env_options()):
        with rasterio.open(candidate.url) as dataset:
            values = _sample_dataset(dataset, frame, transform)
    return pd.Series(values, index=frame.index, dtype="float64", name=candidate.output_id)


def candidate_specs(contract: dict) -> list[RasterCandidate]:
    specs: list[RasterCandidate] = []
    for row in contract["monthly_candidates"]:
        specs.append(
            RasterCandidate(
                output_id=str(row["id"]),
                variable=str(row["chelsa_variable"]),
                construct=str(row["construct"]),
                unit=str(row["unit"]),
                temporal_mode="monthly",
            )
        )
    for row in contract.get("broader_climate_representations", []):
        specs.append(
            RasterCandidate(
                output_id=str(row["id"]),
                variable=str(row["id"]),
                construct=str(row["construct"]),
                unit=str(row.get("unit", "CHELSA native scaled unit")),
                temporal_mode="static",
                url=str(row["url"]),
            )
        )
    ids = [spec.output_id for spec in specs]
    if len(ids) != len(set(ids)):
        raise ValueError("Environment contract contains duplicate candidate IDs")
    return specs


def _numeric_weights(environment: pd.DataFrame, weight_column: str | None) -> pd.Series:
    if weight_column is None:
        return pd.Series(np.ones(len(environment), dtype=float), index=environment.index)
    if weight_column not in environment.columns:
        raise ValueError(f"Missing weight column: {weight_column}")
    weights = pd.to_numeric(environment[weight_column], errors="coerce")
    if weights.isna().any() or (weights <= 0).any() or not np.isfinite(weights.to_numpy(dtype=float)).all():
        raise ValueError("Diagnostic weights must be finite and positive")
    return weights.astype(float)


def coverage_table(environment: pd.DataFrame, variables: Iterable[str], weight_column: str | None = None) -> pd.DataFrame:
    n = len(environment)
    weights = _numeric_weights(environment, weight_column)
    total_weight = float(weights.sum())
    rows = []
    for variable in variables:
        values = pd.to_numeric(environment[variable], errors="coerce")
        finite = np.isfinite(values.to_numpy(dtype=float))
        finite_weights = weights.loc[finite]
        weighted_coverage = float(finite_weights.sum() / total_weight) if total_weight > 0 else np.nan
        rows.append(
            {
                "variable": variable,
                "n_total": n,
                "n_finite": int(finite.sum()),
                "coverage": float(finite.mean()) if n else np.nan,
                "weighted_coverage": weighted_coverage,
                "mean": float(values[finite].mean()) if finite.any() else np.nan,
                "sd": float(values[finite].std(ddof=1)) if finite.sum() > 1 else np.nan,
            }
        )
    return pd.DataFrame(rows)


def _weighted_corr_pair(x: pd.Series, y: pd.Series, weights: pd.Series, rank: bool) -> float:
    frame = pd.DataFrame({"x": pd.to_numeric(x, errors="coerce"), "y": pd.to_numeric(y, errors="coerce"), "w": weights})
    frame = frame.replace([np.inf, -np.inf], np.nan).dropna()
    if len(frame) < 3 or frame["x"].nunique() < 2 or frame["y"].nunique() < 2:
        return np.nan
    if rank:
        frame["x"] = frame["x"].rank(method="average")
        frame["y"] = frame["y"].rank(method="average")
    w = frame["w"].to_numpy(dtype=float)
    xv = frame["x"].to_numpy(dtype=float)
    yv = frame["y"].to_numpy(dtype=float)
    wsum = float(w.sum())
    mx = float(np.sum(w * xv) / wsum)
    my = float(np.sum(w * yv) / wsum)
    dx = xv - mx
    dy = yv - my
    cov = float(np.sum(w * dx * dy) / wsum)
    vx = float(np.sum(w * dx * dx) / wsum)
    vy = float(np.sum(w * dy * dy) / wsum)
    if vx <= 0 or vy <= 0:
        return np.nan
    return cov / math.sqrt(vx * vy)


def correlation_long(
    environment: pd.DataFrame,
    variables: list[str],
    method: str,
    weight_column: str | None = None,
) -> pd.DataFrame:
    if method not in {"pearson", "spearman"}:
        raise ValueError("Correlation method must be pearson or spearman")
    rows = []
    if weight_column is None:
        corr = environment[variables].corr(method=method, min_periods=3)
        for i, left in enumerate(variables):
            for right in variables[i + 1 :]:
                rows.append({"method": method, "variable_a": left, "variable_b": right, "correlation": corr.loc[left, right]})
        return pd.DataFrame(rows)
    weights = _numeric_weights(environment, weight_column)
    label = method + "_equal_taxon_weight" if weight_column == "equal_taxon_weight" else method + "_weighted"
    for i, left in enumerate(variables):
        for right in variables[i + 1 :]:
            rows.append(
                {
                    "method": label,
                    "variable_a": left,
                    "variable_b": right,
                    "correlation": _weighted_corr_pair(environment[left], environment[right], weights, rank=method == "spearman"),
                }
            )
    return pd.DataFrame(rows)


def complete_standardized(
    environment: pd.DataFrame,
    variables: list[str],
    weight_column: str | None = None,
) -> tuple[pd.DataFrame, pd.Series, list[str]]:
    columns = list(variables) + ([weight_column] if weight_column else [])
    data = environment[columns].copy()
    for variable in variables:
        data[variable] = pd.to_numeric(data[variable], errors="coerce")
    if weight_column:
        data[weight_column] = pd.to_numeric(data[weight_column], errors="coerce")
    data = data.replace([np.inf, -np.inf], np.nan).dropna()
    weights = _numeric_weights(data, weight_column)
    retained = [column for column in variables if data[column].nunique(dropna=True) > 1]
    if data.empty or not retained:
        return pd.DataFrame(index=data.index), weights, retained
    standardized = pd.DataFrame(index=data.index)
    final: list[str] = []
    w = weights.to_numpy(dtype=float)
    wsum = float(w.sum())
    for column in retained:
        values = data[column].to_numpy(dtype=float)
        mean = float(np.sum(w * values) / wsum)
        variance = float(np.sum(w * (values - mean) ** 2) / wsum)
        if variance <= 0:
            continue
        standardized[column] = (values - mean) / math.sqrt(variance)
        final.append(column)
    return standardized, weights.loc[standardized.index], final


def vif_table(
    environment: pd.DataFrame,
    variables: list[str],
    weight_column: str | None = None,
) -> tuple[pd.DataFrame, dict]:
    standardized, weights, retained = complete_standardized(environment, variables, weight_column)
    report = {
        "n_complete_rows": int(len(standardized)),
        "n_candidate_variables": len(variables),
        "n_nonconstant_complete_variables": len(retained),
        "matrix_rank": 0,
        "condition_number": None,
        "weight_column": weight_column,
        "complete_row_total_weight": float(weights.sum()) if len(weights) else 0.0,
    }
    if len(retained) == 0 or len(standardized) < 3:
        return pd.DataFrame(columns=["variable", "vif"]), report
    matrix = standardized[retained].to_numpy(dtype=float)
    w = weights.to_numpy(dtype=float)
    sqrt_w = np.sqrt(w / np.mean(w))
    weighted_matrix = matrix * sqrt_w[:, None]
    report["matrix_rank"] = int(np.linalg.matrix_rank(weighted_matrix))
    singular = np.linalg.svd(weighted_matrix, compute_uv=False)
    report["condition_number"] = (
        float(singular[0] / singular[-1]) if singular[-1] > np.finfo(float).eps else "infinite"
    )
    rows = []
    for idx, variable in enumerate(retained):
        y = matrix[:, idx]
        others = np.delete(matrix, idx, axis=1)
        if others.shape[1] == 0:
            vif: float | str = 1.0
        else:
            design = np.column_stack([np.ones(len(others)), others])
            weighted_design = design * sqrt_w[:, None]
            weighted_y = y * sqrt_w
            coef, *_ = np.linalg.lstsq(weighted_design, weighted_y, rcond=None)
            fitted = design @ coef
            mean_y = float(np.sum(w * y) / np.sum(w))
            ss_res = float(np.sum(w * (y - fitted) ** 2))
            ss_tot = float(np.sum(w * (y - mean_y) ** 2))
            r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else 1.0
            vif = "infinite" if 1.0 - r2 <= 1e-12 else 1.0 / (1.0 - r2)
        rows.append({"variable": variable, "vif": vif})
    return pd.DataFrame(rows), report


def redundancy_components(
    correlations: pd.DataFrame,
    variables: list[str],
    threshold: float,
    primary_method: str = "pearson",
) -> list[list[str]]:
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

    primary = correlations.loc[correlations["method"].eq(primary_method)]
    for row in primary.itertuples(index=False):
        if pd.notna(row.correlation) and abs(float(row.correlation)) >= threshold:
            union(str(row.variable_a), str(row.variable_b))
    groups: dict[str, list[str]] = {}
    for variable in variables:
        groups.setdefault(find(variable), []).append(variable)
    return sorted((sorted(values) for values in groups.values()), key=lambda x: (len(x), x), reverse=True)


def diagnose_environment(
    environment: pd.DataFrame,
    variables: list[str],
    threshold: float,
    weight_column: str | None = None,
) -> tuple[dict, dict[str, pd.DataFrame]]:
    if weight_column is None and "equal_taxon_weight" in environment.columns:
        weight_column = "equal_taxon_weight"
    coverage = coverage_table(environment, variables, weight_column)
    unweighted_pearson = correlation_long(environment, variables, "pearson")
    unweighted_spearman = correlation_long(environment, variables, "spearman")
    parts = [unweighted_pearson, unweighted_spearman]
    primary_method = "pearson"
    if weight_column:
        parts.append(correlation_long(environment, variables, "pearson", weight_column))
        parts.append(correlation_long(environment, variables, "spearman", weight_column))
        primary_method = "pearson_equal_taxon_weight" if weight_column == "equal_taxon_weight" else "pearson_weighted"
    correlations = pd.concat(parts, ignore_index=True)
    vif, matrix = vif_table(environment, variables, weight_column)
    clusters = redundancy_components(correlations, variables, threshold, primary_method=primary_method)
    report = {
        "status": "ENVIRONMENT_ONLY_DIAGNOSTICS_COMPLETE",
        "n_rows": int(len(environment)),
        "variables": variables,
        "matrix": matrix,
        "weight_column": weight_column,
        "primary_redundancy_correlation_method": primary_method,
        "absolute_correlation_redundancy_threshold": float(threshold),
        "redundancy_components": clusters,
        "trait_columns_read": 0,
        "representation_frozen": False,
        "note": "Diagnostics flag redundancy only. Final representation must be chosen from biological meaning and environment-only diagnostics before trait fitting.",
    }
    return report, {"coverage": coverage, "correlations": correlations, "vif": vif}


def build_matrix(
    observations: pd.DataFrame,
    contract: dict,
    require_native: bool,
    variables: list[str] | None = None,
) -> tuple[pd.DataFrame, list[str]]:
    frame = validate_observation_frame(observations, require_native=require_native)
    specs = candidate_specs(contract)
    if variables:
        wanted = set(variables)
        known = {spec.output_id for spec in specs}
        unknown = wanted - known
        if unknown:
            raise ValueError("Unknown environment variables: " + ", ".join(sorted(unknown)))
        specs = [spec for spec in specs if spec.output_id in wanted]
    output = frame.copy()
    for spec in specs:
        if spec.temporal_mode == "monthly":
            output[spec.output_id] = sample_monthly_candidate(output, spec, contract)
        elif spec.temporal_mode == "static":
            output[spec.output_id] = sample_static_candidate(output, spec)
        else:
            raise ValueError(f"Unknown temporal mode: {spec.temporal_mode}")
    return output, [spec.output_id for spec in specs]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--observations",
        type=Path,
        required=True,
        help="CSV with obs_id, latitude, longitude, observation_month and optionally native_range_status/equal_taxon_weight",
    )
    parser.add_argument("--contract", type=Path, default=DEFAULT_CONTRACT)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument(
        "--allow-nonnative-pilot",
        action="store_true",
        help="Engineering pilot only; cannot freeze final environment representation",
    )
    parser.add_argument("--variables", nargs="*", default=None)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    contract = load_contract(args.contract)
    source_sha = sha256_file(args.observations)
    observations = pd.read_csv(args.observations, low_memory=False)
    matrix, variables = build_matrix(
        observations,
        contract,
        require_native=not args.allow_nonnative_pilot,
        variables=args.variables,
    )
    threshold = float(contract["redundancy_rule"]["default_absolute_correlation_flag"])
    weight_column = "equal_taxon_weight" if "equal_taxon_weight" in matrix.columns else None
    report, tables = diagnose_environment(matrix, variables, threshold, weight_column=weight_column)
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
    (args.out_dir / "environment_diagnostics.json").write_text(
        json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="utf-8"
    )
    print(
        json.dumps(
            {
                "status": report["status"],
                "n_rows": report["n_rows"],
                "variables": variables,
                "pilot": report["engineering_pilot_only"],
                "weight_column": report["weight_column"],
            }
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
