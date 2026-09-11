#!/usr/bin/env python3
"""Broad-space sensitivity for supported full-27/full-environment atlas rows.

This is a v2-native sensitivity.  It does not import the old nine-endpoint SPDE
result set.  A second-order spherical-coordinate basis removes broad geographic
trend, and residual Moran's I checks whether coarse residual autocorrelation is
still detectable.  Passing this screen is not equivalent to a causal spatial
model or proof of adaptation.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.spatial import cKDTree


PREDICTORS = [
    "chelsa_bio01",
    "chelsa_bio04",
    "chelsa_bio12",
    "chelsa_bio15",
    "chelsa_rsds_mean",
    "chelsa_vpd_mean",
    "chelsa_sfcwind_mean",
    "chelsa_gsp",
    "chelsa_npp",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--analysis-contract", required=True, type=Path)
    parser.add_argument("--traits-long", required=True, type=Path)
    parser.add_argument("--environment", required=True, type=Path)
    parser.add_argument("--atlas-dir", required=True, type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--permutations", type=int, default=999)
    parser.add_argument("--moran-permutations", type=int, default=999)
    parser.add_argument("--moran-maximum-observations", type=int, default=5000)
    return parser.parse_args()


def stable_rng(seed: int, *parts: str) -> np.random.Generator:
    payload = "|".join([str(seed), *map(str, parts)]).encode("utf-8")
    digest = hashlib.sha256(payload).digest()
    return np.random.default_rng(int.from_bytes(digest[:8], "little", signed=False))


def standardize(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    sd = float(np.std(values, ddof=0))
    if not np.isfinite(sd) or sd <= 0:
        raise ValueError("No finite variation")
    return (values - float(np.mean(values))) / sd


def spherical_basis(latitude: np.ndarray, longitude: np.ndarray) -> np.ndarray:
    lat = np.radians(np.asarray(latitude, dtype=float))
    lon = np.radians(np.asarray(longitude, dtype=float))
    x = np.cos(lat) * np.cos(lon)
    y = np.cos(lat) * np.sin(lon)
    z = np.sin(lat)
    return np.column_stack(
        [x, y, z, x * y, x * z, y * z, x * x - y * y, 3.0 * z * z - 1.0]
    )


def demean_by_taxon(values: np.ndarray, taxa: np.ndarray) -> np.ndarray:
    frame = pd.DataFrame(np.asarray(values, dtype=float))
    groups = pd.Series(taxa, index=frame.index)
    return (frame - frame.groupby(groups).transform("mean")).to_numpy(float)


def prepare_design(
    predictor: np.ndarray,
    basis: np.ndarray,
    taxa: np.ndarray | None,
) -> tuple[np.ndarray, np.ndarray]:
    x = np.asarray(predictor, dtype=float)[:, None]
    b = np.asarray(basis, dtype=float)
    if taxa is not None:
        x = demean_by_taxon(x, taxa)
        b = demean_by_taxon(b, taxa)
    x = standardize(x[:, 0])[:, None]
    keep = np.std(b, axis=0, ddof=0) > 1e-10
    b = b[:, keep]
    if b.shape[1]:
        b = np.column_stack([standardize(b[:, column]) for column in range(b.shape[1])])
    return np.column_stack([x, b]), b


def permute_rows_within_groups(
    values: np.ndarray,
    groups: np.ndarray | None,
    rng: np.random.Generator,
) -> np.ndarray:
    if groups is None:
        return values[rng.permutation(len(values))]
    result = values.copy()
    for group in np.unique(groups):
        index = np.flatnonzero(groups == group)
        result[index] = values[rng.permutation(index)]
    return result


def permutation_groups(groups: np.ndarray | None, n_rows: int) -> list[np.ndarray]:
    if groups is None:
        return [np.arange(n_rows)]
    return [np.flatnonzero(groups == group) for group in np.unique(groups)]


def batched_permuted_scores(
    residual: np.ndarray,
    weights: np.ndarray,
    index_groups: list[np.ndarray],
    batch_size: int,
    rng: np.random.Generator,
) -> np.ndarray:
    """Return coefficient numerators for independent row permutations.

    Rows are shuffled independently within each declared group.  Accumulating
    the weighted score group by group avoids constructing a batch by full-data
    permutation cube.
    """
    values = np.asarray(residual, dtype=float)
    vector_response = values.ndim == 1
    scores = (
        np.zeros(batch_size, dtype=float)
        if vector_response
        else np.zeros((batch_size, values.shape[1]), dtype=float)
    )
    for index in index_groups:
        order = np.broadcast_to(np.arange(len(index)), (batch_size, len(index))).copy()
        order = rng.permuted(order, axis=1)
        permuted = values[index][order]
        if vector_response:
            scores += permuted @ weights[index]
        else:
            scores += np.einsum("bnc,n->bc", permuted, weights[index], optimize=True)
    return scores


def multivariate_freedman_lane(
    responses: np.ndarray,
    design: np.ndarray,
    reduced_design: np.ndarray,
    groups: np.ndarray | None,
    permutations: int,
    rng: np.random.Generator,
) -> tuple[np.ndarray, float]:
    full_inverse = np.linalg.pinv(design)
    full_beta = full_inverse @ responses
    observed = float(np.linalg.norm(full_beta[0]))
    if reduced_design.shape[1]:
        reduced_beta = np.linalg.pinv(reduced_design) @ responses
        reduced_fitted = reduced_design @ reduced_beta
    else:
        reduced_fitted = np.zeros_like(responses)
    residual = responses - reduced_fitted
    weights = full_inverse[0]
    constant = weights @ reduced_fitted
    index_groups = permutation_groups(groups, len(residual))
    exceed = 0
    remaining = permutations
    batch_limit = 128
    while remaining:
        size = min(batch_limit, remaining)
        beta = constant[None, :] + batched_permuted_scores(
            residual, weights, index_groups, size, rng
        )
        exceed += int(np.sum(np.linalg.norm(beta, axis=1) >= observed - 1e-15))
        remaining -= size
    return np.asarray(full_beta[0], dtype=float), float((exceed + 1) / (permutations + 1))


def linear_freedman_lane(
    response: np.ndarray,
    design: np.ndarray,
    reduced_design: np.ndarray,
    groups: np.ndarray | None,
    permutations: int,
    rng: np.random.Generator,
) -> tuple[float, float, np.ndarray]:
    full_inverse = np.linalg.pinv(design)
    full_beta = full_inverse @ response
    fitted = design @ full_beta
    if reduced_design.shape[1]:
        reduced_beta = np.linalg.pinv(reduced_design) @ response
        reduced_fitted = reduced_design @ reduced_beta
    else:
        reduced_fitted = np.zeros_like(response)
    reduced_residual = response - reduced_fitted
    weights = full_inverse[0]
    constant = float(weights @ reduced_fitted)
    index_groups = permutation_groups(groups, len(reduced_residual))
    exceed = 0
    remaining = permutations
    batch_limit = 128
    while remaining:
        size = min(batch_limit, remaining)
        beta = constant + batched_permuted_scores(
            reduced_residual, weights, index_groups, size, rng
        )
        exceed += int(np.sum(np.abs(beta) >= abs(full_beta[0]) - 1e-15))
        remaining -= size
    return float(full_beta[0]), float((exceed + 1) / (permutations + 1)), response - fitted


def morans_i_from_neighbors(values: np.ndarray, neighbours: np.ndarray) -> float:
    response = np.asarray(values, dtype=float)
    centred = response - float(np.mean(response))
    denominator = float(np.dot(centred, centred))
    if denominator <= 0 or len(response) <= neighbours.shape[1]:
        return float("nan")
    numerator = float(np.sum(centred[:, None] * centred[neighbours]))
    return float(
        len(response) / (len(response) * neighbours.shape[1]) * numerator / denominator
    )


def morans_i(values: np.ndarray, latitude: np.ndarray, longitude: np.ndarray, k: int = 8) -> float:
    xyz = spherical_basis(latitude, longitude)[:, :3]
    neighbours = cKDTree(xyz).query(xyz, k=k + 1)[1][:, 1:]
    return morans_i_from_neighbors(values, neighbours)


def moran_test(
    residual: np.ndarray,
    latitude: np.ndarray,
    longitude: np.ndarray,
    permutations: int,
    maximum_observations: int,
    rng: np.random.Generator,
) -> tuple[float, float, int]:
    values = np.asarray(residual, dtype=float)
    lat = np.asarray(latitude, dtype=float)
    lon = np.asarray(longitude, dtype=float)
    if len(values) > maximum_observations:
        index = np.sort(rng.choice(len(values), maximum_observations, replace=False))
        values, lat, lon = values[index], lat[index], lon[index]
    xyz = spherical_basis(lat, lon)[:, :3]
    neighbours = cKDTree(xyz).query(xyz, k=9)[1][:, 1:]
    observed = morans_i_from_neighbors(values, neighbours)
    if not np.isfinite(observed):
        return observed, float("nan"), len(values)
    centred = values - float(np.mean(values))
    denominator = float(np.dot(centred, centred))
    scale = 1.0 / (neighbours.shape[1] * denominator)
    exceed = 0
    remaining = permutations
    batch_limit = 32
    while remaining:
        size = min(batch_limit, remaining)
        matrix = np.broadcast_to(centred, (size, len(centred))).copy()
        permuted = rng.permuted(matrix, axis=1)
        numerators = np.sum(
            permuted[:, :, None] * permuted[:, neighbours], axis=(1, 2)
        )
        simulated = numerators * scale
        exceed += int(np.sum(np.abs(simulated) >= abs(observed) - 1e-15))
        remaining -= size
    return observed, float((exceed + 1) / (permutations + 1)), len(values)


def unit_data(
    row: pd.Series,
    traits: pd.DataFrame,
    environment: pd.DataFrame,
) -> pd.DataFrame:
    members = str(row["member_endpoint_ids"]).split("|")
    part = traits[traits["endpoint_id"].isin(members)][
        ["obs_id", "taxon_name", "endpoint_id", "value"]
    ]
    wide = part.pivot(index=["obs_id", "taxon_name"], columns="endpoint_id", values="value").reset_index()
    wide.columns.name = None
    return wide.merge(
        environment.drop(columns="taxon_name"), on="obs_id", how="inner", validate="one_to_one"
    )


def load_inputs(traits_path: Path, environment_path: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    traits = pd.read_csv(
        traits_path,
        usecols=["obs_id", "taxon_name", "endpoint_id", "measurement_available", "value"],
        low_memory=False,
    )
    environment = pd.read_csv(
        environment_path,
        usecols=["obs_id", "taxon_name", "latitude", "longitude", *PREDICTORS],
        low_memory=False,
    )
    traits["obs_id"] = traits["obs_id"].astype(str)
    environment["obs_id"] = environment["obs_id"].astype(str)
    available = traits["measurement_available"].astype(str).str.lower().isin({"true", "1", "yes"})
    traits["value"] = pd.to_numeric(traits["value"], errors="coerce")
    traits = traits[available & traits["value"].notna()].copy()
    return traits, environment


def spatial_within(
    selected: pd.DataFrame,
    traits: pd.DataFrame,
    environment: pd.DataFrame,
    seed: int,
    permutations: int,
    moran_permutations: int,
    moran_maximum: int,
) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for _, atlas_row in selected.iterrows():
        members = str(atlas_row["member_endpoint_ids"]).split("|")
        data = unit_data(atlas_row, traits, environment).dropna(
            subset=[*members, atlas_row["predictor"], "latitude", "longitude"]
        )
        counts = data.groupby("taxon_name").size()
        data = data[data["taxon_name"].isin(counts[counts >= 2].index)].copy()
        taxa = data["taxon_name"].astype(str).to_numpy()
        basis = spherical_basis(data["latitude"], data["longitude"])
        design, reduced = prepare_design(data[atlas_row["predictor"]].to_numpy(float), basis, taxa)
        base = {
            "scale": "within_taxon",
            "unit_id": atlas_row["unit_id"],
            "member_endpoint_ids": atlas_row["member_endpoint_ids"],
            "inferential_unit": atlas_row["inferential_unit"],
            "predictor": atlas_row["predictor"],
            "environment_block": atlas_row["environment_block"],
            "base_beta_std": atlas_row.get("beta_std", np.nan),
            "base_q_fdr_bh_global_family": atlas_row["q_fdr_bh_global_family"],
            "n_observations": int(len(data)),
            "n_taxa": int(data["taxon_name"].nunique()),
        }
        try:
            rng = stable_rng(seed, "spatial_within", atlas_row["unit_id"], atlas_row["predictor"])
            if atlas_row["inferential_unit"] == "linear_endpoint":
                response = demean_by_taxon(data[members[0]].to_numpy(float)[:, None], taxa)[:, 0]
                response = standardize(response)
                beta, p_value, residual = linear_freedman_lane(
                    response, design, reduced, taxa, permutations, rng
                )
                clustered = sm.OLS(response, design).fit(
                    cov_type="cluster", cov_kwds={"groups": taxa}
                )
                fit = {
                    "spatial_beta_std": beta,
                    "spatial_standard_error_clustered": float(clustered.bse[0]),
                    "spatial_confidence_low_95": float(clustered.conf_int()[0, 0]),
                    "spatial_confidence_high_95": float(clustered.conf_int()[0, 1]),
                    "spatial_permutation_p_value": p_value,
                    "same_linear_direction_as_base": bool(
                        np.sign(beta) == np.sign(float(atlas_row["beta_std"]))
                    ),
                }
            else:
                response = demean_by_taxon(data[members].to_numpy(float), taxa)
                response = np.column_stack(
                    [standardize(response[:, column]) for column in range(response.shape[1])]
                )
                beta, p_value = multivariate_freedman_lane(
                    response, design, reduced, taxa, permutations, rng
                )
                residual = response - design @ np.linalg.lstsq(design, response, rcond=None)[0]
                residual = np.linalg.norm(residual, axis=1)
                fit = {
                    "spatial_beta_sine_std": float(beta[0]),
                    "spatial_beta_cosine_std": float(beta[1]),
                    "spatial_effect_magnitude": float(np.linalg.norm(beta)),
                    "spatial_effect_direction_degrees": float(math.degrees(math.atan2(beta[0], beta[1])) % 360),
                    "spatial_permutation_p_value": p_value,
                    "same_linear_direction_as_base": np.nan,
                }
            moran_rng = stable_rng(seed, "moran_within", atlas_row["unit_id"], atlas_row["predictor"])
            statistic, moran_p, moran_n = moran_test(
                residual,
                data["latitude"].to_numpy(float),
                data["longitude"].to_numpy(float),
                moran_permutations,
                moran_maximum,
                moran_rng,
            )
            linear_direction_ok = (
                fit["same_linear_direction_as_base"]
                if atlas_row["inferential_unit"] == "linear_endpoint"
                else True
            )
            rows.append(
                {
                    **base,
                    "status": "ok",
                    **fit,
                    "residual_morans_i": statistic,
                    "residual_morans_p_value": moran_p,
                    "residual_morans_n": moran_n,
                    "broad_spatial_sensitivity_pass": bool(
                        p_value < 0.05
                        and linear_direction_ok
                        and np.isfinite(moran_p)
                        and moran_p >= 0.05
                    ),
                }
            )
        except Exception as error:
            rows.append({**base, "status": f"failed:{type(error).__name__}:{error}"})
    return pd.DataFrame(rows)


def taxon_level_data(
    atlas_row: pd.Series,
    traits: pd.DataFrame,
    environment: pd.DataFrame,
    minimum_trait_observations: int,
) -> pd.DataFrame:
    members = str(atlas_row["member_endpoint_ids"]).split("|")
    part = traits[traits["endpoint_id"].isin(members)][
        ["obs_id", "taxon_name", "endpoint_id", "value"]
    ]
    wide = part.pivot(index=["obs_id", "taxon_name"], columns="endpoint_id", values="value").reset_index()
    wide.columns.name = None
    wide = wide.dropna(subset=members)
    counts = wide.groupby("taxon_name").size().rename("n_trait_observations")
    trait_medians = wide.groupby("taxon_name")[members].median().join(counts)
    trait_medians = trait_medians[trait_medians["n_trait_observations"].ge(minimum_trait_observations)]
    env = environment.copy()
    basis = spherical_basis(env["latitude"], env["longitude"])
    for index in range(basis.shape[1]):
        env[f"spatial_basis_{index}"] = basis[:, index]
    columns = [*PREDICTORS, "latitude", "longitude", *[f"spatial_basis_{i}" for i in range(basis.shape[1])]]
    env_taxon = env.groupby("taxon_name")[columns].median(numeric_only=True)
    return trait_medians.join(env_taxon, how="inner").reset_index()


def spatial_among(
    selected: pd.DataFrame,
    traits: pd.DataFrame,
    environment: pd.DataFrame,
    seed: int,
    permutations: int,
    moran_permutations: int,
    moran_maximum: int,
) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for _, atlas_row in selected.iterrows():
        minimum_trait_observations = int(atlas_row["minimum_trait_observations_per_taxon"])
        members = str(atlas_row["member_endpoint_ids"]).split("|")
        data = taxon_level_data(atlas_row, traits, environment, minimum_trait_observations)
        basis_columns = [f"spatial_basis_{index}" for index in range(8)]
        data = data.dropna(subset=[*members, atlas_row["predictor"], *basis_columns])
        design, reduced = prepare_design(
            data[atlas_row["predictor"]].to_numpy(float), data[basis_columns].to_numpy(float), None
        )
        base = {
            "scale": "among_taxon",
            "scope": atlas_row["scope"],
            "unit_id": atlas_row["unit_id"],
            "member_endpoint_ids": atlas_row["member_endpoint_ids"],
            "inferential_unit": atlas_row["inferential_unit"],
            "predictor": atlas_row["predictor"],
            "environment_block": atlas_row["environment_block"],
            "base_beta_std": atlas_row.get("beta_std", np.nan),
            "base_q_fdr_bh_global_family": atlas_row["q_fdr_bh_global_family"],
            "n_taxa": int(len(data)),
            "minimum_trait_observations_per_taxon": minimum_trait_observations,
        }
        try:
            rng = stable_rng(seed, "spatial_among", atlas_row["unit_id"], atlas_row["predictor"])
            if atlas_row["inferential_unit"] == "linear_endpoint":
                response = standardize(data[members[0]].to_numpy(float))
                beta, p_value, residual = linear_freedman_lane(
                    response, design, reduced, None, permutations, rng
                )
                fit = {
                    "spatial_beta_std": beta,
                    "spatial_permutation_p_value": p_value,
                    "same_linear_direction_as_base": bool(
                        np.sign(beta) == np.sign(float(atlas_row["beta_std"]))
                    ),
                }
            else:
                response = np.column_stack(
                    [standardize(data[member].to_numpy(float)) for member in members]
                )
                beta, p_value = multivariate_freedman_lane(
                    response, design, reduced, None, permutations, rng
                )
                residual = response - design @ np.linalg.lstsq(design, response, rcond=None)[0]
                residual = np.linalg.norm(residual, axis=1)
                fit = {
                    "spatial_beta_sine_std": float(beta[0]),
                    "spatial_beta_cosine_std": float(beta[1]),
                    "spatial_effect_magnitude": float(np.linalg.norm(beta)),
                    "spatial_effect_direction_degrees": float(math.degrees(math.atan2(beta[0], beta[1])) % 360),
                    "spatial_permutation_p_value": p_value,
                    "same_linear_direction_as_base": np.nan,
                }
            moran_rng = stable_rng(seed, "moran_among", atlas_row["unit_id"], atlas_row["predictor"])
            statistic, moran_p, moran_n = moran_test(
                residual,
                data["latitude"].to_numpy(float),
                data["longitude"].to_numpy(float),
                moran_permutations,
                moran_maximum,
                moran_rng,
            )
            linear_direction_ok = (
                fit["same_linear_direction_as_base"]
                if atlas_row["inferential_unit"] == "linear_endpoint"
                else True
            )
            rows.append(
                {
                    **base,
                    "status": "ok",
                    **fit,
                    "residual_morans_i": statistic,
                    "residual_morans_p_value": moran_p,
                    "residual_morans_n": moran_n,
                    "broad_spatial_sensitivity_pass": bool(
                        p_value < 0.05
                        and linear_direction_ok
                        and np.isfinite(moran_p)
                        and moran_p >= 0.05
                    ),
                }
            )
        except Exception as error:
            rows.append({**base, "status": f"failed:{type(error).__name__}:{error}"})
    return pd.DataFrame(rows)


def main() -> int:
    args = parse_args()
    contract = json.loads(args.analysis_contract.read_text(encoding="utf-8"))
    seed = int(contract["seed"])
    traits, environment = load_inputs(args.traits_long, args.environment)
    within_atlas = pd.read_csv(args.atlas_dir / "v2_full27_environment_within.csv", low_memory=False)
    among_atlas = pd.read_csv(args.atlas_dir / "v2_full27_environment_among.csv", low_memory=False)
    selected_within = within_atlas[
        within_atlas["status"].eq("ok")
        & within_atlas["q_fdr_bh_global_family"].lt(0.05)
    ].copy()
    primary_scope = f"among_taxon_min{int(contract['among_taxon']['primary_minimum_trait_observations_per_taxon'])}"
    selected_among = among_atlas[
        among_atlas["scope"].eq(primary_scope)
        & among_atlas["status"].eq("ok")
        & among_atlas["q_fdr_bh_global_family"].lt(0.05)
    ].copy()
    within = spatial_within(
        selected_within,
        traits,
        environment,
        seed,
        args.permutations,
        args.moran_permutations,
        args.moran_maximum_observations,
    )
    among = spatial_among(
        selected_among,
        traits,
        environment,
        seed,
        args.permutations,
        args.moran_permutations,
        args.moran_maximum_observations,
    )
    args.out_dir.mkdir(parents=True, exist_ok=True)
    within.to_csv(args.out_dir / "v2_full27_spatial_within.csv", index=False, encoding="utf-8-sig")
    among.to_csv(args.out_dir / "v2_full27_spatial_among.csv", index=False, encoding="utf-8-sig")
    report = {
        "analysis_id": "geb_v2_full27_broad_spatial_sensitivity_v1",
        "status": "retrospective_selected-row_sensitivity",
        "n_within_rows_entered": int(len(selected_within)),
        "n_within_rows_passed": int(within.get("broad_spatial_sensitivity_pass", pd.Series(dtype=bool)).sum()),
        "n_among_rows_entered": int(len(selected_among)),
        "n_among_rows_passed": int(among.get("broad_spatial_sensitivity_pass", pd.Series(dtype=bool)).sum()),
        "method": contract["spatial_sensitivity"],
        "claim_boundary": "Broad spherical-trend sensitivity plus residual Moran diagnostic; not the older SPDE result, not full removal of spatial process, and not evidence of adaptation by itself.",
    }
    (args.out_dir / "v2_full27_spatial_sensitivity_report.json").write_text(
        json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(report, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
