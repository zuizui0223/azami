#!/usr/bin/env python3
"""Post-hoc full-cohort VIF-step multivariable sensitivity for Chapter 1 v3.

This analysis does not replace the frozen v2 atlas.  It uses the same 46,276
strict-spatial observations and the same nine frozen CHELSA predictors, restores
five fields that were already measured historically but omitted from the old
observation aggregation, then asks whether trait-environment associations remain
when collinearity is reduced by response-blind sequential VIF filtering.

The primary sensitivity threshold is VIF < 10.  VIF < 5 is reported as a stricter
sensitivity.  Predictor selection is based on taxon-median environment values and
never on trait outcomes.  Among-taxon trait models retain the frozen min-5 rule.
Permutation p-values are Freedman-Lane style residual permutations conditional on
all other retained predictors; BH correction is applied within each post-hoc
threshold family.  These are exploratory diagnostics, not a replacement for the
pre-registered marginal v2 family and not causal effect estimates.
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

PREDICTORS = [
    "chelsa_bio01", "chelsa_bio04", "chelsa_bio12", "chelsa_bio15",
    "chelsa_rsds_mean", "chelsa_vpd_mean", "chelsa_sfcwind_mean",
    "chelsa_gsp", "chelsa_npp",
]
RECOVERED = {
    "visible_floret_fraction": "corolla_visible_fraction",
    "corolla_white_pixel_fraction": "corolla_white_fraction",
    "corolla_redmagenta_pixel_fraction": "corolla_redmagenta_fraction",
    "corolla_purple_pixel_fraction": "corolla_purple_fraction",
    "corolla_yellow_pixel_fraction": "corolla_yellow_fraction",
}
EXPECTED_TRAIT_SHA = "d775794f2bce2dfd0c1f63c5c8e01778c518f6eeb327bf0d9944045143a02344"
EXPECTED_ENV_SHA = "e242aa7ce69d12b11937c1335e84b9638799c50b42ef36b95725e77190df98e7"
EXPECTED_ENV_ROWS = 46276
EXPECTED_ENV_TAXA = 259


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def stable_rng(seed: int, *parts: str) -> np.random.Generator:
    payload = "|".join([str(seed), *parts]).encode()
    digest = hashlib.sha256(payload).digest()
    return np.random.default_rng(int.from_bytes(digest[:8], "little", signed=False))


def zscore(a: np.ndarray) -> np.ndarray:
    a = np.asarray(a, float)
    sd = np.std(a, ddof=0)
    if not np.isfinite(a).all() or not np.isfinite(sd) or sd <= 0:
        raise ValueError("nonfinite or constant vector")
    return (a - np.mean(a)) / sd


def bh_adjust(values: pd.Series) -> pd.Series:
    p = np.asarray(values, float)
    n = len(p)
    if n == 0:
        return pd.Series(dtype=float, index=values.index)
    order = np.argsort(p)
    ranked = p[order]
    adj = ranked * n / np.arange(1, n + 1)
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    adj = np.clip(adj, 0, 1)
    out = np.empty(n, float)
    out[order] = adj
    return pd.Series(out, index=values.index)


def vif_table(frame: pd.DataFrame, columns: list[str]) -> pd.DataFrame:
    z = np.column_stack([zscore(frame[c].to_numpy(float)) for c in columns])
    rows = []
    for j, name in enumerate(columns):
        y = z[:, j]
        other = np.delete(z, j, axis=1)
        if other.shape[1] == 0:
            vif = 1.0
            r2 = 0.0
        else:
            x = np.column_stack([np.ones(len(y)), other])
            fit = x @ np.linalg.lstsq(x, y, rcond=None)[0]
            sse = float(np.sum((y - fit) ** 2))
            sst = float(np.sum((y - np.mean(y)) ** 2))
            r2 = 1.0 - sse / sst
            vif = float("inf") if 1 - r2 <= 1e-12 else 1.0 / (1.0 - r2)
        rows.append({"predictor": name, "vif": float(vif), "r2_against_others": float(r2)})
    return pd.DataFrame(rows).sort_values(["vif", "predictor"], ascending=[False, True]).reset_index(drop=True)


def vifstep(frame: pd.DataFrame, threshold: float) -> tuple[list[str], list[dict[str, Any]], pd.DataFrame]:
    selected = list(PREDICTORS)
    history: list[dict[str, Any]] = []
    step = 0
    while True:
        table = vif_table(frame, selected)
        maximum = float(table.iloc[0].vif)
        history.append({
            "step": step,
            "n_predictors": len(selected),
            "max_vif": maximum,
            "max_vif_predictor": str(table.iloc[0].predictor),
            "predictors": list(selected),
        })
        if maximum < threshold or len(selected) <= 1:
            return selected, history, table
        drop = str(table.iloc[0].predictor)
        selected.remove(drop)
        step += 1


def inferential_units(contract: pd.DataFrame) -> list[dict[str, Any]]:
    work = contract.copy()
    work["circular_group"] = work["circular_group"].fillna("").astype(str).str.strip()
    units: list[dict[str, Any]] = []
    for row in work[work.circular_group.eq("")].to_dict("records"):
        units.append({
            "unit_id": row["endpoint_id"],
            "members": [row["endpoint_id"]],
            "module": row["module"],
            "analysis_tier": row["analysis_tier"],
            "validation_status": row["validation_status"],
            "inferential_unit": "linear_endpoint",
        })
    for group, part in work[work.circular_group.ne("")].groupby("circular_group", sort=True):
        members = part.endpoint_id.tolist()
        if len(members) != 2 or not any("sin" in x for x in members) or not any("cos" in x for x in members):
            raise ValueError(f"bad circular group {group}: {members}")
        units.append({
            "unit_id": str(group),
            "members": members,
            "module": str(part.iloc[0].module),
            "analysis_tier": str(part.iloc[0].analysis_tier),
            "validation_status": str(part.iloc[0].validation_status),
            "inferential_unit": "circular_joint",
        })
    return sorted(units, key=lambda x: (x["module"], x["unit_id"]))


def restore_historical_fields(
    traits: pd.DataFrame, environment: pd.DataFrame, contract: pd.DataFrame, recovered_path: Path
) -> tuple[pd.DataFrame, dict[str, Any]]:
    recovered = pd.read_csv(recovered_path, low_memory=False)
    recovered["obs_id"] = recovered.obs_id.astype(str)
    if recovered.obs_id.duplicated().any():
        raise ValueError("recovered display table not unique by obs_id")
    env_taxon = environment.set_index("obs_id").taxon_name.astype(str)
    meta = contract.set_index("endpoint_id")
    report: dict[str, Any] = {}
    for endpoint, column in RECOVERED.items():
        if endpoint not in meta.index or column not in recovered.columns:
            raise ValueError(f"missing recovery mapping {endpoint} <- {column}")
        values = pd.to_numeric(recovered[column], errors="coerce")
        mapping = dict(zip(recovered.obs_id, values))
        existing = traits.endpoint_id.eq(endpoint)
        existing_ids = set(traits.loc[existing, "obs_id"])
        mapped = traits.loc[existing, "obs_id"].map(mapping)
        usable = mapped.notna()
        update_idx = traits.loc[existing].index[usable]
        traits.loc[update_idx, "value"] = mapped[usable].to_numpy(float)
        traits.loc[update_idx, "measurement_available"] = True
        add_ids = [obs for obs, value in mapping.items() if pd.notna(value) and obs not in existing_ids and obs in env_taxon.index]
        if add_ids:
            m = meta.loc[endpoint]
            add = pd.DataFrame({
                "obs_id": add_ids,
                "taxon_name": [env_taxon.loc[x] for x in add_ids],
                "endpoint_id": endpoint,
                "module": str(m.module),
                "analysis_tier": str(m.analysis_tier),
                "measurement_available": True,
                "value": [float(mapping[x]) for x in add_ids],
            })
            traits = pd.concat([traits, add], ignore_index=True, sort=False)
        report[endpoint] = {
            "updated_existing_rows": int(len(update_idx)),
            "appended_rows": int(len(add_ids)),
            "finite_rows_after_restore": int(pd.to_numeric(
                traits.loc[traits.endpoint_id.eq(endpoint), "value"], errors="coerce").notna().sum()),
        }
    if traits.duplicated(["obs_id", "endpoint_id"]).any():
        raise ValueError("duplicate obs_id/endpoint_id after restoration")
    return traits, report


def load_inputs(args: argparse.Namespace):
    if sha256(args.traits_long) != EXPECTED_TRAIT_SHA:
        raise ValueError("frozen trait-universe hash mismatch")
    if sha256(args.environment) != EXPECTED_ENV_SHA:
        raise ValueError("nine-predictor environment hash mismatch")
    contract = pd.read_csv(args.contract, dtype=str, keep_default_na=False)
    traits = pd.read_csv(
        args.traits_long,
        usecols=["obs_id", "taxon_name", "endpoint_id", "module", "analysis_tier", "measurement_available", "value"],
        low_memory=False,
    )
    environment = pd.read_csv(args.environment, low_memory=False)
    traits["obs_id"] = traits.obs_id.astype(str)
    traits["taxon_name"] = traits.taxon_name.astype(str)
    environment["obs_id"] = environment.obs_id.astype(str)
    environment["taxon_name"] = environment.taxon_name.astype(str)
    if len(environment) != EXPECTED_ENV_ROWS or environment.obs_id.nunique() != EXPECTED_ENV_ROWS:
        raise ValueError("wrong environment cohort")
    if environment.taxon_name.nunique() != EXPECTED_ENV_TAXA:
        raise ValueError("wrong environment taxon count")
    if traits.duplicated(["obs_id", "endpoint_id"]).any():
        raise ValueError("duplicate trait rows")
    for p in PREDICTORS:
        environment[p] = pd.to_numeric(environment[p], errors="coerce")
        if environment[p].notna().mean() < .98:
            raise ValueError(f"low coverage {p}")
    traits["value"] = pd.to_numeric(traits.value, errors="coerce")
    available = traits.measurement_available.astype(str).str.strip().str.lower().isin({"true", "1", "yes"})
    traits = traits.loc[available & traits.value.notna()].copy()
    traits, recovery = restore_historical_fields(traits, environment, contract, args.recovered_display)
    return contract, traits, environment, recovery


def taxon_environment(environment: pd.DataFrame) -> pd.DataFrame:
    out = environment.groupby("taxon_name", sort=True)[PREDICTORS].median()
    out = out.dropna(subset=PREDICTORS)
    if len(out) < 200:
        raise ValueError("unexpectedly small taxon environment table")
    return out


def linear_taxon_table(traits: pd.DataFrame, env_taxon: pd.DataFrame, endpoint: str, minimum: int) -> pd.DataFrame:
    part = traits.loc[traits.endpoint_id.eq(endpoint), ["taxon_name", "value"]]
    g = part.groupby("taxon_name").value.agg(["median", "count"])
    g = g.loc[g["count"] >= minimum]
    return g.join(env_taxon, how="inner").dropna()


def circular_taxon_table(traits: pd.DataFrame, env_taxon: pd.DataFrame, members: list[str], minimum: int) -> pd.DataFrame:
    sin_id = next(x for x in members if "sin" in x)
    cos_id = next(x for x in members if "cos" in x)
    frames = []
    for endpoint, name in [(sin_id, "sine"), (cos_id, "cosine")]:
        part = traits.loc[traits.endpoint_id.eq(endpoint), ["taxon_name", "value"]]
        g = part.groupby("taxon_name").value.agg(["median", "count"])
        g = g.loc[g["count"] >= minimum][["median"]].rename(columns={"median": name})
        frames.append(g)
    return frames[0].join(frames[1], how="inner").join(env_taxon, how="inner").dropna()


def reduced_residual(y: np.ndarray, x_other: np.ndarray) -> np.ndarray:
    if x_other.shape[1] == 0:
        design = np.ones((len(y), 1))
    else:
        design = np.column_stack([np.ones(len(y)), x_other])
    return y - design @ np.linalg.lstsq(design, y, rcond=None)[0]


def residualized_x(x: np.ndarray, x_other: np.ndarray) -> np.ndarray:
    return reduced_residual(x, x_other)


def permutation_linear(y: np.ndarray, x: np.ndarray, other: np.ndarray, permutations: int, rng: np.random.Generator) -> float:
    yres = reduced_residual(y, other)
    xres = residualized_x(x, other)
    denom = float(np.dot(xres, xres))
    if denom <= 1e-12:
        raise ValueError("zero partial predictor variance")
    observed = float(np.dot(xres, yres) / denom)
    exceed = 0
    remaining = permutations
    while remaining:
        size = min(256, remaining)
        matrix = np.broadcast_to(yres, (size, len(yres))).copy()
        permuted = rng.permuted(matrix, axis=1)
        simulated = permuted @ xres / denom
        exceed += int(np.sum(np.abs(simulated) >= abs(observed) - 1e-15))
        remaining -= size
    return float((exceed + 1) / (permutations + 1))


def permutation_circular(sine: np.ndarray, cosine: np.ndarray, x: np.ndarray, other: np.ndarray,
                         permutations: int, rng: np.random.Generator) -> float:
    sres = reduced_residual(sine, other)
    cres = reduced_residual(cosine, other)
    xres = residualized_x(x, other)
    denom = float(np.dot(xres, xres))
    if denom <= 1e-12:
        raise ValueError("zero partial predictor variance")
    bs = float(np.dot(xres, sres) / denom)
    bc = float(np.dot(xres, cres) / denom)
    observed = float(math.hypot(bs, bc))
    exceed = 0
    remaining = permutations
    pair = np.column_stack([sres, cres])
    while remaining:
        size = min(256, remaining)
        # Keep sine/cosine residual pairs together under each permutation.
        sims = np.empty((size, 2), float)
        for i in range(size):
            idx = rng.permutation(len(xres))
            sims[i, 0] = float(np.dot(xres, pair[idx, 0]) / denom)
            sims[i, 1] = float(np.dot(xres, pair[idx, 1]) / denom)
        exceed += int(np.sum(np.hypot(sims[:, 0], sims[:, 1]) >= observed - 1e-15))
        remaining -= size
    return float((exceed + 1) / (permutations + 1))


def fit_linear(frame: pd.DataFrame, selected: list[str], unit: dict[str, Any], threshold: float,
               permutations: int, seed: int) -> list[dict[str, Any]]:
    y = zscore(frame["median"].to_numpy(float))
    xz = np.column_stack([zscore(frame[p].to_numpy(float)) for p in selected])
    design = sm.add_constant(xz, has_constant="add")
    model = sm.OLS(y, design).fit().get_robustcov_results(cov_type="HC3", use_t=False)
    rows = []
    for j, predictor in enumerate(selected):
        other = np.delete(xz, j, axis=1)
        p_perm = permutation_linear(
            y, xz[:, j], other, permutations,
            stable_rng(seed, str(threshold), unit["unit_id"], predictor, "linear"),
        )
        rows.append({
            "vif_threshold": threshold,
            "unit_id": unit["unit_id"], "module": unit["module"],
            "analysis_tier": unit["analysis_tier"], "validation_status": unit["validation_status"],
            "inferential_unit": unit["inferential_unit"], "predictor": predictor,
            "n_taxa": len(frame), "beta_std_partial": float(model.params[j + 1]),
            "hc3_se": float(model.bse[j + 1]), "hc3_ci_low": float(model.conf_int()[j + 1, 0]),
            "hc3_ci_high": float(model.conf_int()[j + 1, 1]), "p_hc3": float(model.pvalues[j + 1]),
            "p_perm": p_perm, "beta_sine_std_partial": np.nan, "beta_cosine_std_partial": np.nan,
            "effect_magnitude_partial": np.nan, "effect_direction_degrees_partial": np.nan,
        })
    return rows


def fit_circular(frame: pd.DataFrame, selected: list[str], unit: dict[str, Any], threshold: float,
                 permutations: int, seed: int) -> list[dict[str, Any]]:
    sine = zscore(frame["sine"].to_numpy(float))
    cosine = zscore(frame["cosine"].to_numpy(float))
    xz = np.column_stack([zscore(frame[p].to_numpy(float)) for p in selected])
    design = sm.add_constant(xz, has_constant="add")
    fit_s = sm.OLS(sine, design).fit().get_robustcov_results(cov_type="HC3", use_t=False)
    fit_c = sm.OLS(cosine, design).fit().get_robustcov_results(cov_type="HC3", use_t=False)
    rows = []
    for j, predictor in enumerate(selected):
        other = np.delete(xz, j, axis=1)
        p_perm = permutation_circular(
            sine, cosine, xz[:, j], other, permutations,
            stable_rng(seed, str(threshold), unit["unit_id"], predictor, "circular"),
        )
        bs, bc = float(fit_s.params[j + 1]), float(fit_c.params[j + 1])
        rows.append({
            "vif_threshold": threshold,
            "unit_id": unit["unit_id"], "module": unit["module"],
            "analysis_tier": unit["analysis_tier"], "validation_status": unit["validation_status"],
            "inferential_unit": unit["inferential_unit"], "predictor": predictor,
            "n_taxa": len(frame), "beta_std_partial": np.nan, "hc3_se": np.nan,
            "hc3_ci_low": np.nan, "hc3_ci_high": np.nan, "p_hc3": np.nan,
            "p_perm": p_perm, "beta_sine_std_partial": bs, "beta_cosine_std_partial": bc,
            "effect_magnitude_partial": float(math.hypot(bs, bc)),
            "effect_direction_degrees_partial": float(math.degrees(math.atan2(bs, bc)) % 360.0),
        })
    return rows


def verify_frozen_univariate(units: list[dict[str, Any]], traits: pd.DataFrame, env_taxon: pd.DataFrame,
                             frozen_path: Path, minimum: int) -> dict[str, Any]:
    frozen = pd.read_csv(frozen_path, low_memory=False)
    frozen = frozen.loc[frozen.scope.eq(f"among_taxon_min{minimum}") & frozen.status.eq("ok")].copy()
    checks = []
    max_error = 0.0
    for unit in units:
        if unit["inferential_unit"] == "linear_endpoint":
            table = linear_taxon_table(traits, env_taxon, unit["members"][0], minimum)
            if len(table) < 3:
                continue
            for predictor in PREDICTORS:
                row = frozen.loc[(frozen.unit_id.eq(unit["unit_id"])) & frozen.predictor.eq(predictor)]
                if row.empty:
                    continue
                beta = float(np.corrcoef(zscore(table["median"].to_numpy(float)), zscore(table[predictor].to_numpy(float)))[0, 1])
                expected = float(row.iloc[0].beta_std)
                err = abs(beta - expected)
                max_error = max(max_error, err)
                checks.append((unit["unit_id"], predictor, err))
        else:
            table = circular_taxon_table(traits, env_taxon, unit["members"], minimum)
            if len(table) < 3:
                continue
            for predictor in PREDICTORS:
                row = frozen.loc[(frozen.unit_id.eq(unit["unit_id"])) & frozen.predictor.eq(predictor)]
                if row.empty:
                    continue
                x = zscore(table[predictor].to_numpy(float))
                bs = float(np.dot(x, zscore(table.sine.to_numpy(float))) / np.dot(x, x))
                bc = float(np.dot(x, zscore(table.cosine.to_numpy(float))) / np.dot(x, x))
                err = max(abs(bs - float(row.iloc[0].beta_sine_std_among)), abs(bc - float(row.iloc[0].beta_cosine_std_among)))
                max_error = max(max_error, err)
                checks.append((unit["unit_id"], predictor, err))
    # Restored endpoints were not present numerically in the frozen v2 output; all comparable rows must match.
    if max_error > 1e-10:
        worst = sorted(checks, key=lambda x: x[2], reverse=True)[:10]
        raise ValueError(f"taxon aggregation does not reproduce frozen univariate coefficients: {worst}")
    return {"n_comparable_rows": len(checks), "maximum_absolute_coefficient_error": max_error}


def run_threshold(threshold: float, units: list[dict[str, Any]], traits: pd.DataFrame,
                  env_taxon: pd.DataFrame, global_selected: list[str], minimum: int,
                  permutations: int, seed: int) -> tuple[pd.DataFrame, pd.DataFrame]:
    results: list[dict[str, Any]] = []
    selections = []
    for unit in units:
        if unit["inferential_unit"] == "linear_endpoint":
            frame = linear_taxon_table(traits, env_taxon, unit["members"][0], minimum)
        else:
            frame = circular_taxon_table(traits, env_taxon, unit["members"], minimum)
        if len(frame) <= len(global_selected) + 3:
            selections.append({"vif_threshold": threshold, "unit_id": unit["unit_id"], "n_taxa": len(frame),
                               "status": "insufficient_taxa", "selected_predictors": "", "max_vif": np.nan})
            continue
        selected = list(global_selected)
        # Fail-safe: subset-specific missingness can recreate collinearity. Remove further predictors only if needed.
        while len(selected) > 1:
            vt = vif_table(frame, selected)
            if float(vt.iloc[0].vif) < threshold:
                break
            selected.remove(str(vt.iloc[0].predictor))
        vt = vif_table(frame, selected)
        selections.append({
            "vif_threshold": threshold, "unit_id": unit["unit_id"], "n_taxa": len(frame), "status": "ok",
            "selected_predictors": ";".join(selected), "max_vif": float(vt.vif.max()),
        })
        if unit["inferential_unit"] == "linear_endpoint":
            results.extend(fit_linear(frame, selected, unit, threshold, permutations, seed))
        else:
            results.extend(fit_circular(frame, selected, unit, threshold, permutations, seed))
    result = pd.DataFrame(results)
    if len(result):
        result["q_perm_bh_posthoc_family"] = bh_adjust(result.p_perm)
        result["posthoc_fdr_significant_0_05"] = result.q_perm_bh_posthoc_family.lt(.05)
    return result, pd.DataFrame(selections)


def compare_original_hits(frozen_path: Path, result10: pd.DataFrame, selections10: pd.DataFrame) -> pd.DataFrame:
    frozen = pd.read_csv(frozen_path, low_memory=False)
    original = frozen.loc[
        frozen.scope.eq("among_taxon_min5") & frozen.status.eq("ok") &
        pd.to_numeric(frozen.q_fdr_bh_global_family, errors="coerce").lt(.05)
    ].copy()
    keep = ["unit_id", "module", "inferential_unit", "predictor", "beta_std", "p_value",
            "beta_sine_std_among", "beta_cosine_std_among", "effect_magnitude", "effect_direction_degrees",
            "q_fdr_bh_global_family"]
    original = original[keep].rename(columns={c: "frozen_" + c for c in keep if c not in {"unit_id", "module", "inferential_unit", "predictor"}})
    merged = original.merge(result10, on=["unit_id", "module", "inferential_unit", "predictor"], how="left")
    selection_map = selections10.set_index("unit_id").selected_predictors.to_dict()
    merged["vif10_predictor_retained"] = [
        predictor in str(selection_map.get(unit, "")).split(";")
        for unit, predictor in zip(merged.unit_id, merged.predictor)
    ]
    merged["vif10_result_status"] = np.where(
        merged.vif10_predictor_retained,
        np.where(merged.p_perm.notna(), "refit", "retained_but_no_result"),
        "removed_by_vifstep",
    )
    return merged


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--contract", type=Path, required=True)
    p.add_argument("--traits-long", type=Path, required=True)
    p.add_argument("--environment", type=Path, required=True)
    p.add_argument("--recovered-display", type=Path, required=True)
    p.add_argument("--frozen-among", type=Path, required=True)
    p.add_argument("--out-dir", type=Path, required=True)
    p.add_argument("--minimum", type=int, default=5)
    p.add_argument("--permutations", type=int, default=9999)
    p.add_argument("--seed", type=int, default=20260909)
    args = p.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    contract, traits, environment, recovery = load_inputs(args)
    env_taxon = taxon_environment(environment)
    units = inferential_units(contract)
    verification = verify_frozen_univariate(units, traits, env_taxon, args.frozen_among, args.minimum)

    global_selection = {}
    for threshold in (10.0, 5.0):
        selected, history, final_vif = vifstep(env_taxon, threshold)
        global_selection[str(int(threshold))] = {
            "selected_predictors": selected,
            "removed_predictors": [row["max_vif_predictor"] for row in history[:-1]],
            "history": history,
            "final_vif": final_vif.to_dict("records"),
        }
        final_vif.assign(vif_threshold=threshold).to_csv(args.out_dir / f"global_vif_final_{int(threshold)}.csv", index=False)

    result10, select10 = run_threshold(10.0, units, traits, env_taxon,
                                       global_selection["10"]["selected_predictors"], args.minimum,
                                       args.permutations, args.seed)
    result5, select5 = run_threshold(5.0, units, traits, env_taxon,
                                     global_selection["5"]["selected_predictors"], args.minimum,
                                     args.permutations, args.seed)
    result10.to_csv(args.out_dir / "vif10_full_multivariable_results.csv", index=False)
    result5.to_csv(args.out_dir / "vif5_full_multivariable_results.csv", index=False)
    select10.to_csv(args.out_dir / "vif10_endpoint_predictor_sets.csv", index=False)
    select5.to_csv(args.out_dir / "vif5_endpoint_predictor_sets.csv", index=False)

    hits10 = result10.loc[result10.posthoc_fdr_significant_0_05].sort_values(["q_perm_bh_posthoc_family", "unit_id", "predictor"])
    hits5 = result5.loc[result5.posthoc_fdr_significant_0_05].sort_values(["q_perm_bh_posthoc_family", "unit_id", "predictor"])
    hits10.to_csv(args.out_dir / "vif10_posthoc_fdr_hits.csv", index=False)
    hits5.to_csv(args.out_dir / "vif5_posthoc_fdr_hits.csv", index=False)

    original_compare = compare_original_hits(args.frozen_among, result10, select10)
    original_compare.to_csv(args.out_dir / "original_fdr_hits_vif10_comparison.csv", index=False)

    report = {
        "status": "FULL_COHORT_VIFSTEP_POSTHOC_COMPLETE",
        "claim_boundary": [
            "Uses the full frozen 46,276-observation / 259-taxon strict-spatial v2 cohort; no native-range filter is used.",
            "VIF filtering is response-blind but retrospective and does not replace the frozen marginal v2 family.",
            "Partial coefficients remain observational and cannot be interpreted as causal environmental effects.",
            "Post-hoc BH q-values belong only to this VIF-filtered sensitivity family.",
            "Restored display/composition fields were already measured historically; no new image operation was performed.",
            "Four colour-composition fractions are a closed composition and are not four independent biological discoveries.",
        ],
        "input_sha256": {
            "traits_long": sha256(args.traits_long), "environment": sha256(args.environment),
            "recovered_display": sha256(args.recovered_display), "frozen_among": sha256(args.frozen_among),
        },
        "cohort": {"observations": len(environment), "taxa": int(environment.taxon_name.nunique()),
                   "taxon_environment_complete": len(env_taxon), "minimum_trait_observations_per_taxon": args.minimum},
        "recovery": recovery,
        "frozen_univariate_reproduction": verification,
        "global_vif_selection": global_selection,
        "permutations": args.permutations,
        "vif10": {
            "successful_tests": int(len(result10)), "posthoc_fdr_hits": int(len(hits10)),
            "significant_rows": hits10[["unit_id", "module", "inferential_unit", "predictor", "n_taxa",
                                        "beta_std_partial", "beta_sine_std_partial", "beta_cosine_std_partial",
                                        "effect_magnitude_partial", "effect_direction_degrees_partial", "p_perm",
                                        "q_perm_bh_posthoc_family"]].to_dict("records"),
        },
        "vif5": {
            "successful_tests": int(len(result5)), "posthoc_fdr_hits": int(len(hits5)),
            "significant_rows": hits5[["unit_id", "module", "inferential_unit", "predictor", "n_taxa",
                                       "beta_std_partial", "beta_sine_std_partial", "beta_cosine_std_partial",
                                       "effect_magnitude_partial", "effect_direction_degrees_partial", "p_perm",
                                       "q_perm_bh_posthoc_family"]].to_dict("records"),
        },
        "original_fdr_hits_vif10": original_compare[[
            "unit_id", "module", "inferential_unit", "predictor", "frozen_q_fdr_bh_global_family",
            "vif10_predictor_retained", "vif10_result_status", "n_taxa", "beta_std_partial",
            "effect_magnitude_partial", "p_perm", "q_perm_bh_posthoc_family"
        ]].replace({np.nan: None}).to_dict("records"),
    }
    (args.out_dir / "full_vifstep_report.json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")

    compact = {
        "status": report["status"],
        "cohort": report["cohort"],
        "vif10_selected": global_selection["10"]["selected_predictors"],
        "vif10_removed": global_selection["10"]["removed_predictors"],
        "vif10_final_max_vif": max(row["vif"] for row in global_selection["10"]["final_vif"]),
        "vif10_posthoc_fdr_hits": report["vif10"]["posthoc_fdr_hits"],
        "vif10_significant": report["vif10"]["significant_rows"],
        "vif5_selected": global_selection["5"]["selected_predictors"],
        "vif5_removed": global_selection["5"]["removed_predictors"],
        "vif5_posthoc_fdr_hits": report["vif5"]["posthoc_fdr_hits"],
        "original_fdr_hits_vif10": report["original_fdr_hits_vif10"],
        "univariate_reproduction": verification,
    }
    print("FULL_VIFSTEP_COMPACT_JSON=" + json.dumps(compact, separators=(",", ":"), allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
