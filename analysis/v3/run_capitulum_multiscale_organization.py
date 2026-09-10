#!/usr/bin/env python3
"""Compare whole-capitulum phenotype organization within and among species.

This v3 lane treats the capitulum as an integrated multivariate phenotype without
assuming that it is one syndrome or imposing a fixed number of syndromes.
Primary + candidate image-derived endpoints are analysed jointly.  The recovered
visible-floret fraction is added to the historical complete-18 phenotype, giving
19 endpoint dimensions / 18 inferential units because hue sine/cosine are one
circular unit.

Questions:
1. How strongly are traits organized within species versus among species?
2. Is the pairwise organization geometry shared across biological scales?
3. Do predeclared environmental blocks organize the whole phenotype differently
   within versus among species?

The analysis is observational.  Correlation/modularity patterns are not genetic
syndromes, and environmental blocks are contexts rather than causal mechanisms.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from analysis import run_capitulum_environment_blocks as envmod
from analysis import run_capitulum_functional_space as spacemod

HUE_COMPONENTS = ["corolla_hue_sin", "corolla_hue_cos"]
HUE_UNIT = "corolla_hue"
RECOVERED_ENDPOINT = "visible_floret_fraction"
RECOVERED_SOURCE = "corolla_visible_fraction"


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--traits-long", required=True, type=Path)
    p.add_argument("--endpoint-contract", required=True, type=Path)
    p.add_argument("--recovered-display", required=True, type=Path)
    p.add_argument("--environment", required=True, type=Path)
    p.add_argument("--block-contract", required=True, type=Path)
    p.add_argument("--out-dir", required=True, type=Path)
    p.add_argument("--thresholds", nargs="+", type=int, default=[5, 2])
    p.add_argument("--matrix-permutations", type=int, default=9999)
    p.add_argument("--module-permutations", type=int, default=9999)
    p.add_argument("--environment-permutations", type=int, default=9999)
    p.add_argument("--structure-bootstrap", type=int, default=1000)
    p.add_argument("--environment-bootstrap", type=int, default=500)
    p.add_argument("--seed", type=int, default=20260910)
    return p.parse_args()


def as_bool(series: pd.Series) -> pd.Series:
    return series.astype(str).str.strip().str.lower().isin({"true", "1", "yes"})


def selected_contract(path: Path) -> tuple[pd.DataFrame, list[str], dict[str, str]]:
    contract = pd.read_csv(path, dtype=str, keep_default_na=False, low_memory=False)
    selected = contract[contract["analysis_tier"].isin(["primary", "candidate"])].copy()
    endpoint_ids = sorted(selected["endpoint_id"].tolist())
    if len(endpoint_ids) != 19:
        raise ValueError(f"Expected 19 primary+candidate endpoint dimensions, found {len(endpoint_ids)}")
    if not set(HUE_COMPONENTS).issubset(endpoint_ids):
        raise ValueError("Hue sine/cosine components are missing from the selected contract")
    if RECOVERED_ENDPOINT not in endpoint_ids:
        raise ValueError("Recovered visible-floret endpoint is not in the primary+candidate contract")
    module_map = selected.set_index("endpoint_id")["module"].to_dict()
    return selected, endpoint_ids, module_map


def materialize_complete19(
    traits_path: Path,
    recovered_path: Path,
    environment: pd.DataFrame,
    endpoint_ids: list[str],
    module_map: dict[str, str],
) -> tuple[pd.DataFrame, dict[str, Any]]:
    columns = [
        "obs_id", "taxon_name", "endpoint_id", "module", "analysis_tier",
        "measurement_available", "value",
    ]
    traits = pd.read_csv(traits_path, usecols=columns, low_memory=False)
    traits["obs_id"] = traits["obs_id"].astype(str)
    traits["taxon_name"] = traits["taxon_name"].astype(str)
    traits["value"] = pd.to_numeric(traits["value"], errors="coerce")
    traits = traits[
        traits["endpoint_id"].isin(endpoint_ids)
        & as_bool(traits["measurement_available"])
        & traits["value"].notna()
    ].copy()

    environment_ids = set(environment["obs_id"].astype(str))
    traits = traits[traits["obs_id"].isin(environment_ids)].copy()
    if traits.duplicated(["obs_id", "endpoint_id"]).any():
        raise ValueError("Historical trait rows are not unique by obs_id/endpoint_id")

    recovered = pd.read_csv(recovered_path, low_memory=False)
    recovered["obs_id"] = recovered["obs_id"].astype(str)
    if RECOVERED_SOURCE not in recovered.columns:
        raise ValueError(f"Recovered table lacks {RECOVERED_SOURCE}")
    recovered[RECOVERED_SOURCE] = pd.to_numeric(recovered[RECOVERED_SOURCE], errors="coerce")
    recovered = recovered[
        recovered["obs_id"].isin(environment_ids) & recovered[RECOVERED_SOURCE].notna()
    ][["obs_id", RECOVERED_SOURCE]].copy()
    if recovered["obs_id"].duplicated().any():
        raise ValueError("Recovered display table is not unique by obs_id")

    existing = set(traits.loc[traits["endpoint_id"].eq(RECOVERED_ENDPOINT), "obs_id"])
    recovered = recovered[~recovered["obs_id"].isin(existing)].copy()
    taxon_lookup = environment.assign(obs_id=environment["obs_id"].astype(str)).set_index("obs_id")["taxon_name"]
    add = pd.DataFrame({
        "obs_id": recovered["obs_id"],
        "taxon_name": recovered["obs_id"].map(taxon_lookup),
        "endpoint_id": RECOVERED_ENDPOINT,
        "module": module_map[RECOVERED_ENDPOINT],
        "analysis_tier": "candidate",
        "measurement_available": True,
        "value": recovered[RECOVERED_SOURCE].to_numpy(float),
    })
    add = add[add["taxon_name"].notna()].copy()
    traits = pd.concat([traits, add[columns]], ignore_index=True, sort=False)
    if traits.duplicated(["obs_id", "endpoint_id"]).any():
        raise ValueError("Materialized trait rows are not unique by obs_id/endpoint_id")

    wide = traits.pivot(index=["obs_id", "taxon_name"], columns="endpoint_id", values="value")
    present = [x for x in endpoint_ids if x in wide.columns]
    if present != endpoint_ids:
        missing = sorted(set(endpoint_ids) - set(present))
        raise ValueError(f"Selected endpoint dimensions absent after materialization: {missing}")
    complete = wide[endpoint_ids].dropna().reset_index()
    report = {
        "environment_rows_supplied": int(len(environment)),
        "environment_taxa_supplied": int(environment["taxon_name"].nunique()),
        "selected_endpoint_dimensions": len(endpoint_ids),
        "recovered_visible_floret_rows_added": int(len(add)),
        "complete19_observations": int(len(complete)),
        "complete19_taxa": int(complete["taxon_name"].nunique()),
    }
    return complete, report


def linear_endpoint_ids(endpoint_ids: list[str]) -> list[str]:
    return [x for x in endpoint_ids if x not in HUE_COMPONENTS]


def signed_upper(corr: pd.DataFrame, names: list[str]) -> np.ndarray:
    a = corr.loc[names, names].to_numpy(float)
    return a[np.triu_indices_from(a, k=1)]


def rms_matrix_difference(a: np.ndarray, b: np.ndarray) -> float:
    return float(np.sqrt(np.mean((np.asarray(a, float) - np.asarray(b, float)) ** 2)))


def sign_agreement(a: np.ndarray, b: np.ndarray) -> float:
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    usable = np.isfinite(a) & np.isfinite(b) & (np.abs(a) > 1e-12) & (np.abs(b) > 1e-12)
    if not usable.any():
        return float("nan")
    return float(np.mean(np.sign(a[usable]) == np.sign(b[usable])))


def effective_dimensionality(corr: pd.DataFrame) -> tuple[float, float, float]:
    values = np.linalg.eigvalsh(corr.to_numpy(float))[::-1]
    values = np.clip(values, 0.0, None)
    total = float(values.sum())
    if total <= 0:
        return float("nan"), float("nan"), float("nan")
    frac = values / total
    participation = float(1.0 / np.sum(frac ** 2))
    return participation, float(frac[0]), float(np.cumsum(frac)[min(2, len(frac) - 1)])


def signed_matrix_label_permutation_p(
    within_corr: pd.DataFrame,
    among_corr: pd.DataFrame,
    names: list[str],
    permutations: int,
    rng: np.random.Generator,
) -> float:
    a = within_corr.loc[names, names].to_numpy(float)
    b = among_corr.loc[names, names].to_numpy(float)
    idx = np.triu_indices_from(a, k=1)
    observed = spacemod.spearman(a[idx], b[idx])
    exceed = 0
    for _ in range(permutations):
        perm = rng.permutation(len(names))
        bp = b[np.ix_(perm, perm)]
        if spacemod.spearman(a[idx], bp[idx]) >= observed - 1e-15:
            exceed += 1
    return float((exceed + 1) / (permutations + 1))


def bootstrap_structure(
    table: pd.DataFrame,
    endpoint_ids: list[str],
    units: list[str],
    module_by_unit: dict[str, str],
    n_boot: int,
    seed: int,
    scope: str,
) -> dict[str, Any]:
    taxa = table["taxon_name"].drop_duplicates().tolist()
    rng = spacemod.stable_rng(seed, scope, "v3_multiscale_bootstrap")
    linear = linear_endpoint_ids(endpoint_ids)
    strength_similarity: list[float] = []
    signed_similarity: list[float] = []
    module_difference: list[float] = []
    rms_difference: list[float] = []
    for _ in range(n_boot):
        sampled = rng.choice(taxa, size=len(taxa), replace=True)
        chunks = []
        for j, taxon in enumerate(sampled):
            g = table[table["taxon_name"].eq(taxon)].copy()
            g["taxon_name"] = f"boot_{j:04d}"
            chunks.append(g)
        b = pd.concat(chunks, ignore_index=True)
        try:
            wc, ac, wu, au = spacemod.matrices_for_table(b, endpoint_ids, units)
            strength_similarity.append(spacemod.spearman(spacemod.upper_values(wu), spacemod.upper_values(au)))
            ws = signed_upper(wc, linear)
            ass = signed_upper(ac, linear)
            signed_similarity.append(spacemod.spearman(ws, ass))
            module_difference.append(
                spacemod.module_contrast(wu, module_by_unit)[0]
                - spacemod.module_contrast(au, module_by_unit)[0]
            )
            rms_difference.append(rms_matrix_difference(ws, ass))
        except Exception:
            continue

    def ci(values: list[float]) -> list[float]:
        arr = np.asarray(values, float)
        if len(arr) == 0:
            return [float("nan"), float("nan")]
        return [float(np.quantile(arr, 0.025)), float(np.quantile(arr, 0.975))]

    return {
        "successful_replicates": int(len(strength_similarity)),
        "strength_similarity_ci95": ci(strength_similarity),
        "signed_similarity_ci95": ci(signed_similarity),
        "within_minus_among_module_contrast_ci95": ci(module_difference),
        "signed_matrix_rms_difference_ci95": ci(rms_difference),
    }


def structure_analysis(
    complete: pd.DataFrame,
    endpoint_ids: list[str],
    module_map: dict[str, str],
    thresholds: list[int],
    args: argparse.Namespace,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    units = spacemod.unit_names(endpoint_ids)
    if len(units) != 18:
        raise ValueError(f"Expected 18 inferential units from 19 endpoint dimensions, found {len(units)}")
    module_by_unit = {u: spacemod.unit_module(u, module_map) for u in units}
    linear = linear_endpoint_ids(endpoint_ids)

    summaries: list[dict[str, Any]] = []
    module_rows: list[dict[str, Any]] = []
    endpoint_mats: list[pd.DataFrame] = []
    unit_mats: list[pd.DataFrame] = []
    eigens: list[pd.DataFrame] = []

    for threshold in thresholds:
        counts = complete.groupby("taxon_name").size()
        keep = counts[counts >= threshold].index
        table = complete[complete["taxon_name"].isin(keep)].copy()
        scope = f"complete19_min{threshold}"
        wc, ac, wu, au = spacemod.matrices_for_table(table, endpoint_ids, units)

        wmc, w_in, w_between = spacemod.module_contrast(wu, module_by_unit)
        amc, a_in, a_between = spacemod.module_contrast(au, module_by_unit)
        strength_similarity = spacemod.spearman(spacemod.upper_values(wu), spacemod.upper_values(au))
        strength_p = spacemod.matrix_permutation_p(
            wu, au, args.matrix_permutations,
            spacemod.stable_rng(args.seed, scope, "strength_matrix_label_perm"),
        )
        ws = signed_upper(wc, linear)
        ass = signed_upper(ac, linear)
        signed_similarity = spacemod.spearman(ws, ass)
        signed_p = signed_matrix_label_permutation_p(
            wc, ac, linear, args.matrix_permutations,
            spacemod.stable_rng(args.seed, scope, "signed_matrix_label_perm"),
        )
        rms_diff = rms_matrix_difference(ws, ass)
        sign_same = sign_agreement(ws, ass)
        wdim, wpc1, wpc3 = effective_dimensionality(wc)
        adim, apc1, apc3 = effective_dimensionality(ac)
        boot = bootstrap_structure(
            table, endpoint_ids, units, module_by_unit,
            args.structure_bootstrap, args.seed, scope,
        )

        summaries.append({
            "scope": scope,
            "min_complete_observations_per_taxon": threshold,
            "n_complete_observations": int(len(table)),
            "n_taxa": int(table["taxon_name"].nunique()),
            "n_endpoint_dimensions": len(endpoint_ids),
            "n_inferential_units": len(units),
            "within_module_contrast": wmc,
            "among_module_contrast": amc,
            "within_minus_among_module_contrast": wmc - amc,
            "association_strength_matrix_spearman": strength_similarity,
            "association_strength_trait_label_permutation_p": strength_p,
            "signed_linear_matrix_spearman": signed_similarity,
            "signed_linear_trait_label_permutation_p": signed_p,
            "signed_linear_pair_sign_agreement": sign_same,
            "signed_linear_matrix_rms_difference": rms_diff,
            "within_effective_dimensionality": wdim,
            "among_effective_dimensionality": adim,
            "within_pc1_variance_fraction": wpc1,
            "among_pc1_variance_fraction": apc1,
            "within_pc1_to_pc3_cumulative_fraction": wpc3,
            "among_pc1_to_pc3_cumulative_fraction": apc3,
            "strength_similarity_ci95_low": boot["strength_similarity_ci95"][0],
            "strength_similarity_ci95_high": boot["strength_similarity_ci95"][1],
            "signed_similarity_ci95_low": boot["signed_similarity_ci95"][0],
            "signed_similarity_ci95_high": boot["signed_similarity_ci95"][1],
            "within_minus_among_module_contrast_ci95_low": boot["within_minus_among_module_contrast_ci95"][0],
            "within_minus_among_module_contrast_ci95_high": boot["within_minus_among_module_contrast_ci95"][1],
            "signed_matrix_rms_difference_ci95_low": boot["signed_matrix_rms_difference_ci95"][0],
            "signed_matrix_rms_difference_ci95_high": boot["signed_matrix_rms_difference_ci95"][1],
            "bootstrap_successful_replicates": boot["successful_replicates"],
        })

        for scale, matrix, contrast, mean_in, mean_between in [
            ("within_taxon", wu, wmc, w_in, w_between),
            ("among_taxon", au, amc, a_in, a_between),
        ]:
            module_rows.append({
                "scope": scope,
                "scale": scale,
                "module_contrast": contrast,
                "mean_within_module_strength": mean_in,
                "mean_between_module_strength": mean_between,
                "module_label_permutation_p": spacemod.module_permutation_p(
                    matrix, module_by_unit, args.module_permutations,
                    spacemod.stable_rng(args.seed, scope, scale, "module_label_perm"),
                ),
            })

        endpoint_mats += [
            spacemod.long_matrix(wc, scope, "within_taxon", "signed_endpoint_correlation"),
            spacemod.long_matrix(ac, scope, "among_taxon", "signed_endpoint_correlation"),
        ]
        unit_mats += [
            spacemod.long_matrix(wu, scope, "within_taxon", "inferential_unit_association_strength"),
            spacemod.long_matrix(au, scope, "among_taxon", "inferential_unit_association_strength"),
        ]
        eigens += [
            spacemod.eigenspectrum(wc, scope, "within_taxon"),
            spacemod.eigenspectrum(ac, scope, "among_taxon"),
        ]

    return (
        pd.DataFrame(summaries),
        pd.DataFrame(module_rows),
        pd.concat(endpoint_mats, ignore_index=True),
        pd.concat(unit_mats, ignore_index=True),
        pd.concat(eigens, ignore_index=True),
    )


def environment_analysis(
    complete: pd.DataFrame,
    environment: pd.DataFrame,
    endpoint_ids: list[str],
    module_map: dict[str, str],
    block_contract: dict[str, Any],
    thresholds: list[int],
    args: argparse.Namespace,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    env = environment.copy()
    env["obs_id"] = env["obs_id"].astype(str)
    merged = complete.merge(env.drop(columns=["taxon_name"], errors="ignore"), on="obs_id", how="left", validate="one_to_one")
    tests: list[dict[str, Any]] = []
    coeffs: list[dict[str, Any]] = []
    energies: list[dict[str, Any]] = []
    geometry: list[dict[str, Any]] = []

    for threshold in thresholds:
        counts = merged.groupby("taxon_name").size()
        keep = counts[counts >= threshold].index
        base = merged[merged["taxon_name"].isin(keep)].copy()
        scope = f"complete19_env_min{threshold}"
        for block in block_contract["blocks"]:
            predictors = list(block["predictors"])
            missing = [x for x in predictors if x not in base.columns]
            if missing:
                continue
            work = base.dropna(subset=predictors).copy()
            coverage = len(work) / len(base) if len(base) else 0.0
            if coverage < float(block_contract.get("minimum_environment_coverage", 0.98)):
                continue

            yw, xw, weights, groups = envmod.within_arrays(work, endpoint_ids, predictors)
            rw, bw = envmod.fit_weighted_multivariate(yw, xw, weights)
            pw = envmod.permutation_p_within(
                yw, xw, weights, groups, rw, args.environment_permutations,
                envmod.stable_rng(args.seed, scope, block["block_id"], "within_perm"),
            )
            ya, xa, taxa = envmod.among_arrays(work, endpoint_ids, predictors)
            ra, ba = envmod.fit_unweighted_multivariate(ya, xa)
            pa = envmod.permutation_p_among(
                ya, xa, ra, args.environment_permutations,
                envmod.stable_rng(args.seed, scope, block["block_id"], "among_perm"),
            )

            for scale, nobs, r2, pval, beta in [
                ("within_taxon", len(work), rw, pw, bw),
                ("among_taxon", len(taxa), ra, pa, ba),
            ]:
                tests.append({
                    "scope": scope,
                    "block_id": block["block_id"],
                    "tier": block["tier"],
                    "construct": block["construct"],
                    "scale": scale,
                    "n_observations_or_taxa": int(nobs),
                    "n_taxa": int(work["taxon_name"].nunique()),
                    "n_predictors": len(predictors),
                    "n_response_endpoint_dimensions": len(endpoint_ids),
                    "multivariate_r2": r2,
                    "permutation_p": pval,
                })
                coeffs += envmod.coefficient_rows(
                    scope, block, scale, predictors, endpoint_ids, module_map, beta
                )
                energies += envmod.module_energy_rows(
                    scope, block, scale, endpoint_ids, module_map, beta
                )

            cosine = envmod.cosine(envmod.flatten_beta(bw), envmod.flatten_beta(ba))
            low, high, nboot = envmod.bootstrap_cosine(
                work, endpoint_ids, predictors, args.environment_bootstrap,
                args.seed, scope, block["block_id"],
            )
            geometry.append({
                "scope": scope,
                "block_id": block["block_id"],
                "tier": block["tier"],
                "coefficient_matrix_cosine_within_vs_among": cosine,
                "bootstrap_ci95_low": low,
                "bootstrap_ci95_high": high,
                "bootstrap_successful_replicates": int(nboot),
                "within_r2": rw,
                "among_r2": ra,
                "among_minus_within_r2": ra - rw,
            })

    testdf = pd.DataFrame(tests)
    if not testdf.empty:
        testdf["q_bh_across_available_blocks"] = np.nan
        for (_, scale), idx in testdf.groupby(["scope", "scale"]).groups.items():
            testdf.loc[idx, "q_bh_across_available_blocks"] = envmod.bh_adjust(
                testdf.loc[idx, "permutation_p"].astype(float)
            )
        testdf["fdr_supported_0_05"] = testdf["q_bh_across_available_blocks"].lt(0.05)
    return testdf, pd.DataFrame(coeffs), pd.DataFrame(energies), pd.DataFrame(geometry)


def main() -> int:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    selected, endpoint_ids, module_map = selected_contract(args.endpoint_contract)
    environment = pd.read_csv(args.environment, low_memory=False)
    environment["obs_id"] = environment["obs_id"].astype(str)
    environment["taxon_name"] = environment["taxon_name"].astype(str)
    if environment["obs_id"].duplicated().any():
        raise ValueError("Environment input must be unique by obs_id")

    complete, cohort_report = materialize_complete19(
        args.traits_long, args.recovered_display, environment,
        endpoint_ids, module_map,
    )
    if len(complete) < 1000:
        raise ValueError(f"Unexpectedly small complete19 cohort: {len(complete)}")

    structure, module_integration, endpoint_matrix, unit_matrix, eigens = structure_analysis(
        complete, endpoint_ids, module_map, args.thresholds, args
    )
    block_contract = json.loads(args.block_contract.read_text(encoding="utf-8"))
    env_tests, env_coeffs, env_energy, env_geometry = environment_analysis(
        complete, environment, endpoint_ids, module_map,
        block_contract, args.thresholds, args,
    )

    structure.to_csv(args.out_dir / "capitulum_multiscale_structure_summary.csv", index=False)
    module_integration.to_csv(args.out_dir / "capitulum_multiscale_module_integration.csv", index=False)
    endpoint_matrix.to_csv(args.out_dir / "capitulum_multiscale_endpoint_correlations.csv", index=False)
    unit_matrix.to_csv(args.out_dir / "capitulum_multiscale_unit_strengths.csv", index=False)
    eigens.to_csv(args.out_dir / "capitulum_multiscale_eigenspectra.csv", index=False)
    env_tests.to_csv(args.out_dir / "capitulum_multiscale_environment_blocks.csv", index=False)
    env_coeffs.to_csv(args.out_dir / "capitulum_multiscale_environment_coefficients.csv", index=False)
    env_energy.to_csv(args.out_dir / "capitulum_multiscale_environment_module_energy.csv", index=False)
    env_geometry.to_csv(args.out_dir / "capitulum_multiscale_environment_cross_scale_geometry.csv", index=False)

    report = {
        "analysis_id": "ch1_capitulum_multiscale_organization_v1",
        "status": "whole_capitulum_multilevel_pattern_first_reanalysis",
        "cohort": cohort_report,
        "endpoint_dimensions": endpoint_ids,
        "inferential_units": spacemod.unit_names(endpoint_ids),
        "selection_rule": "primary + candidate endpoint dimensions only; descriptive_only and validation_only excluded",
        "syndrome_assumption": "none; no clustering or fixed syndrome count imposed",
        "within_rule": "taxon-centred endpoint values with inverse taxon sample-size weights; equal total taxon weight",
        "among_rule": "taxon medians",
        "environment_rule": "six frozen 2026-08-27 environmental blocks tested separately at within- and among-taxon scales",
        "claim_boundary": "phenotypic organization and environmental association geometry only; not genetic syndrome structure, adaptation, or causal mechanism",
        "thresholds": args.thresholds,
        "matrix_permutations": args.matrix_permutations,
        "module_permutations": args.module_permutations,
        "environment_permutations": args.environment_permutations,
        "structure_bootstrap": args.structure_bootstrap,
        "environment_bootstrap": args.environment_bootstrap,
        "seed": args.seed,
    }
    (args.out_dir / "capitulum_multiscale_report.json").write_text(
        json.dumps(report, indent=2), encoding="utf-8"
    )
    print(json.dumps({
        "cohort": cohort_report,
        "structure": structure.to_dict("records"),
        "environment_supported": env_tests[env_tests.get("fdr_supported_0_05", False).eq(True)].to_dict("records") if not env_tests.empty else [],
        "environment_geometry": env_geometry.to_dict("records"),
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
