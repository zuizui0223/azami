#!/usr/bin/env python3
"""Validate the frozen Chapter 1 submission output contract.

The validator is intentionally independent of GitHub Actions artifact layout. It
searches a bundle recursively for required filenames, checks dimensions and
scientific invariants, and rejects accidental inclusion of obsolete AI-state
columns or model-weight files.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any

import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--bundle-root", required=True)
    parser.add_argument(
        "--config",
        default=str(Path(__file__).with_name("submission_config.json")),
    )
    parser.add_argument("--out-json", required=True)
    return parser.parse_args()


def find_unique(root: Path, filename: str) -> Path:
    matches = sorted(path for path in root.rglob(filename) if path.is_file())
    if len(matches) != 1:
        raise ValueError(
            f"Expected exactly one {filename!r} below {root}, found {len(matches)}: "
            f"{[str(path) for path in matches[:10]]}"
        )
    return matches[0]


def read_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"Expected JSON object: {path}")
    return value


def bool_series(values: pd.Series) -> pd.Series:
    if values.dtype == bool:
        return values
    return values.astype(str).str.strip().str.lower().isin({"true", "1", "yes"})


def validate_bundle(root: Path, config: dict[str, Any]) -> dict[str, Any]:
    root = root.resolve()
    if not root.is_dir():
        raise ValueError(f"Bundle root is not a directory: {root}")

    required_paths = {
        filename: find_unique(root, filename)
        for filename in config["required_files"]
    }
    checks: list[dict[str, Any]] = []

    def check(name: str, passed: bool, detail: Any) -> None:
        checks.append({"check": name, "passed": bool(passed), "detail": detail})
        if not passed:
            raise ValueError(f"Submission contract failed: {name}: {detail}")

    forbidden_suffixes = {
        suffix.lower() for suffix in config.get("forbidden_submission_suffixes", [])
    }
    forbidden_files = [
        str(path.relative_to(root))
        for path in root.rglob("*")
        if path.is_file() and path.suffix.lower() in forbidden_suffixes
    ]
    check("no_model_weights_in_result_bundle", not forbidden_files, forbidden_files)

    qc = read_json(required_paths["qc_bias_and_join_report.json"])
    observation = pd.read_csv(
        required_paths["integrated_primary_auxiliary_observation.csv"],
        low_memory=False,
    )
    species = pd.read_csv(
        required_paths["integrated_primary_auxiliary_species.csv"],
        low_memory=False,
    )
    registry = pd.read_csv(
        required_paths["integrated_continuous_endpoint_registry.csv"],
        dtype=str,
        keep_default_na=False,
    )
    prephylo = pd.read_csv(
        required_paths["continuous_environment_model_coefficients.csv"],
        low_memory=False,
    )
    tree = read_json(required_paths["historical_tree_audit_report.json"])
    model = read_json(required_paths["historical_constraint_model_report.json"])
    historical = pd.read_csv(
        required_paths["historical_constraint_model_coefficients_clean.csv"],
        low_memory=False,
    )
    random_summary = pd.read_csv(
        required_paths["historical_constraint_random_tree_summary.csv"],
        low_memory=False,
    )
    evidence = pd.read_csv(
        required_paths["ploidy_hybrid_evidence_queue.csv"],
        dtype=str,
        keep_default_na=False,
    )

    expected = config["expected"]
    check("head_count", int(qc.get("n_heads", -1)) == expected["n_heads"], qc.get("n_heads"))
    check("observation_count", len(observation) == expected["n_observations"], len(observation))
    check("taxon_count", len(species) == expected["n_taxa"], len(species))
    check(
        "unique_species_taxa",
        "taxon_name" in species.columns
        and species["taxon_name"].notna().all()
        and species["taxon_name"].is_unique,
        species.get("taxon_name", pd.Series(dtype=str)).nunique(),
    )
    check(
        "observation_ids_unique",
        "obs_id" in observation.columns and observation["obs_id"].astype(str).is_unique,
        observation.get("obs_id", pd.Series(dtype=str)).nunique(),
    )

    tier_counts = registry["analysis_tier"].value_counts().to_dict()
    check(
        "primary_endpoint_count",
        int(tier_counts.get("main", 0)) == expected["n_primary_endpoints"],
        tier_counts,
    )
    check(
        "auxiliary_endpoint_count",
        int(tier_counts.get("auxiliary", 0)) == expected["n_auxiliary_endpoints"],
        tier_counts,
    )

    for row in registry.to_dict("records"):
        observation_column = row["observation_variable"]
        species_column = row["species_variable"]
        check(
            f"observation_endpoint_present:{observation_column}",
            observation_column in observation.columns,
            observation_column,
        )
        check(
            f"species_endpoint_present:{species_column}",
            species_column in species.columns,
            species_column,
        )
        for label, frame, column in (
            ("observation", observation, observation_column),
            ("species", species, species_column),
        ):
            values = pd.to_numeric(frame[column], errors="coerce")
            finite = values.dropna().map(math.isfinite).all()
            check(f"finite_{label}_endpoint:{column}", bool(finite), int(values.notna().sum()))

    forbidden_columns = set(config.get("forbidden_analysis_columns", []))
    leaked = sorted(
        forbidden_columns.intersection(observation.columns)
        | forbidden_columns.intersection(species.columns)
    )
    check("no_legacy_ai_or_holdout_columns", not leaked, leaked)

    required_prephylo = {
        "analysis_tier",
        "model_scale",
        "cohort",
        "estimate_standardized",
        "p_fdr_within_tier_scale_cohort",
    }
    check(
        "prephylogenetic_schema",
        required_prephylo.issubset(prephylo.columns),
        sorted(required_prephylo.difference(prephylo.columns)),
    )
    strict_within = prephylo.loc[
        prephylo["analysis_tier"].eq("main")
        & prephylo["model_scale"].eq("observation_within_species")
        & prephylo["cohort"].eq("positional_accuracy_le_10km")
    ]
    check("strict_within_species_rows_present", len(strict_within) == 36, len(strict_within))

    gbotb = tree.get("gbotb_lcvp", {})
    check("tree_taxon_count", tree.get("n_input_taxa") == expected["n_taxa"], tree.get("n_input_taxa"))
    check(
        "direct_backbone_tip_count",
        gbotb.get("n_direct_backbone_tips") == expected["n_direct_backbone_tips"],
        gbotb.get("n_direct_backbone_tips"),
    )
    check(
        "random_tree_count_tree_report",
        gbotb.get("n_randomized_scenario2_trees") == expected["n_random_trees"],
        gbotb.get("n_randomized_scenario2_trees"),
    )
    check(
        "historical_models_all_succeeded",
        model.get("n_successful_model_fits") == expected["n_historical_model_fits"]
        and model.get("n_failed_model_fits") == 0,
        {
            "success": model.get("n_successful_model_fits"),
            "failed": model.get("n_failed_model_fits"),
        },
    )
    check(
        "historical_endpoint_counts",
        model.get("n_main_endpoints") == expected["n_primary_endpoints"]
        and model.get("n_auxiliary_endpoints") == expected["n_auxiliary_endpoints"],
        {
            "main": model.get("n_main_endpoints"),
            "auxiliary": model.get("n_auxiliary_endpoints"),
        },
    )

    lambda_rows = historical.loc[historical["covariance_model"].eq("PGLS_Pagel_lambda")]
    lambdas = pd.to_numeric(lambda_rows["pagel_lambda"], errors="coerce")
    check(
        "pagel_lambda_bounded",
        lambdas.notna().all() and lambdas.between(0, 1).all(),
        {"min": float(lambdas.min()), "max": float(lambdas.max()), "n": int(len(lambdas))},
    )
    check(
        "random_summary_row_count",
        len(random_summary) == expected["n_random_tree_endpoint_predictor_rows"],
        len(random_summary),
    )
    successful_random = pd.to_numeric(
        random_summary["n_successful_random_trees"], errors="coerce"
    )
    expected_random = pd.to_numeric(
        random_summary["n_expected_random_trees"], errors="coerce"
    )
    complete_support = bool_series(random_summary["complete_random_tree_support"])
    check(
        "all_random_tree_models_complete",
        successful_random.eq(expected["n_random_trees"]).all()
        and expected_random.eq(expected["n_random_trees"]).all()
        and complete_support.all(),
        {
            "minimum_successful": int(successful_random.min()),
            "maximum_successful": int(successful_random.max()),
        },
    )

    check("evidence_queue_taxon_count", len(evidence) == expected["n_taxa"], len(evidence))
    check(
        "evidence_queue_unique_taxa",
        "accepted_taxon" in evidence.columns and evidence["accepted_taxon"].is_unique,
        evidence.get("accepted_taxon", pd.Series(dtype=str)).nunique(),
    )
    ploidy_values = set(evidence.get("ploidy_class", pd.Series(dtype=str)).astype(str))
    hybrid_values = set(evidence.get("hybrid_status", pd.Series(dtype=str)).astype(str))
    check(
        "allowed_ploidy_classes",
        ploidy_values.issubset(set(config["allowed_ploidy_classes"])),
        sorted(ploidy_values),
    )
    check(
        "allowed_hybrid_classes",
        hybrid_values.issubset(set(config["allowed_hybrid_classes"])),
        sorted(hybrid_values),
    )

    return {
        "analysis_id": config["analysis_id"],
        "analysis_version": config["analysis_version"],
        "bundle_root": str(root),
        "status": "pass",
        "n_checks": len(checks),
        "checks": checks,
        "resolved_files": {
            name: str(path.relative_to(root)) for name, path in required_paths.items()
        },
    }


def main() -> None:
    args = parse_args()
    config = read_json(Path(args.config))
    report = validate_bundle(Path(args.bundle_root), config)
    target = Path(args.out_json)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps({
        "status": report["status"],
        "analysis_version": report["analysis_version"],
        "n_checks": report["n_checks"],
        "out_json": str(target),
    }, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
