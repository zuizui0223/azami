#!/usr/bin/env python3
"""Run a conservative climate-association atlas for automated Cirsium image traits.

This is deliberately an association and validation workflow, not a causal or
phylogenetic-comparative analysis. Each focal trait state is modelled against
all other analyzable states of the same trait. Results are emitted for both the
all-ensemble measurement mode and the strict-consensus mode.

Model families:
  1. climate_only: standardized CHELSA predictors only; evaluated by 10-degree
     spatial-block folds and by held-out-species folds.
  2. species_adjusted_within_species: CHELSA predictors plus species dummies,
     restricted to species with observed within-species outcome variation;
     evaluated only with spatial-block folds and only on test rows from
     training-represented variable species.

No p-values, causal language, or random-split headline metrics are produced.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import balanced_accuracy_score, roc_auc_score
from sklearn.model_selection import GroupKFold
from sklearn.preprocessing import StandardScaler


REQUIRED_COLUMNS = {
    "obs_id",
    "trait_id",
    "taxon_name",
    "environment_primary_complete",
    "spatial_block_10deg",
    "coordinate_precision_tier",
    "observation_ai_all_state",
    "observation_ai_conservative_state",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run spatial and species-aware climate association atlas for automated trait states.")
    parser.add_argument("--environment-long", required=True)
    parser.add_argument("--environment-provenance", required=True)
    parser.add_argument("--plan", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--n-folds", type=int, default=5)
    parser.add_argument("--seed", type=int, default=20260706)
    return parser.parse_args()


def text(value: Any) -> str:
    return "" if pd.isna(value) else str(value).strip()


def as_bool(value: Any) -> bool:
    return text(value).lower() in {"true", "1", "yes"}


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def require_columns(frame: pd.DataFrame, columns: Iterable[str], label: str) -> None:
    missing = sorted(set(columns).difference(frame.columns))
    if missing:
        raise ValueError(f"{label} lacks required columns: {missing}")


def finite_frame(frame: pd.DataFrame, predictors: list[str]) -> pd.DataFrame:
    work = frame.copy()
    for predictor in predictors:
        work[predictor] = pd.to_numeric(work[predictor], errors="coerce")
    valid = np.isfinite(work[predictors].to_numpy(dtype=float)).all(axis=1)
    return work.loc[valid].copy()


def fold_iterator(groups: pd.Series, n_requested: int):
    unique = pd.Index(groups.astype(str).unique())
    n_splits = min(int(n_requested), len(unique))
    if n_splits < 2:
        return None, 0
    return GroupKFold(n_splits=n_splits), n_splits


def make_classifier(C: float, class_weight: str | None, seed: int) -> LogisticRegression:
    return LogisticRegression(
        penalty="l2",
        C=float(C),
        class_weight=class_weight,
        solver="lbfgs",
        max_iter=5000,
        random_state=int(seed),
    )


def full_climate_coefficients(
    frame: pd.DataFrame,
    predictors: list[str],
    y: np.ndarray,
    C: float,
    class_weight: str | None,
    seed: int,
    common: dict[str, Any],
) -> tuple[list[dict[str, Any]], str]:
    if len(np.unique(y)) < 2:
        return [], "outcome_has_one_class"
    try:
        scaler = StandardScaler()
        x = scaler.fit_transform(frame[predictors].to_numpy(dtype=float))
        classifier = make_classifier(C, class_weight, seed)
        classifier.fit(x, y)
    except Exception as error:  # explicit output status is more useful than a silent missing model
        return [], f"fit_error:{type(error).__name__}:{error}"
    rows: list[dict[str, Any]] = []
    for predictor, coefficient in zip(predictors, classifier.coef_[0].tolist()):
        rows.append({
            **common,
            "model_family": "climate_only",
            "predictor": predictor,
            "coefficient_standardized": float(coefficient),
            "odds_ratio_per_one_sd_penalized": float(math.exp(float(coefficient))),
            "n_observations_fit": int(len(frame)),
            "n_species_fit": int(frame["taxon_name"].nunique()),
            "n_spatial_blocks_fit": int(frame["spatial_block_10deg"].nunique()),
            "interpretation": "Penalized descriptive association coefficient; not a causal effect estimate.",
        })
    return rows, "success"


def within_species_subset(frame: pd.DataFrame, y: np.ndarray) -> tuple[pd.DataFrame, np.ndarray, list[str]]:
    work = frame.copy()
    work["_outcome"] = y.astype(int)
    variation = work.groupby("taxon_name", observed=True)["_outcome"].nunique()
    species = sorted(variation.loc[variation >= 2].index.astype(str).tolist())
    restricted = work.loc[work["taxon_name"].astype(str).isin(species)].copy()
    return restricted.drop(columns=["_outcome"]), restricted["_outcome"].to_numpy(dtype=int), species


def full_species_adjusted_coefficients(
    frame: pd.DataFrame,
    predictors: list[str],
    y: np.ndarray,
    C: float,
    class_weight: str | None,
    seed: int,
    common: dict[str, Any],
) -> tuple[list[dict[str, Any]], str, int]:
    restricted, y_restricted, variable_species = within_species_subset(frame, y)
    if len(variable_species) < 2 or len(np.unique(y_restricted)) < 2:
        return [], "insufficient_within_species_variation", len(variable_species)
    try:
        scaler = StandardScaler()
        climate = scaler.fit_transform(restricted[predictors].to_numpy(dtype=float))
        dummies = pd.get_dummies(restricted["taxon_name"].astype(str), drop_first=True, dtype=float)
        x = np.column_stack([climate, dummies.to_numpy(dtype=float)])
        classifier = make_classifier(C, class_weight, seed)
        classifier.fit(x, y_restricted)
    except Exception as error:
        return [], f"fit_error:{type(error).__name__}:{error}", len(variable_species)
    rows: list[dict[str, Any]] = []
    for predictor, coefficient in zip(predictors, classifier.coef_[0][: len(predictors)].tolist()):
        rows.append({
            **common,
            "model_family": "species_adjusted_within_species",
            "predictor": predictor,
            "coefficient_standardized": float(coefficient),
            "odds_ratio_per_one_sd_penalized": float(math.exp(float(coefficient))),
            "n_observations_fit": int(len(restricted)),
            "n_species_fit": int(len(variable_species)),
            "n_spatial_blocks_fit": int(restricted["spatial_block_10deg"].nunique()),
            "interpretation": "Species-dummy-adjusted within-observed-species association diagnostic; not a phylogenetic comparative estimate.",
        })
    return rows, "success", len(variable_species)


def cv_climate_only(
    frame: pd.DataFrame,
    predictors: list[str],
    y: np.ndarray,
    groups: pd.Series,
    group_scheme: str,
    n_requested: int,
    C: float,
    class_weight: str | None,
    seed: int,
    common: dict[str, Any],
) -> list[dict[str, Any]]:
    splitter, n_splits = fold_iterator(groups, n_requested)
    if splitter is None:
        return [{**common, "model_family": "climate_only", "validation_scheme": group_scheme, "fold": "", "status": "insufficient_groups", "n_groups": int(groups.astype(str).nunique())}]
    rows: list[dict[str, Any]] = []
    x_all = frame[predictors].to_numpy(dtype=float)
    groups_array = groups.astype(str).to_numpy()
    for fold, (train_idx, test_idx) in enumerate(splitter.split(x_all, y, groups_array), start=1):
        y_train, y_test = y[train_idx], y[test_idx]
        base = {
            **common,
            "model_family": "climate_only",
            "validation_scheme": group_scheme,
            "fold": int(fold),
            "n_folds": int(n_splits),
            "n_train": int(len(train_idx)),
            "n_test": int(len(test_idx)),
            "n_train_groups": int(pd.Series(groups_array[train_idx]).nunique()),
            "n_test_groups": int(pd.Series(groups_array[test_idx]).nunique()),
            "n_test_positive": int(y_test.sum()),
            "positive_prevalence_test": float(y_test.mean()) if len(y_test) else np.nan,
        }
        if len(np.unique(y_train)) < 2:
            rows.append({**base, "status": "unscorable_train_one_class", "roc_auc": np.nan, "balanced_accuracy": np.nan})
            continue
        if len(np.unique(y_test)) < 2:
            rows.append({**base, "status": "unscorable_test_one_class", "roc_auc": np.nan, "balanced_accuracy": np.nan})
            continue
        try:
            scaler = StandardScaler()
            x_train = scaler.fit_transform(x_all[train_idx])
            x_test = scaler.transform(x_all[test_idx])
            classifier = make_classifier(C, class_weight, seed + fold)
            classifier.fit(x_train, y_train)
            probability = classifier.predict_proba(x_test)[:, 1]
            prediction = (probability >= 0.5).astype(int)
            rows.append({
                **base,
                "status": "success",
                "roc_auc": float(roc_auc_score(y_test, probability)),
                "balanced_accuracy": float(balanced_accuracy_score(y_test, prediction)),
            })
        except Exception as error:
            rows.append({**base, "status": f"fit_error:{type(error).__name__}:{error}", "roc_auc": np.nan, "balanced_accuracy": np.nan})
    return rows


def cv_species_adjusted_spatial(
    frame: pd.DataFrame,
    predictors: list[str],
    y: np.ndarray,
    groups: pd.Series,
    n_requested: int,
    C: float,
    class_weight: str | None,
    seed: int,
    common: dict[str, Any],
) -> list[dict[str, Any]]:
    splitter, n_splits = fold_iterator(groups, n_requested)
    if splitter is None:
        return [{**common, "model_family": "species_adjusted_within_species", "validation_scheme": "spatial_block_10deg_within_species", "fold": "", "status": "insufficient_groups", "n_groups": int(groups.astype(str).nunique())}]
    rows: list[dict[str, Any]] = []
    groups_array = groups.astype(str).to_numpy()
    x_all = frame[predictors].to_numpy(dtype=float)
    taxa = frame["taxon_name"].astype(str).to_numpy()
    for fold, (train_idx, test_idx) in enumerate(splitter.split(x_all, y, groups_array), start=1):
        train = frame.iloc[train_idx].copy()
        train_y = y[train_idx]
        test = frame.iloc[test_idx].copy()
        test_y = y[test_idx]
        train["_outcome"] = train_y
        variation = train.groupby("taxon_name", observed=True)["_outcome"].nunique()
        variable_species = sorted(variation.loc[variation >= 2].index.astype(str).tolist())
        train_keep = np.isin(taxa[train_idx], variable_species)
        test_keep = np.isin(taxa[test_idx], variable_species)
        base = {
            **common,
            "model_family": "species_adjusted_within_species",
            "validation_scheme": "spatial_block_10deg_within_species",
            "fold": int(fold),
            "n_folds": int(n_splits),
            "n_train": int(train_keep.sum()),
            "n_test": int(test_keep.sum()),
            "n_train_groups": int(pd.Series(groups_array[train_idx][train_keep]).nunique()) if train_keep.any() else 0,
            "n_test_groups": int(pd.Series(groups_array[test_idx][test_keep]).nunique()) if test_keep.any() else 0,
            "n_variable_species_train": int(len(variable_species)),
        }
        if train_keep.sum() < 2 or test_keep.sum() < 2:
            rows.append({**base, "status": "insufficient_variable_species_overlap", "n_test_positive": np.nan, "positive_prevalence_test": np.nan, "roc_auc": np.nan, "balanced_accuracy": np.nan})
            continue
        y_train = train_y[train_keep]
        y_test = test_y[test_keep]
        if len(np.unique(y_train)) < 2:
            rows.append({**base, "status": "unscorable_train_one_class", "n_test_positive": int(y_test.sum()), "positive_prevalence_test": float(y_test.mean()), "roc_auc": np.nan, "balanced_accuracy": np.nan})
            continue
        if len(np.unique(y_test)) < 2:
            rows.append({**base, "status": "unscorable_test_one_class", "n_test_positive": int(y_test.sum()), "positive_prevalence_test": float(y_test.mean()), "roc_auc": np.nan, "balanced_accuracy": np.nan})
            continue
        try:
            scaler = StandardScaler()
            climate_train = scaler.fit_transform(x_all[train_idx][train_keep])
            climate_test = scaler.transform(x_all[test_idx][test_keep])
            train_taxa = pd.Series(taxa[train_idx][train_keep], name="taxon_name")
            test_taxa = pd.Series(taxa[test_idx][test_keep], name="taxon_name")
            train_dummy = pd.get_dummies(train_taxa, drop_first=True, dtype=float)
            test_dummy = pd.get_dummies(test_taxa, drop_first=True, dtype=float).reindex(columns=train_dummy.columns, fill_value=0.0)
            x_train = np.column_stack([climate_train, train_dummy.to_numpy(dtype=float)])
            x_test = np.column_stack([climate_test, test_dummy.to_numpy(dtype=float)])
            classifier = make_classifier(C, class_weight, seed + fold)
            classifier.fit(x_train, y_train)
            probability = classifier.predict_proba(x_test)[:, 1]
            prediction = (probability >= 0.5).astype(int)
            rows.append({
                **base,
                "status": "success",
                "n_test_positive": int(y_test.sum()),
                "positive_prevalence_test": float(y_test.mean()),
                "roc_auc": float(roc_auc_score(y_test, probability)),
                "balanced_accuracy": float(balanced_accuracy_score(y_test, prediction)),
            })
        except Exception as error:
            rows.append({**base, "status": f"fit_error:{type(error).__name__}:{error}", "n_test_positive": int(y_test.sum()), "positive_prevalence_test": float(y_test.mean()), "roc_auc": np.nan, "balanced_accuracy": np.nan})
    return rows


def tier_for_trait(trait_id: str, plan: dict[str, Any]) -> str:
    for tier, traits in plan["trait_tiers"].items():
        if trait_id in traits:
            return tier
    return "not_predeclared"


def outcome_summary_markdown(
    statuses: pd.DataFrame,
    coefficients: pd.DataFrame,
    validation: pd.DataFrame,
    predictors: list[str],
    plan: dict[str, Any],
) -> str:
    eligible = statuses.loc[statuses["eligibility_status"].eq("eligible")].copy() if not statuses.empty else statuses
    lines = [
        "# Ch.1 trait–environment association atlas",
        "",
        plan["semantic_status"],
        "",
        "## Inputs",
        f"- Predictors: {', '.join(predictors)}",
        "- Measurement modes: all-ensemble winner and conservative strict-consensus.",
        "- Validation: grouped 10-degree spatial blocks, plus held-out species for climate-only models.",
        "- No random-split headline metrics, p-values, phylogenetic claims, or causal interpretation are included.",
        "",
        "## Eligible state outcomes",
        f"- Number eligible for modelling: {len(eligible)}",
    ]
    if not eligible.empty:
        lines.append(eligible[["measurement_mode", "trait_id", "trait_state", "trait_tier", "n_positive", "n_negative", "n_positive_species", "n_positive_spatial_blocks_10deg"]].to_markdown(index=False))
    lines.extend(["", "## Validation summary", ""])
    if not validation.empty:
        successful = validation.loc[validation["status"].eq("success")].copy()
        if successful.empty:
            lines.append("No grouped validation fold was jointly scorable for any eligible state outcome.")
        else:
            aggregate = successful.groupby(["measurement_mode", "trait_id", "trait_state", "trait_tier", "model_family", "validation_scheme"], as_index=False).agg(
                n_scorable_folds=("fold", "count"),
                mean_roc_auc=("roc_auc", "mean"),
                mean_balanced_accuracy=("balanced_accuracy", "mean"),
            )
            lines.append(aggregate.to_markdown(index=False, floatfmt=".3f"))
    lines.extend(["", "## Coefficient output", ""])
    lines.append(f"- Standardized penalized coefficients are provided for {len(coefficients)} predictor/outcome/model rows.")
    lines.append("- Coefficient sign and magnitude are descriptive associations conditional on the stated model family; they are not adaptive-effect estimates.")
    return "\n".join(lines) + "\n"


def main() -> None:
    args = parse_args()
    if args.n_folds < 2:
        raise ValueError("--n-folds must be at least 2")
    environment_path = Path(args.environment_long)
    provenance_path = Path(args.environment_provenance)
    plan_path = Path(args.plan)
    for path in (environment_path, provenance_path, plan_path):
        if not path.is_file():
            raise FileNotFoundError(path)

    long = pd.read_csv(environment_path, dtype=str, keep_default_na=False)
    require_columns(long, REQUIRED_COLUMNS, "environment-long table")
    if long.duplicated(["obs_id", "trait_id"]).any():
        raise ValueError("environment-long table has duplicate obs_id/trait_id rows")
    environment_provenance = json.loads(provenance_path.read_text(encoding="utf-8"))
    if not text(environment_provenance.get("semantic_status", "")).startswith("Environmental covariates joined"):
        raise ValueError("environment provenance does not state fully automated trait/environment join semantics")
    plan = json.loads(plan_path.read_text(encoding="utf-8"))
    predictors = [str(value) for value in plan["input_policy"]["required_predictors"]]
    require_columns(long, predictors, "environment-long table")
    for mode in plan["input_policy"]["measurement_modes"]:
        if mode["state_column"] not in long.columns:
            raise ValueError(f"Missing measurement-state column: {mode['state_column']}")

    complete = long["environment_primary_complete"].map(as_bool)
    work = long.loc[complete].copy()
    work = finite_frame(work, predictors)
    work["taxon_name"] = work["taxon_name"].map(text)
    work["spatial_block_10deg"] = work["spatial_block_10deg"].map(text)
    work = work.loc[work["taxon_name"].ne("") & work["spatial_block_10deg"].ne("")].copy()
    if work.empty:
        raise ValueError("No environment-complete trait rows remain after numeric and grouping checks")

    eligibility = plan["eligibility"]
    exclusion = {text(value) for value in plan["input_policy"]["exclude_states"]}
    excluded_traits = {text(value) for value in plan["input_policy"]["exclude_traits"]}
    C = float(plan["models"]["climate_only"]["regularization_C"])
    class_weight = plan["models"]["climate_only"].get("class_weight")

    feasibility_rows: list[dict[str, Any]] = []
    status_rows: list[dict[str, Any]] = []
    coefficient_rows: list[dict[str, Any]] = []
    validation_rows: list[dict[str, Any]] = []

    for mode in plan["input_policy"]["measurement_modes"]:
        mode_id = str(mode["id"])
        state_column = str(mode["state_column"])
        mode_data = work.copy()
        mode_data["trait_state"] = mode_data[state_column].map(text)
        mode_data = mode_data.loc[~mode_data["trait_state"].isin(exclusion)].copy()
        for trait_id, trait_data in mode_data.groupby("trait_id", sort=True):
            trait_id = text(trait_id)
            if trait_id in excluded_traits:
                continue
            trait_tier = tier_for_trait(trait_id, plan)
            states = sorted(state for state in trait_data["trait_state"].unique().tolist() if text(state) not in exclusion)
            for trait_state in states:
                frame = trait_data.copy()
                y = frame["trait_state"].eq(trait_state).astype(int).to_numpy(dtype=int)
                n_positive = int(y.sum())
                n_negative = int(len(y) - n_positive)
                positive = frame.loc[y.astype(bool)]
                common = {
                    "measurement_mode": mode_id,
                    "trait_id": trait_id,
                    "trait_state": trait_state,
                    "trait_tier": trait_tier,
                }
                eligible = (
                    n_positive >= int(eligibility["minimum_positive_observations"])
                    and n_negative >= int(eligibility["minimum_negative_observations"])
                    and int(positive["taxon_name"].nunique()) >= int(eligibility["minimum_positive_species"])
                    and int(positive["spatial_block_10deg"].nunique()) >= int(eligibility["minimum_positive_spatial_blocks_10deg"])
                    and int(frame["spatial_block_10deg"].nunique()) >= int(eligibility["minimum_total_spatial_blocks_10deg"])
                )
                feasibility_rows.append({
                    **common,
                    "n_observations": int(len(frame)),
                    "n_positive": n_positive,
                    "n_negative": n_negative,
                    "n_species": int(frame["taxon_name"].nunique()),
                    "n_positive_species": int(positive["taxon_name"].nunique()),
                    "n_spatial_blocks_10deg": int(frame["spatial_block_10deg"].nunique()),
                    "n_positive_spatial_blocks_10deg": int(positive["spatial_block_10deg"].nunique()),
                    "n_precision_le_10km": int(frame["coordinate_precision_tier"].isin(["high_le_1km", "moderate_1_to_10km"]).sum()),
                    "eligibility_status": "eligible" if eligible else "insufficient_support",
                    "eligibility_rule": "positive>=100; negative>=100; positive species>=10; positive 10-degree blocks>=10; total 10-degree blocks>=10",
                })
                if not eligible:
                    status_rows.append({**common, "model_family": "all", "status": "not_run_insufficient_support"})
                    continue

                climate_coefficients, climate_status = full_climate_coefficients(
                    frame, predictors, y, C, class_weight, args.seed, common
                )
                coefficient_rows.extend(climate_coefficients)
                status_rows.append({**common, "model_family": "climate_only", "status": climate_status})
                validation_rows.extend(
                    cv_climate_only(
                        frame, predictors, y, frame["spatial_block_10deg"], "spatial_block_10deg", args.n_folds,
                        C, class_weight, args.seed, common,
                    )
                )
                validation_rows.extend(
                    cv_climate_only(
                        frame, predictors, y, frame["taxon_name"], "held_out_species", args.n_folds,
                        C, class_weight, args.seed + 1000, common,
                    )
                )

                species_coefficients, species_status, n_variable_species = full_species_adjusted_coefficients(
                    frame, predictors, y, C, class_weight, args.seed + 2000, common
                )
                coefficient_rows.extend(species_coefficients)
                status_rows.append({
                    **common,
                    "model_family": "species_adjusted_within_species",
                    "status": species_status,
                    "n_species_within_state_variation": n_variable_species,
                })
                if n_variable_species >= int(eligibility["minimum_species_within_state_variation"]):
                    validation_rows.extend(
                        cv_species_adjusted_spatial(
                            frame, predictors, y, frame["spatial_block_10deg"], args.n_folds,
                            C, class_weight, args.seed + 3000, common,
                        )
                    )
                else:
                    validation_rows.append({
                        **common,
                        "model_family": "species_adjusted_within_species",
                        "validation_scheme": "spatial_block_10deg_within_species",
                        "fold": "",
                        "status": "not_run_insufficient_species_within_state_variation",
                        "n_variable_species_train": n_variable_species,
                    })

    feasibility = pd.DataFrame(feasibility_rows)
    statuses = pd.DataFrame(status_rows)
    coefficients = pd.DataFrame(coefficient_rows)
    validation = pd.DataFrame(validation_rows)
    if feasibility.empty:
        raise ValueError("No analyzable trait states were found after state-exclusion rules")

    output = Path(args.out_dir)
    output.mkdir(parents=True, exist_ok=True)
    feasibility.to_csv(output / "trait_environment_state_feasibility.csv", index=False, encoding="utf-8-sig")
    statuses.to_csv(output / "trait_environment_model_status.csv", index=False, encoding="utf-8-sig")
    coefficients.to_csv(output / "trait_environment_model_coefficients.csv", index=False, encoding="utf-8-sig")
    validation.to_csv(output / "trait_environment_grouped_validation.csv", index=False, encoding="utf-8-sig")
    outcome_summary_markdown(statuses, coefficients, validation, predictors, plan).replace("\r\n", "\n")
    (output / "trait_environment_association_report.md").write_text(
        outcome_summary_markdown(statuses, coefficients, validation, predictors, plan), encoding="utf-8"
    )

    summary = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "semantic_status": plan["semantic_status"],
        "n_environment_complete_trait_rows": int(len(work)),
        "n_environment_complete_observations": int(work["obs_id"].nunique()),
        "n_trait_state_outcomes_screened": int(len(feasibility)),
        "n_trait_state_outcomes_eligible": int(feasibility["eligibility_status"].eq("eligible").sum()),
        "n_coefficient_rows": int(len(coefficients)),
        "n_validation_rows": int(len(validation)),
        "predictors": predictors,
        "analysis_plan_version": plan.get("plan_version", ""),
        "environment_provenance_sha256": sha256_file(provenance_path),
        "input_sha256": {
            "environment_long": sha256_file(environment_path),
            "analysis_plan": sha256_file(plan_path),
        },
        "measurement_modes": [item["id"] for item in plan["input_policy"]["measurement_modes"]],
        "validation": plan["validation"],
        "sensitivity": plan["sensitivity"],
    }
    (output / "trait_environment_association_provenance.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
