#!/usr/bin/env python3
"""Hardened entry point for the Ch.1 trait–environment association atlas.

This wrapper loads the base implementation in 43_run_trait_environment_association_atlas.py
and replaces only the species-adjusted spatial-CV function. The replacement uses
one category schema fixed from the training fold for both train and test species
indicators, avoiding fold-specific dummy baselines.
"""

from __future__ import annotations

import importlib.util
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from sklearn.metrics import balanced_accuracy_score, roc_auc_score
from sklearn.preprocessing import StandardScaler


BASE_PATH = Path(__file__).with_name("43_run_trait_environment_association_atlas.py")
spec = importlib.util.spec_from_file_location("ch1_trait_environment_base", BASE_PATH)
if spec is None or spec.loader is None:
    raise RuntimeError(f"Unable to load base association atlas: {BASE_PATH}")
base = importlib.util.module_from_spec(spec)
spec.loader.exec_module(base)


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
    splitter, n_splits = base.fold_iterator(groups, n_requested)
    if splitter is None:
        return [{
            **common,
            "model_family": "species_adjusted_within_species",
            "validation_scheme": "spatial_block_10deg_within_species",
            "fold": "",
            "status": "insufficient_groups",
            "n_groups": int(groups.astype(str).nunique()),
        }]

    rows: list[dict[str, Any]] = []
    groups_array = groups.astype(str).to_numpy()
    x_all = frame[predictors].to_numpy(dtype=float)
    taxa = frame["taxon_name"].astype(str).to_numpy()

    for fold, (train_idx, test_idx) in enumerate(splitter.split(x_all, y, groups_array), start=1):
        train_y_all = y[train_idx]
        test_y_all = y[test_idx]
        train_taxa_all = taxa[train_idx]
        test_taxa_all = taxa[test_idx]
        variation = pd.DataFrame({"taxon_name": train_taxa_all, "outcome": train_y_all}).groupby("taxon_name", observed=True)["outcome"].nunique()
        variable_species = sorted(variation.loc[variation >= 2].index.astype(str).tolist())
        train_keep = np.isin(train_taxa_all, variable_species)
        test_keep = np.isin(test_taxa_all, variable_species)
        base_row = {
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
            rows.append({**base_row, "status": "insufficient_variable_species_overlap", "n_test_positive": np.nan, "positive_prevalence_test": np.nan, "roc_auc": np.nan, "balanced_accuracy": np.nan})
            continue

        y_train = train_y_all[train_keep]
        y_test = test_y_all[test_keep]
        if len(np.unique(y_train)) < 2:
            rows.append({**base_row, "status": "unscorable_train_one_class", "n_test_positive": int(y_test.sum()), "positive_prevalence_test": float(y_test.mean()), "roc_auc": np.nan, "balanced_accuracy": np.nan})
            continue
        if len(np.unique(y_test)) < 2:
            rows.append({**base_row, "status": "unscorable_test_one_class", "n_test_positive": int(y_test.sum()), "positive_prevalence_test": float(y_test.mean()), "roc_auc": np.nan, "balanced_accuracy": np.nan})
            continue

        try:
            scaler = StandardScaler()
            climate_train = scaler.fit_transform(x_all[train_idx][train_keep])
            climate_test = scaler.transform(x_all[test_idx][test_keep])
            categories = sorted(set(train_taxa_all[train_keep].tolist()))
            train_categories = pd.Categorical(train_taxa_all[train_keep], categories=categories)
            test_categories = pd.Categorical(test_taxa_all[test_keep], categories=categories)
            train_dummy = pd.get_dummies(train_categories, drop_first=True, dtype=float)
            test_dummy = pd.get_dummies(test_categories, drop_first=True, dtype=float)
            # Both categorical vectors share the training-fold category universe.
            test_dummy = test_dummy.reindex(columns=train_dummy.columns, fill_value=0.0)
            x_train = np.column_stack([climate_train, train_dummy.to_numpy(dtype=float)])
            x_test = np.column_stack([climate_test, test_dummy.to_numpy(dtype=float)])
            classifier = base.make_classifier(C, class_weight, seed + fold)
            classifier.fit(x_train, y_train)
            probability = classifier.predict_proba(x_test)[:, 1]
            prediction = (probability >= 0.5).astype(int)
            rows.append({
                **base_row,
                "status": "success",
                "n_test_positive": int(y_test.sum()),
                "positive_prevalence_test": float(y_test.mean()),
                "roc_auc": float(roc_auc_score(y_test, probability)),
                "balanced_accuracy": float(balanced_accuracy_score(y_test, prediction)),
            })
        except Exception as error:
            rows.append({
                **base_row,
                "status": f"fit_error:{type(error).__name__}:{error}",
                "n_test_positive": int(y_test.sum()),
                "positive_prevalence_test": float(y_test.mean()),
                "roc_auc": np.nan,
                "balanced_accuracy": np.nan,
            })
    return rows


base.cv_species_adjusted_spatial = cv_species_adjusted_spatial

if __name__ == "__main__":
    base.main()
