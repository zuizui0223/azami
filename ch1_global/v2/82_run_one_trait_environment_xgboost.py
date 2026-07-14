#!/usr/bin/env python3
from __future__ import annotations

import argparse
import importlib.util
import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.inspection import permutation_importance
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
from sklearn.model_selection import GroupKFold
from xgboost import XGBRegressor


def load_base():
    path = Path(__file__).with_name("76_run_exhaustive_within_species_environment_rf.py")
    spec = importlib.util.spec_from_file_location("env76", path)
    mod = importlib.util.module_from_spec(spec)
    assert spec and spec.loader
    spec.loader.exec_module(mod)
    return mod


def weighted_metrics(y, pred, weights):
    return (
        float(r2_score(y, pred, sample_weight=weights)),
        float(mean_squared_error(y, pred, sample_weight=weights) ** 0.5),
        float(mean_absolute_error(y, pred, sample_weight=weights)),
    )


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--environment", required=True)
    p.add_argument("--trait", required=True)
    p.add_argument("--out-dir", required=True)
    p.add_argument("--n-splits", type=int, default=5)
    p.add_argument("--n-estimators", type=int, default=1000)
    p.add_argument("--learning-rate", type=float, default=0.03)
    p.add_argument("--max-depth", type=int, default=6)
    p.add_argument("--min-child-weight", type=float, default=5.0)
    p.add_argument("--subsample", type=float, default=0.8)
    p.add_argument("--colsample-bytree", type=float, default=0.8)
    p.add_argument("--reg-lambda", type=float, default=1.0)
    p.add_argument("--reg-alpha", type=float, default=0.0)
    p.add_argument("--random-state", type=int, default=20260714)
    p.add_argument("--permutation-max-n", type=int, default=2000)
    p.add_argument("--permutation-repeats", type=int, default=1)
    args = p.parse_args()

    mod = load_base()
    if args.trait not in mod.TRAITS:
        raise SystemExit(f"Unknown trait: {args.trait}")

    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(args.environment, low_memory=False)
    climate = list(mod.CLIMATE)
    topo = list(mod.TOPOGRAPHY)
    soil = [f"soil_{prop}_0_30cm" for prop in mod.SOIL_PROPERTIES]
    groups = {
        "climate": climate,
        "climate_topography": climate + topo,
        "climate_soil": climate + soil,
        "climate_topography_soil": climate + topo + soil,
    }

    fold_rows = []
    importance_rows = []
    for model_group, predictors in groups.items():
        work = mod.prepare_within(df, args.trait, predictors)
        if len(work) < 100 or work["taxon_name"].nunique() < 3 or work["spatial_block_10deg"].nunique() < args.n_splits:
            fold_rows.append({
                "trait": args.trait,
                "model_group": model_group,
                "status": "insufficient",
                "n_observations": len(work),
                "n_species": work["taxon_name"].nunique(),
                "n_spatial_blocks": work["spatial_block_10deg"].nunique(),
            })
            continue

        xcols = [f"x__{p}" for p in predictors]
        X = work[xcols].to_numpy(float)
        y = work["y_within"].to_numpy(float)
        weights = work["species_equal_weight"].to_numpy(float)
        spatial_groups = work["spatial_block_10deg"].astype(str).to_numpy()
        splitter = GroupKFold(n_splits=args.n_splits)

        for fold, (train_idx, test_idx) in enumerate(splitter.split(X, y, spatial_groups), start=1):
            model = XGBRegressor(
                objective="reg:squarederror",
                tree_method="hist",
                n_estimators=args.n_estimators,
                learning_rate=args.learning_rate,
                max_depth=args.max_depth,
                min_child_weight=args.min_child_weight,
                subsample=args.subsample,
                colsample_bytree=args.colsample_bytree,
                reg_lambda=args.reg_lambda,
                reg_alpha=args.reg_alpha,
                random_state=args.random_state + fold,
                n_jobs=-1,
            )
            model.fit(X[train_idx], y[train_idx], sample_weight=weights[train_idx], verbose=False)
            pred = model.predict(X[test_idx])
            r2, rmse, mae = weighted_metrics(y[test_idx], pred, weights[test_idx])
            fold_rows.append({
                "trait": args.trait,
                "model_group": model_group,
                "status": "ok",
                "fold": fold,
                "n_train": len(train_idx),
                "n_test": len(test_idx),
                "n_species_total": int(work["taxon_name"].nunique()),
                "n_test_spatial_blocks": int(pd.Series(spatial_groups[test_idx]).nunique()),
                "weighted_r2": r2,
                "weighted_rmse": rmse,
                "weighted_mae": mae,
            })

            if model_group == "climate_topography_soil":
                rng = np.random.default_rng(args.random_state + fold)
                use = test_idx
                if len(use) > args.permutation_max_n:
                    use = np.sort(rng.choice(use, size=args.permutation_max_n, replace=False))
                result = permutation_importance(
                    model,
                    X[use],
                    y[use],
                    scoring="neg_mean_squared_error",
                    n_repeats=args.permutation_repeats,
                    random_state=args.random_state + fold,
                    n_jobs=-1,
                )
                for predictor, mean_imp, sd_imp in zip(predictors, result.importances_mean, result.importances_std):
                    importance_rows.append({
                        "trait": args.trait,
                        "fold": fold,
                        "predictor": predictor,
                        "permutation_importance_mean": float(mean_imp),
                        "permutation_importance_sd": float(sd_imp),
                    })

    folds = pd.DataFrame(fold_rows)
    folds.to_csv(out / "xgb_spatial_cv_folds.csv", index=False, encoding="utf-8-sig")
    ok = folds.loc[folds.get("status", pd.Series(index=folds.index, dtype=str)).eq("ok")].copy()
    summary = (
        ok.groupby(["trait", "model_group"], as_index=False)
        .agg(
            n_folds=("fold", "nunique"),
            mean_weighted_r2=("weighted_r2", "mean"),
            sd_weighted_r2=("weighted_r2", "std"),
            mean_weighted_rmse=("weighted_rmse", "mean"),
            mean_weighted_mae=("weighted_mae", "mean"),
        )
        if not ok.empty else pd.DataFrame()
    )
    summary.to_csv(out / "xgb_spatial_cv_summary.csv", index=False, encoding="utf-8-sig")
    pd.DataFrame(importance_rows).to_csv(out / "xgb_full_model_permutation_importance.csv", index=False, encoding="utf-8-sig")

    if not summary.empty:
        wide = summary.set_index("model_group")["mean_weighted_r2"]
        base = wide.get("climate", np.nan)
        ablation = pd.DataFrame([{
            "trait": args.trait,
            "r2_climate": base,
            "r2_climate_topography": wide.get("climate_topography", np.nan),
            "r2_climate_soil": wide.get("climate_soil", np.nan),
            "r2_full": wide.get("climate_topography_soil", np.nan),
            "delta_topography_over_climate": wide.get("climate_topography", np.nan) - base,
            "delta_soil_over_climate": wide.get("climate_soil", np.nan) - base,
            "delta_full_over_climate": wide.get("climate_topography_soil", np.nan) - base,
        }])
    else:
        ablation = pd.DataFrame()
    ablation.to_csv(out / "xgb_environment_group_ablation.csv", index=False, encoding="utf-8-sig")
    (out / "trait_xgb_report.json").write_text(json.dumps({
        "trait": args.trait,
        "model": "XGBRegressor",
        "tree_method": "hist",
        "purpose": "nonlinear predictive benchmark; SPDE-INLA remains primary inferential analysis",
    }, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
