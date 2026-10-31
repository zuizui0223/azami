"""Supplement-only, source-verified S1 Random Forest trait predictability screen.

This exploratory analysis does NOT refit/replace the frozen S1/S2 Holm tests or H3c.
No original observation, taxon, uploader, coordinate, cell or prediction rows leave
this runner. The only generated file is an aggregate JSON with source digests.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

from analysis.v4.native_integration import FEATURES
from analysis.v4.s1_s2_extreme_lego import ALL_ENV, BASE, EXTRA

SOURCE_N = 313350
SEED = 20261031
CORE4 = ("chelsa_bio18", "chelsa_bio05", "chelsa_rsds_mean", "topo_elevation")
FULL12 = tuple(ALL_ENV)
PREDICTOR_SETS = {"fixed_four_predictor_rf": CORE4, "expanded_twelve_predictor_rf": FULL12}
GROUPS = {
    "temperature": ("chelsa_bio01", "chelsa_bio04", "chelsa_bio05", "chelsa_bio06"),
    "precipitation": ("chelsa_bio12", "chelsa_bio13", "chelsa_bio14", "chelsa_bio15", "chelsa_bio17", "chelsa_bio18"),
    "radiation": ("chelsa_rsds_mean",),
    "elevation": ("topo_elevation",),
}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def _require_columns(frame: pd.DataFrame, required, label: str):
    missing = set(required).difference(frame.columns)
    if missing:
        raise ValueError(label + "_MISSING_COLUMNS: " + ",".join(sorted(missing)))


def verified_source(native_path: Path, receipt_path: Path,
                    extra_path: Path, extra_sha256: str) -> tuple[pd.DataFrame, dict]:
    """Original S1 grain, decode/bytes, and one-to-one external 7-axis join."""
    receipt = json.loads(Path(receipt_path).read_text(encoding="utf-8"))
    native_hash = sha256(native_path)
    if (receipt.get("prepared_table_sha256") != native_hash or
            receipt.get("native_observations") != SOURCE_N or
            receipt.get("environment_decoding_status") !=
            "physical_units_gsp_invalid_encoding_masked_v1"):
        raise ValueError("S1_ORIGINAL_SOURCE_OR_DECODE_NOT_VERIFIED")
    if len(extra_sha256) != 64 or sha256(extra_path) != extra_sha256:
        raise ValueError("EXTERNAL_SEVEN_AXIS_SOURCE_SHA_MISMATCH")
    native = pd.read_csv(native_path, dtype={"obs_id": str}, low_memory=False)
    extra = pd.read_csv(extra_path, dtype={"obs_id": str}, low_memory=False)
    end_cols = [c for names in FEATURES.values() for c in names]
    _require_columns(native, ["obs_id", "accepted_key", "latitude", "longitude",
                              *BASE, *end_cols, *(v + "__eligible" for v in end_cols)], "S1")
    _require_columns(extra, ["obs_id", *EXTRA], "EXTRA")
    if (len(native) != SOURCE_N or native.obs_id.isna().any() or
            not native.obs_id.is_unique or extra.obs_id.isna().any() or
            not extra.obs_id.is_unique or set(EXTRA).intersection(native.columns)):
        raise ValueError("ORIGINAL_S1_IDENTITY_OR_EXTERNAL_JOIN_CONFLICT")
    joined = native.merge(extra[["obs_id", *EXTRA]], on="obs_id", how="left",
                          validate="one_to_one")
    if len(joined) != SOURCE_N:
        raise ValueError("SOURCE_JOIN_CHANGED_OBSERVATION_GRAIN")
    return joined, {"source_s1_sha256": native_hash, "extra_environment_sha256": extra_sha256,
                    "s1_observations": SOURCE_N, "source_verified": True,
                    "original_individual_records_exported": False}


def _valid_measurement(frame: pd.DataFrame, endpoints: list[str]) -> np.ndarray:
    mask = np.ones(len(frame), dtype=bool)
    for name in endpoints:
        mask &= frame[name + "__eligible"].eq(True).to_numpy(bool)
        mask &= np.isfinite(pd.to_numeric(frame[name], errors="coerce").to_numpy(float))
    return mask


def aggregate_endpoint(frame: pd.DataFrame, construct: str, *, min_taxa: int = 25,
                       min_cells: int = 8, min_photos: int = 40,
                       min_per_cell: int = 3, degrees: int = 5,
                       predictors: tuple[str, ...] = FULL12) -> pd.DataFrame:
    """One source-photo median per taxon×5° cell; outcome-specific QC.

    Identical complete required-predictor cells underlie the paired RF fits.
    The original function default remains the verified, frozen-source 12 axes.
    No response values imputed. Hue sine and cosine use a common photo set.
    """
    if construct not in FEATURES or degrees != 5:
        raise ValueError("UNKNOWN_CONSTRUCT_OR_UNDECLARED_SPATIAL_GRAIN")
    names = list(FEATURES[construct])
    if not set(FULL12).issubset(predictors) or len(set(predictors)) != len(predictors):
        raise ValueError("EXPANDED_PREDICTOR_UNIVERSE_MUST_INCLUDE_VERIFIED_TWELVE")
    _require_columns(frame, ["accepted_key", "latitude", "longitude", *predictors,
                              *names, *(v + "__eligible" for v in names)], "TRAIT")
    columns = list(predictors) + ["latitude", "longitude"] + names
    z = frame[columns].apply(pd.to_numeric, errors="coerce")
    coords = z[["latitude", "longitude"]]
    finite = np.isfinite(z).all(axis=1).to_numpy(dtype=bool, copy=True)
    finite &= coords.latitude.between(-90, 90).to_numpy(bool)
    finite &= coords.longitude.between(-180, 180).to_numpy(bool)
    finite &= frame.accepted_key.notna().to_numpy(bool)
    finite &= _valid_measurement(frame, names)
    w = z.loc[finite, list(predictors) + names].copy()
    w["taxon"] = frame.loc[finite, "accepted_key"].astype(str).to_numpy()
    lat = coords.loc[finite, "latitude"].to_numpy(float)
    lon = coords.loc[finite, "longitude"].to_numpy(float)
    w["lat_cell"] = np.floor((lat + 90) / degrees).astype(int)
    w["lon_cell"] = np.floor((lon + 180) / degrees).astype(int)
    w["region"] = np.select((lon < -30, lon < 60), (0, 1), default=2).astype(int)
    g = w.groupby(["taxon", "lat_cell", "lon_cell"], sort=True, observed=True)
    out = g[list(predictors) + names].median()
    out["photos"] = g.size()
    out["region"] = g.region.agg(lambda s: int(s.mode().iloc[0]))
    out = out.reset_index()
    out = out.loc[out.photos >= min_per_cell].copy()
    qual = out.groupby("taxon").agg(n_cells=("photos", "size"), n_photos=("photos", "sum"))
    keep = qual.index[(qual.n_cells >= min_cells) & (qual.n_photos >= min_photos)]
    out = out.loc[out.taxon.isin(keep)].reset_index(drop=True)
    if construct == "hue":
        # The median of sin/cos coordinates can be close to zero when the
        # local photographed hues conflict; direction is then undefined.
        norm = np.linalg.norm(out[names].to_numpy(float), axis=1)
        out = out.loc[norm >= 0.20].reset_index(drop=True)
    if construct == "hue":
        # Re-apply the taxon support gate after dropping ambiguous cell hue.
        counts = out.groupby("taxon").agg(n_cells=("photos", "size"),
                                           n_photos=("photos", "sum"))
        keep_hue = counts.index[(counts.n_cells >= min_cells) &
                                (counts.n_photos >= min_photos)]
        out = out.loc[out.taxon.isin(keep_hue)].reset_index(drop=True)
    if out.taxon.nunique() < min_taxa or out.region.nunique() < 3:
        raise ValueError("INSUFFICIENT_QUALIFIED_TAXA_OR_THREE_REGIONS_" + construct)
    if (out[list(predictors) + names].isna().any().any() or
            not np.isfinite(out[list(predictors) + names]).all().all()):
        raise ValueError("NONFINITE_AGGREGATED_VALUES")
    return out


def _normalise_hue(v: np.ndarray) -> np.ndarray:
    v = np.asarray(v, dtype=float)
    norm = np.linalg.norm(v, axis=1, keepdims=True)
    return v / np.maximum(norm, 1e-12)


def _loss(y: np.ndarray, prediction: np.ndarray, weight: np.ndarray,
          *, circular: bool) -> float:
    if circular:
        y = _normalise_hue(y)
        prediction = _normalise_hue(prediction)
    diff = y - prediction
    return float(np.sum(weight * np.sum(diff * diff, axis=1)))


def grouped_folds(cells: pd.DataFrame, kind: str, *, n_taxon_folds: int = 5):
    """Completely held-out taxa OR broad longitude sectors, never random photos."""
    if kind == "unseen_taxa":
        from sklearn.model_selection import GroupKFold
        groups = cells.taxon.to_numpy()
        n = min(n_taxon_folds, len(np.unique(groups)))
        if n < 3:
            raise ValueError("TOO_FEW_TAXA_FOR_GROUPED_EVALUATION")
        return list(GroupKFold(n_splits=n).split(cells, groups=groups))
    if kind == "unseen_taxa_balanced":
        # Posthoc robustness: equal numbers of complete taxa in each fold;
        # no species appears on both sides. The original GroupKFold remains
        # the untouched primary exploratory analysis.
        species = np.unique(cells.taxon.astype(str).to_numpy())
        n = min(n_taxon_folds, len(species))
        if n < 3:
            raise ValueError("TOO_FEW_TAXA_FOR_BALANCED_GROUP_EVALUATION")
        shuffled = np.random.default_rng(SEED).permutation(species)
        folds = []
        for species_holdout in np.array_split(shuffled, n):
            is_test = cells.taxon.isin(species_holdout).to_numpy(bool)
            train, test = np.flatnonzero(~is_test), np.flatnonzero(is_test)
            if not len(train) or not len(test):
                raise ValueError("EMPTY_BALANCED_TAXON_FOLD")
            folds.append((train, test))
        return folds
    if kind == "unseen_longitude_sectors":
        from sklearn.model_selection import LeaveOneGroupOut
        groups = cells.region.to_numpy()
        if len(np.unique(groups)) != 3:
            raise ValueError("THREE_WIDE_SECTORS_REQUIRED")
        return list(LeaveOneGroupOut().split(cells, groups=groups))
    raise ValueError("UNKNOWN_GENERALISATION_TARGET")


def score_models(cells: pd.DataFrame, construct: str, scheme: str, *,
                 trees: int = 160, repeats: int = 3, seed: int = SEED,
                 predictor_sets: dict[str, tuple[str, ...]] | None = None,
                 predictor_groups: dict[str, tuple[str, ...]] | None = None) -> dict:
    """Paired held-out evaluation; grouped permutation on test folds only.

    Skill is 1 - weighted SSE(model) / weighted SSE(training-fold constant).
    A negative value means worse than the training-fold constant prediction.
    """
    from sklearn.ensemble import RandomForestRegressor
    predictor_sets = PREDICTOR_SETS if predictor_sets is None else predictor_sets
    predictor_groups = GROUPS if predictor_groups is None else predictor_groups
    if len(predictor_sets) != 2:
        raise ValueError("EXACTLY_TWO_PAIRED_RF_PREDICTOR_SETS_REQUIRED")
    endpoints = FEATURES[construct]
    circular = construct == "hue"
    groups = grouped_folds(cells, scheme)
    species_counts = cells.groupby("taxon").size()
    weights = (1.0 / cells.taxon.map(species_counts)).to_numpy(float)
    all_results = {}
    for model_name, predictors in predictor_sets.items():
        columns = list(predictors)
        x = cells[columns].to_numpy(float)
        y = cells[endpoints].to_numpy(float)
        base_loss = model_loss = 0.0
        perm = {k: 0.0 for k, vars_ in predictor_groups.items() if set(vars_) & set(columns)}
        folds = []
        for fold, (train, test) in enumerate(groups):
            if len(train) < 25 or len(test) < 5:
                raise ValueError("INSUFFICIENT_GROUPED_TRAIN_OR_TEST_CELLS")
            est = RandomForestRegressor(n_estimators=trees, max_features=0.75,
                                        min_samples_leaf=5, n_jobs=1,
                                        random_state=seed + fold)
            est.fit(x[train], y[train].ravel() if len(endpoints) == 1 else y[train],
                    sample_weight=weights[train])
            truth = y[test]
            estimate = np.asarray(est.predict(x[test])).reshape(-1, len(endpoints))
            reference = np.average(y[train], axis=0, weights=weights[train])[None, :]
            const = np.repeat(reference, len(test), axis=0)
            fold_base = _loss(truth, const, weights[test], circular=circular)
            fold_model = _loss(truth, estimate, weights[test], circular=circular)
            if fold_base <= 1e-12:
                raise ValueError("DEGENERATE_REFERENCE_LOSS")
            base_loss += fold_base
            model_loss += fold_model
            folds.append({"holdout": fold, "test_cells": int(len(test)),
                          "test_taxa": int(cells.iloc[test].taxon.nunique()),
                          "skill": float(1 - fold_model / fold_base)})
            rng = np.random.default_rng(seed + 100 * fold)
            for category, members in predictor_groups.items():
                cols = [columns.index(name) for name in members if name in columns]
                if not cols:
                    continue
                for rep in range(repeats):
                    shuffled = x[test].copy()
                    order = rng.permutation(len(test))
                    shuffled[:, cols] = shuffled[order][:, cols]
                    prediction = np.asarray(est.predict(shuffled)).reshape(-1, len(endpoints))
                    perm[category] += (_loss(truth, prediction, weights[test], circular=circular)
                                      - fold_model) / repeats
        all_results[model_name] = {
            "predictors": columns, "n_folds": len(folds),
            "out_of_fold_skill_vs_training_mean": float(1 - model_loss / base_loss),
            "baseline_weighted_sse": float(base_loss),
            "model_weighted_sse": float(model_loss),
            "predictor_group_skill_loss_permuted_on_heldout_cells":
                {k: float(v / base_loss) for k, v in perm.items()},
            "fold_diagnostics_aggregate_only": folds,
        }
    key0, key1 = list(predictor_sets)
    gain = float(all_results[key1]["out_of_fold_skill_vs_training_mean"] -
                 all_results[key0]["out_of_fold_skill_vs_training_mean"])
    result = {"construct": construct, "generalisation": scheme,
            "n_cells": int(len(cells)), "n_taxa": int(cells.taxon.nunique()),
            "photo_contributions": int(cells.photos.sum()),
            "hue_is_joint_circular_response": bool(circular),
            "models": all_results,
            "paired_skill_gain_upper_minus_lower": gain,
            "paired_model_order": [key0, key1]}
    if key0 == "fixed_four_predictor_rf" and key1 == "expanded_twelve_predictor_rf":
        result["twelve_minus_four_out_of_fold_skill"] = gain
    return result


def run(native: Path, receipt: Path, extra: Path, extra_sha256: str,
        output: Path, *, trees: int = 160, repeats: int = 3) -> dict:
    frame, provenance = verified_source(native, receipt, extra, extra_sha256)
    results = []
    for name in FEATURES:
        cells = aggregate_endpoint(frame, name)
        for split in ("unseen_taxa", "unseen_longitude_sectors"):
            results.append(score_models(cells, name, split, trees=trees, repeats=repeats))
    report = {
        "schema": "azami_v4_exploratory_rf_predictability_v1",
        "status": "source_fitted_exploratory_not_original_frozen_holm",
        "source": provenance,
        "seed": SEED, "n_estimators": trees, "permutation_repeats": repeats,
        "cell_grain": "median per accepted taxon x 5 degree lat-lon cell (>=3 source photographs)",
        "taxon_gate": ">=40 source photographs and >=8 eligible cells per construct",
        "models": {k: list(v) for k, v in PREDICTOR_SETS.items()},
        "evaluation": "paired unseen-taxon 5-fold and 3 broad held-out longitude sectors",
        "environment_scope": "12 existing, source-verified CHELSA/RSDS/elevation variables only; not soil or pollinators",
        "scope": "out-of-sample predictive association, not causal variable importance or adaptation",
        "original_S1_S2_Holm_families_unchanged": True,
        "original_H3c_analogue_results_unchanged": True,
        "source_taxon_observer_photo_coordinate_and_predictions_not_exported": True,
        "results": results,
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    return report


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--native-table", type=Path, required=True)
    p.add_argument("--native-receipt", type=Path, required=True)
    p.add_argument("--extra-table", type=Path, required=True)
    p.add_argument("--extra-sha256", type=str, required=True)
    p.add_argument("--output", type=Path, required=True)
    args = p.parse_args()
    run(args.native_table, args.native_receipt, args.extra_table, args.extra_sha256,
        args.output)


if __name__ == "__main__":
    main()