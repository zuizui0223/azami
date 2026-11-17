"""Exploratory within-taxon 1-degree floral deviations under held-out 5-degree blocks.

A different estimand from the original unseen-taxon RF: predict *local deviations*
for species already observed in training. The taxon baseline and predictor means
use TRAINING blocks only. Held-out 5-degree block identity is never an RF feature.
Only anonymised aggregate skill leaves the authorised original-S1 runner.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from analysis.v4.native_integration import FEATURES
from analysis.v4.rf_archived_atmosphere import verify_original_1981_2010_airfields
from analysis.v4.rf_extended_environment import BIO9, TERRAIN5, verify_extended_source
from analysis.v4.rf_multitrait_predictability import (
    aggregate_endpoint, sha256, verified_source, _loss,
)
from analysis.v4.rf_worldclim_terrain import MODEL_SET

SCOPED_CONSTRUCTS = ("angle", "chroma", "aspect")
MIN_TRAIN_CELLS_PER_TAXON = 3
FOLD_SEED = 20261031


def heldout_spatial_folds(cells: pd.DataFrame, n_splits: int = 5):
    """Intact 5-degree geographical blocks; no spatial block is shared."""
    from sklearn.model_selection import GroupKFold
    if not {"lat_cell", "lon_cell", "taxon"}.issubset(cells):
        raise ValueError("GEO_CELL_KEYS_REQUIRED_FOR_HELDOUT_BLOCKS")
    if ((cells.lat_cell < 0).any() or (cells.lat_cell >= 180).any() or
            (cells.lon_cell < 0).any() or (cells.lon_cell >= 360).any()):
        raise ValueError("ONE_DEGREE_LAT_LON_CELL_KEYS_REQUIRED")
    blocks = (cells.lat_cell.to_numpy(int) // 5) * 72 + cells.lon_cell.to_numpy(int) // 5
    if np.unique(blocks).size < n_splits:
        raise ValueError("INSUFFICIENT_DISTINCT_GEOGRAPHIC_BLOCKS")
    splits = list(GroupKFold(n_splits=n_splits).split(cells, groups=blocks))
    for train, test in splits:
        if set(blocks[train]).intersection(blocks[test]):
            raise AssertionError("GEOGRAPHIC_BLOCK_LEAKAGE")
    return splits


def within_taxon_heldout(cells: pd.DataFrame, construct: str, *,
                         trees: int = 80, min_train_cells: int = MIN_TRAIN_CELLS_PER_TAXON,
                         n_splits: int = 5) -> dict:
    """Taxon-train-mean baseline versus RF on train-centred local anomalies.

    Same source observations, same target and fold masks for paired 23 and 28.
    Test species lacking train data are never assigned an invented intercept.
    """
    from sklearn.ensemble import RandomForestRegressor
    if (construct not in SCOPED_CONSTRUCTS or trees < 8 or min_train_cells < 3 or
            set(cells.columns).issuperset({"photo_id", "photographer_id"})):
        raise ValueError("UNREGISTERED_LOCAL_TAXON_PREDICTION_DESIGN")
    target_cols = FEATURES[construct]
    if len(target_cols) != 1:
        raise ValueError("LOCAL_PILOT_LINEAR_PHOTOGRAPHIC_TRAITS_ONLY")
    cols28 = MODEL_SET["chelsa_vpd_wind23_plus_worldclim_terrain28"]
    missing = set(cols28).difference(cells.columns)
    if missing:
        raise ValueError("MISSING_28_SOURCE_VERIFIED_PREDICTORS")
    x = cells[list(cols28)].to_numpy(float)
    y = cells[list(target_cols)].to_numpy(float)
    if not np.isfinite(x).all() or not np.isfinite(y).all():
        raise ValueError("NONFINITE_OR_IMPUTED_SOURCE_VALUES")
    folds = heldout_spatial_folds(cells, n_splits=n_splits)
    kept = []
    stats = {name: {"model_sse": 0.0, "baseline_sse": 0.0, "fold_skills": []}
             for name in MODEL_SET}
    for index, (tr, te) in enumerate(folds):
        train = cells.iloc[tr]; test = cells.iloc[te]
        train_counts = train.groupby("taxon").size()
        allowed = set(train_counts.index[train_counts >= min_train_cells])
        keep = test.taxon.isin(allowed).to_numpy(bool)
        te_used = te[keep]
        if len(te_used) < 5 or not allowed:
            raise ValueError("GEOGRAPHIC_FOLD_NO_TRAIN_COVERED_TEST_TAXA")
        test_species = cells.iloc[te_used].taxon.to_numpy(str)
        train_species = train.taxon.to_numpy(str)
        train_y = pd.DataFrame(y[tr], index=train_species).groupby(level=0).mean()
        train_x = pd.DataFrame(x[tr], index=train_species).groupby(level=0).mean()
        original_train_baseline = train_y.loc[train_species].to_numpy(float)
        original_test_baseline = train_y.loc[test_species].to_numpy(float)
        original_train_x = train_x.loc[train_species].to_numpy(float)
        original_test_x = train_x.loc[test_species].to_numpy(float)
        residual_y = y[tr] - original_train_baseline
        residual_x = x[tr] - original_train_x
        test_residual_x = x[te_used] - original_test_x
        fit_w = (1 / train.taxon.map(train_counts)).to_numpy(float)
        n_test_by_species = cells.iloc[te_used].groupby("taxon").size()
        test_w = (1 / cells.iloc[te_used].taxon.map(n_test_by_species)).to_numpy(float)
        baseline_loss = _loss(y[te_used], original_test_baseline, test_w,
                              circular=False)
        if baseline_loss <= 1e-10:
            raise ValueError("ZERO_WITHIN_TAXON_HELDOUT_BASELINE_LOSS")
        record = {"fold": index, "test_cells_with_known_taxon": len(te_used),
                  "test_taxa_with_training_baseline": len(set(test_species)),
                  "heldout_test_cells_without_known_taxon": int(len(te) - len(te_used))}
        for name, prednames in MODEL_SET.items():
            sel = [list(cols28).index(v) for v in prednames]
            model = RandomForestRegressor(n_estimators=trees, min_samples_leaf=5,
                                           max_features=0.75, n_jobs=1,
                                           random_state=FOLD_SEED + index)
            model.fit(residual_x[:, sel], residual_y.ravel(), sample_weight=fit_w)
            fitted = original_test_baseline + model.predict(test_residual_x[:, sel]).reshape(-1, 1)
            error = _loss(y[te_used], fitted, test_w, circular=False)
            stats[name]["model_sse"] += error
            stats[name]["baseline_sse"] += baseline_loss
            stats[name]["fold_skills"].append(1 - error / baseline_loss)
        kept.append(record)
    models = {name:{"skill_over_species_training_means": float(1 - v["model_sse"] / v["baseline_sse"]),
                    "model_weighted_sse": float(v["model_sse"]),
                    "baseline_weighted_sse": float(v["baseline_sse"]),
                    "fold_skill_values": [float(x) for x in v["fold_skills"]]}
              for name,v in stats.items()}
    names = list(MODEL_SET)
    return {
        "construct": construct, "status": "source_matched_posthoc_within_taxon_geoblock_RF",
        "n_taxa": int(cells.taxon.nunique()), "n_taxon_cells": int(len(cells)),
        "min_train_cells_per_known_taxon": min_train_cells,
        "heldout_block_width_degrees": 5, "cell_aggregation_degrees": 1,
        "folds_anonymous_counts": kept,
        "train_species_means_never_use_heldout_trait_values": True,
        "heldout_unknown_taxa_excluded_not_predicted": True,
        "models": models,
        "paired_skill_gain_28_minus_23": float(
            models[names[1]]["skill_over_species_training_means"] -
            models[names[0]]["skill_over_species_training_means"]),
    }


def run(native: Path, receipt: Path, archived_seven: Path, seven_sha: str,
        bio_nine: Path, nine_sha: str, nine_manifest: Path,
        terrain_table: Path, terrain_manifest: Path, five_deg_verified: Path,
        original_metadata: Path, output: Path, *, trees: int = 80):
    """Authenticate same S1/BIO21/DEM source as completed terrain analysis."""
    old = json.loads(five_deg_verified.read_text())
    prov = old.get("terrain_provenance", {})
    if (old.get("status") != "original_source_terrain_derived_and_evaluated" or
            old.get("source",{}).get("s1_observations") != 313350 or
            prov.get("mismatch_n") != 0 or
            sha256(terrain_table) != prov.get("terrain_private_table_sha256") or
            sha256(terrain_manifest) != prov.get("terrain_private_manifest_sha256")):
        raise ValueError("ORIGINAL_DEM_S1_VERIFIED_SOURCE_REQUIRED")
    frame, base = verified_source(native, receipt, archived_seven, seven_sha)
    frame,_ = verify_extended_source(frame,bio_nine,nine_sha,nine_manifest,
                                      required_columns=BIO9)
    verify_original_1981_2010_airfields(frame,receipt,original_metadata)
    frame,_ = verify_extended_source(frame,terrain_table,sha256(terrain_table),
                                      terrain_manifest,required_columns=TERRAIN5)
    if old["source"]["source_s1_sha256"] != base["source_s1_sha256"]:
        raise ValueError("DIFFERENT_SOURCE_S1_THAN_TERRAIN_REFERENCE")
    results=[]
    for name in SCOPED_CONSTRUCTS:
        cells=aggregate_endpoint(frame,name,
            predictors=MODEL_SET["chelsa_vpd_wind23_plus_worldclim_terrain28"],
            min_per_cell=3,degrees=1)
        results.append(within_taxon_heldout(cells,name,trees=trees))
    report={"schema":"azami_s1_1degree_within_taxon_5degree_spatial_holdout_RF_v1",
            "status":"source_fitted_posthoc_distinct_estimand",
            "source_s1_sha256":base["source_s1_sha256"],
            "n_original_s1_observations":base["s1_observations"],
            "verified_original_dem_mismatch_n":prov["mismatch_n"],
            "registered_predictor_universes":{k:list(v) for k,v in MODEL_SET.items()},
            "estimand":"Within-known-taxon local photographic deviations, with taxon intercept and predictor centring estimated ONLY from training 5degree blocks",
            "original_unseen_taxon_and_sector_RF_not_replaced":True,
            "original_GEB_S1_S2_H3c_unmodified":True,
            "no_source_identifiers_locations_per_photo_predictions_or_model_weights_exported":True,
            "results":results}
    output.parent.mkdir(parents=True,exist_ok=True)
    output.write_text(json.dumps(report,ensure_ascii=False,indent=2,allow_nan=False)+"\n",encoding="utf-8")
    return report


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for arg in ("native","receipt","archived-seven","seven-sha",
                "bio-nine","nine-sha","nine-manifest",
                "terrain-table","terrain-manifest","five-deg-verified",
                "original-metadata","output"):
        p.add_argument("--"+arg,required=True,
            type=str if arg.endswith("sha") else Path)
    p.add_argument("--trees",type=int,default=80)
    a=p.parse_args()
    run(a.native,a.receipt,a.archived_seven,a.seven_sha,a.bio_nine,a.nine_sha,
        a.nine_manifest,a.terrain_table,a.terrain_manifest,a.five_deg_verified,
        a.original_metadata,a.output,trees=a.trees)

if __name__=="__main__":
    main()