"""Synthetic-only 44-predictor registration, provenance, and holdout tests."""
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

pytest.importorskip("sklearn")

from analysis.v4.native_integration import FEATURES
from analysis.v4.rf_multitrait_predictability import (
    CORE4, FULL12, aggregate_endpoint, score_models,
)
import analysis.v4.rf_extended_environment as extended


def synthetic(n_taxa=30, n_cells=8, n_photos=6):
    rng = np.random.default_rng(42)
    rows = []
    for i in range(n_taxa):
        region = i % 3
        for j in range(n_cells):
            lat = -55 + 5 * j + 0.02
            lon = (-100, 0, 85)[region] + 5 * j + .02
            for r in range(n_photos):
                xx = rng.normal(size=44) + region*.2 + j*.1
                row = {"obs_id": f"obs_{i}_{j}_{r}",
                       "accepted_key": f"taxon_{i}", "latitude": lat,
                       "longitude": lon,
                       **dict(zip(extended.TIERS["vegetation44"], xx))}
                row.update({
                    "orientation_image_vertical_angle": 90 + xx[4] + xx[16],
                    "corolla_lab_lightness": 55 + xx[1],
                    "corolla_lab_chroma": 12 + xx[6],
                    "corolla_hue_sin": np.sin(xx[5]),
                    "corolla_hue_cos": np.cos(xx[5]),
                    "capitulum_outline_aspect_ratio": 1.5 + xx[0]*.1,
                    "capitulum_outline_circularity": .7 + xx[2]*.01,
                    "capitulum_outline_solidity": .8 + xx[3]*.01,
                    "capitulum_width_profile_cv": .4 + xx[6]*.01,
                })
                row.update({v + "__eligible": True for a in FEATURES.values() for v in a})
                rows.append(row)
    return pd.DataFrame(rows)


def write_external(tmp_path: Path, frame: pd.DataFrame, *, corrupt=None):
    data = frame[["obs_id", *extended.NEW32]].copy()
    if corrupt == "missing_id":
        data.loc[0, "obs_id"] = "not_original"
    if corrupt == "duplicate_id":
        data.loc[0, "obs_id"] = data.loc[1, "obs_id"]
    table = tmp_path / "extra32.csv"
    data.to_csv(table, index=False)
    digest = extended.sha256(table)
    manifest = {
        "schema": "azami_rf_extra32_source_manifest_v1",
        "table_sha256": digest,
        "layers": {
            v: {"dataset": extended.DATASET_FAMILIES[v],
                "temporal_support": extended.EXPECTED_PERIODS[v],
                "source_bundle_sha256": "a"*64,
                "unit": "source unit", "derivation": "documented mean or layer statistic",
                "sampling_method": "original coordinate pixel sample",
                "source_url": "https://www.example.org/source/reference",
                **({"depth_interval": "0-5cm", "statistic": "Q0.50"}
                   if v in extended.SOIL6 else {})}
            for v in extended.NEW32
        }
    }
    manifest_path = tmp_path / "manifest.json"
    manifest_path.write_text(json.dumps(manifest))
    return table, digest, manifest_path


def test_tier_coverage_and_separate_predictor_blocks():
    assert len(extended.NEW32) == 32
    assert len(extended.TIERS["vegetation44"]) == 44
    assert list(map(lambda x:len(x[1]), extended.STAGES)) == [4, 12, 21, 26, 28, 33, 39, 44]
    assert len(extended.PAIRS) == 7
    assert set().union(*map(set, extended.GROUPS44.values())) == set(extended.TIERS["vegetation44"])
    assert set(CORE4) <= set(FULL12)
    assert not {"latitude", "longitude", "accepted_key", "obs_id"} & set(extended.TIERS["vegetation44"])


def test_44_tier_paired_same_cohort_and_holdouts():
    frame = synthetic()
    cells = aggregate_endpoint(frame, "angle", predictors=extended.TIERS["vegetation44"])
    assert cells.taxon.nunique() == 30
    assert len(cells) == 240
    result = score_models(cells, "angle", "unseen_taxa", trees=8, repeats=1,
                          predictor_sets={"soil39": extended.TIERS["soil39"],
                                          "vegetation44": extended.TIERS["vegetation44"]},
                          predictor_groups=extended.GROUPS44)
    assert result["paired_model_order"] == ["soil39", "vegetation44"]
    assert np.isfinite(result["paired_skill_gain_upper_minus_lower"])
    assert set(result["models"]["vegetation44"]["predictor_group_skill_loss_permuted_on_heldout_cells"]) == set(extended.GROUPS44)
    assert "taxon_" not in str(result)


def test_21_and_44_require_different_complete_case_support():
    frame = synthetic()
    frame.loc[0:10, "soil_phh2o_0_5cm"] = np.nan
    lower = aggregate_endpoint(frame, "angle", predictors=extended.TIERS["complete_bioclim21"])
    upper = aggregate_endpoint(frame, "angle", predictors=extended.TIERS["vegetation44"])
    assert lower.photos.sum() == 30 * 8 * 6
    assert upper.photos.sum() < lower.photos.sum()
    # Whole taxon_0 loses its eight-cell support when its first two cells
    # are incomplete; the enlarged candidate cohort must not be invented.
    assert len(upper) == len(lower) - 8
    assert "taxon_0" not in set(upper.taxon)
    assert not upper[list(extended.TIERS["vegetation44"])].isna().any().any()


def test_verification_manifest_and_observation_identifiers(tmp_path):
    frame = synthetic(n_taxa=1, n_cells=1, n_photos=4)
    original = frame.drop(columns=list(extended.NEW32)).copy()
    table, digest, manifest = write_external(tmp_path, frame)
    joined, receipt = extended.verify_extended_source(original, table, digest, manifest)
    assert joined.obs_id.tolist() == original.obs_id.tolist()
    assert set(extended.NEW32) <= set(joined)
    assert receipt["source_identifiers_not_exported"]
    assert not any("obs_" in str(v) for v in receipt.values())
    assert receipt["declaration_only_raw_raster_bytes_not_rechecked_here"]
    with pytest.raises(ValueError, match="SHA_MISMATCH"):
        extended.verify_extended_source(original, table, "0"*64, manifest)
    m = json.loads(manifest.read_text())
    m["layers"][extended.BIO9[0]]["dataset"] = "invented_dataset"
    manifest.write_text(json.dumps(m))
    with pytest.raises(ValueError, match="DECLARATION"):
        extended.verify_extended_source(original, table, digest, manifest)


def test_identifier_alignment_and_duplicate_rejected(tmp_path):
    frame = synthetic(n_taxa=1, n_cells=1, n_photos=4)
    original = frame.drop(columns=list(extended.NEW32)).copy()
    table, digest, manifest = write_external(tmp_path, frame, corrupt="missing_id")
    with pytest.raises(ValueError, match="MISSING_SOURCE_OBSERVATIONS"):
        extended.verify_extended_source(original, table, digest, manifest)
    table, digest, manifest = write_external(tmp_path, frame, corrupt="duplicate_id")
    with pytest.raises(ValueError, match="IDENTITY_INVALID"):
        extended.verify_extended_source(original, table, digest, manifest)


def test_invalid_source_stratum_does_not_trigger_silent_drop(tmp_path):
    frame = synthetic(n_taxa=1, n_cells=1, n_photos=4)
    original = frame.drop(columns=list(extended.NEW32)).copy()
    table, digest, manifest = write_external(tmp_path, frame)
    # Additional columns are prohibited; they can carry sensitive data or a
    # response-derived predictor into the RF without audit.
    bad = frame[["obs_id", *extended.NEW32, "orientation_image_vertical_angle"]]
    bad.to_csv(table, index=False)
    m = json.loads(manifest.read_text())
    m["table_sha256"] = extended.sha256(table)
    manifest.write_text(json.dumps(m))
    with pytest.raises(ValueError, match="IDENTITY_INVALID"):
        extended.verify_extended_source(original, table, extended.sha256(table), manifest)


def test_runner_requires_declared_tier_before_reading_source(tmp_path):
    with pytest.raises(ValueError, match="UNKNOWN_OR_BASE_TIER"):
        extended.run_extended(tmp_path / "missing", tmp_path / "missing", tmp_path / "missing",
                              "x", tmp_path / "missing", "x", tmp_path / "missing",
                              tmp_path / "out.json", upper_tier="not_registered")


def test_partial_21_tier_needs_only_first_nine_nonmissing_layers(tmp_path):
    frame = synthetic(n_taxa=1, n_cells=1, n_photos=4)
    original = frame.drop(columns=list(extended.NEW32)).copy()
    nine = extended.added_columns("complete_bioclim21")
    table, digest, manifest = write_external(tmp_path, frame)
    # Only the 9 CHELSA data columns have been acquired; no invented soil/NDVI.
    frame[["obs_id", *nine]].to_csv(table, index=False)
    doc = json.loads(manifest.read_text())
    doc["layers"] = {k: doc["layers"][k] for k in nine}
    doc["table_sha256"] = extended.sha256(table)
    manifest.write_text(json.dumps(doc))
    joined, metadata = extended.verify_extended_source(
        original, table, extended.sha256(table), manifest, required_columns=nine)
    assert len(nine) == 9 and len(joined) == len(original)
    assert set(nine).issubset(joined)
    assert not set(extended.SOIL6) & set(joined)
    assert len(metadata["nonmissing_counts_by_axis"]) == 9