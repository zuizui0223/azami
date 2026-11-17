"""Tests use synthetic measurements only; no scientific source results implied."""
import numpy as np
import pandas as pd
import pytest

pytest.importorskip("sklearn")

from analysis.v4.native_integration import FEATURES
from analysis.v4.rf_multitrait_predictability import (
    CORE4, FULL12, aggregate_endpoint, grouped_folds, score_models,
    _loss, _normalise_hue,
)


def tiny_synthetic(n_taxa=30, n_cells=8, reps=6):
    rows = []
    rng = np.random.default_rng(24)
    for taxon in range(n_taxa):
        region = taxon % 3
        for cell in range(n_cells):
            lat = -50 + 5 * cell + 0.01
            lon = (-100, 5, 90)[region] + 5 * cell + 0.01
            for repeat in range(reps):
                x = rng.normal(size=len(FULL12)) + cell * 0.01 + region
                row = {"accepted_key": "taxon_" + str(taxon),
                       "latitude": lat, "longitude": lon,
                       **dict(zip(FULL12, x))}
                row.update({"orientation_image_vertical_angle": 50 + x[0] * 8,
                            "corolla_lab_lightness": 50 + x[1] * 3,
                            "corolla_lab_chroma": 15 + x[2],
                            "corolla_hue_sin": np.sin(x[3]),
                            "corolla_hue_cos": np.cos(x[3]),
                            "capitulum_outline_aspect_ratio": 1.5 + x[4] / 10,
                            "capitulum_outline_circularity": .5 + x[5] / 100,
                            "capitulum_outline_solidity": .8 + x[6] / 100,
                            "capitulum_width_profile_cv": .3 + x[7] / 100})
                row.update({v + "__eligible": True for members in FEATURES.values()
                            for v in members})
                rows.append(row)
    return pd.DataFrame(rows)


def test_four_is_subset_of_source_twelve():
    assert len(FULL12) == 12
    assert len(CORE4) == 4
    assert set(CORE4) <= set(FULL12)
    assert set().union(*[set(x) for x in FEATURES.values()]).isdisjoint(FULL12)


def test_endpoint_qc_and_equal_eligible_base_for_both_models():
    data = tiny_synthetic()
    data.loc[0, "corolla_lab_chroma__eligible"] = False
    data.loc[1, "chelsa_bio18"] = np.nan
    a = aggregate_endpoint(data, "angle")
    b = aggregate_endpoint(data, "chroma")
    assert len(a) == len(b) == 30 * 8
    assert a.taxon.nunique() == 30
    assert a[list(FULL12)].notna().all().all()
    assert a.photos.min() >= 3


def test_no_species_or_sector_appears_on_both_sides():
    cells = aggregate_endpoint(tiny_synthetic(), "angle")
    for train, test in grouped_folds(cells, "unseen_taxa"):
        assert set(cells.iloc[train].taxon).isdisjoint(cells.iloc[test].taxon)
    for train, test in grouped_folds(cells, "unseen_longitude_sectors"):
        assert set(cells.iloc[train].region).isdisjoint(cells.iloc[test].region)


def test_run_is_deterministic_and_does_not_return_source_keys():
    cells = aggregate_endpoint(tiny_synthetic(), "angle")
    a = score_models(cells, "angle", "unseen_taxa", trees=8, repeats=1)
    b = score_models(cells, "angle", "unseen_taxa", trees=8, repeats=1)
    assert a == b
    assert set(a["models"]) == {"fixed_four_predictor_rf", "expanded_twelve_predictor_rf"}
    assert "taxon_" not in str(a)
    for model in a["models"].values():
        assert np.isfinite(model["out_of_fold_skill_vs_training_mean"])
        assert model["n_folds"] == 5
        assert all(np.isfinite(list(model["predictor_group_skill_loss_permuted_on_heldout_cells"].values())))


def test_circular_hue_joint_chord_loss():
    y = np.array([[0., 1.], [1., 0.]])
    assert _loss(y, y, np.ones(2), circular=True) == 0
    flip = np.array([[0., -1.], [-1., 0.]])
    assert _loss(y, flip, np.ones(2), circular=True) == pytest.approx(8.)
    assert np.allclose(np.linalg.norm(_normalise_hue(y), axis=1), 1.)


def test_rejects_missing_traits_not_silently_imputed():
    frame = tiny_synthetic()
    frame = frame.drop(columns=["corolla_lab_lightness__eligible"])
    with pytest.raises(ValueError, match="TRAIT_MISSING_COLUMNS"):
        aggregate_endpoint(frame, "lightness")


def test_provenance_gate_and_immutable_join(tmp_path, monkeypatch):
    import hashlib
    import json
    import analysis.v4.rf_multitrait_predictability as mod
    # No private source data: make a miniature structurally equivalent receipt.
    monkeypatch.setattr(mod, "SOURCE_N", 4)
    tiny = tiny_synthetic(n_taxa=1, n_cells=1, reps=4)
    tiny.insert(0, "obs_id", ["a", "b", "c", "d"])
    native = tiny.drop(columns=list(mod.EXTRA))
    extra = tiny[["obs_id", *mod.EXTRA]]
    native_file = tmp_path / "native.csv"
    extra_file = tmp_path / "extra.csv"
    receipt_file = tmp_path / "receipt.json"
    native.to_csv(native_file, index=False)
    extra.to_csv(extra_file, index=False)
    digest = mod.sha256(extra_file)
    receipt_file.write_text(json.dumps({
        "prepared_table_sha256": mod.sha256(native_file),
        "native_observations": 4,
        "environment_decoding_status": "physical_units_gsp_invalid_encoding_masked_v1"
    }))
    joined, record = mod.verified_source(native_file, receipt_file, extra_file, digest)
    assert joined.obs_id.tolist() == ["a", "b", "c", "d"]
    assert record["source_verified"] is True
    assert "obs_id" not in str(record)
    with pytest.raises(ValueError, match="SHA_MISMATCH"):
        mod.verified_source(native_file, receipt_file, extra_file, "0" * 64)
    bad = pd.concat([extra, extra.iloc[[0]]], ignore_index=True)
    bad.to_csv(extra_file, index=False)
    with pytest.raises(ValueError, match="DUPLICATE|CONFLICT"):
        mod.verified_source(native_file, receipt_file, extra_file, mod.sha256(extra_file))
