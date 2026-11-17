"""Source receipts, 21/23 predictor split, no raw output tested on synthetic data."""
import json
import numpy as np
import pytest
from analysis.v4.rf_archived_atmosphere import (
    ATMOSPHERE2, BASE21, ATMOSPHERE23, MODEL_PAIR,
    PREDICTOR_GROUPS, verify_original_1981_2010_airfields,
)
from analysis.v4.rf_multitrait_predictability import aggregate_endpoint, score_models
from test_v4_rf_multitrait_predictability import tiny_synthetic


def test_archived_chelsa_atmosphere_not_already_in_21():
    assert len(BASE21)==21
    assert len(ATMOSPHERE23)==23
    assert not set(ATMOSPHERE2)&set(BASE21)
    assert set(MODEL_PAIR)=={"complete_bioclim21","bio21_plus_chelsa_vpd_wind23"}
    assert set(ATMOSPHERE2) == set(PREDICTOR_GROUPS["atmospheric_demand_CHELSA"] +
                                    PREDICTOR_GROUPS["wind"][-1:])


def test_archived_fields_reject_unpinned_source(tmp_path):
    frame=tiny_synthetic()
    for k in ATMOSPHERE2:
        frame[k]=np.linspace(0.1,2.1,len(frame))
    receipt={"environment_decoding_status":"physical_units_gsp_invalid_encoding_masked_v1",
             "environment_decoding":{k:{"source_sha256":"a"*64,"scale":.1,"offset":0}
                                     for k in ATMOSPHERE2}}
    metadata={k:{"source_sha256":"a"*64,"scale":.1,"offset":0} for k in ATMOSPHERE2}
    rec=tmp_path/"rec.json";meta=tmp_path/"meta.json"
    rec.write_text(json.dumps(receipt));meta.write_text(json.dumps(metadata))
    result=verify_original_1981_2010_airfields(frame,rec,meta)
    assert all(v["finite_s1_records"]==len(frame) for v in result.values())
    receipt["environment_decoding"][ATMOSPHERE2[0]]["scale"] = 10
    rec.write_text(json.dumps(receipt))
    with pytest.raises(ValueError,match="SOURCE_SHA_SCALE_DRIFT"):
        verify_original_1981_2010_airfields(frame,rec,meta)


def test_paired_21_v_23_uses_same_eligible_photo_cell_rows():
    frame=tiny_synthetic()
    rng=np.random.default_rng(42)
    for name in (n for n in BASE21 if n not in frame.columns):
        frame[name] = rng.normal(size=len(frame))
    for k in ATMOSPHERE2:
        frame[k]=rng.normal(size=len(frame))
    cells=aggregate_endpoint(frame,"angle",predictors=ATMOSPHERE23)
    x=score_models(cells,"angle","unseen_taxa_balanced",trees=8,repeats=1,
                   predictor_sets=MODEL_PAIR,predictor_groups=PREDICTOR_GROUPS)
    assert x["n_taxa"]==30
    assert x["models"]["complete_bioclim21"]["baseline_weighted_sse"] == \
           x["models"]["bio21_plus_chelsa_vpd_wind23"]["baseline_weighted_sse"]
    assert np.isfinite(x["paired_skill_gain_upper_minus_lower"])
    assert set(ATMOSPHERE2) <= set(x["models"]["bio21_plus_chelsa_vpd_wind23"]["predictors"])