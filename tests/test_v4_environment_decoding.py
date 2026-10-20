import numpy as np
import pytest
import hashlib
import json
import pandas as pd

from analysis.v4 import decode_environment

from analysis.v4.decode_environment import decode_column, GSP_SOURCE_SHA


def test_gsp_declared_and_undeclared_sentinels_mask_before_scaling():
    metadata = {"nodata": 2147483647, "dtype": "uint32", "scale": .1,
                "offset": 0, "source_sha256": GSP_SOURCE_SHA}
    result, extra = decode_column(
        [0, 6690, 2147483647, 4294967295, np.nan], metadata, "chelsa_gsp")
    np.testing.assert_allclose(result[:2], [0,669])
    assert np.isnan(result[2:]).all()
    assert extra == 1
    with pytest.raises(ValueError, match="source changed"):
        decode_column([6690], {**metadata, "source_sha256": "other"}, "chelsa_gsp")


def test_decoding_preserves_valid_zero_and_applies_offset_once():
    metadata = {"nodata": -9999, "scale": .1, "offset": -273.15}
    result, extra = decode_column([0, 2827, -9999], metadata, "chelsa_bio01")
    np.testing.assert_allclose(result[:2], [-273.15,9.55])
    assert np.isnan(result[-1])
    assert extra == 0


def test_archive_decoding_keeps_traits_and_rows_but_changes_environment_support(tmp_path, monkeypatch):
    columns = decode_environment.FOUR + [
        "chelsa_bio04", "chelsa_bio15", "chelsa_sfcwind_mean", "chelsa_gsp", "chelsa_npp"]
    frame = pd.DataFrame({c:[100.,200.,300.] for c in columns})
    frame["chelsa_gsp"] = [100.,4294967295.,300.]
    frame["obs_id"] = [1,2,3]
    frame["trait"] = [1.2,np.nan,3.4]
    frame["trait__eligible"] = [True,False,True]
    frame["chelsa_bio01__between"] = [123.,123.,123.]
    table = tmp_path/"raw.csv.gz"
    frame.to_csv(table,index=False)
    digest = hashlib.sha256(table.read_bytes()).hexdigest()
    monkeypatch.setattr(decode_environment,"RAW_INPUT_SHA",digest)
    receipt = tmp_path/"receipt.json"
    receipt.write_text(json.dumps({"prepared_table_sha256":digest,"native_observations":3,
                                  "endpoints":{"trait":{"finite":2,"eligible":2,
                                                       "eligible_with_four_climate":2}}}))
    spec = {"nodata":2147483647,"dtype":"uint32","scale":.1,"offset":0,
            "source_sha256":GSP_SOURCE_SHA}
    metadata = tmp_path/"metadata.json"
    metadata.write_text(json.dumps({c:spec for c in columns}))
    out = tmp_path/"decoded"
    decode_environment.execute(table,receipt,metadata,out)
    decoded = pd.read_csv(out/"observations.csv.gz")
    np.testing.assert_allclose(decoded.trait,frame.trait,equal_nan=True)
    assert decoded.obs_id.tolist() == [1,2,3]
    assert decoded.trait__eligible.tolist() == [True,False,True]
    assert "chelsa_bio01__between" not in decoded
    r = json.loads((out/"measurement_receipt.json").read_text())
    assert r["complete_nine_climate"] == 2
    assert r["native_observations"] == 3
    with pytest.raises(ValueError,match="refuse repeat"):
        decode_environment.execute(out/"observations.csv.gz",out/"measurement_receipt.json",
                                   metadata,tmp_path/"twice")
