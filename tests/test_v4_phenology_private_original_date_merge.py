"""Original metadata and private flowering event dates, never camera-derived."""
import gzip
import io
import json
import hashlib
from pathlib import Path
import zipfile

import pandas as pd
import pytest
# AESGCM is an optional private-data dependency; keep generic suite importable.

from analysis.v4 import phenology_private_original_date_merge as m


def generate(tmp_path,monkeypatch,*,meta_date=True,native_date=False,conflict=False):
    AESGCM=pytest.importorskip('cryptography.hazmat.primitives.ciphers.aead').AESGCM
    dataset=pd.DataFrame({
        "obs_id":["o1","o1","o2","o3","o4"],
        "user_id":["u1","u1","u2","u3","u3"],
        "quality_grade":["research"]*5,
    })
    if meta_date:
        dataset["observed_on"]=["2023-03-04","2023-03-04","2023-04-11",
                                "2024-03-08","2024-04-12"]
    zipbuf=io.BytesIO()
    with zipfile.ZipFile(zipbuf,"w",zipfile.ZIP_DEFLATED) as z:
        z.writestr("photo_metadata_merged.csv",dataset.to_csv(index=False))
    zipbytes=zipbuf.getvalue()
    key=b"R"*32
    nonce=b"N"*12
    secret=nonce+AESGCM(key).encrypt(nonce,zipbytes,m.AAD)
    archive=tmp_path/"metadata.aesgcm"
    archive.write_bytes(secret)
    monkeypatch.setenv("AZAMI_METADATA_KEY",key.hex())
    monkeypatch.setattr(m,"SEALED_SHA256",hashlib.sha256(secret).hexdigest())
    monkeypatch.setattr(m,"SOURCE_ZIP_SHA256",hashlib.sha256(zipbytes).hexdigest())
    monkeypatch.setattr(m,"EXPECTED_SOURCE_ROWS",5)
    monkeypatch.setattr(m,"EXPECTED_SOURCE_OBSERVATIONS",4)
    monkeypatch.setattr(m,"EXPECTED_NATIVE_N",3)
    monkeypatch.setattr(m,"EXPECTED_NATIVE_OBSERVERS",3)
    raw=pd.DataFrame({"obs_id":["o1","o2","o3"]})
    if native_date:
        raw["observed_on"]=["2023-03-04","2023-04-11",
                            "2024-04-08" if conflict else "2024-03-08"]
    source=tmp_path/"native.csv.gz"
    raw.to_csv(source,index=False,compression="gzip")
    return source,archive,tmp_path/"private.csv.gz"


def test_original_observation_date_from_metadata_only(tmp_path,monkeypatch):
    s,e,o=generate(tmp_path,monkeypatch)
    results,path=m.original_dates(s,e,o)
    assert path.exists()
    assert results["status"]=="ORIGINAL_OBSERVATION_DATE_SOURCE_VERIFIED_PRIVATE_ONLY"
    assert results["date_valid_native_observations"]==3
    assert results["original_observer_count"]==3
    assert results["month_distribution_original_events"]["3"]==2
    assert results["month_distribution_original_events"]["4"]==1
    private=pd.read_csv(o,compression="gzip")
    assert list(private)==["obs_id","user_id","observed_on"]
    assert not results["flowering_phenology_true_anthesis_observed"]
    assert "u1" not in json.dumps(results) and "o1" not in json.dumps(results)


def test_strict_agreement_between_native_and_archive_dates(tmp_path,monkeypatch):
    s,e,o=generate(tmp_path,monkeypatch,meta_date=True,native_date=True)
    results,_=m.original_dates(s,e,o)
    assert results["date_valid_native_observations"]==3


def test_conflicting_original_observation_dates_not_silently_patched(tmp_path,monkeypatch):
    s,e,o=generate(tmp_path,monkeypatch,meta_date=True,native_date=True,conflict=True)
    with pytest.raises(ValueError,match="disagree"):
        m.original_dates(s,e,o)


def test_no_original_source_date_is_not_filled_from_photo_upload_time(tmp_path,monkeypatch):
    s,e,o=generate(tmp_path,monkeypatch,meta_date=False,native_date=False)
    results,out=m.original_dates(s,e,o)
    assert out is None
    assert results["status"].startswith("STOP_")
    assert not o.exists()


def test_archive_source_hash_is_immutable(tmp_path,monkeypatch):
    s,e,o=generate(tmp_path,monkeypatch)
    e.write_bytes(e.read_bytes()+b"corrupt")
    with pytest.raises(ValueError,match="encrypted"):
        m.original_dates(s,e,o)
