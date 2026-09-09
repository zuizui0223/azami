from __future__ import annotations

import pytest

from analysis.v3.reconcile_wave_c_inventory import reconcile


def batch():
    chunks = {f"c{i:06d}": {"observations": 1, "photo_jobs": 1, "requests": 1} for i in range(18, 146)}
    return {
        "schema_version": 1,
        "status": "bounded_native_raw_measurement_batch_no_ecology",
        "plan_id": "1" * 64,
        "chunks": chunks,
        "aggregate": {"observations": 128, "photo_jobs": 128, "requests": 128},
    }


def asset(chunk: str, *, final=True, aid=1):
    suffix = ".zip" if final else "-incomplete-123-1.zip"
    return {
        "id": aid,
        "name": f"v3-raw-{'1'*16}-{chunk}{suffix}",
        "state": "uploaded",
        "size": 100,
        "digest": "sha256:" + "a" * 64,
    }


def test_complete_inventory_requires_all_128_unique_final_assets():
    assets = [asset(f"c{i:06d}", aid=i) for i in range(18, 146)]
    report = reconcile(batch(), assets)
    assert report["status"] == "WAVE_C_PROTECTED_FINAL_INVENTORY_COMPLETE_128_OF_128_NO_RESUME"
    assert report["protected_final_chunks"] == 128
    assert report["missing_chunks"] == []
    assert report["resume_chunks"] == []
    assert report["trait_values_read"] == 0
    assert report["ecological_models_executed"] == 0


def test_missing_final_is_resume_target_and_checkpoint_is_only_context():
    assets = [asset(f"c{i:06d}", aid=i) for i in range(18, 146) if i not in {22, 111}]
    assets.append(asset("c000022", final=False, aid=999))
    report = reconcile(batch(), assets)
    assert report["status"] == "WAVE_C_PROTECTED_FINAL_INVENTORY_INCOMPLETE_RESUME_ONLY_MISSING"
    assert report["missing_chunks"] == ["c000022", "c000111"]
    assert report["missing_with_protected_checkpoint"] == ["c000022"]
    assert report["missing_without_protected_checkpoint"] == ["c000111"]
    assert report["resume_chunks"] == report["missing_chunks"]


def test_duplicate_final_asset_stops_instead_of_choosing():
    assets = [asset(f"c{i:06d}", aid=i) for i in range(18, 146)]
    assets.append(asset("c000022", aid=1000))
    with pytest.raises(ValueError, match="Duplicate protected final assets"):
        reconcile(batch(), assets)
