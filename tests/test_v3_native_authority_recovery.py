import base64
import hashlib
import json
from pathlib import Path
from unittest.mock import patch

import pandas as pd
import pytest

from analysis.v3 import recover_native_source_authority as recovery
from analysis.v3.workflow import ROOT, digest


def response(payload):
    class Reply:
        status_code = 200
        content = json.dumps(payload, separators=(",", ":")).encode()

        def raise_for_status(self):
            pass
    return Reply()


def test_cache_reuses_exact_response_and_offline_pin(tmp_path):
    cache = recovery.ResponseCache(tmp_path / "cache")
    with patch.object(recovery.requests, "get", return_value=response({"results": []})) as get:
        original = cache.get("https://example.test/source")
        assert cache.get("https://example.test/source") == original
        assert get.call_count == 1
    pin = cache.seal()
    offline = recovery.ResponseCache(tmp_path / "cache", offline=True, expected_manifest_sha256=pin)
    with patch.object(recovery.requests, "get", side_effect=AssertionError("network forbidden")):
        assert offline.get("https://example.test/source") == original
        assert offline.seal() == pin


def test_offline_requires_manifest_pin_and_cache_rejects_tampering(tmp_path):
    cache = recovery.ResponseCache(tmp_path / "cache")
    with patch.object(recovery.requests, "get", return_value=response({"results": []})):
        cache.get("https://example.test/source")
    pin = cache.seal()
    with pytest.raises(ValueError, match="caller-pinned"):
        recovery.ResponseCache(tmp_path / "cache", offline=True)
    item = next(p for p in (tmp_path / "cache").glob("*.json") if p.name != "authority_cache_manifest.json")
    item.write_text("{}", encoding="utf-8")
    offline = recovery.ResponseCache(tmp_path / "cache", offline=True, expected_manifest_sha256=pin)
    with pytest.raises(ValueError, match="missing or changed"):
        offline.get("https://example.test/source")


def test_unsealed_cache_also_checks_response_identity(tmp_path):
    cache = recovery.ResponseCache(tmp_path)
    with patch.object(recovery.requests, "get", return_value=response({"results": []})):
        cache.get("https://example.test/source")
    item = next(tmp_path.glob("*.json"))
    data = json.loads(item.read_text())
    data["content_base64"] = base64.b64encode(b'{"changed":true}').decode()
    item.write_text(json.dumps(data))
    with pytest.raises(ValueError, match="identity mismatch"):
        cache.get("https://example.test/source")


def test_frozen_cache_cannot_expand_or_omit_requests(tmp_path):
    cache = recovery.ResponseCache(tmp_path)
    with patch.object(recovery.requests, "get", return_value=response({"results": []})):
        cache.get("https://example.test/a")
        cache.get("https://example.test/b")
    pin = cache.seal()
    offline = recovery.ResponseCache(tmp_path, offline=True, expected_manifest_sha256=pin)
    offline.get("https://example.test/a")
    with pytest.raises(ValueError, match="exact frozen request set"):
        offline.seal()
    with pytest.raises(ValueError, match="missing or changed"):
        offline.get("https://example.test/c")


def test_lf_recovery_preserves_original_and_requires_exact_pin(tmp_path):
    source = tmp_path / "windows.csv"
    source.write_bytes(b"a,b\r\n1,2\r\n")
    target = tmp_path / "lf.csv"
    with pytest.raises(ValueError, match="historical source hash"):
        recovery.recover_lf_source(source, target, "0" * 64)
    assert not target.exists()
    result = recovery.recover_lf_source(source, target, hashlib.sha256(b"a,b\n1,2\n").hexdigest())
    assert result["crlf_pairs_replaced"] == 2
    assert source.read_bytes() == b"a,b\r\n1,2\r\n"
    with pytest.raises(FileExistsError):
        recovery.recover_lf_source(source, target, result["output_sha256"])


@pytest.mark.parametrize("path", [ROOT, ROOT / "reproducibility/private", ROOT / "analysis/private"])
def test_private_output_cannot_enter_public_tree(path):
    with pytest.raises(ValueError, match="must use local_data"):
        recovery.private_directory(path)


def test_isolated_helpers_preserve_original_resolution_and_distribution_rules(tmp_path):
    import analysis.rebuild_frozen_native_status as legacy
    cache = recovery.ResponseCache(tmp_path)
    helpers = recovery.historical_helpers(cache)
    assert legacy.request_json is not helpers.request_json
    payload = {"results": [{"canonicalName": "Cirsium test", "key": 12,
                           "taxonomicStatus": "ACCEPTED", "scientificName": "Cirsium test"}]}
    assert helpers.resolve_exact_name("Cirsium test", payload) == legacy.resolve_exact_name("Cirsium test", payload)
    distributions = {"results": [
        {"locationId": "TDWG:AAA"},
        {"locationId": "TDWG:BBB", "establishmentMeans": "INTRODUCED"},
        {"locationId": "TDWG:CCC", "occurrenceStatus": "ABSENT"},
        {"locationId": "other:DDD", "establishmentMeans": "NATIVE"}], "endOfRecords": True}
    frame = pd.DataFrame({"accepted_key": [12]})
    with patch.object(legacy, "request_json", return_value=distributions):
        expected = legacy.fetch_distributions(frame, 1, 1, 0)
    with patch.object(recovery.requests, "get", return_value=response(distributions)):
        actual = helpers.fetch_distributions(frame, 1, 1, 0)
    pd.testing.assert_frame_equal(actual, expected)
    assert actual["classified_status"].tolist() == ["native", "introduced", "uncertain", "uncertain"]


def synthetic_run_inputs(tmp_path):
    contract = json.loads((ROOT / "analysis/v3/native_range_join_contract.json").read_text())
    geometry = {"type": "FeatureCollection", "features": [{"type": "Feature",
        "properties": {"LEVEL3_COD": "AAA"}, "geometry": {"type": "Polygon",
        "coordinates": [[[0, 0], [2, 0], [2, 2], [0, 2], [0, 0]]]}}]}
    contract["tdwg_level3"]["expected_sha256"] = hashlib.sha256(response(geometry).content).hexdigest()
    contract_path = tmp_path / "contract.json"
    contract_path.write_text(json.dumps(contract))
    frame = pd.DataFrame({"obs_id": [1, 2], "source_taxon_name": ["Cirsium test"] * 2,
        "source_taxon_rank": ["species"] * 2, "analysis_latitude": [1., 5.],
        "analysis_longitude": [1., 5.], "coordinate_status": ["public_location_present_precision_not_gated"] * 2,
        "captive_state": ["false"] * 2, "date_status": ["exact_day"] * 2,
        "observation_month": [5, 6], "position_accuracy_status": ["present"] * 2,
        "position_accuracy_m": [10., 10.]})
    source_path = tmp_path / "source.csv"
    frame.to_csv(source_path, index=False)

    def request(url, **kwargs):
        if "/dataset/" in url:
            spec = contract["source_taxonomy_and_distribution"]
            return response({"doi": spec["dataset_doi"], "modified": spec["dataset_modified"]})
        if "/species/search?" in url:
            return response({"results": [{"canonicalName": "Cirsium test", "scientificName": "Cirsium test",
                        "key": 12, "taxonomicStatus": "ACCEPTED"}]})
        if "/distributions?" in url:
            return response({"results": [{"locationId": "TDWG:AAA"}], "endOfRecords": True})
        if "raw.githubusercontent.com" in url:
            return response(geometry)
        raise AssertionError("Unexpected request")
    return source_path, contract_path, request


def test_complete_recovery_offline_replay_and_no_ecological_authorization(tmp_path):
    source, contract, request = synthetic_run_inputs(tmp_path)
    options = dict(source_csv=source, expected_source_sha256=digest(source),
                   cache_dir=tmp_path / "cache", contract_path=contract)
    with patch.object(recovery.requests, "get", side_effect=request):
        acquired = recovery.run(out_dir=tmp_path / "acquired", **options)
    assert acquired["n_observations"] == 2
    assert acquired["status_counts"] == {"native": 1, "unmapped_or_unusable_location": 1}
    assert acquired["primary_native_wild_public_exact_date_resolved_rows"] == 1
    assert not acquired["historical_cohort_identity_verified"]
    assert not acquired["ecological_fitting_authorized"]
    with patch.object(recovery.requests, "get", side_effect=AssertionError("network forbidden")):
        replay = recovery.run(out_dir=tmp_path / "replay", offline=True,
            expected_cache_manifest_sha256=acquired["authority_cache_manifest_sha256"], **options)
    assert replay["status"] == "NATIVE_AUTHORITY_REPLAYED"
    assert acquired["outputs_sha256"] == replay["outputs_sha256"]


def test_source_or_authority_drift_stops_before_output(tmp_path):
    source, contract, request = synthetic_run_inputs(tmp_path)
    kwargs = dict(source_csv=source, out_dir=tmp_path / "out", cache_dir=tmp_path / "cache",
                  contract_path=contract, expected_source_sha256="0" * 64)
    with patch.object(recovery.requests, "get", side_effect=AssertionError("network forbidden")):
        with pytest.raises(ValueError, match="input identity"):
            recovery.run(**kwargs)
    kwargs["expected_source_sha256"] = digest(source)
    with patch.object(recovery.requests, "get", return_value=response({"doi": "changed"})):
        with pytest.raises(ValueError, match="DOI differs"):
            recovery.run(**kwargs)
    assert not (tmp_path / "out").exists()
