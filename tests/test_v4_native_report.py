import hashlib
import json
import pytest
from analysis.v4.render_native_report import read_sources


def sources(root):
    payload = json.dumps({"input_sha256": "same-input"})
    (root/"native_trait_models.json").write_text(payload)
    digest = hashlib.sha256(payload.encode()).hexdigest()
    for name in ("native_integration", "native_spatial"):
        (root/(name+".json")).write_text(json.dumps(
            {"input_sha256": "same-input", "model_sha256": digest}))
    (root/"native_integration_audit.json").write_text(payload)


def test_report_checks_fitted_model_identity(tmp_path):
    sources(tmp_path)
    assert len(read_sources(tmp_path)) == 4
    path = tmp_path/"native_trait_models.json"
    path.write_text(path.read_text()+" ")
    with pytest.raises(ValueError, match="fitted models"):
        read_sources(tmp_path)


def test_report_rejects_different_population_audit(tmp_path):
    sources(tmp_path)
    (tmp_path/"native_integration_audit.json").write_text(
        json.dumps({"input_sha256": "another-input"}))
    with pytest.raises(ValueError, match="different inputs"):
        read_sources(tmp_path)
