import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from analysis.v3.production_environment import CONTRACT
from analysis.v3.select_production_environment import select
from analysis.v3.workflow import digest


def fixture(tmp_path):
    contract = json.loads(CONTRACT.read_text())
    contract["source"].update(observations=120, taxa=3)
    specification = tmp_path / "contract.json"
    specification.write_text(json.dumps(contract))
    frame = pd.DataFrame(np.random.default_rng(20260908).normal(size=(120, 15)), columns=contract["acquire"])
    frame["obs_id"] = [str(i) for i in range(120)]
    frame["accepted_key"] = np.repeat(["a", "b", "c"], 40)
    frame["forbidden_trait"] = "not a number"
    matrix = tmp_path / "matrix.csv"
    frame.to_csv(matrix, index=False)
    return matrix, specification


def test_candidate_selection_reads_only_environment_and_does_not_promote(tmp_path):
    matrix, contract = fixture(tmp_path)
    out = tmp_path / "result.json"
    result = select(matrix, digest(matrix), out, contract)
    assert result["complete_rows"] == 120 and result["complete_taxa"] == 3
    assert result["candidate_only"] is True
    assert result["trait_values_read"] == result["ecological_models_executed"] == 0
    assert result["active_ecological_contract_changed"] is False
    with pytest.raises(ValueError, match="Preserve"):
        select(matrix, digest(matrix), out, contract)


def test_changed_matrix_is_rejected_before_output(tmp_path):
    matrix, contract = fixture(tmp_path)
    with pytest.raises(ValueError, match="identity"):
        select(matrix, "0" * 64, tmp_path / "result.json", contract)
    assert not (tmp_path / "result.json").exists()


def test_duplicate_observations_not_silently_deduplicated(tmp_path):
    matrix, contract = fixture(tmp_path)
    frame = pd.read_csv(matrix)
    frame.loc[1, "obs_id"] = frame.loc[0, "obs_id"]
    frame.to_csv(matrix, index=False)
    with pytest.raises(ValueError, match="counts or identities"):
        select(matrix, digest(matrix), tmp_path / "result.json", contract)


def test_unmasked_gsp_integer_maximum_blocks_selection_without_repair(tmp_path):
    matrix, contract = fixture(tmp_path)
    frame = pd.read_csv(matrix)
    frame.loc[0, "GSP"] = 429496729.5
    frame.to_csv(matrix, index=False)
    before = digest(matrix)
    result = select(matrix, before, tmp_path / "result.json", contract)
    assert result["status"] == "ENVIRONMENT_SELECTION_BLOCKED_BY_INVALID_SOURCE_VALUES"
    assert result["invalid_GSP_uint32_max_rows"] == 1
    assert result["selected"] == result["trace"] == []
    assert digest(matrix) == before
