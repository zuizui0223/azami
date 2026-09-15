import json
from pathlib import Path

import numpy as np
import pandas as pd

from analysis.v3 import run_construct_scale_upgrade as upgrade
from analysis.v3 import run_rv_estimator_validity as validity


ROOT = Path(__file__).resolve().parents[1]


def synthetic_common() -> pd.DataFrame:
    rows = []
    for taxon_index in range(6):
        for obs_index in range(4):
            row = {
                "obs_id": f"o{taxon_index}_{obs_index}",
                "taxon_name": f"t{taxon_index}",
            }
            for feature_index, construct in enumerate(validity.CORE):
                row[f"{construct}__0"] = (
                    (taxon_index + 1) * (feature_index + 1)
                    + (obs_index - 1.5) * (1 + feature_index % 3)
                    + 0.01 * obs_index * feature_index
                )
            rows.append(row)
    frame = pd.DataFrame(rows)
    frame.attrs["feature_map"] = {c: [f"{c}__0"] for c in validity.CORE}
    return frame


def test_contract_declares_posthoc_status_and_fixed_rules():
    contract = json.loads(
        (ROOT / "analysis/ch1/RV_ESTIMATOR_VALIDITY_CONTRACT_20260915.json").read_text()
    )
    assert contract["status"] == "posthoc_estimator_validity_sensitivity"
    assert contract["frozen_input_scope"]["common_complete_observations"] == 1734
    assert contract["frozen_input_scope"]["common_taxa"] == 42
    assert contract["diagnostics"]["equal_n_within_resampling"]["replicates"] == 1000
    assert contract["diagnostics"]["permutation_null_centering"]["permutations_per_pair_per_scale"] == 499


def test_equal_n_resampling_uses_one_row_per_taxon():
    frame = synthetic_common()
    _, _, pairs = upgrade.matrices_from_common(frame)
    result = validity.equal_n_within_resampling(frame, pairs, replicates=7, seed=123)
    assert len(result) == 7
    assert set(result.n_rows_within) == {6}
    assert set(result.n_taxa) == {6}
    assert result.relations_stronger_among.between(0, 36).all()
    assert np.isfinite(result.among_minus_within_median_rv).all()


def test_coordinate_standardization_sets_taxon_median_sd_to_one():
    frame = synthetic_common()
    standardized = validity.standardize_coordinates_by_taxon_medians(frame)
    for construct, cols in standardized.attrs["feature_map"].items():
        for col in cols:
            medians = standardized.groupby("taxon_name")[col].median()
            assert np.isclose(medians.std(ddof=0), 1.0)
