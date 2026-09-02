from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "analysis" / "run_geb_v2_full27_environment_atlas.py"
CONTRACT = ROOT / "ch1_global" / "v2" / "ontology" / "ch1_continuous_trait_contract.csv"
ANALYSIS_CONTRACT = ROOT / "analysis" / "ch1" / "v2_full27_environment_atlas_contract.json"


def load_module():
    spec = importlib.util.spec_from_file_location("v2_full27_environment_atlas", SCRIPT)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_contract_starts_from_27_and_yields_26_joint_units():
    module = load_module()
    contract, analysis_contract = module.load_contracts(CONTRACT, ANALYSIS_CONTRACT)
    units = module.inferential_units(contract)
    assert len(contract) == 27
    assert len(units) == 26
    hue = [unit for unit in units if unit["unit_id"] == "corolla_hue"]
    assert len(hue) == 1
    assert set(hue[0]["member_endpoint_ids"]) == {"corolla_hue_sin", "corolla_hue_cos"}
    assert analysis_contract["trait_universe"]["old_v1_nine_endpoint_subset_used"] is False


def test_all_nine_predictors_are_covered_once_by_six_blocks():
    module = load_module()
    analysis_contract = json.loads(ANALYSIS_CONTRACT.read_text(encoding="utf-8"))
    mapping = module.environment_block_map(analysis_contract)
    assert list(mapping) == module.EXPECTED_PREDICTORS
    assert len(set(mapping.values())) == 6
    assert mapping["chelsa_sfcwind_mean"] == "mechanical"
    assert mapping["chelsa_npp"] == "resource_productivity"


def test_inventory_retains_unexecuted_and_validation_only_endpoints():
    module = load_module()
    contract = pd.read_csv(CONTRACT, dtype=str, keep_default_na=False)
    measured_ids = [
        "orientation_image_vertical_angle",
        "involucre_surface_edge_density",
    ]
    traits = pd.DataFrame(
        {
            "obs_id": ["1", "1"],
            "taxon_name": ["Cirsium one", "Cirsium one"],
            "endpoint_id": measured_ids,
            "value": [1.0, 2.0],
        }
    )
    environment = pd.DataFrame({"obs_id": ["1"], "taxon_name": ["Cirsium one"]})
    inventory = module.build_inventory(contract, traits, environment)
    assert len(inventory) == 27
    status = inventory.set_index("endpoint_id")["measurement_status"]
    assert status["involucre_surface_edge_density"] == "measured"
    assert status["visible_floret_fraction"] == "unexecuted_no_measurement"
    assert inventory["registered_in_v2_starting_27"].all()


def test_variance_decomposition_reports_information_below_taxon_means():
    module = load_module()
    contract = pd.read_csv(CONTRACT, dtype=str, keep_default_na=False).head(1)
    endpoint = contract.iloc[0]["endpoint_id"]
    traits = pd.DataFrame(
        {
            "endpoint_id": [endpoint] * 4,
            "taxon_name": ["a", "a", "b", "b"],
            "value": [0.0, 2.0, 10.0, 12.0],
        }
    )
    result = module.build_variance_decomposition(contract, traits, seed=1, bootstrap_repeats=10)
    assert result.iloc[0]["status"] == "ok"
    assert np.isclose(result.iloc[0]["fraction_below_taxon_means"], 4.0 / 104.0)
    assert np.isclose(result.iloc[0]["fraction_retained_between_taxa"], 100.0 / 104.0)


def test_fdr_family_does_not_split_or_remove_rows_by_old_tier():
    module = load_module()
    frame = pd.DataFrame(
        {
            "status": ["ok", "ok", "unexecuted_no_measurement"],
            "p_value": [0.001, 0.04, np.nan],
            "analysis_tier_metadata_only": ["primary", "validation_only", "descriptive_only"],
        }
    )
    result = module.add_global_fdr(frame)
    assert len(result) == 3
    assert result["result_retained_regardless_of_q"].all()
    assert np.isclose(result.loc[0, "q_fdr_bh_global_family"], 0.002)
    assert np.isclose(result.loc[1, "q_fdr_bh_global_family"], 0.04)
    assert pd.isna(result.loc[2, "q_fdr_bh_global_family"])


def test_batched_within_taxon_circular_permutation_preserves_null_scheme():
    module = load_module()
    x = np.array([-1.0, 1.0, -1.0, 1.0])
    sine = x.copy()
    cosine = np.zeros_like(x)
    taxa = np.array(["a", "a", "b", "b"])
    result = module.permutation_p_circular_within(
        sine,
        cosine,
        x,
        taxa,
        permutations=1999,
        rng=np.random.default_rng(7),
    )
    # Each two-member taxon can independently retain or reverse its order.
    # Half of all joint permutations reproduce the observed absolute slope.
    assert np.isclose(result["effect_magnitude"], 1.0)
    assert 0.45 <= result["p_value"] <= 0.55
    assert result["p_value_method"] == "within_taxon_predictor_permutation_1999"
