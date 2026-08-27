from __future__ import annotations

import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]


def load_module():
    path = ROOT / "analysis" / "run_geb_v2_full27_sampling_composition_sensitivity.py"
    spec = importlib.util.spec_from_file_location("sampling_composition", path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_stable_top_taxa_uses_count_then_alphabetical_tie_break():
    module = load_module()
    frame = pd.DataFrame(
        {
            "taxon_name": ["b", "b", "a", "a", "c"],
            "obs_id": [1, 2, 3, 4, 5],
        }
    )
    assert module.stable_top_taxa(frame, 3) == ["a", "b", "c"]


def test_equal_taxon_weights_give_each_taxon_unit_total_weight():
    module = load_module()
    taxa = pd.Series(["abundant", "abundant", "abundant", "rare"])
    weights = module.equal_taxon_weights(taxa)
    totals = pd.DataFrame({"taxon": taxa, "weight": weights}).groupby("taxon").weight.sum()
    assert np.allclose(totals.to_numpy(), 1.0)


def test_equal_taxon_weighted_within_fit_differs_from_observation_weighted_fit():
    module = load_module()
    frame = pd.DataFrame(
        {
            "taxon_name": ["a"] * 6 + ["b"] * 2,
            "trait": [-3, -2, -1, 1, 2, 3, -1, 1],
            "env": [-3, -2, -1, 1, 2, 3, 1, -1],
        }
    )
    ordinary = module.fit_within(
        frame,
        ["trait"],
        "env",
        "linear_endpoint",
        minimum_observations=1,
        minimum_taxa=1,
        equal_weight=False,
    )
    equal = module.fit_within(
        frame,
        ["trait"],
        "env",
        "linear_endpoint",
        minimum_observations=1,
        minimum_taxa=1,
        equal_weight=True,
    )
    assert ordinary["status"] == "ok"
    assert equal["status"] == "ok"
    assert ordinary["beta_std_sensitivity"] > 0
    assert abs(equal["beta_std_sensitivity"]) < abs(ordinary["beta_std_sensitivity"])


def test_circular_alignment_requires_less_than_ninety_degree_rotation():
    module = load_module()
    assert np.isclose(module.circular_alignment(1, 0, 1, 0), 1.0)
    assert module.circular_alignment(1, 0, -1, 0) < 0
    assert np.isclose(module.circular_alignment(1, 0, 0, 1), 0.0)


def test_scenario_grid_contains_declared_families_without_dropping_raw_rows():
    module = load_module()
    top = [f"taxon_{index}" for index in range(10)]
    regions = ["Asia", "Europe"]
    within = module.scenario_definitions(top, regions, "within_taxon")
    among = module.scenario_definitions(top, regions, "among_taxon")
    assert len(within) == 1 + 10 + 1 + 2 + 1
    assert len(among) == 10 + 1 + 2 + 1
    assert {row["sensitivity_family"] for row in within} == {
        "equal_taxon_weight",
        "dominant_taxon_omission",
        "leave_one_broad_region_out",
        "native_only",
    }
    assert all(row["scenario"] != "baseline_replacement" for row in within + among)


def test_native_and_region_filters_apply_to_trait_and_environment_scopes():
    module = load_module()
    unit = pd.DataFrame(
        {
            "obs_id": ["1", "2", "3"],
            "taxon_name": ["a", "a", "b"],
            "broad_region": ["Asia", "Europe", "Europe"],
            "native_range_status": ["native", "introduced", "native"],
        }
    )
    metadata = unit.copy()
    scenario = {
        "excluded_taxa": [],
        "excluded_region": "Asia",
        "native_only": True,
    }
    trait_scope, metadata_scope = module.apply_scenario(unit, metadata, scenario)
    assert trait_scope["obs_id"].tolist() == ["3"]
    assert metadata_scope["obs_id"].tolist() == ["3"]
