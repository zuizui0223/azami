from __future__ import annotations

import importlib.util
import unittest
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]


def load(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


PRIMARY = load(ROOT / "analysis" / "run_bias_control_reanalysis.py", "bias_control_primary")
INVOLUCRE = load(ROOT / "analysis" / "run_involucre_resolution_audit.py", "bias_control_involucre")
HEMISPHERE = load(
    ROOT / "analysis" / "run_hemisphere_season_sensitivity.py",
    "bias_control_hemisphere",
)
NATIVE = load(
    ROOT / "analysis" / "run_native_range_sensitivity.py",
    "bias_control_native_range",
)


class TestBiasControlReanalysis(unittest.TestCase):
    def test_cyclic_residualization_removes_taxon_specific_seasonality(self) -> None:
        rows = []
        rng = np.random.default_rng(20260826)
        for taxon_index, taxon in enumerate(("Cirsium alpha", "Cirsium beta", "Cirsium gamma")):
            for day in np.linspace(1, 365, 60, endpoint=False):
                phase = 2 * np.pi * (day - 1) / 365
                seasonal = np.sin(phase + taxon_index * 0.3)
                rows.append({
                    "taxon_name": taxon,
                    "observed_on": pd.Timestamp("2025-01-01") + pd.Timedelta(days=float(day - 1)),
                    "trait": seasonal + rng.normal(0, 0.03),
                    "climate": seasonal + rng.normal(0, 0.03),
                })
        frame = PRIMARY.add_cyclic_date(pd.DataFrame(rows))
        baseline = PRIMARY.fit_model(frame, "trait", "climate", 5, "taxon_mean_only")
        adjusted = PRIMARY.fit_model(
            frame, "trait", "climate", 5, "taxon_specific_cyclic_doy"
        )
        self.assertGreater(abs(float(baseline["beta_std_within"])), 0.9)
        self.assertLess(abs(float(adjusted["beta_std_within"])), 0.2)

    def test_omission_scenarios_are_count_ranked_and_outcome_blind(self) -> None:
        frame = pd.DataFrame({
            "taxon_name": ["b"] * 3 + ["a"] * 5 + ["c"] * 2,
            "unused_outcome": [1000, -1, 5, 3, 2, 1, 0, 999, 998, 997],
        })
        concentration, scenarios = PRIMARY.omission_scenarios(frame, 2)
        self.assertEqual(concentration.iloc[0]["taxon_name"], "a")
        self.assertEqual(concentration.iloc[1]["taxon_name"], "b")
        self.assertEqual(scenarios[-1]["omitted_taxa"], ["a", "b"])

    def test_southern_dates_are_shifted_by_half_cycle(self) -> None:
        frame = pd.DataFrame({
            "taxon_name": ["Cirsium alpha", "Cirsium alpha"],
            "observed_on": ["2025-01-01", "2025-01-01"],
            "latitude": [35.0, -35.0],
        })
        result = HEMISPHERE.add_hemisphere_season_terms(frame)
        self.assertAlmostEqual(
            float(result.iloc[0]["aligned_doy_sin"]),
            -float(result.iloc[1]["aligned_doy_sin"]),
            places=12,
        )
        self.assertAlmostEqual(
            float(result.iloc[0]["aligned_doy_cos"]),
            -float(result.iloc[1]["aligned_doy_cos"]),
            places=12,
        )

    def test_hemisphere_specific_curves_remove_opposite_calendar_seasons(self) -> None:
        rows = []
        rng = np.random.default_rng(17)
        for taxon in ("Cirsium alpha", "Cirsium beta", "Cirsium gamma"):
            for south, latitude in ((0, 40.0), (1, -40.0)):
                for day in np.linspace(1, 365, 36, endpoint=False):
                    phase = 2 * np.pi * (day - 1) / 365
                    seasonal = np.sin(phase + south * np.pi)
                    rows.append({
                        "taxon_name": taxon,
                        "observed_on": pd.Timestamp("2025-01-01") + pd.Timedelta(days=float(day - 1)),
                        "latitude": latitude,
                        "trait": seasonal + rng.normal(0, 0.03),
                        "climate": seasonal + rng.normal(0, 0.03),
                    })
        frame = HEMISPHERE.add_hemisphere_season_terms(pd.DataFrame(rows))
        adjusted = HEMISPHERE.fit_model(
            frame, "trait", "climate", 5, "taxon_cyclic_doy_by_hemisphere"
        )
        self.assertLess(abs(float(adjusted["beta_std_within"])), 0.2)

    def test_wcvp_name_resolution_fails_closed_on_multiple_targets(self) -> None:
        payload = {"results": [
            {"canonicalName": "Cirsium alpha", "taxonomicStatus": "ACCEPTED", "key": 1},
            {"canonicalName": "Cirsium alpha", "taxonomicStatus": "SYNONYM", "acceptedKey": 2},
        ]}
        result = NATIVE.resolve_exact_name("Cirsium alpha", payload)
        self.assertEqual(result["resolution_status"], "unresolved")

    def test_wcvp_name_resolution_accepts_convergent_exact_records(self) -> None:
        payload = {"results": [
            {"canonicalName": "Cirsium alpha", "taxonomicStatus": "ACCEPTED", "key": 7, "scientificName": "Cirsium alpha A."},
            {"canonicalName": "Cirsium alpha", "taxonomicStatus": "SYNONYM", "acceptedKey": 7, "accepted": "Cirsium alpha A."},
        ]}
        result = NATIVE.resolve_exact_name("Cirsium alpha", payload)
        self.assertEqual(result["resolution_status"], "resolved_unique_accepted_key")
        self.assertEqual(result["accepted_key"], 7)

    def test_involucre_preparation_retains_quality_covariates(self) -> None:
        head = pd.DataFrame([
            {
                "obs_id": "1", "taxon_name": "Cirsium alpha", "involucre_status": "usable",
                "original_min_dimension": 180, "original_sharpness": 90,
                "involucre_projection_roughness": 0.1, "involucre_spread_fraction": 0.2,
                "spine_relative_length_max_proxy": 0.3,
            },
            {
                "obs_id": "1", "taxon_name": "Cirsium alpha", "involucre_status": "usable",
                "original_min_dimension": 220, "original_sharpness": 110,
                "involucre_projection_roughness": 0.2, "involucre_spread_fraction": 0.3,
                "spine_relative_length_max_proxy": 0.4,
            },
        ])
        environment = pd.DataFrame([{
            "obs_id": "1", "taxon_name": "Cirsium alpha",
            "coordinate_precision_tier": "high_le_1km",
            "env_chelsa_bio01_native": 1, "env_chelsa_bio04_native": 2,
            "env_chelsa_bio12_native": 3, "env_chelsa_bio15_native": 4,
        }])
        result = INVOLUCRE.prepare_observations(head, environment)
        self.assertEqual(len(result), 1)
        self.assertAlmostEqual(float(result.iloc[0]["min_dimension_median"]), 200)
        self.assertEqual(result.iloc[0]["resolution_stratum"], "200_299")
        self.assertTrue(np.isfinite(float(result.iloc[0]["log1p_sharpness"])))


if __name__ == "__main__":
    unittest.main()
