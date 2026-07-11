from __future__ import annotations

import argparse
import importlib.util
import unittest
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "ch1_global" / "v2" / "71_build_within_species_expansion_queue.py"
spec = importlib.util.spec_from_file_location("within_expansion", SCRIPT)
assert spec is not None and spec.loader is not None
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


class TestWithinSpeciesExpansion(unittest.TestCase):
    def fixture(self) -> pd.DataFrame:
        rows = []
        for species, n, blocks in (("Cirsium alpha", 80, 10), ("Cirsium beta", 30, 6)):
            for i in range(n):
                block = i % blocks
                rows.append({
                    "obs_id": f"{species}_{i}",
                    "photo_id": f"p_{species}_{i}",
                    "photo_index": "0",
                    "taxon_name": species,
                    "taxon_rank": "species",
                    "quality_grade": "research",
                    "captive": "false",
                    "latitude": str(-10 + block * 2.1),
                    "longitude": str(20 + block * 2.1),
                    "coordinate_usable_for_environment": "true",
                    "positional_accuracy": "1000",
                    "geoprivacy": "open",
                    "obscured": "false",
                    "observed_on": "2026-01-01",
                    "small_image_url": "s",
                    "medium_image_url": "m",
                    "large_image_url": "l",
                    "source_file_stem": f"stem_{species}_{i}",
                })
        return pd.DataFrame(rows)

    def args(self) -> argparse.Namespace:
        return argparse.Namespace(
            quality_grades="research",
            max_positional_accuracy_m=10_000,
            min_observations_per_species=60,
            min_spatial_blocks_per_species=8,
            max_observations_per_species=70,
            max_per_species_spatial_block=10,
            spatial_block_deg=2.0,
            target_species="",
            seed=1,
        )

    def test_only_well_sampled_species_is_selected(self) -> None:
        queue, audit = module.build_queue(self.fixture(), self.args())
        self.assertEqual(set(queue["taxon_name"]), {"Cirsium alpha"})
        self.assertEqual(len(queue), 70)
        self.assertGreaterEqual(queue["spatial_block"].nunique(), 8)
        beta = audit.loc[audit["taxon_name"].eq("Cirsium beta")].iloc[0]
        self.assertFalse(bool(beta["eligible_for_expansion"]))

    def test_selection_is_reproducible(self) -> None:
        first, _ = module.build_queue(self.fixture(), self.args())
        second, _ = module.build_queue(self.fixture(), self.args())
        self.assertEqual(list(first["obs_id"]), list(second["obs_id"]))

    def test_trait_columns_do_not_affect_selection(self) -> None:
        data = self.fixture()
        data["orientation_angle_degrees"] = range(len(data))
        selected, _ = module.build_queue(data, self.args())
        self.assertNotIn("orientation_angle_degrees", selected.columns)


if __name__ == "__main__":
    unittest.main()
