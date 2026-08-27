from __future__ import annotations

import importlib.util
from pathlib import Path
import unittest

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]


def load_module(name: str, relative: str):
    spec = importlib.util.spec_from_file_location(name, ROOT / relative)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


INC = load_module("cap_env_incremental", "analysis/run_capitulum_environment_incremental.py")
COHORT = load_module("cap_env_cohort", "analysis/build_capitulum_environment_cohort.py")


class CapitulumEnvironmentIncrementalTests(unittest.TestCase):
    def test_complete18_rule_uses_completeness_not_trait_magnitude(self):
        endpoints = [f"e{i:02d}" for i in range(18)]
        rows = []
        for obs_id, scale in [("a", 1.0), ("b", 1000000.0), ("c", -1000000.0)]:
            for j, endpoint in enumerate(endpoints):
                if obs_id == "c" and j == 17:
                    continue
                rows.append({
                    "obs_id": obs_id,
                    "endpoint_id": endpoint,
                    "analysis_tier": "primary" if j < 9 else "candidate",
                    "analysis_eligible": True,
                    "value": scale * (j + 1),
                })
        traits = pd.DataFrame(rows)
        ids, observed_endpoints = COHORT.complete18_ids(traits)
        self.assertEqual(observed_endpoints, endpoints)
        self.assertEqual(set(ids), {"a", "b"})

    def test_nested_effect_detects_extension_beyond_core(self):
        rng = np.random.default_rng(20260827)
        n = 160
        core = rng.normal(size=(n, 4))
        extension = rng.normal(size=(n, 1))
        signal = extension @ np.linspace(0.4, 1.2, 18)[None, :]
        y = signal + rng.normal(scale=0.35, size=(n, 18)) + core @ rng.normal(scale=0.05, size=(4, 18))
        yz = INC.standardize(y)
        xc = INC.standardize(core)
        xe = INC.standardize(extension)
        r2c, r2f, delta, partial, _, _ = INC.nested_effect_ols(yz, xc, xe)
        self.assertGreater(r2f, r2c)
        self.assertGreater(delta, 0.25)
        self.assertGreater(partial, 0.25)

    def test_within_permutation_never_moves_rows_between_taxa(self):
        values = np.arange(24).reshape(12, 2)
        groups = np.array(["a"] * 4 + ["b"] * 4 + ["c"] * 4)
        shuffled = INC.permute_rows_within(values, groups, np.random.default_rng(5))
        for group in np.unique(groups):
            idx = np.flatnonzero(groups == group)
            self.assertEqual(
                {tuple(x) for x in shuffled[idx]},
                {tuple(x) for x in values[idx]},
            )


if __name__ == "__main__":
    unittest.main()
