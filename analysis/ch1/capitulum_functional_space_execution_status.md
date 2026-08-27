# Capitulum functional-space execution status — 2026-08-27

## Frozen before process-extension outcomes

The following files were committed before any results from the new CHELSA process variables were inspected:

- `manuscript/CAPITULUM_FUNCTIONAL_SPACE_HYPOTHESES_2026-08-27.md`;
- `analysis/ch1/capitulum_environment_blocks_contract.json`;
- `analysis/ch1/chelsa_process_environment_sources.json`;
- `analysis/ch1/capitulum_function_annotation.csv`;
- `analysis/run_capitulum_functional_space.py`;
- `analysis/run_capitulum_environment_blocks.py`;
- `analysis/sample_chelsa_process_environment.py`.

The new process variables are mean shortwave radiation, mean VPD, mean near-surface wind, growing-season precipitation and potential NPP. They enter only through the predeclared environmental blocks. No all-BIOCLIM fishing screen is permitted.

## Already available baseline data

The immutable final GEB-v2 artifact contains the 18 measured inferential endpoints and the frozen four CHELSA predictors (BIO1, BIO4, BIO12, BIO15). Local dry runs against that artifact reproduced:

- 1,874 observations complete for all 18 endpoints from 124 taxa;
- main >=5 complete observations/taxon scope: 1,734 observations / 42 taxa;
- >=2 sensitivity: 1,825 observations / 75 taxa;
- main within-module association contrast: 0.164502;
- main among-module association contrast: 0.088475;
- main within-vs-among association-matrix Spearman: 0.366299;
- main whole-space thermal-block within-taxon multivariate R2: 0.008510.

These local values are provisional execution checks until a GitHub Actions artifact records the exact committed code, environment, contracts and source artifact ID.

## Process-extension outcome status

**Not yet interpreted.** The first remote artifact-backed run must complete before any shortwave-radiation, VPD, wind, growing-season-precipitation or NPP result is added to the manuscript or EAzami handoff.

If a frozen CHELSA source is inaccessible or fails the minimum-coverage gate, the workflow must fail. The predictor may not be replaced after observing trait outcomes.
