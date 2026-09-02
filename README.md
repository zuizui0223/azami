# azami analysis and reproducibility repository

This repository retains the code and frozen machine-readable products needed to reproduce and audit the image-phenomics analyses.

Prepublication manuscript prose and journal-submission packages are intentionally kept outside GitHub.

## Tracked scientific material

- `ch1_global/v2/` — image collection, screening, measurement, QC and historical-analysis implementation.
- `analysis/` — statistical analyses, sensitivity analyses and figure-generation code.
- `analysis/ch1/` — frozen analysis contracts and registries needed to define the estimands.
- `analysis_outputs/` — frozen machine-readable result tables and validation reports.
- `reproducibility/figures/` — generated result figures together with the source/provenance material required to rebuild them.
- `reproducibility/contracts/` — compact trait definitions needed to reproduce measurement/analysis decisions.
- `tests/` — regression and invariant tests.
- `src/azami_ch1/` — reusable utilities.

