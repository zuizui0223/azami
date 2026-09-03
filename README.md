# azami analysis and reproducibility repository

This repository retains the code and frozen machine-readable products needed to reproduce and audit the image-phenomics analyses.

Prepublication manuscript prose and journal-submission packages are intentionally kept outside GitHub.

## Recovery status

The frozen Chapter 1 v2 numerical analysis is fully restored on `main`: canonical code, frozen outputs, direct archived inputs, sampling auxiliaries, historical sensitivity, figures and reproducibility CI are all present and fail-closed. The machine-readable completeness record is `reproducibility/recovery_inventory.json`.

The complete pre-cleanup scientific tree remains permanently recoverable from immutable tag `azami-ch1-v2-2026-08-27` and is also browseable at branch `archive/ch1-precleanup-20260827`. Superseded v1/reviewer/submission-era files are intentionally not reintroduced into the active code-only tree.

A small set of upstream GitHub Actions artifacts is not needed for the frozen v2 numerical rerun but is still useful for full upstream rollback and independent-audit preservation. Their retention deadlines and digests are tracked in issue #85 and `reproducibility/actions_artifact_catalog.json`.

## Tracked scientific material

- `ch1_global/v2/` — image collection, screening, measurement, QC and historical-analysis implementation.
- `analysis/` — statistical analyses, sensitivity analyses and figure-generation code.
- `analysis/ch1/` — frozen analysis contracts and registries needed to define the estimands.
- `analysis_outputs/` — frozen machine-readable result tables and validation reports.
- `reproducibility/figures/` — generated result figures together with the source/provenance material required to rebuild them.
- `reproducibility/contracts/` — compact trait definitions needed to reproduce measurement/analysis decisions.
- `reproducibility/recovery_inventory.json` — active-v2 completeness, historical recovery locations and external artifact gates.
- `tests/` — regression and invariant tests.
- `src/azami_ch1/` — reusable utilities.
