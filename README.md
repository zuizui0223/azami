# azami analysis and reproducibility repository

This repository retains the code and frozen machine-readable products needed to audit the Chapter 1 v2 image-phenomics analysis.

Prepublication manuscript prose and journal-submission packages are intentionally kept outside GitHub.

## Reproducibility status

The frozen Chapter 1 v2 analysis itself is fully recovered and internally reproducible by the owner: canonical code, frozen outputs, SHA-locked analysis inputs, sampling auxiliaries, historical sensitivity, figures, and integrity CI are preserved.

**Independent third-party full numerical reproduction is not yet complete.** The required analysis-input ZIPs are preserved in owner-controlled storage, but they have not yet been deposited in a public credential-free data archive. Therefore a reader can currently audit the public code, frozen results, validation reports, and figures, but cannot perform the complete numerical rerun without obtaining the unpublished input package.

The public-release gate is machine-readable in `reproducibility/public_release_manifest.json`. The canonical third-party procedure is `reproducibility/README.md`. Do not claim full independent reproducibility until the public release manifest contains a stable public DOI/URL and immutable code ref and the released files match the recorded SHA-256 values.

The complete pre-cleanup scientific tree remains recoverable from immutable tag `azami-ch1-v2-2026-08-27` and branch `archive/ch1-precleanup-20260827`.

## Tracked scientific material

- `ch1_global/v2/` — image collection, screening, measurement, QC and historical-analysis implementation.
- `analysis/` — statistical analyses, sensitivity analyses and figure-generation code.
- `analysis/ch1/` — frozen analysis contracts and registries.
- `analysis_outputs/` — frozen machine-readable result tables and validation reports.
- `reproducibility/README.md` — independent third-party reproduction procedure.
- `reproducibility/public_release_manifest.json` — public-data release gate and required file hashes.
- `reproducibility/material_availability.json` — audit of public GitHub material, owner preservation material, Git-history inputs and public external sources.
- `reproducibility/figures/` — frozen figures and tracked figure source material.
- `reproducibility/recovery_inventory.json` — owner-side recovery completeness and historical recovery locations.
- `reproducibility/durable_archive_manifest.json` — owner preservation archive identities; not a public distribution endpoint.
- `tests/` — regression and invariant tests.
- `src/azami_ch1/` — reusable utilities.
