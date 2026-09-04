# azami analysis and reproducibility repository

This repository retains the code and frozen machine-readable products needed to reproduce and audit the image-phenomics analyses.

Prepublication manuscript prose and journal-submission packages are intentionally kept outside GitHub.

## Third-party reproducibility status

The frozen Chapter 1 v2 analysis code, contracts, frozen outputs, validation reports and figure sources are public and auditable in this repository. The independent-reader procedure is [`reproducibility/README.md`](reproducibility/README.md).

A public-release staging bundle containing the exact four minimum numerical input archives has now been built and SHA-verified. Its identity is recorded in [`reproducibility/public_release_manifest.json`](reproducibility/public_release_manifest.json):

```text
azami_ch1_v2_reproduction_inputs_2026-09-04.zip
56,942,044 bytes
SHA-256 50ec15b1280d4660839ca4bf0d55c970a5f49b4d4feaabb7a073b73500253677
```

The immutable code ref paired to that bundle is `584af97b050d15701f26ce1facea212d5b648d4d`. Prepared public-deposition metadata is in [`reproducibility/zenodo_metadata.json`](reproducibility/zenodo_metadata.json), and the author-side publication/clean-session verification procedure is [`reproducibility/PUBLIC_DATA_RELEASE.md`](reproducibility/PUBLIC_DATA_RELEASE.md).

**Full independent third-party numerical reproduction is not yet claimed.** The remaining publication gate is to place that prepared bundle in a credential-free public repository with a stable URL or DOI, independently re-download it, verify the bundle and embedded SHA-256 values, and rerun the frozen validations. This gate is tracked in GitHub issue #88.

The five frozen CHELSA v2.1 process rasters may optionally be included in that public deposit for a fully offline rerun. If not, their exact public URLs remain frozen in the source registry and the reconstructed nine-predictor environment must match its recorded SHA-256.

## Owner-side recovery status

Separately from public distribution, the source artifacts are durably preserved in owner-controlled storage, and the complete pre-cleanup scientific tree remains permanently recoverable from immutable tag `azami-ch1-v2-2026-08-27` and branch `archive/ch1-precleanup-20260827`. These preservation locations are recovery infrastructure, not the public third-party distribution endpoint.

## Tracked scientific material

- `ch1_global/v2/` — image collection, screening, measurement, QC and historical-analysis implementation.
- `analysis/` — statistical analyses, sensitivity analyses and figure-generation code.
- `analysis/ch1/` — frozen analysis contracts and registries needed to define the estimands.
- `analysis_outputs/` — frozen machine-readable result tables and validation reports.
- `reproducibility/figures/` — generated result figures together with source/provenance material required to rebuild them.
- `reproducibility/contracts/` — compact trait definitions needed to reproduce measurement/analysis decisions.
- `reproducibility/public_release_manifest.json` — public third-party release gate, minimum input hashes, prepared bundle identity and code ref.
- `reproducibility/material_availability.json` — public/private/external material boundary.
- `reproducibility/recovery_inventory.json` — owner recovery completeness and separate public-release status.
- `reproducibility/durable_archive_manifest.json` — owner-side preservation identities and reassembly contract.
- `tests/` — regression and invariant tests.
- `src/azami_ch1/` — reusable utilities.
