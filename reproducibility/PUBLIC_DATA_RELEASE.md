# Public data release procedure

This document records the author-side publication procedure that converts the preserved Chapter 1 v2 inputs into a credential-free third-party reproducibility package.

## Published Zenodo v1

The prepared bundle has been deposited and published as:

```text
DOI:      10.5281/zenodo.22295791
record:   https://zenodo.org/records/22295791
version:  v1
bundle:   azami_ch1_v2_reproduction_inputs_2026-09-04.zip
size:     56,942,044 bytes
SHA-256:  50ec15b1280d4660839ca4bf0d55c970a5f49b4d4feaabb7a073b73500253677
```

The bundle contains the exact four minimum binary inputs recorded in `reproducibility/public_release_manifest.json`, plus a SHA ledger, public-bundle manifest, release README and Zenodo metadata. Each embedded binary was re-read from the finished bundle and verified against its frozen SHA-256 before publication.

The immutable code ref paired with this bundle is:

```text
584af97b050d15701f26ce1facea212d5b648d4d
```

## Publication status

The DOI publication step is complete. The release must **not** yet be marked fully third-party reproducible solely because the DOI exists.

Current state:

```text
published_pending_anonymous_redownload_verification
```

## Independent post-publication verification

From a clean session that is not authenticated to the author account:

1. open `https://doi.org/10.5281/zenodo.22295791`;
2. download `azami_ch1_v2_reproduction_inputs_2026-09-04.zip`;
3. require bundle size `56,942,044` bytes and SHA-256 `50ec15b1280d4660839ca4bf0d55c970a5f49b4d4feaabb7a073b73500253677`;
4. extract it;
5. clone `https://github.com/zuizui0223/azami` and checkout `584af97b050d15701f26ce1facea212d5b648d4d`;
6. run `python reproducibility/verify_local_materials.py --archive-dir <extracted-directory> --scope numerical`;
7. require all four embedded input SHA-256 values to match `reproducibility/public_release_manifest.json`;
8. follow `reproducibility/README.md` through the canonical rebuild;
9. require `PASS 24/24`, `PASS 7/7`, `PASS 11/11`, and figure build `status: ok` with 5 main and 7 supporting figures.

## Close the publication gate

Only after the anonymous credential-free re-download and validation succeed:

- set `public_data_release.publicly_downloadable_without_owner_credentials` to `true`;
- set `public_data_release.anonymous_redownload_verified` to `true`;
- change the public-release and third-party statuses to ready;
- update `reproducibility/material_availability.json` and `reproducibility/recovery_inventory.json` consistently;
- close GitHub issue #88;
- run `Reproducibility integrity` on `main` and require success.

The five frozen CHELSA v2.1 process rasters may optionally be added to a future Zenodo version for a fully offline rerun. If omitted, their exact public URLs remain frozen in `analysis/ch1/chelsa_process_environment_sources.json`, and the reconstructed nine-predictor environment is still required to match SHA-256 `e242aa7ce69d12b11937c1335e84b9638799c50b42ef36b95725e77190df98e7`.
