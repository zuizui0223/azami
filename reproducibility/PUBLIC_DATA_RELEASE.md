# Public data release procedure

This is the author-side publication procedure that converts the preserved Chapter 1 v2 inputs into a credential-free third-party reproducibility package.

## Prepared bundle

The prepared staging file is:

```text
azami_ch1_v2_reproduction_inputs_2026-09-04.zip
```

Expected properties:

```text
size:    56,942,044 bytes
SHA-256: 50ec15b1280d4660839ca4bf0d55c970a5f49b4d4feaabb7a073b73500253677
```

It contains the exact four minimum binary inputs recorded in `reproducibility/public_release_manifest.json`, plus a SHA ledger, public-bundle manifest, release README and Zenodo metadata. Each embedded binary was re-read from the finished bundle and verified against its frozen SHA-256 before this contract was committed.

The immutable code ref paired with this bundle is:

```text
584af97b050d15701f26ce1facea212d5b648d4d
```

## Deposit

1. Create a dataset deposit in a credential-free public research-data repository that provides a stable URL or DOI. A DOI-bearing repository such as Zenodo, OSF or Figshare is preferred.
2. Upload `azami_ch1_v2_reproduction_inputs_2026-09-04.zip` without modifying or recompressing it.
3. Use `reproducibility/zenodo_metadata.json` as the prepared metadata template. Do not invent an ORCID or a blanket data license; add those only after the author has explicitly selected them and confirmed compatibility with the underlying source-data terms.
4. Publish the deposit.

## Independent post-publication verification

Do not mark third-party full numerical reproducibility ready immediately after upload. From a clean session that is not authenticated to the author account:

1. download the public bundle;
2. verify the bundle SHA-256 above;
3. extract it;
4. clone `https://github.com/zuizui0223/azami` and checkout `584af97b050d15701f26ce1facea212d5b648d4d`;
5. run `python reproducibility/verify_local_materials.py --archive-dir <extracted-directory> --scope numerical`;
6. follow `reproducibility/README.md` through the canonical rebuild;
7. require `PASS 24/24`, `PASS 7/7`, `PASS 11/11`, and figure build `status: ok` with 5 main and 7 supporting figures.

## Close the publication gate

Only after the credential-free re-download and validation succeed:

- set `public_data_release.doi` and `public_data_release.url` in `reproducibility/public_release_manifest.json`;
- set `public_data_release.publicly_downloadable_without_owner_credentials` to `true`;
- change the public-release and third-party statuses from blocked/prepared to ready;
- update `reproducibility/material_availability.json` and `reproducibility/recovery_inventory.json` consistently;
- close GitHub issue #88;
- run `Reproducibility integrity` on `main` and require success.

The five frozen CHELSA v2.1 process rasters may optionally be added to the deposit for a fully offline rerun. If omitted, their exact public URLs remain frozen in `analysis/ch1/chelsa_process_environment_sources.json`, and the reconstructed nine-predictor environment is still required to match SHA-256 `e242aa7ce69d12b11937c1335e84b9638799c50b42ef36b95725e77190df98e7`.
