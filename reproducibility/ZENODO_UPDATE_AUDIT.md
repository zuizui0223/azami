# Zenodo update audit

## Decision — 2026-09-27

Zenodo is the **durable analysis-input store** for Chapter 1. It is not a mirror of the GitHub repository.

The current v3 Zenodo package must contain only the exact processed numerical inputs needed by the GitHub replay:

| Input | Frozen identity | Zenodo treatment |
|---|---|---|
| Continuous image-derived measurements | Actions artifact `9612943217` | preserve exact artifact ZIP |
| Nine-predictor environment | Actions artifact `9633419268` | preserve exact artifact ZIP |
| Broad-region lookup | Actions artifact `8983877726` | preserve exact artifact ZIP |
| Historical-placement trees | Actions artifact `8227254443` | preserve exact artifact ZIP |
| Observation native status | exact SHA `c01eeb9ff245d7f73da1a12fa4eede904dd9770467655f20e3d85de2ac8dd84a` | preserve exact normalized CSV recovered from checksum-verified artifact `10292140117` |

A small `README.txt`, `input_contract.json`, `release_manifest.json`, `checksums.json` and final Zenodo metadata are allowed because they describe and verify the data bytes.

## Explicit exclusions

The Zenodo data package does **not** contain:

- analysis code or a Git repository snapshot;
- Main manuscript, Supporting Information, title page or cover letter;
- figures or figure provenance exports;
- the 16 current reference outputs;
- replay receipts;
- fitted-model outputs;
- original third-party photographs;
- detector training images, annotations or weights.

Those materials remain versioned or documented in GitHub where appropriate. The Zenodo manifest pins the GitHub repository and exact code commit used with the inputs.

## Existing public v2 record

The existing public record remains unchanged:

- record DOI: `10.5281/zenodo.22295791`;
- concept DOI: `10.5281/zenodo.22295790`;
- title: *Azami Chapter 1 v2 reproducibility input package*;
- file: `azami_ch1_v2_reproduction_inputs_2026-09-04.zip`;
- size: 56,942,044 bytes;
- SHA-256: `50ec15b1280d4660839ca4bf0d55c970a5f49b4d4feaabb7a073b73500253677`;
- creator metadata: `ZHANG, Ruiqi`.

The v2 ZIP already demonstrates the intended archive role: preserve exact numerical inputs, while executable analysis history lives in GitHub.

## Why a new version is needed

The existing v2 package does not contain the current nine-predictor environment input and does not expose the current frozen native-status CSV as a standalone replay input. The new v3 version therefore updates the **input bytes**, not the code or manuscript.

The current builder `reproducibility/build_current_release_bundle.py` verifies the exact artifact identities defined by `reproducibility.run_current_analysis.INPUTS`, extracts and verifies native status, and builds a deterministic data-only ZIP.

## Durable private recovery state

All current-only material needed to recover the five inputs has already been copied away from expiring Actions storage and checksum recorded in `CURRENT_RELEASE_STAGING_20260912.json`. Additional result/figure artifacts may remain in the private owner archive for project durability, but they are **not part of the Zenodo v3 public package**.

## Reproducibility model

The public reproduction path is intentionally split:

```text
Zenodo DOI
  └─ exact processed numerical inputs
          +
GitHub pinned commit
  └─ analysis code + runbook + reference outputs
          ↓
8-stage numerical replay
          ↓
16-file reference comparison
```

This is sufficient for the supported numerical-reproduction claim. It does not claim reconstruction from original photographs.

## Release metadata

Prepared metadata: `reproducibility/zenodo_release_metadata.prepared.json`.

The dataset creator is inherited from the existing v2 concept as `ZHANG, Ruiqi`. This is a dataset-deposition creator field and does not define manuscript authorship.

Prepared title:

> Azami Chapter 1 v3 numerical analysis input package

The prepared metadata remains `release_approved=false` until the final data-only wording, exact GitHub code commit and licensing statement are approved.

## Remaining public-release gate

Only one pre-build gate remains:

- `release_metadata`.

After approval:

1. bind metadata to the exact final GitHub commit;
2. build the final data-only ZIP with `release_ready=true`;
3. publish it as a new version under concept DOI `10.5281/zenodo.22295790`;
4. download the published package without owner credentials;
5. verify every checksum;
6. run the GitHub replay against those downloaded inputs and confirm the 16-file reference comparison passes.

A private Drive copy or Actions artifact is durability evidence, not a public-release substitute.
