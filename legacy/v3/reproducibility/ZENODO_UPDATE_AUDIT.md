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

## Verified v3 data-only staging bundle — 2026-09-27

The data-only builder passed on main commit `810637d9dfdde9c4b810c2472f43f254786ce614`.

- workflow run: `36287521925`;
- Actions artifact: `10920813215`;
- Actions wrapper SHA-256: `0afaed34da6654819bf979e676f50f877a816beac1d281fe92b6735eee5f6faa`;
- inner archive: `azami_ch1_v3_analysis_inputs.zip`;
- inner archive size: 59,461,974 bytes;
- inner archive SHA-256: `c4ec876206e6f8ef0cd69d126fa31b2b71aa1b1919172aeb371f75bbabcd6d10`;
- durable owner copy: Drive file `1Bqj_5pLAd7DQmC5OLs9x426i8XwHsgZh`.

The inner archive contains exactly the five analysis inputs plus `README.txt`, `input_contract.json`, `release_manifest.json` and `checksums.json`. Inspection confirmed:

- `code_included=false`;
- `manuscript_files_included=false`;
- `figures_included=false`;
- `reference_outputs_included=false`;
- `replay_receipts_included=false`.

The only remaining release gap reported by the builder is `release_metadata`.

## Release-ready v3 data-only candidate — 2026-09-27

The final data-only builder passed on main commit `7483ae8e20fa5df2439ad524312d3ac9a0fbf71f`.

- workflow run: `36288310738`;
- Actions artifact: `10921221306`;
- wrapper SHA-256: `941c70768664ae9df82cd18f522dda9f4fd352ad723e08ad1e9a8d9e23f0996e`;
- final upload file: `azami_ch1_v3_analysis_inputs.zip`;
- final upload SHA-256: `be57e9ba80773c97afdc200f9a3982eaa85044cfc6791b62187ccfcd41f54438`;
- durable owner copy: Drive file `16pildr5nOL08HAb6ATEHepUPei5oEn0t`.

The inner archive was unpacked independently after CI. All embedded checksum rows matched. It contains exactly ten files: the five analysis inputs plus `README.txt`, `input_contract.json`, `checksums.json`, `release_manifest.json`, and `zenodo_release_metadata.json`.

The final manifest reports:

- `release_ready=true`;
- `release_gaps=[]`;
- `analysis_input_count=5`;
- `code_included=false`;
- `manuscript_files_included=false`;
- `figures_included=false`;
- `reference_outputs_included=false`;
- `replay_receipts_included=false`.

The generated final metadata binds the archive to GitHub commit `7483ae8e20fa5df2439ad524312d3ac9a0fbf71f`, creator `ZHANG, Ruiqi`, archive strategy `new_version_existing_concept`, and record license `CC-BY-4.0`.

Therefore all repository-side preparation gates are closed. The archive is **release-ready but not yet public**.

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

Repository-side preparation is complete. The only remaining external step is to publish the verified file `azami_ch1_v3_analysis_inputs.zip` as a new version under concept DOI `10.5281/zenodo.22295790`.

After publication:

1. download the published ZIP without owner credentials;
2. verify outer SHA-256 `be57e9ba80773c97afdc200f9a3982eaa85044cfc6791b62187ccfcd41f54438`;
3. verify the embedded checksum manifest;
4. run the GitHub replay at commit `7483ae8e20fa5df2439ad524312d3ac9a0fbf71f` against those downloaded inputs;
5. confirm the 16-file reference comparison passes;
6. record the new Zenodo record DOI and public-verification receipt.

A private Drive copy or Actions artifact is durability evidence, not a public-release substitute.
