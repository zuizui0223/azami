# Zenodo update audit

Checked on 2026-09-11 against the live [record API](https://zenodo.org/api/records/22295791), the downloaded ZIP directory and its SHA-256. **A new version is required for the current construct-level analysis. The existing v2 deposition remains unchanged.**

## Existing public record

- Record DOI: [10.5281/zenodo.22295791](https://doi.org/10.5281/zenodo.22295791).
- Concept DOI: `10.5281/zenodo.22295790`.
- Title: *Azami Chapter 1 v2 reproducibility input package*.
- Last update reported by the API at the 2026-09-11 audit: `2026-09-04T16:04:51.328319+09:00`.
- File: `azami_ch1_v2_reproduction_inputs_2026-09-04.zip`, 56,942,044 bytes.
- SHA-256: `50ec15b1280d4660839ca4bf0d55c970a5f49b4d4feaabb7a073b73500253677` (local downloaded copy verified).
- Published metadata license: CC BY 4.0. This is not the new code MIT license.

The ZIP contains its manifest/readme/checksums/metadata and exactly four embedded input archives: continuous measurements `9612943217`, v2 multilevel output `9632715852`, trees `8227254443`, and region lookup `8983877726`. It does not contain the full current analysis package.

## 2026-09-12 durable staging update

The public Zenodo record has **not** been changed. However, the principal current-only artifacts that were still dependent on expiring GitHub Actions storage have now been recovered by exact artifact identity, re-verified against the frozen SHA-256 values, and copied into the existing non-public durable Drive archive.

The machine-readable receipt is [`CURRENT_RELEASE_STAGING_20260912.json`](CURRENT_RELEASE_STAGING_20260912.json).

New durable copies:

| Role | Actions artifact | Verified ZIP SHA-256 | Durable status |
|---|---:|---|---|
| Nine-predictor process environment | `9633419268` | `d7c0c466f55b67695d06ae46c21a6452dbe6cfd92a52db8042caa200429e97f4` | copied to durable Drive archive |
| Biological axes + sensitivity outputs | `10130210432` | `7568ec98b0709f6cda027b928298d58c7e6aa42af0f6c0947e1e91c3d1abbab2` | copied to durable Drive archive |
| Construct-scale integration outputs | `10131007603` | `a483e116c46211023df95bef9388aa89deaef3441d85701fc93792867884514a` | copied to durable Drive archive |
| Complete-construct upgrade + direct scale contrast | `10136229131` | `7b4c5703b048ff108cd794921394ae22d1423461c9aaf330b47029c5d1f4908f` | copied to durable Drive archive |
| Assessability + technical-stress outputs | `10135679053` | `fa80f6079d676186d48665a574b5e68f906037a466a9e3101a57d87ee02f5ba7` | copied to durable Drive archive |

Together, the four current-output archives above are the source archives for the 15 checksum-verified current reference files indexed in `current_reference/manifest.json`. Continuous measurements `9612943217`, broad-region lookup `8983877726`, and the 52-tree historical input `8227254443` were already recorded in `durable_archive_manifest.json`.

This closes an **ephemeral-storage risk**, not the public-release gate. A private durable copy is not a Zenodo publication and must not be described as credential-free reproducibility.

## Required contents of a new version

| Material | Current staging status | Action for public current release |
|---|---|---|
| Continuous measurements | exact original artifact already durable and present in v2 release | Retain original bytes and checksum |
| Full nine-predictor process environment, artifact `9633419268` | checksum-verified and durably staged on 2026-09-12 | Include exact staged archive; v2 multilevel output is not this input |
| Broad-region lookup and 52 placement trees | exact artifacts already durable and present in v2 release | Retain original bytes and checksums |
| Native-status table used in sensitivity, not primary filtering | recoverable from immutable Git tag; not yet a dedicated staged release file | Include exact checksum-verified table for offline reproduction |
| Current code and numerical dependencies | repository code exists; no final release archive has yet been frozen | Include code archive pinned to the final cleaned main commit and `requirements-current.txt` or equivalent pinned environment |
| Current aggregate reference outputs | all four source ZIPs for the 15-file manifest are checksum-verified and durably staged | Include the 15 checksum-verified files indexed in `current_reference/manifest.json` |
| Current figures and provenance | current figure recipes exist, but final manuscript figure package is not yet frozen | Add final current figure source code, input map and export files after document/figure QA |
| Full replay validation | local seven-stage replay and 15-file comparison receipts exist | Include current input-verification and comparison receipts, runtime versions and execution command |
| License metadata | existing data record says CC BY 4.0; repository code is MIT | Identify MIT software separately; retain and verify third-party/data terms |

Current numerical input identities are executable in `run_current_analysis.py`; current reference output identities are recorded in `current_reference/manifest.json`. The existing v2 release and its DOI must remain intact. Use a new version under the existing concept, or linked code/data records if different licensing requires separation; do not silently replace the v2 package.

## Remaining release gate

The material-recovery problem is now narrower. The remaining work is to freeze/package the exact native-status input, final code/dependency revision, replay receipts, and final figure/provenance surface; resolve release metadata/licensing; publish the approved new version; then perform the anonymous redownload test.

Before claiming current credential-free reproducibility: publish the approved new version, download it without owner credentials, verify every checksum, unpack into a clean directory, and complete the numerical run and reference comparison. Retain unsupported rows and sensitivity failures. GitHub Actions success, a private Drive copy, or local file presence alone does not close this public-archive requirement.

Original photographs and upstream detector training are outside the minimum numerical replay. If claiming reproduction from original photos rather than from frozen measurements, separately audit image access/licenses, model-weight availability, all upstream inputs and software versions. The present v2 minimum bundle and this staging update do not establish that stronger claim.
