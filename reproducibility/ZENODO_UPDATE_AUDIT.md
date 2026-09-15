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

## Durable staging updates

The public Zenodo record has **not** been changed. Current-only artifacts that would otherwise depend on expiring GitHub Actions storage have been recovered by exact artifact identity, re-verified against SHA-256 values and copied into the existing non-public durable Drive archive. The machine-readable receipt is [`CURRENT_RELEASE_STAGING_20260912.json`](CURRENT_RELEASE_STAGING_20260912.json); it has been extended as new submission-readiness analyses closed.

Durable copies now include:

| Role | Actions artifact | Verified ZIP SHA-256 | Durable status |
|---|---:|---|---|
| Nine-predictor process environment | `9633419268` | `d7c0c466f55b67695d06ae46c21a6452dbe6cfd92a52db8042caa200429e97f4` | copied to durable Drive archive |
| Biological axes + sensitivity outputs | `10130210432` | `7568ec98b0709f6cda027b928298d58c7e6aa42af0f6c0947e1e91c3d1abbab2` | copied to durable Drive archive |
| Construct-scale integration outputs | `10131007603` | `a483e116c46211023df95bef9388aa89deaef3441d85701fc93792867884514a` | copied to durable Drive archive |
| Complete-construct upgrade + direct scale contrast | `10136229131` | `7b4c5703b048ff108cd794921394ae22d1423461c9aaf330b47029c5d1f4908f` | copied to durable Drive archive |
| Assessability + technical-stress outputs | `10135679053` | `fa80f6079d676186d48665a574b5e68f906037a466a9e3101a57d87ee02f5ba7` | copied to durable Drive archive |
| WCVP accepted-name sensitivity | `10292140117` | `dfb6eec3001e3a984662d5aba06cda5fa80e144b36ccb4af9fdf45973854edc5` | copied to durable Drive archive; contains exact native-status input |
| RV estimator-validity sensitivity | `10382387052` | `345866a7e333f78677ad3797811e2cbe82d3e09ee061f7642df7f4c4d5ec008e` | copied to durable Drive archive |
| Pre-estimator scale-integration figure | `10291656193` | `2fed9448c2210af4ded7a4ccc5cbe6f543b64e8f19ba870f5200fa095290e766` | durable historical figure artifact |
| Estimator-validity scale-integration figure | `10382578412` | `d538c428f3acb3eab6a203c0a0a435ce0ae51e0da9b4f23b7d42e92fcc0e3af7` | current manuscript figure candidate; visual QA passed |

The current reference manifest now contains **16 checksum-verified files**. The six pre-estimator current-output archives plus the estimator-validity artifact supply those reference files; continuous measurements `9612943217`, broad-region lookup `8983877726`, and the 52-tree historical input `8227254443` were already recorded in `durable_archive_manifest.json`.

This closes an **ephemeral-storage risk**, not the public-release gate. A private durable copy is not a Zenodo publication and must not be described as credential-free reproducibility.

## Required contents of a new version

| Material | Current staging status | Action for public current release |
|---|---|---|
| Continuous measurements | exact original artifact already durable and present in v2 release | Retain original bytes and checksum |
| Full nine-predictor process environment, artifact `9633419268` | checksum-verified and durably staged | Include exact staged archive; v2 multilevel output is not this input |
| Broad-region lookup and 52 placement trees | exact artifacts already durable and present in v2 release | Retain original bytes and checksums |
| Native-status table used in sensitivity, not primary filtering | exact input is preserved inside the durable WCVP sensitivity artifact and has frozen SHA identity | Extract the exact member into the final release input set and reverify the permitted newline-normalized SHA |
| Current code and numerical dependencies | repository code exists; no final release archive has yet been frozen | Include code archive pinned to the final cleaned main commit and `requirements-current.txt` or equivalent pinned environment |
| Current aggregate reference outputs | 16-file manifest is checksum-verified; estimator-validity summary is included and its source artifact is durable | Include all 16 files indexed in `current_reference/manifest.json` |
| Current figures and provenance | current Figure 3 estimator-validity candidate is checksum-verified, durably staged and visually QA'd; the complete final manuscript figure package is not yet frozen | Add final current figure source code, input map and all submitted exports after document/figure QA |
| Full replay validation | earlier seven-stage / 15-file receipts exist, but the current runner now has eight stages and a 16-file reference surface | Produce a new eight-stage replay and 16-file comparison receipt from the final code surface |
| License metadata | existing data record says CC BY 4.0; repository code is MIT | Identify MIT software separately; retain and verify third-party/data terms |

Current numerical input identities are executable in `run_current_analysis.py`; current reference output identities are recorded in `current_reference/manifest.json`. The existing v2 release and its DOI must remain intact. Use a new version under the existing concept, or linked code/data records if different licensing requires separation; do not silently replace the v2 package.

## Remaining release gate

The material-recovery problem is now narrow. The remaining work is to freeze/package the final code/dependency revision, produce the updated eight-stage replay/16-file comparison receipts, freeze the complete submitted figure/provenance surface, resolve release metadata/licensing, publish the approved new version, then perform the anonymous redownload test.

Before claiming current credential-free reproducibility: publish the approved new version, download it without owner credentials, verify every checksum, unpack into a clean directory, and complete the numerical run and reference comparison. Retain unsupported rows and sensitivity failures. GitHub Actions success, a private Drive copy, or local file presence alone does not close this public-archive requirement.

Original photographs and upstream detector training are outside the minimum numerical replay. If claiming reproduction from original photos rather than from frozen measurements, separately audit image access/licenses, model-weight availability, all upstream inputs and software versions. The present v2 minimum bundle and this staging update do not establish that stronger claim.
