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
| Current eight-stage replay + 16-file validation | `10382954095` | `2064b63366681aa3bd708fef00d2e7ac50d20a7fad771849b793c0d2b6957a0e` | copied to durable Drive archive; validation PASS |
| Pre-estimator scale-integration figure | `10291656193` | `2fed9448c2210af4ded7a4ccc5cbe6f543b64e8f19ba870f5200fa095290e766` | durable historical figure artifact |
| Estimator-validity scale-integration figure | `10382578412` | `d538c428f3acb3eab6a203c0a0a435ce0ae51e0da9b4f23b7d42e92fcc0e3af7` | current manuscript figure candidate; visual QA passed |

The current reference manifest now contains **16 checksum-verified files**. The six pre-estimator current-output archives plus the estimator-validity artifact supply those reference files; continuous measurements `9612943217`, broad-region lookup `8983877726`, and the 52-tree historical input `8227254443` were already recorded in `durable_archive_manifest.json`.

This closes an **ephemeral-storage risk**, not the public-release gate. A private durable copy is not a Zenodo publication and must not be described as credential-free reproducibility.

## Assembled current v3 Zenodo staging bundle — 2026-09-26

The current material surface has now been assembled end-to-end from current main commit `5c0ae2af9a92bb79295fd11d77ad5bbc25f39757`.

- workflow run: `36239203959`;
- Actions artifact: `10905302325`;
- Actions artifact SHA-256: `8b96fe69f7ae539d4118a29550684c7c1f473959a83659b0532ef67aee4870be`;
- inner candidate archive: `azami_ch1_v3_zenodo_staging.zip`;
- inner archive SHA-256: `b9866cdd7a5a8980d93d60c39ddb7ac65d0143561e8802e556d64fc264397c5b`;
- durable owner copy: Drive file `1texEEfD1aT0NvUKgrhNQuxYHyHqNN3P-`.

The staging archive contains the four frozen numerical archives, checksum-recovered native status, the current 16-file reference surface, current eight-stage replay receipts, a `git archive` snapshot of the packaged checkout, and the checksum-verified current GEB figure surface.

The builder reports exactly two unresolved release gaps:

1. `release_metadata`;
2. `figure_document_qa`.

Therefore the material-assembly problem is no longer open. The archive is intentionally labelled **staging**, has `release_ready=false`, and must not be uploaded as the final Zenodo version without closing those two gates and rebuilding at the exact final Git head.

## Required contents of a new version

| Material | Current staging status | Action for public current release |
|---|---|---|
| Continuous measurements | exact original artifact already durable and present in v2 release | Retain original bytes and checksum |
| Full nine-predictor process environment, artifact `9633419268` | checksum-verified and durably staged | Include exact staged archive; v2 multilevel output is not this input |
| Broad-region lookup and 52 placement trees | exact artifacts already durable and present in v2 release | Retain original bytes and checksums |
| Native-status table used in sensitivity, not primary filtering | exact input is preserved inside durable WCVP sensitivity artifact `10292140117`; release builder now knows its archive SHA, member path and normalized native-status SHA | **Packaging route closed:** final builder auto-extracts and verifies the exact member when the taxonomy artifact is supplied; the resulting CSV is written into the release `inputs/` surface |
| Current code and numerical dependencies | repository code exists; direct numerical requirements now match the completed replay environment, including Biopython 1.88; no final release archive has yet been frozen | Include a `git archive` snapshot pinned to the final cleaned main commit together with the exact `requirements-current.txt` |
| Current aggregate reference outputs | 16-file manifest is checksum-verified; estimator-validity summary is included and its source artifact is durable | Include all 16 files indexed in `current_reference/manifest.json` |
| Current figures and provenance | **render surface complete and durably staged**: artifact `10903835882` contains 5 Main figures, Figures S1.1–S1.7 and S2.1 as PNG/PDF plus checksum/provenance files; visual QA passed outside Word | Synchronize the Word manuscript to this surface, complete caption/numbering/pagination QA, then promote the same hashes in a QA-finalized figure manifest |
| Full replay validation | **complete**: the 2026-09-15 eight-stage replay matched all 16 current aggregate reference files; artifact `10382954095` is durably staged | Include `current_replay_execution.json`, `current_replay_validation.json` and the frozen replay environment in the final release; rerun after publication as the credential-free verification |
| License metadata | existing data record says CC BY 4.0; repository code is MIT | Identify MIT software separately; retain and verify third-party/data terms |

Current numerical input identities are executable in `run_current_analysis.py`; current reference output identities are recorded in `current_reference/manifest.json`. The existing v2 release and its DOI must remain intact. Use a new version under the existing concept, or linked code/data records if different licensing requires separation; do not silently replace the v2 package.

## Remaining release gate

The material-recovery and archive-assembly problems are now closed for a staging build. The native-status recovery route is checksum-gated, the direct numerical dependencies match the completed replay, the eight-stage replay/16-file comparison has passed, the 5-Main + 8-Supporting figure surface is packaged, and a full current-v3 Zenodo staging ZIP has been assembled and durably copied. The remaining work is document-level figure/pagination QA, author-owned release metadata/licensing, one final rebuild at the exact final Git head with `release_ready=true`, publication of the new Zenodo version, and the credential-free redownload/clean replay.

Before claiming current credential-free reproducibility: publish the approved new version, download it without owner credentials, verify every checksum, unpack into a clean directory, and complete the numerical run and reference comparison. Retain unsupported rows and sensitivity failures. GitHub Actions success, a private Drive copy, or local file presence alone does not close this public-archive requirement.

Original photographs and upstream detector training are outside the minimum numerical replay. If claiming reproduction from original photos rather than from frozen measurements, separately audit image access/licenses, model-weight availability, all upstream inputs and software versions. The present v2 minimum bundle and this staging update do not establish that stronger claim. See [UPSTREAM_DATA_BOUNDARY_20260926.md](UPSTREAM_DATA_BOUNDARY_20260926.md) for the frozen distinction between source provenance, image-to-trait production and current numerical replay.
