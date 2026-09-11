# Zenodo update audit

Checked on 2026-09-11 against the live [record API](https://zenodo.org/api/records/22295791), the downloaded ZIP directory and its SHA-256. **A new version is required for the current construct-level analysis. No deposition was changed or published during this audit.**

## Existing public record

- Record DOI: [10.5281/zenodo.22295791](https://doi.org/10.5281/zenodo.22295791).
- Concept DOI: `10.5281/zenodo.22295790`.
- Title: *Azami Chapter 1 v2 reproducibility input package*.
- Last update reported by the API: `2026-09-04T16:04:51.328319+09:00`.
- File: `azami_ch1_v2_reproduction_inputs_2026-09-04.zip`, 56,942,044 bytes.
- SHA-256: `50ec15b1280d4660839ca4bf0d55c970a5f49b4d4feaabb7a073b73500253677` (local downloaded copy verified).
- Published metadata license: CC BY 4.0. This is not the new code MIT license.

The ZIP contains its manifest/readme/checksums/metadata and exactly four embedded input archives: continuous measurements `9612943217`, v2 multilevel output `9632715852`, trees `8227254443`, and region lookup `8983877726`. It does not contain the full current analysis package.

## Required contents of a new version

| Material | Existing v2 ZIP | Action for current release |
|---|---|---|
| Continuous measurements | Present, exact original archive | Retain original bytes and checksum |
| Full nine-predictor process environment, artifact `9633419268` | Absent | Add verified archive; v2 multilevel output is not this input |
| Broad-region lookup and 52 placement trees | Present | Retain original bytes and checksum |
| Native-status table used in sensitivity, not primary filtering | Recoverable from immutable Git history only | Include exact checksum-verified table for offline reproduction |
| Current code and numerical dependencies | Not the current code revision | Include code archive pinned to the final cleaned main commit and `requirements-current.txt` |
| Current aggregate reference outputs | Absent | Include the 15 checksum-verified files indexed in `current_reference/manifest.json` and any further figure source tables needed by the final paper |
| Current figures and provenance | v2 package is not a current figure package | Add final current figure source code, input map and export files after document/figure QA |
| Full replay validation | Historical v2 checks only | Include current input-verification and 15-file numerical-comparison receipts, runtime versions and execution command |
| License metadata | Existing data record says CC BY 4.0 | Identify MIT software separately; retain and verify third-party/data terms |

Current numerical input identities are executable in `run_current_analysis.py`; current reference output identities are recorded in `current_reference/manifest.json`. The existing v2 release and its DOI must remain intact. Use a new version under the existing concept, or linked code/data records if different licensing requires separation; do not silently replace the v2 package.

## Publication check

Before claiming current credential-free reproducibility: publish the approved new version, download it without owner credentials, verify every checksum, unpack into a clean directory, and complete the numerical run and reference comparison. Retain unsupported rows and sensitivity failures. GitHub Actions success or local file presence alone does not close this public-archive requirement.

Original photographs and upstream detector training are outside the minimum numerical replay. If claiming reproduction from original photos rather than from frozen measurements, separately audit image access/licenses, model-weight availability, all upstream inputs and software versions. The present v2 minimum bundle and this cleanup do not establish that stronger claim.
