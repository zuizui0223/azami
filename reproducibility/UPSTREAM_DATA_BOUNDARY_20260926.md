# Chapter 1 upstream-data and reproducibility boundary — 2026-09-26

This document separates **source provenance**, **image-to-trait production**, and **current manuscript numerical reproduction**. They are related but they are not the same reproducibility claim.

## Five-layer evidence chain

| Layer | Material | Current role | Current reproducibility status |
|---|---|---|---|
| U0 — external source material | Public observation photographs and metadata; CHELSA climate sources; Natural Earth geography; WCVP/POWO taxonomic material; phylogenetic backbones and third-party model dependencies | Provenance and upstream source context | Source identities/terms are documented where used, but the minimum numerical package does not redistribute every original upstream byte |
| U1 — image-to-trait production | Capitulum localization, image QC and deterministic continuous-trait measurement; historical detector training/provenance | Produces the frozen image-derived measurements | Historical implementations and provenance are retained under `legacy/` and reproducibility ledgers, but the current v3 runner does not redownload photographs, retrain YOLO or rerun the full image-to-trait pipeline |
| N0 — frozen numerical inputs | Continuous measurements; nine-predictor environment table; broad-region lookup; 52 historical-placement trees; frozen native-status table | **Starting point of the current manuscript replay** | Exact identities and SHA-256 contracts are fixed and checked by `reproducibility.run_current_analysis` |
| N1 — current analysis code | `analysis/v3/`, `reproducibility/run_current_analysis.py`, current validator, plus the two retained spatial/historical helpers under `legacy/v2/analysis/` | Reproduces the statistical analyses reported by the current manuscript | Current runner executes **8 numerical stages** |
| O0 — current reference outputs | Frozen aggregate tables/reports used to validate the replay | Numerical target surface for the manuscript | `reproducibility/current_reference/manifest.json` contains **16 checksum-locked files**; the 2026-09-15 replay matched all 16 |

## What the current reproduction claim means

The supported claim is:

> The current Chapter 1 statistical results can be replayed from frozen image-derived measurements and fixed numerical auxiliary inputs using the current analysis code, with the replay checked against 16 frozen aggregate reference files.

The current reproduction claim does **not** mean:

- every original photograph is redistributed in the numerical archive;
- every third-party raster, taxonomy source, map or phylogenetic backbone is relicensed by this repository;
- the detector is retrained from scratch during the current replay;
- the current replay reconstructs the continuous measurements from the original photographs;
- image-derived values are independently validated physical traits merely because the numerical replay passes.

Photograph rights and other third-party terms remain separate from the repository MIT software license; see `NOTICE.md`.

## Exact current numerical starting surface

The current runner verifies these numerical materials before fitting any manuscript model:

1. continuous image-derived measurements — Actions artifact `9612943217`;
2. nine-predictor process environment — Actions artifact `9633419268`;
3. broad-region lookup — Actions artifact `8983877726`;
4. historical-placement trees — Actions artifact `8227254443`;
5. frozen native-status table — SHA-256 `c01eeb9ff245d7f73da1a12fa4eede904dd9770467655f20e3d85de2ac8dd84a` after the permitted newline normalization.

These materials are the **upstream inputs for the current statistical analysis**, even though some of them are themselves derived from still-earlier public or third-party sources.

## Current code and replay state

The current code is present in the repository. The canonical entry point is:

```bash
python -m reproducibility.run_current_analysis --download
```

The 2026-09-15 replay receipt records:

- 8 completed numerical stages;
- 16 aggregate files compared;
- validation status `PASS`;
- no new photographs downloaded.

See `reproducibility/current_replay_execution.json` and `reproducibility/current_replay_validation.json`.

## Public-archive status

The public Zenodo record `10.5281/zenodo.22295791` is the historical **v2 minimum numerical package**. It does not contain the complete current v3 construct-level release.

Current v3 materials are checksum-verified and durably staged, and the current eight-stage replay has already passed. A final current public archive still requires the final code/dependency freeze, final submitted figure/provenance manifest, approved release metadata/licensing, publication of the new Zenodo version (or linked records if licensing requires separation), and a credential-free post-publication redownload plus clean replay.

## Double-anonymous review boundary

The anonymous peer-review bundle builder packages the current processed numerical inputs, analysis code and 16-file reference surface without public repository/author/permanent-DOI identifiers, and self-replays the eight-stage analysis before producing the ZIP.

That bundle is a **numerical review package**, not a raw-photograph archive. It still needs to be placed behind an approved stable anonymous reviewer-access route before the blinded manuscript can cite it.

## Language to use in the manuscript

Prefer:

> Numerical reproduction begins from frozen image-derived continuous measurements and fixed environmental, spatial and historical auxiliary inputs. The archived analysis code reproduces the manuscript's statistical analyses and validates them against checksum-locked aggregate reference outputs. Reproduction from original third-party photographs is a separate upstream workflow and is not claimed by the minimum numerical archive.

Avoid:

> All raw source data and the complete image-analysis pipeline are reproduced by the manuscript archive.

unless a separate upstream archive has actually been assembled, licensed, published and independently tested.
