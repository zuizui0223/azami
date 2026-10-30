# Azami | Global photographic phenomics of *Cirsium*

> **Current working paper:** *Global photographic phenomics reveals climate-associated capitulum traits and local floral resemblance in Cirsium* — a **GEB Research Article candidate, not submitted or published**. [PR #159](https://github.com/zuizui0223/azami/pull/159) was **merged into `main`** on 2026-10-30 (merge commit `e2229a20045044977afe235255355ca3892186fa`). This is the canonical working manuscript, **not** a submitted or published paper.
>
> **現在地:** 最新の論文はResults **3.1 降水と向き、3.2 日射・標高と花色、3.3 似た気候での多形質の類似**。古い「形質組合せ／向き／花色」の3部構成を現在の本文と混同しないでください。

## Read the current research

| Purpose | Canonical entry |
| --- | --- |
| **Latest manuscript and figures** | [Research Article](manuscripts/v4_s1_s2_frozen_20261025/GEB_RESEARCH_ARTICLE.md) · [manuscript guide](manuscripts/README.md) |
| Results, cohort definitions and pending decisions | [GEB preparation checklist](manuscripts/v4_s1_s2_frozen_20261025/GEB_SUBMISSION_CHECKLIST.md) · [claim and submission gates](manuscripts/v4_s1_s2_frozen_20261025/CLAIM_LEDGER_AND_SUBMISSION_GATES.md) |
| Figures 1–5, provenance and licenses | [figure directory](manuscripts/v4_s1_s2_frozen_20261025/figures/) · [asset audit](manuscripts/v4_s1_s2_frozen_20261025/GEB_SUBMISSION_ASSET_AUDIT_20261030_JA.md) |
| Frozen scientific evidence | [V4 analysis](analysis/v4/README.md) · [24 October freeze](analysis/v4/FINAL_V4_BIO18_HOLM_LEGO_FREEZE_20261024_JA.md) |
| Repeat an allowed verification / prepare journal materials | [reproduction and submission runbook](docs/REPRODUCIBILITY_AND_SUBMISSION_RUNBOOK.md) |
| Prepare a V4 Zenodo data-only deposit (not yet published) | [V4 Zenodo handoff](docs/ZENODO_V4_RELEASE_HANDOFF.md) · [CC BY 4.0 data-license notice](manuscripts/v4_s1_s2_frozen_20261025/reproducibility_v4/LICENSE_DATA_CC_BY_4_0.md) · [release-candidate receipt](manuscripts/v4_s1_s2_frozen_20261025/reproducibility_v4/CCBY4_RELEASE_CANDIDATE_RECEIPT_20261030.json) |
| Find old work without treating it as a current result | [repository map](docs/REPOSITORY_MAP.md) · [archive guide](legacy/README.md) |
| Branch/PR cleanup without losing provenance | [branch triage](docs/BRANCH_ARCHIVE_TRIAGE_20261030.md) |

## Three results in the **current** article

1. **3.1 — Precipitation and photographed capitulum orientation:** among 178 native taxa (168,907 angle-eligible records), warmest-quarter precipitation (BIO18) is positively associated with image-frame head angle (standardised β=+0.045447; Holm P=0.00920635 in the original 80-test S1 family). This is *not* a measurement of gravity-based nodding.
2. **3.2 — Floral colour along radiation and elevation:** separate within-range models for native/introduced *C. arvense* and *C. vulgare* yield three introduced-range associations surviving Holm correction in the original S2 128-test family. The four original colour-cohort sizes are **22,936/26,171** (*C. arvense*, native/introduced) and **17,335/27,200** (*C. vulgare*, native/introduced). The analysis did **not** formally test range-by-climate interaction or evolution after introduction.
3. **3.3 — Multivariate resemblance under analogous local climates:** 63 taxa in the primary 1° analysis (7,117 taxon-by-cell values; 153,843 contributing photographic records *before* aggregation) show closer joint orientation–chroma–outline-aspect phenotypes at climate-matched cells than under a within-taxon whole-vector permutation null. Distance advantages: **+0.04096** at 1° and **+0.09433** at 5° (exploratory Holm-2 P=0.010 at each grain). The 1° signal depends on trait/region choices, documented in SI.

The shared **156,683 photographs across 174 taxa** form a separate eight-construct covariance cohort; they are **not** a universal denominator for results 3.1–3.3. S1 and S2 are original source cohorts, not three independent sets of images. The H3c exploratory correction is separate from their 80/128-test families. See the [manuscript guide](manuscripts/README.md) before reusing any sample size.

## What the work contributes — and what it does not claim

**Methodological:** observable, continuous measures of capitulum presentation, colour and two-dimensional outline can remain attached to the same photographic record and its location, enabling analyses beyond a single species mean.

**Ecological:** capitulum appearance has climatic associations among native taxa, within the ranges of two species, and across different taxa at locally analogous climates. Rain avoidance, photoprotection and thermal absorption are **functional hypotheses**, not demonstrated mechanisms. Image top is not gravity vertical; visible JPEG colour is not calibrated pigment chemistry. This study does not demonstrate adaptive convergence or fitness advantages.

## Repository boundaries

- [`manuscripts/`](manuscripts/README.md): the current article, its SI, figure scripts, provenance and submission builders. Generated Word/PDF/ZIP review packages are produced outside Git, not scientific inputs.
- [`analysis/`](analysis/README.md): V4 source contracts, scripts, checked aggregate outputs, and historical Chapter 1 code kept separate.
- [`legacy/`](legacy/README.md): earlier V2/V3 and deferred experiments; preserved for provenance, not current claims.
- [`tests/`](tests/): independent tests and frozen numerical/claim contracts.
- [`.github/workflows/`](.github/workflows/): six current V4 verification or source-recovery workflows. These include protected-data operations; do not run or redistribute their private intermediates casually.

**Science freeze:** the original S1/S2 statistical families were fixed on **2026-10-24**. H3c was added as a separate exploratory analysis on **2026-10-29**. The active manuscript keeps both histories explicit. The older [Japanese three-part working outline](manuscripts/v4_s1_s2_frozen_20261025/GEB_THREE_PART_STRUCTURE_JA.md) is historical planning context, **not** the current Results index.

## Publication, licenses and privacy

**Data licensing choice (2026-10-30):** CC BY 4.0 is approved for original author-controlled, legally redistributable aggregate numerical outputs; original software remains MIT, and third-party rights require separate clearance. The [new CC BY 4.0 source-preserving candidate](manuscripts/v4_s1_s2_frozen_20261025/reproducibility_v4/AGGREGATE_RELEASE_CANDIDATE_README.md) passed nine offline replay checks. No V4 data DOI or final public/anonymous reviewer deposit is established. An Actions artifact can expire; a successful local check against an existing aggregate cannot recreate restricted image-level data. Raw observation identifiers, precise locations, protected source photographs and other third-party material must follow their original access restrictions. Code's MIT license does **not** override image/data licenses: [LICENSE](LICENSE), [NOTICE](NOTICE.md), [Figure 1 credits](manuscripts/v4_s1_s2_frozen_20261025/FIGURE1_PHOTO_CREDITS_20261030.json). The historical Zenodo DOI for V2 is not a V4 dataset DOI.

**Maintainership policy:** edit the dated GEB paper via reviewable PRs onto `main`, keep source-locked scientific numbers stable, and record new verification separately. Do not bulk-move paths or delete branches before checking import paths, GitHub Actions, existing PR bases and archived artifact receipts. See [branch triage](docs/BRANCH_ARCHIVE_TRIAGE_20261030.md).
