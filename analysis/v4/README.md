# Canonical frozen V4 S1–S2 (GEB Research Article)

This directory contains a **curated copy of the pre-established V4 S1/S2 analysis**, not a new fit. The final scientific decision and all four Holm-supported effects were frozen on **24 October 2026** at commit `7e25e409d39a7e8bb8a1ce1043d6bbfd41b7d159`.

- [Immutable scientific conclusions](FINAL_V4_BIO18_HOLM_LEGO_FREEZE_20261024_JA.md)
- [Numerical status and endpoint model comparisons](V4_S1_S2_FOUR_MODEL_CORRECTED_UNCORRECTED_COMPARISON_20261024_JA.md)
- [Primary and exploratory environment contract](V4_S1_S2_EXTREME_LEGO_CONTRACT_20261024.md)
- [Original frozen artifact auditor](final_v4_bio18_holm_freeze.py) (stdlib only)
- [GEB manuscript](../../manuscripts/v4_s1_s2_frozen_20261025/GEB_RESEARCH_ARTICLE.md) and [claim ledger](../../manuscripts/v4_s1_s2_frozen_20261025/CLAIM_LEDGER_AND_SUBMISSION_GATES.md)

Frozen primary exposures: CHELSA BIO18, BIO05, RSDS and WorldClim elevation. Source-run Actions [40062126006](https://github.com/zuizui0223/azami/actions/runs/40062126006), aggregate artifact 12547876704, source ZIP SHA-256 `6e081dc081bb8ee2311a05ac5b7fcddd3f6309be2fa909e87d432789f4e0e205`. Check via `python -m analysis.v4.final_v4_bio18_holm_freeze --artifact-zip <verified.zip> --output <aggregate.json>`.

All original private source-reconstruction dependencies and cached assets remain reproducible at the **pinned source execution commit** `d2c45b00018a0c04e7c69a29857207dc1ac06d0c`. The scoped V4 code here reproduces the source calculation logic subject to the original private inputs; do not claim that `main` contains protected photos, decrypted metadata or new validated measurements. T2/T3 live outside this canonical S1/S2 manuscript scope.

V3 remains at [`legacy/v3/`](../../legacy/v3/README.md), numerically separate.

## Verified elevation-free analyses, separated from the October frozen model

The **historical frozen S1/S2 model** on 2026-10-24 includes elevation. The user's **new analysis direction** (31 October) excludes elevation from current climate-only alternatives; this does not relabel earlier four-predictor coefficients as three-predictor results.

- **S1/S2 three-predictor ablation:** [PR #161](https://github.com/zuizui0223/azami/pull/161) compares BIO18/BIO05/RSDS against the frozen four-predictor baseline on identical data. The source-verified [31 October result](https://github.com/zuizui0223/azami/blob/analysis/v4-elevation-ablation-20261031/analysis/v4/ELEVATION_ABLATION_SOURCE_VERIFIED_RESULT_20261031_JA.md) preserves the original coefficients. Full eight-construct Holm-adjusted re-estimation is **not yet complete**; do not claim final climate-only manuscript effects.
- **H3c three- and five-climate-axis comparison:** [PR #162](https://github.com/zuizui0223/azami/pull/162) tests BIO12/BIO01/RSDS; BIO18/BIO05/RSDS; BIO18/BIO01/RSDS; and BIO18/BIO01/RSDS/BIO04/BIO15, all **without elevation**. [Source-verified biological interpretation](results/20261031_H3C_FOUR_ELEVATION_FREE_CLIMATE_AXES_SOURCE_VERIFIED_JA.md) and [machine-readable aggregate](results/20261031_H3C_FOUR_CLIMATE_AXES_PAIRED_SUMMARY.json) document both the positive five-axis resemblance and the **lack of a supported increment from adding seasonality**. This is a post-result exploratory comparison; it does not replace the original H3c Holm-2 endpoint.

## GEB scope boundary (25 October 2026)

For the present manuscript, S2 compares **photographic phenotype–environment coefficients and visible-trait covariance within four species×range strata**; it does **not** compare the distributions of BIO18/climatic niches between native and introduced ranges. The subsequent source-real niche results have been deferred to [`legacy/v4_s2_precipitation_niche_20261025/`](../../legacy/v4_s2_precipitation_niche_20261025/README.md). **Frozen 24 October models, coefficients, Holm families and original artifacts remain unchanged.**
