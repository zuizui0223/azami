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
