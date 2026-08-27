# Chapter 1 manuscript workspace

Use this directory for submission assembly. Do not copy numbers from legacy notes, old figures or workflow logs.

## Canonical submission entry point

1. `SUBMISSION_MANUSCRIPT.md` — current manuscript order and scientific status;
2. `AZAMI_CH1_SPATIAL_SCALE_CANONICAL_2026-08-27.md` — trait x gradient x spatial-scale story and EAzami temporal handoff;
3. `results/GEB_V2_FULL27_FULL_ENVIRONMENT_RESULTS_2026-08-27.md` — canonical all-27/all-nine results, provenance and claim ceiling;
4. `COHORT_FLOW_AND_ANALYSIS_LEDGER.md` — immutable cohort names, counts and analysis permissions;
5. `final_claims.json` — reviewer-revised machine-readable claims;
6. `BIAS_CONTROL_REANALYSIS_PROTOCOL.md` — locked raw-calendar and hemisphere-aware seasonal, dominant-taxon, native-range, image-quality and image-remeasurement rules;
7. `results/reviewer_precision_summary.json` — legacy precision audit and replacement results;
8. `SUBMISSION_FIGURE_MAP.md` — active and excluded figure definitions;
9. `DETECTOR_INDEPENDENT_AUDIT_PROTOCOL.md` — leakage-free detector validation design and stop rule;
10. `EXTERNAL_COMPLETION_GATES.md` — ordered measurement, taxonomy, spatial, niche and archive gates.

**Current status: submission hold for external validation.** The full 27-endpoint / nine-predictor atlas, retrospective sampling-composition audit and sequential spatial/historical gates are complete and locally validated. Sixteen of 26 selected within-taxon and all ten selected among-taxon directions are stable throughout the sampling audit. Five within-taxon rows pass the broad-space screen, but only four are also sampling-composition stable; chroma-radiation and orientation-annual-precipitation pass all sampling scenarios and all 52 placement trees and remain adaptive-pattern candidates under current controls. Earlier A/B/C rows remain provenance rather than entry rules for this lane. Image remeasurement and independent validation gates remain open.

The frozen Azami Chapter 1 v2 method/ecology defence, including the bounded anthocyanin and rain-shielding hypotheses, is `AZAMI_CH1_V2_DEFENSE_2026-08-27.md`.

## Current statistical correction

The former environmental-responsiveness index was the RMS of unpooled absolute species slopes. It was strongly confounded with sample size and standard error, so the rho = -0.333 result and all median-split quadrants are withdrawn.

The revised analysis contains 101 taxa with all seven linear endpoints and four climate predictors, gives equal weight to orientation, colour and shape, and uses every archived slope standard error. The primary hierarchical variance meta-regression detects no common coupling between visible variation and latent slope magnitude.

## Current detector gate

The production YOLO11n detector was trained on open-vocabulary pseudo-labels. Its earlier 1,000-image development packet is therefore not independent validation. The new workflow excludes every prior proposal/training photo and observation, generates a new blinded 1,000-image packet, assigns 25% for double annotation and keeps frozen predictions hidden.

Detector precision, recall and F1 remain unreported until the human boxes are complete and adjudicated. The production threshold is fixed at 0.25; the audit uses IoU 0.50 and retains predictions down to 0.01 only for descriptive threshold diagnostics.

## Drafting rules

- Write Results from `final_claims.json` and the current result ledgers, not from memory.
- Name the exact cohort whenever reporting an FDR count.
- Start the canonical endpoint atlas from all 27 registered endpoints and nine predictors; never use the old nine-endpoint or A/B/C results to select rows.
- Keep the 6,626-head image atlas, the 406,582-observation exhaustive layer and the 46,276-observation spatially thinned primary cohort separate.
- Treat hue as a joint circular endpoint.
- Use *visible dispersion*, *within-species spatial environment–trait association* and *among-species environmental sorting*.
- Do not use *climate tracking* or *environmental responsiveness* as temporal or experimental claims.
- Do not restore legacy quadrants or raw absolute-slope RMS rankings.
- Do not report detector accuracy from pseudo-label training metrics or production detections.
- Mark biological mechanisms as hypotheses.
- Lead with continuous trait x environmental gradient x spatial scale; use whole-capitulum organization only as a secondary synthesis.
- Keep EAzami state thresholds and transition models predeclared; do not discretize Azami endpoints post hoc.
- Do not restore the withdrawn involucral BIO4 headline.
- Do not restore chroma–BIO12 to the biological headline: it failed native-only BH correction; the paired colour negative control remains required for broader colour-specific interpretation.
- Do not call day-of-year adjustment developmental-stage control.

## Reproducibility

The revised lability analysis is implemented in `analysis/reanalyze_lability_precision.py`. The detector gate is implemented by `analysis/build_independent_detector_audit_packet.py`, `analysis/detector_audit_annotation_app.py`, `analysis/finalize_detector_audit_annotations.py` and `analysis/evaluate_independent_detector_audit.py`. Executable stages are registered in `analysis/ch1/pipeline.json`.

The obsolete supervisor-review entry documents and invalid legacy figures have been removed from the active branch. Historical strategy files are superseded wherever they conflict with this README, `SUBMISSION_MANUSCRIPT.md` or `final_claims.json`.
