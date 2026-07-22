# Chapter 1 manuscript workspace

Use this directory for submission assembly. Do not copy numbers from legacy notes, old figures or workflow logs.

## Canonical submission entry point

1. `SUBMISSION_MANUSCRIPT.md` — current manuscript order and scientific status;
2. `COHORT_FLOW_AND_ANALYSIS_LEDGER.md` — immutable cohort names, counts and analysis permissions;
3. `final_claims.json` — reviewer-revised machine-readable claims;
4. `results/reviewer_precision_summary.json` — legacy precision audit and replacement results;
5. `SUBMISSION_FIGURE_MAP.md` — active and excluded figure definitions;
6. `DETECTOR_INDEPENDENT_AUDIT_PROTOCOL.md` — leakage-free detector validation design and stop rule;
7. `EXTERNAL_COMPLETION_GATES.md` — ordered measurement, taxonomy, spatial, niche and archive gates.

## Current statistical correction

The former environmental-responsiveness index was the RMS of unpooled absolute species slopes. It was strongly confounded with sample size and standard error, so the rho = -0.333 result and all median-split quadrants are withdrawn.

The revised analysis contains 101 taxa with all seven linear endpoints and four climate predictors, gives equal weight to orientation, colour and shape, and uses every archived slope standard error. The primary hierarchical variance meta-regression detects no common coupling between visible variation and latent slope magnitude.

## Current detector gate

The production YOLO11n detector was trained on open-vocabulary pseudo-labels. Its earlier 1,000-image development packet is therefore not independent validation. The new workflow excludes every prior proposal/training photo and observation, generates a new blinded 1,000-image packet, assigns 25% for double annotation and keeps frozen predictions hidden.

Detector precision, recall and F1 remain unreported until the human boxes are complete and adjudicated. The production threshold is fixed at 0.25; the audit uses IoU 0.50 and retains predictions down to 0.01 only for descriptive threshold diagnostics.

## Drafting rules

- Write Results from `final_claims.json` and the current result ledgers, not from memory.
- Name the exact cohort whenever reporting an FDR count.
- Keep the 6,626-head image atlas, the 406,582-observation exhaustive layer and the 46,276-observation spatially thinned primary cohort separate.
- Treat hue as a joint circular endpoint.
- Use *visible dispersion*, *within-species spatial environment–trait association* and *among-species environmental sorting*.
- Do not use *climate tracking* or *environmental responsiveness* as temporal or experimental claims.
- Do not restore legacy quadrants or raw absolute-slope RMS rankings.
- Do not report detector accuracy from pseudo-label training metrics or production detections.
- Mark biological mechanisms as hypotheses.

## Reproducibility

The revised lability analysis is implemented in `analysis/reanalyze_lability_precision.py`. The detector gate is implemented by `analysis/build_independent_detector_audit_packet.py`, `analysis/detector_audit_annotation_app.py`, `analysis/finalize_detector_audit_annotations.py` and `analysis/evaluate_independent_detector_audit.py`. Executable stages are registered in `analysis/ch1/pipeline.json`.

The obsolete supervisor-review entry documents and invalid legacy figures have been removed from the active branch. Historical strategy files are superseded wherever they conflict with this README, `SUBMISSION_MANUSCRIPT.md` or `final_claims.json`.
