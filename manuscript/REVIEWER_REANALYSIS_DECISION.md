# Reviewer-driven reanalysis decision

## Decision

The original species-level environmental-responsiveness index and all conclusions based on its raw RMS magnitude are withdrawn from the submission manuscript.

## Reason

The legacy score used absolute unpooled species-specific slope estimates. Across 102 taxa it was strongly associated with median slope sample size (Spearman rho = -0.919) and median slope standard error (rho = 0.914). The original visible-variation association (rho = -0.333) became 0.056 after controlling median sample size (P = 0.577). Stability across minimum-n thresholds did not remove this bias.

## Replacement

A revised 101-taxon cohort requires all seven linear endpoints, all four climate predictors and at least 10 observations for every slope. Orientation, colour and shape receive equal weight. The replacement evidence consists of:

1. a transparent sampling-noise-adjusted association-energy estimator using `beta^2 - SE^2`;
2. a hierarchical variance meta-regression over all 2,828 slope estimates and their standard errors.

Neither detects a common coupling between visible variation and latent slope magnitude.

## Submission consequences

- Withdraw rho = -0.333 as a biological result.
- Withdraw all four median-split quadrants and taxon classifications based on them.
- Exclude raw absolute-slope RMS module comparisons from submission figures.
- Exclude minimum-n rank stability of the biased score as supporting evidence.
- Replace *environmental responsiveness* and *climate tracking* with *within-species spatial environment–trait association*.
- Keep the conceptual distinction among visible dispersion, within-species spatial association and among-species environmental sorting.

## Superseded documents

Where they conflict with the revised claim registry, the following historical strategy documents are superseded:

- `FINAL_MANUSCRIPT_STRATEGY.md`
- `RESULT_LED_STORY_REVIEW.md`
- `CIRSIUM_EVOLUTIONARY_ECOLOGY_AND_DISCUSSION_STRATEGY.md`
- `MANUSCRIPT_OUTLINE.md`
- `FINAL_STORY_AND_ANALYSIS.md`
- earlier figure maps and source ledgers referring to the legacy two-axis quadrants.

The authoritative sources are now `SUBMISSION_MANUSCRIPT.md`, `COHORT_FLOW_AND_ANALYSIS_LEDGER.md`, `SUBMISSION_FIGURE_MAP.md`, `final_claims.json` and `results/reviewer_precision_summary.json`.
