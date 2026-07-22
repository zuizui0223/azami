# Chapter 1 manuscript workspace

Use this directory for submission assembly. Do not copy numbers from legacy notes, old figures or workflow logs.

## Canonical submission entry point

1. `SUBMISSION_MANUSCRIPT.md` — current manuscript order and scientific status;
2. `COHORT_FLOW_AND_ANALYSIS_LEDGER.md` — immutable cohort names, counts and analysis permissions;
3. `final_claims.json` — reviewer-revised machine-readable claims;
4. `results/reviewer_precision_summary.json` — legacy precision audit and replacement results;
5. `SUBMISSION_FIGURE_MAP.md` — active and excluded figure definitions;
6. `EXTERNAL_COMPLETION_GATES.md` — independent detector/measurement, taxonomy, spatial and archive gates.

## Current correction

The former environmental-responsiveness index was the RMS of unpooled absolute species slopes. It was strongly confounded with sample size and standard error, so the rho = -0.333 result and all median-split quadrants are withdrawn.

The revised analysis contains 101 taxa with all seven linear endpoints and four climate predictors, gives equal weight to orientation, colour and shape, and uses every archived slope standard error. The primary hierarchical variance meta-regression detects no common coupling between visible variation and latent slope magnitude.

## Drafting rules

- Write Results from `final_claims.json` and the current result ledgers, not from memory.
- Name the exact cohort whenever reporting an FDR count.
- Keep the 6,626-head image atlas, the 406,582-observation exhaustive layer and the 46,276-observation spatially thinned primary cohort separate.
- Treat hue as a joint circular endpoint.
- Use *visible dispersion*, *within-species spatial environment–trait association* and *among-species environmental sorting*.
- Do not use *climate tracking* or *environmental responsiveness* as temporal or experimental claims.
- Do not restore legacy quadrants or raw absolute-slope RMS rankings.
- Mark biological mechanisms as hypotheses.

## Reproducibility

The revised analysis is implemented in `analysis/reanalyze_lability_precision.py` and `.github/workflows/ch1-reviewer-precision-reanalysis.yml`. All headline claims are validated by `analysis/validate_final_claims.py`.

The obsolete supervisor-review entry documents and invalid legacy figures have been removed from the active branch. Historical strategy files may still discuss the former analysis and are superseded wherever they conflict with this README, `SUBMISSION_MANUSCRIPT.md` or `final_claims.json`.
