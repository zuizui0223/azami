# Chapter 1 continuous-trait reanalysis v1

## Status

This is a compact, review-facing release of the local PR #69 × continuous-trait integration run. It is **not** a replacement for the frozen 46,276-observation primary family and creates no new submission claim.

- Contract: 27 category-free continuous endpoints across eight modules.
- Executed locally: 17 endpoints; 10 endpoints remain unexecuted and are not imputed.
- Observation grain: 3,725 atlas observations and 216 source-assigned taxa.
- Species-summary grain: 214 taxa after the upstream minimum-observation rule.
- Diagnostic environmental table: 935 observations at positional accuracy ≤10 km.
- Model-complete high-resolution candidate cohort: 904 observations and 165 taxa.
- Generic screen: 48 successful endpoint–predictor models; 0 primary and 4 candidate rows passed tier-wise BH.
- PR #69 gate comparison: all four candidate screen rows have a matching quality-adjusted model; 0 passed the locked adjusted BH rule and 0 are submission-eligible.
- Locked involucre reproduction: 72 coefficient rows reproduced within `1e-9`; 0/3 endpoints passed the full resolution/sharpness retention rule.

## Why the generic signals are not results

The generic registry screen fits one predictor at a time and is only a discovery diagnostic. PR #69 predeclared a stricter model for high-resolution involucre endpoints: all four CHELSA predictors, log minimum head dimension and log1p sharpness, followed by BH across the 12 climate rows and a BIO4 sign-consistency check across fixed resolution strata. Independent continuous botanical reference validation remains pending.

`candidate_screen_vs_pr69_gates.csv` records the failure reason for every generic candidate signal. None may be promoted from the generic q value alone.

## Files

- `continuous_trait_contract_frozen.csv`: endpoint definitions, formula identifiers, units, tiers, validation status and prohibited interpretations.
- `continuous_trait_universe_report.json`: execution and coverage audit.
- `diagnostic_climate_screen_coefficients.csv`: registry-wide diagnostic coefficients, explicitly marked screening-only.
- `diagnostic_climate_screen_coverage.csv`: endpoint–predictor support counts.
- `diagnostic_climate_screen_report.json`: model-family summary and required gates.
- `pr69_bridge_report.json`: key coverage, aliases, source hashes and exact-model reproduction audit.
- `pr69_canonical_involucre_adjusted_models.csv`: locked quality-adjusted models under canonical endpoint names.
- `pr69_canonical_involucre_gate.csv`: endpoint-level PR #69 retention decisions.
- `candidate_screen_vs_pr69_gates.csv`: generic-to-locked gate comparison.
- `candidate_screen_vs_pr69_gates_report.json`: comparison summary.

The full observation-level long table remains a local generated artifact because it is reproducible from the frozen source artifact and would add about 33 MB to Git. The executed notebook is stored at `analysis/notebooks/ch1_pr69_continuous_reanalysis.ipynb`.
