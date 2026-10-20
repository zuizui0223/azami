# V4 native integration implementation plan

Goal: retain V3's question about within- versus between-species capitulum
integration, using the full native source population and conditional environmental
relationships. User approved this direction on 2026-10-20. PR143 stays closed.
This supersedes the earlier "no new V4 execution" handoff restriction for this
branch; it does not change frozen V3 results.

## Scientific contract

- Native source: 313,350 records from 665,115 source observations, not 313,350
  complete trait measurements. Preserve endpoint eligibility and missingness.
- Main modules: presentation angle; whole floral colour (L*, C*, joint circular
  hue); four whole-head outline measures. Bract contour is SI.
- Questions: how do trait associations differ within and between species; how
  much of those associations is shared with environmental gradients; do rainfall
  and presentation, and colour and radiation, correspond under simultaneous
  adjustment? Photographic phenotypes do not directly establish adaptation.
- No blanket thinning. Native membership and spatial dependence are separate.
- Use corrected source decoding, not the old finite uint32 sentinel values.
- Reuse measurements; no new acquisition, hypervolumes, molecular analyses or
  remeasurement project in this workflow.
- Historical native-model results have already been seen. This is a revised,
  explicitly documented analysis, not an independent preregistered confirmation.

## Gates and implementation

1. **Input audit (current executable stage).** Restore the authenticated encrypted
   native table from run 39408603488. Reuse the checksum-locked raster decoder
   from commit a7b870d. `analysis/v4/native_integration_audit.py` reports all
   endpoint counts, pairwise/common-module overlap, availability of sampling
   fields, and environment correlations at raw, species-mean and within levels.
   Upload aggregate results only. Tests check overlap eligibility and that
   species means use the population, not a phenotype-complete subset.
2. **Environment contract.** Before new trait fitting, assign candidate roles,
   inspect the audit and freeze selection/tie-breaking rules. Retain BIO12 as the
   rainfall hypothesis exposure and radiation for colour. Evaluate direct
   climate candidates for temperature, water, seasonality and radiation; GSP,
   NPP and wind require an explicit additional rationale rather than automatic
   inclusion. Correlation/VIF diagnose redundancy, not causal confounding.
   Do not choose exclusions by trait significance. Check decomposed designs,
   not just raw environments. If an anchor remains poorly separable, report that
   rather than repeatedly altering models.
3. **Common model engine.** Reuse and verify observation-level random-intercept
   fitting with within/between environmental decomposition. Responses stay at
   observation level; scaling parameters are stored and effects back-transformed.
   Use each endpoint's eligible rows. Freeze inferential families before fitting.
4. **Integration.** Estimate raw and environment-adjusted within/between
   covariance on explicitly reported overlap samples. Do not conflate covariance
   of species intercepts, species medians and observation residuals. Select the
   covariance estimator and its uncertainty method before execution; account for
   species precision, unequal sample sizes and circular hue. Pair-specific
   analyses retain data, with a matched core-module cohort for comparisons across
   pairs. Never treat changing cohorts as changes in biological integration.
   Existing V3 RV estimator-validity results constrain interpretation.
5. **Space and sampling.** Audit available coordinates/observer identities first.
   Freeze spatial block definition and resampling/uncertainty method before
   fitting; assess residual dependence. Compare all-data slopes with a fixed
   spatially balanced sensitivity sample. Do not select spatial settings by
   significance, or infer independence from a non-significant Moran test.
6. **Report.** One evidence chain: source/QC, integration at two scales,
   environmental correspondence, spatial robustness, limits. Show gradients and
   intervals in meaningful units; preserve unsupported hypotheses and separate
   image association from mechanism. Archive input/code/run/output identities.

## Execution and acceptance

One manual GitHub workflow, `v4-native-integration.yml`, initially exposes only
the input audit. Later stages extend that entry point after the preceding gates,
not parallel competing workflows. No numerical fitting/testing runs locally.
Tests run in Actions. Audit success alone is not scientific completion.
Full completion requires validated model fits, matched integration comparisons,
spatial evaluation, rendered figures and source-backed interpretation.
