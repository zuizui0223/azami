# Consolidated V4 environment contract

Status: selection rule fixed after the input-only audit, before new trait fits.
Earlier exploratory model results had already been inspected; this is a revised analysis, not an
independent preregistration. All selected exposures enter every main module.

## Candidate roles

Candidate set: BIO1 (thermal level), BIO4 (thermal seasonality), BIO12 (annual
rainfall), BIO15 (rainfall seasonality), RSDS (radiation), VPD (atmospheric demand).
BIO12 and RSDS correspond to the existing V3 anchors. BIO1 remains important for
interpreting flower-colour climate correspondence; seasonality cannot substitute
for mean temperature.

GSP is an alternative precipitation timescale, not an additional required
adjustment for annual rainfall. NPP is climate-derived potential productivity,
not an independent measured productivity exposure. Neither belongs in this
direct-climate main candidate set. Wind is a possible separate mechanical
question, but is not required for the present rainfall/colour/integration
questions. No automatic wind addition or separate wind-hypothesis expansion.
These exclusions are based on the estimand, not trait significance.

## Deterministic diagnostic rule

Use the full native table without any trait-outcome or trait-eligibility filter.
Use the same six-candidate complete rows for candidate comparisons. Define
species means using all available native environments before that complete-row
filter. Check raw environments, equal-species means and within-species deviations.

1. Protect BIO12, RSDS and BIO1 as the predefined core exposures.
2. For absolute pairwise correlation above 0.8 at any of the three levels, remove
   an unprotected member by a fixed priority (lowest priority removed first:
   VPD, BIO4, BIO15). Repeat. A protected/protected conflict is reported as
   unresolved; do not silently delete the hypothesis exposure.
3. Compute VIF for the remaining set at each level. Threshold 5 is a diagnostic
   decision rule, not a proof of identifiability or absence of confounding.
   If any VIF exceeds 5, remove the unprotected variable with the largest maximum
   VIF across levels, resolving ties with the same priority, and repeat.
   If only protected variables remain with excessive VIF, stop selection for
   scientific review and report the poorly separable gradients.
4. Record every exclusion, all matrices, complete-row counts and the final set.
   This is a declared multi-level selection algorithm, not a claim to have run
   the R `usdm::vifstep` function. It extends the requested VIF approach to the
   actual within/between question.
5. Check the resulting endpoint-specific decomposed model designs after QC.
   Do not reselect by endpoint or by outcome P value. Report any remaining
   weak support, high VIF or rank deficiency explicitly.

VPD is a valid environmental descriptor; removing it for redundancy does not
show that atmospheric drying is biologically irrelevant. Likewise, retaining
temperature does not identify a temperature mechanism independently of VPD.

## Scaling and interpretation

Store fixed native-population means and standard deviations for retained
environments. Use the same denominator for each variable's between and within
components, and return effects to physical increments for display. Do not
standardize responses separately within species and erase between-species
variation. Colour hue remains a circular pair, not two independent functions.

Adjusted gradients are conditional observational associations. VIF cannot decide
which variable is a confounder or mediator. Conditioning on climate does not
establish protection from rain or anthocyanin physiology.
