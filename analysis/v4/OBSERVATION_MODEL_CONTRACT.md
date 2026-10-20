# Observation-level primary analysis, 2026-10-20

This replaces the species/cell-mean exploratory models. The already inspected
predictor correlations are known; these are retrospective hypothesis evaluations,
not preregistered confirmatory tests.

Use encrypted input artifact 12257794094 from run 39408603488 and its aggregate
receipt 12257893990. Keep every eligible observation separately for each response.
Retain species with fewer than five records. Do not require all colour coordinates
or other trait modules to be present for a response's fitted marginal model.

Three modules:

- Whole visible colour: L*, C*, archived hue sine and cosine jointly tested.
- Presentation: image-vertical orientation angle (higher means more downward in
  the image reference, not an independently measured gravity angle).
- Bract contour projection: relative 95th-percentile projection and outward
  spread fraction jointly tested. These are existing image descriptors, not
  anatomical phyllary opening, botanical spine identity, or measured defence.

All modules use temperature, radiation, VPD and precipitation simultaneously.
For each predictor include the species mean over the native source population
and the observation's deviation from that mean. Standardize both by the same
native-population predictor SD. Do not average trait responses.

Fit a Gaussian random species intercept for each response by REML. Sufficient
crossproducts reproduce the observation-level likelihood; they do not change
the observational unit. Stack species-cluster scores across responses for joint
module Wald F tests, including cross-response covariance and endpoint-specific
missingness. This is a marginal working-covariance approach, not a full
unrestricted multivariate Gaussian random-effects model. Inference is approximate
and conditional on estimated variance components; species-specific slopes and
residual spatial dependence are not fitted here.

Primary multiplicity family: 3 modules x 4 climate predictors x 2 scales = 24
joint tests, Holm adjustment. Component intervals are pointwise descriptions of
direction, not additional significance screens. Unidentifiable or nonconverged
fits remain explicit and invalidate completion of the corresponding inference.

Sensitivity: repeat the identical model after removing the two and five taxa
with most records in the native source, fixed before fitting. Do not select
exclusions to improve a coefficient or P value. These are influence diagnostics.
Soil, broad-space and phylogenetic checks remain additional stages.

GEB positioning: broad-scale, general ecological hypotheses and continuous
within-species variation beyond categorical species descriptions. Relate observed
patterns to functional predictions; do not equate associations with adaptation,
pigment concentration or physical defence. Journal scope:
https://onlinelibrary.wiley.com/page/journal/14668238/homepage/aims.htm

## Hypothesis-led interpretation and attribution audit

The scientific framing is hypothesis-led observational evaluation. Existing
outcomes have already been inspected: neither the 24 tests nor the following
diagnostics are retroactively preregistered confirmation.

- Rain protection predicts a more downward presentation in wetter environments.
  Evaluate between- and within-species coefficients separately. Annual rainfall
  is exposure context, not measured floral wetting; image vertical is not gravity.
- Visible colour can reflect several functional pathways. Test the whole colour
  response and describe lightness, chroma and circular hue jointly. Do not impose
  a universal temperature/pigment sign, or identify anthocyanin from photographs.
- Protective-envelope function requires validated botanical geometry. Archived
  contour projection is a limited proxy, not a measurement of antagonist access.
  Its 1,303 observations are inherited coverage, not a full-source remeasurement.
- Trait covariation can distinguish phenotypic coupling across scales. A direct
  comparison is required; significant in one scale and nonsignificant in another
  is not evidence of a difference. Phenotypic covariance cannot establish
  genetic LEGO, developmentally independent modules or selection.

Before further spatial/phylogenetic extensions, audit the exact endpoint-specific
eight-column design. Retain the four-predictor model. Supplement it with joint
temperature-plus-VPD tests (3 modules x 2 scales, a separate six-test Holm family)
and same-observation omission diagnostics for temperature and VPD. Omission fits
show attribution sensitivity, not a criterion for choosing a favourable model.
Do not call design correlations causal confounding estimates. If independent
attribution is unstable, report joint environmental association and uncertainty,
rather than claim temperature or VPD is the unique mechanism.

Source comparison inspected on 2026-10-20:

- [EAzami V10](https://github.com/zuizui0223/EAzami/blob/c9870e91604ed1514344ddc8f765eae5237dbf27/docs/chapter2/MANUSCRIPT_JEB_V10_FUNCTIONAL_ASSEMBLY_ORIGIN_DECOUPLING.md):
  repeated functional assembly and present ecological association do not identify
  a historical selective trigger; a post-result orientation regime is labelled.
- [aza3 higher-order LEGO](https://github.com/zuizui0223/aza3/blob/2dabd7750a80686f002849873ab941b180135131/docs/ANTO_HIGHER_ORDER_LEGO_MEASUREMENT_PLAN_V2.md):
  same-individual functional/genomic resources and fitness evidence are additional
  requirements, beyond photographic covariance.
- [aza3 empirical head comparison](https://github.com/zuizui0223/aza3/blob/2dabd7750a80686f002849873ab941b180135131/docs/CIRSIUM_HEAD_CORE_VS_BRACT_TIP_SPREAD_EMPIRICAL_20261020.md):
  reproductive dimensions and relative external spread can differ in coupling;
  the existing comparison is exploratory and does not establish genetic modules.
