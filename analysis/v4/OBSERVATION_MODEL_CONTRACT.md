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
