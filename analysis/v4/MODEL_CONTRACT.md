# Shared model component

Fixed before fitting this workflow's selected five-environment specification.
This is a revised analysis informed by previous work, not independent replication.

All response models retain endpoint-specific eligible observations. Species
means of environments are calculated on the native population before response
selection. Each response has a species random intercept and both within and
between terms for every selected environment. Existing REML crossproduct code
is reused and checked against dense GLS and statsmodels MixedLM in Actions.
It does not replace observations with species means.

Stacked species-score covariance supports joint tests for colour and outline.
These are marginal component models, not an estimated full multivariate
random-intercept covariance model. Do not interpret covariance of coefficient
estimates as biological trait covariance.

Main environmental family: three modules x five exposures x two scales = 30
joint tests, Holm adjusted. SI bract: 10 tests, separately Holm adjusted and
explicitly supplementary. Component intervals are pointwise, not simultaneous.
Report all slopes and uncertainty, especially the V3 rainfall/orientation and
radiation/colour anchors, regardless of adjusted significance.

Colour uses L*, C* and a joint sine/cosine hue embedding. These do not measure
anthocyanin concentration. No synthetic Lab a*/b* from mixed colour systems.

The environment component is only one stage of the V3-question workflow.
Raw and environmentally adjusted integration, spatial uncertainty and sampling
assessment remain required before consolidated conclusions.
