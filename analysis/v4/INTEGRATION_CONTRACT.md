# Native phenotype integration contract

Purpose: compare within- and between-species phenotype associations, before and
after removal of the shared model's fixed environmental component. This is
phenotypic covariance, not genetic covariance or a fitted multivariate random
effect covariance matrix.

## Features and samples

Eight named constructs: angle, lightness, chroma, circular hue (sine and cosine
together), aspect ratio, circularity, solidity and width-profile CV. These are
nine measured coordinates in three main biological modules. No phenotype PCA,
invented Lab coordinates or averaging of multiple outline traits is required.
Bract remains separate SI and never constrains the main integration cohort.

Report two scopes: pair-specific eligible observations, and the intersection of
all nine main coordinates for comparing the complete association geometry.
Use selected-environment-complete observations in both raw and adjusted versions.
For each scope retain species with at least five paired/complete observations;
require at least 20 species. Keep samples identical before and after adjustment.
These requirements apply to covariance estimation, not the full-data gradient
models. Report the resulting counts separately.

Species means of environmental exposures and scaling stay those of the native
population. Subtract only fixed between/within environmental contributions using
the saved endpoint model coefficients. Do not subtract fitted species intercepts:
doing so would erase the between-species variation being studied.

## Estimands and uncertainty

Within: species-centred observation covariance with equal total species weight.
Between: equally weighted species means of the same paired observations.
Between means define a descriptive scale; the gradient models still fit raw
responses. These means carry finite within-species sampling variation and are
not inferred latent genetic or species random-effect values.

Scalar pairs: signed correlations as well as RV (= squared correlation for scalar
pairs). Hue pairs: RV of the two-dimensional circular embedding, never independent
sine/cosine ecological tests. Fix coordinate scaling from the full eligible
native trait distribution, using a shared scale for hue sine/cosine so rotation
of hue zero does not alter its geometry.

Report 1,000 species-cluster bootstrap intervals. Raw/adjusted and within/between
comparisons use the same resampled species. Adjusted intervals are explicitly
conditional on the saved environmental coefficients; they do not propagate
coefficient estimation uncertainty. Thus raw-minus-adjusted contrasts describe
the fitted decomposition, not a formal claim that environment causally explains
a percentage of biological integration.

RV has a dimension- and sample-size-dependent baseline. In the common cohort,
repeat V3's 1,000 equal-n comparisons: one centred observation per species versus
one species mean per species, with the same draw for raw and adjusted data.
Those replicate ranges are sampling diagnostics, not confidence intervals.
Do not infer stronger between-species integration solely from uncorrected RV.

No pairwise discovery P-value family is added in this stage. Display all pairs,
their intervals and paired contrasts. Environmental joint tests retain the
previously fixed main/SI families. Spatial assessment remains a later required
stage; species bootstrap alone does not resolve spatial dependence.
