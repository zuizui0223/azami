# V5 photographed-angle flag prevalence audit (2026-11-20)

**Source-verified diagnostic, not an estimated misclassification rate.** The question is: how many originally analysed head-image measurements are potentially difficult to interpret as true upward/downward orientation, given the projection and principal-axis algorithm?

## Completed original S2 79,421-measurement-record audit

- **Successful private source GitHub Actions run:** https://github.com/zuizui0223/azami/actions/runs/44666062718
- Original WCVP native/explicit-introduced labels, three original image-measured head traits, exact original five CHELSA v2.1 climate variables and filtering; exactly 79,421 records matched original V5 analysis strata.
- Original measurements, IDs, user identifiers and private photo coordinates remain in Actions runner. Only aggregate counts are printed/uploaded.
- Comparison unit is an **observation/photo-derived measurement record**, not a distinct JPEG or detected head.

| Diagnostic candidate (not an error label) | N / 79,421 | Fraction |
|---|---:|---:|
| Rotated projected outline aspect near 1 (1/1.10 ≤ aspect ≤ 1.10) | 50,254 | 63.28% |
| Image-based angle > 90° | 25,167 | 31.69% |
| Both near-aspect and angle > 90° | 15,702 | 19.77% |
| Both near-aspect and angle ≥ 110° | 8,786 | 11.06% |
| Aspect within 1/1.05–1.05 | 29,908 | 37.66% |
| More liberal near 1/1.20–1.20 | 69,655 | 87.70% |

The fraction of angle >90° was **31.25% among near-aspect records** and **32.45% among others**. A near-one shape index alone is not evidence that the source angle is wrong; these near-aspect measures are not the PCA eigenvalue anisotropy in the original angle estimator.

### Four verified source strata

| Species × WCVP status | Original V5 N | Near aspect 1.10 N | Both near and angle>90 N (% of stratum) |
|---|---:|---:|---:|
| *C. arvense* × native | 18,505 | 11,707 | 3,942 (21.30%) |
| *C. arvense* × introduced | 22,269 | 14,432 | 5,123 (23.01%) |
| *C. vulgare* × native | 13,115 | 8,229 | 2,231 (17.01%) |
| *C. vulgare* × introduced | 25,532 | 15,886 | 4,406 (17.26%) |

Across equal-record annual-rainfall quintiles in these two species, the joint near-aspect+angle>90 diagnostic fraction ranged 19.27–20.36%; no clear monotonic increase appears in this *unadjusted source-record diagnostic*. This does not exclude photographer, taxon, phenology or spatial confounding, and is not the same population as Part I.

## Interpretation and next checks

The source continuous angle is computed from PCA of foreground projected head coordinates, with the axis signed by the visible floral-pixel centroid and compared with **EXIF-oriented image-up**, not gravity. A longitudinal PCA axis may be uncertain if eigenvalues are similar or the floral-end sign is weak, but the registered **rotated bounding-box aspect ratio** in these records is **not** the eigenvalue ratio and cannot establish this instability on its own.

The *C. pannonicum* (0.732°), *C. mohavense* (111.148°) and *C. subcoriaceum* (113.284°) archived poster examples are illustration-only and are not a representative validation sample. Counting cases near 1 as actual mistakes would overstate the observed evidence.

**Next:** obtain the exact S1 annual V5 157,822-record flag prevalence and its two grid-specific subsets; then independently blind-annotate at least a geographically/taxonomically stratified subset of image orientations, including near-circular and elongated phenotypes. Refit the BIO12–image-angle gradient and multivariate analogue test after reweighting or excluding *independently validated* uncertain angles, keeping the same geography/taxon eligibility for the sensitivity when possible.

S1 source issue: the legacy receipt artifact `v4-all665115-native-climate-aggregate` was not returned from historic run 39408603488 despite the encrypted input being available. The original complete WCVP metadata and measurement Release assets, rather than that missing receipt, provide an independent source-authenticated recovery route. Pending its exact cohort-count and grid-count gates, **the S2 79,421 rates cannot be generalized to all 157,822 native Part-I records**.
