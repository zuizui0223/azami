# Supporting Information: post-freeze S1 BIO18 spatial-model sensitivity (29 October 2026)

This is **an ecological-model specification sensitivity, not evidence that climate has no association with capitulum presentation**. The biological interpretation in GEB Results 3.1 and 3.2 remains the original frozen source-cohort association and the original introduced-stratum colour associations; this analysis does not add tests to either original Holm-adjusted family.

## Spatial-block resampling without latitude/longitude regression covariates (new source-verified result)

An additional exact-source sensitivity used **the original four environmental predictors only**, with no geographic terms added to the regression design. Each of the original 178 taxa was assigned to one dominant **20° geographic block according to the number of occupied 5° cells** (not the number of photographs). The frozen among-taxon coefficient was recovered exactly (**β=+0.045446843**); bootstrap resampling drew whole taxon blocks, retaining each taxon's full multivariate climatic and photographic mean.

- **26** source-observed blocks; **399/399** predeclared block draws identifiable.
- Spatial-block bootstrap percentile **95% interval [+0.014569,+0.084483]**, mean slope +0.044476.
- **99.50%** of block bootstrap samples retained the positive original BIO18 sign.
- Dropping one at a time the five blocks with the largest **source-only taxon contributions** retained a positive coefficient in **5/5** fits; the five range from **+0.03374 to +0.05048**.
- The earlier independent 5° cell-equalisation without geographic terms gave β **+0.048516**, HC3 interval [+0.022397,+0.074635], consistent in direction with the photo-weighted source.
- **No climate-versus-space causal partition is implied**: these blocks are not phylogenetic independent replicates, and climate and geography co-vary. The bootstrap interval is **exploratory**, not the 80-test-family Holm-adjusted CI.

Source-verifiable artifact: [Actions #40855350461](https://github.com/zuizui0223/azami/actions/runs/40855350461), aggregate artifact **12889980119**, [stable machine-readable block receipt](../../analysis/v4/results/20261029_S1_BIO18_SPATIAL_BLOCK_SOURCE_VERIFIED.json). Together with the separate geographic-model sensitivity below, this documents **two different meanings of sensitivity**: the association's positive direction is stable to reweighting/resampling geographic units, whereas its partial coefficient is sensitive to conditioning on spatially aligned geographic basis terms. Neither cancels or proves the cause of the other. The positive original 3.1 result remains a global observational association, not a geographically universal or adaptive rainfall response.

## Source and estimand

- Exact original V4 angle-eligible native S1: **168,907 observation records and 178 taxa**, obtained from the private, checksum-audited photograph/climate source.
- Original model: image-referenced head angle ~ CHELSA BIO18, BIO05, RSDS and WorldClim elevation, separating within- and between-taxon components. Original between-taxon standardised β **+0.045446843** and original Holm P **0.00920635058** were reproduced independently in [Actions #40837107465](https://github.com/zuizui0223/azami/actions/runs/40837107465). The model used for the added taxon-level comparisons has the same coefficient but a different HC3 standard error.
- New specifications: (1) unadjusted geography, (2) three response-blind longitude-sector proportions, (3) eight spherical quadratic spatial terms. Each with original photo-weighted taxon means and independently with occupied 5° cells equalised within taxa. The same fixed covariates enter every model.
- The new **HC3 confidence intervals are pointwise and exploratory**; 3 longitude sectors are not biogeographic realms or phylogenetic sister-clade replicates. No new multiple-testing-family conclusion is drawn.

## Source-reproduced findings

| Within-taxon sampling scheme | Geographic model | Standardised BIO18 β | Pointwise HC3 95% CI |
|---|---|---:|---|
| Photographs | No spatial covariates | +0.04545 | [+0.02156, +0.06934] |
| Photographs | 3 longitude-sector proportions | +0.02003 | [−0.01083, +0.05090] |
| Photographs | Spherical second-order terms | −0.00329 | [−0.06300, +0.05641] |
| Equal 5° geographic cells | No spatial covariates | +0.04852 | [+0.02240, +0.07464] |
| Equal 5° geographic cells | 3 longitude-sector proportions | +0.03093 | [−0.00122, +0.06307] |
| Equal 5° geographic cells | Spherical second-order terms | +0.00927 | [−0.05102, +0.06955] |

Leaving out 38 east-sector dominant taxa yielded **−0.01037** (HC3 CI −0.06846, +0.04772; 140 retained) under photograph weighting; equal-cell weighting left 138 taxa and yielded **−0.02565** (CI −0.09033, +0.03903). All source-real aggregate coefficients and code are preserved in [the 2026-10-29 archived JSON](../../analysis/v4/results/20261029_S1_BIO18_GEOGRAPHY_SOURCE_VERIFIED.json) and the [interpretive audit](../../analysis/v4/results/20261029_S1_BIO18_GEOGRAPHY_SOURCE_VERIFIED_JA.md).

## Interpretation boundary

Spatial predictors share information with climate because geographic variation in precipitation is itself spatially structured; including a broad geographic basis can remove both confounding and the part of the *real* environmental gradient that overlaps geography. Consequently the original geographic association **is model-dependent**, but no causal partition of a "spatial effect" and an "environmental effect" follows from the difference between coefficients. Nor does failure to exclude zero in new pointwise sensitivity intervals establish the absence of a causal role of rainfall. The data do not measure head inclination to gravity, direct rainfall shielding, pollen survival, or evolution of posture.

This analysis remains **SI only**. Further mechanistic discrimination would need phylogenetic contrasts, overlapping climate gradients within regions and calibrated photographs/field measurements, rather than declaring the original BIO18 result false.
