# Supplementary Section S4 — BIO18 × image angle: geographic sensitivity after the original V4 freeze

This analysis is exploratory and **postdates** the fixed S1 80-test Holm family. It was calculated on the **original, cryptographically authenticated private native S1 observation source** with independently resampled CHELSA BIO18/BIO05 and WorldClim elevation pixels (see [source-verification report](../../analysis/v4/results/20261029_BIO18_GEOGRAPHIC_EMPIRICAL_RESULT_JA.md)). The frozen V4 between-taxon coefficient, original 80 Holm decisions and all other original results **remain unchanged**.

## Sampling and analysis estimand

- Original S1 source records: **313,350**. Five-climate angle eligibility: **170,066**.
- New four-predictor source-level completeness: **169,020** angle-and-environment eligible records. With original name-conflict and ≥5-photo gates, the equal-taxon geographic models include **171 source taxa**. The original publication's primary model used 168,907 records/178 taxa from a different all-environment completeness definition.
- All fitted models use BIO18 + BIO05 + shortwave radiation + elevation; outcome is head orientation **relative to the image frame**, never a gravity-based nodding angle.
- Unit: **degrees of image angle per 100 mm/month BIO18**. Equal taxon weights, fixed numeric scaling and HC3 pointwise 95% intervals. These taxon-mean regressions do not equal the originally frozen observational decomposition of within-/between-taxon coefficients.
- Geographic alternatives: no geographic term; two indicators of source-photo **longitude-sector composition**; fixed eight-variable second-order spherical coordinate basis. Within-taxon geographic sampling is (i) original photo mean and (ii) equal 5° cell-mean. Fixed regions are geographic bins, not independent phylogenetic radiations.

## Table S4.1 — Geographic adjustments

| Within-taxon photo weighting | No geographic controls | Fixed longitude-sector proportions | Spherical second-order basis |
| --- | --- | --- | --- |
| Original photo-weighted | **+0.934** [+0.380,+1.489] | +0.437 [−0.267,+1.141] | +0.060 [−1.421,+1.541] |
| Within-taxon 5° cells equally weighted | **+0.983** [+0.381,+1.585] | +0.678 [−0.050,+1.406] | +0.566 [−0.932,+2.064] |

All six models use **171 taxa**. The four geography-adjusted pointwise intervals all intersect zero, whereas the two no-geography intervals do not. These are several exploratory analyses, not six independent multiple-testing-adjusted discoveries.

## Table S4.2 — Leave one fixed longitude sector out

| Excluded sector of the dominant observed source | No-geo photo-weighted | No-geo within-taxon 5°-cell balanced |
| --- | --- | --- |
| 180°W–30°W | +1.216 [+0.429,+2.002]; n=87 taxa retained | +1.287 [+0.451,+2.123]; n=87 |
| 30°W–60°E | +0.785 [+0.227,+1.342]; n=117 | +0.853 [+0.243,+1.462]; n=119 |
| **60°E–180°E** | **−0.186 [−1.495,+1.123]; n=138** | **−0.558 [−2.012,+0.897]; n=136** |

Original photo-weighted band counts are 84,54,33 taxa; 5° balanced counts 84,52,35. The east-band removal sharply reduces the positive global coefficient, suggesting regional taxon composition is important to this observational association.

## Table S4.3 — Within fixed longitude sector

| Region | Photo-weighted n taxa | Photo-weighted slope and pointwise 95% CI | 5°-balanced n | 5°-balanced slope and CI |
| --- | ---: | --- | ---: | --- |
| 180°W–30°W | 84 | +0.003 [−1.178,+1.184] | 84 | −0.314 [−1.705,+1.077] |
| 30°W–60°E | 54 | −0.228 [−6.730,+6.275] | 52 | −2.054 [−8.588,+4.481] |
| 60°E–180°E | 33 | +0.792 [−0.470,+2.052] | 35 | +1.178 [−0.008,+2.365] |

No fixed longitude sector individually has a positive 95% pointwise interval excluding zero under these four-predictor, taxon-level diagnostics. This is **not proof of no regional effect**: the CIs are wide and geographic sampling is uneven.

## Source verification and interpretation

- Original decoded source SHA-256 `342d78e9a8a960e8d9bdcf0bc2cdacd13bee48214a8b962419d5a71d5aef46ed`.
- Original rasters sampled for this follow-up (source SHA-256): BIO18 `68b0b1cb905cdf6972364e947e137f714a68d287f38863b1516161a2a4bbd7ec`; BIO05 `90078d99eeb0b08dd1a115a5ffa60b7c19e5849686854db950283cbac9d3817d`; WorldClim elevation `8b9ac1eddbb1e51ea4aa860f637a5ec30096e1668a85303febe0d91abf901017`.
- Source-real GitHub Actions [#40829203411](https://github.com/zuizui0223/azami/actions/runs/40829203411), aggregated artifact **12877544186**. Content: `S1_BIO18_GEO_SPATIAL_AGGREGATES.json`; original aggregate JSON SHA-256 `1e9563df0abfa0bda8b77453c4bf45d13e78c1f63b0d1c767147320c03d795c9`.
- **Interpretation:** Both rainfall and presentation are geographically patterned among taxa. The global positive association could reflect (at least in part) the geographic turnover of taxon assemblages. Strong spatial covariates and small within-region support can also make a true effect imprecisely estimated. Without phylogenetically structured replication, photographs independently calibrated to gravity, and direct rain-exposure observations, neither a global repeated adaptive rain-avoidance response nor its negation is established.
- The previous BIO12 analysis used **annual** rainfall and a **five-predictor** model, and so is useful context rather than an independent replication of this BIO18/four-variable result.
- Original individual photographs, coordinates, observer IDs and source files were not included in the aggregate artifact or this SI; rights-compliant anonymous long-term data/code deposition remains a submission prerequisite.
