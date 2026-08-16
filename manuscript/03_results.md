# Results

## Global image phenomics recovered repeated continuous capitulum phenotypes

The balanced atlas contained 6,626 detected heads from 3,725 photographs and 216 source-assigned taxa, whereas the exhaustive post-detection stream contained 406,582 detector-positive observations from 286 taxa. After positional-accuracy restriction and taxon × 0.25° spatial thinning, 46,276 observations from 259 taxa formed the primary within-taxon cohort. Trait-specific QC retained 87.2% of atlas heads for colour, 80.4% for outline and 69.2% for orientation. The high-resolution involucral layer retained 1,443 heads from 1,292 observations and 210 taxa.

The leakage-free detector-audit packet contained 1,000 source images spanning 323 species and 85 ten-degree spatial blocks, with 250 images assigned for blinded double annotation and zero photo or observation overlap with detector-development data. Human adjudication is incomplete; therefore no independent detector precision or recall is reported, and detector development metrics are not used as evidence for biological accuracy.

The production workflow therefore generated repeated **numerical phenotype observations** rather than one categorical trait label per taxon. Every assessable head contributed values such as an orientation angle, colour coordinates or contour statistics, and those values were retained as distributions across observations within taxa.

## Most visible primary-trait variation occurred below taxon means

Across the nine primary endpoints, the combined fraction of visible image sums of squares below assigned taxon means ranged from 0.589 for hue sine to 0.931 for width-profile variation. The among-photograph component within taxa was 0.440–0.691, while variation among heads within the same photograph contributed 0.143–0.379.

Selecting one deterministic head per photograph retained below-taxon fractions of 0.582–0.899. In equal-replication sensitivities using 10 photographs per eligible taxon, median fractions remained 0.528–0.879. Thus, the large below-taxon component was not explained solely by multiple detections from one photograph or by unequal taxon replication.

Taxon-level means nevertheless retained substantial multivariate structure. PC1 explained 32.9% of variance, PCs 1–2 explained 56.1% and PCs 1–3 explained 69.3%, with all nine primary endpoints contributing to the PCA. Variation was therefore substantial at both within- and among-taxon scales rather than reducible to one taxon-level trait axis.

## Within-taxon environmental associations were trait specific

In the exhaustive spatially thinned primary cohort, eight of 36 endpoint-component rows passed BH correction.

| Endpoint | Predictor | Standardized β | 95% CI | BH q |
|---|---|---:|---|---:|
| Orientation angle | BIO1 annual mean temperature | +0.0171 | 0.0046 to 0.0295 | 0.0363 |
| Corolla chroma | BIO12 annual precipitation | +0.0393 | 0.0102 to 0.0684 | 0.0363 |
| Hue sine | BIO4 temperature seasonality | +0.0534 | 0.0424 to 0.0643 | <0.0001 |
| Hue sine | BIO12 annual precipitation | -0.0440 | -0.0608 to -0.0272 | <0.0001 |
| Hue cosine | BIO4 temperature seasonality | -0.0283 | -0.0453 to -0.0113 | 0.0132 |
| Hue cosine | BIO15 precipitation seasonality | +0.0201 | 0.0069 to 0.0334 | 0.0207 |
| Outline aspect ratio | BIO4 temperature seasonality | +0.0137 | 0.0048 to 0.0226 | 0.0207 |
| Outline aspect ratio | BIO12 annual precipitation | -0.0108 | -0.0187 to -0.0028 | 0.0363 |

Higher annual mean temperature was associated with a larger image-referenced orientation angle, corresponding to a shift from upward toward more horizontal or downward heads relative to image vertical. Corolla chroma increased with annual precipitation. Outline aspect ratio increased with temperature seasonality and decreased with annual precipitation. Circular hue components were associated with temperature seasonality, annual precipitation and precipitation seasonality, but the biological colour direction cannot be assigned from one sine or cosine coefficient alone.

## Grouped SPDE-INLA concentrated robust support in orientation and visible colour

All 36 grouped SPDE-INLA fits completed with no CPO failures. The most stable spatial-model effects were:

- orientation angle: positive BIO1 association;
- corolla lightness: positive BIO1 association and negative soil-pH association;
- corolla chroma: negative BIO1 association and positive BIO12 association;
- hue sine: positive BIO4 association;
- hue cosine: negative BIO1 association and positive soil-pH association.

The BIO1 effect on hue sine was q-supported only in the climate model and changed sign across model groups, so it is retained as model-dependent rather than a stable direction. No topographic predictor and no SPDE effect for aspect ratio, circularity, solidity or width-profile CV passed global BH correction. Thus, the aspect-ratio climate associations are specific to the exhaustive standardized within-taxon analysis, whereas soil-pH colour associations emerge only when wider spatial predictor sets are considered.

## High-resolution involucral architecture tracked temperature seasonality

In the ≤10 km auxiliary cohort of 904 observations from 165 taxa, all three final auxiliary contour proxies increased with temperature seasonality and survived BH correction across the 12-test auxiliary family:

- involucre projection roughness: β = +0.0975, 95% CI 0.0227–0.1723, q = 0.0424;
- outward spread fraction: β = +0.0937, 95% CI 0.0271–0.1603, q = 0.0424;
- maximum spine-like projection: β = +0.0911, 95% CI 0.0222–0.1599, q = 0.0424.

The corresponding all-coordinate slopes had the same positive direction but did not survive FDR correction. These associations concern outward contour architecture and do not identify whether the underlying structure is a longer bract tip, greater spreading or recurvature.

## Among-taxon environmental sorting was concentrated in orientation and colour

Among 148 taxa complete for all primary traits and environmental variables, 10,000 trait-label permutations showed non-random environmental sorting mainly for orientation and visible colour. Orientation, chroma and both hue components retained BH-supported centroid separation and lower-than-null niche overlap. Width-profile CV retained centroid-distance support only. Lightness, aspect ratio, circularity and solidity were unsupported in both all-taxa metrics.

In the 49-taxon complete direct-backbone sensitivity, support narrowed. Hue cosine retained both centroid and overlap support; lightness and hue sine retained overlap support only. Thus, the strongest all-taxon environmental sorting was not a generic artifact of quartile splitting, but its historical robustness varied among traits.

## Some among-taxon climate associations persisted across alternative phylogenetic placements

Taxon-level PGLS identified six associations that were BH-supported in 100% of the 50 randomized scenario-2 trees:

| Endpoint | Predictor | Median standardized β | 2.5–97.5% across trees | FDR-supported trees |
|---|---|---:|---|---:|
| Orientation angle | BIO12 annual precipitation | +0.304 | +0.304 to +0.313 | 50/50 |
| Hue sine | BIO1 annual mean temperature | -0.410 | -0.410 to -0.410 | 50/50 |
| Hue sine | BIO4 temperature seasonality | -0.274 | -0.274 to -0.274 | 50/50 |
| Hue sine | BIO12 annual precipitation | -0.211 | -0.211 to -0.211 | 50/50 |
| Solidity | BIO4 temperature seasonality | -0.357 | -0.357 to -0.347 | 50/50 |
| Width-profile CV | BIO1 annual mean temperature | +0.248 | +0.234 to +0.269 | 50/50 |

Pagel-λ was estimated at zero in the randomized-tree fits for these rows, so coefficient estimates were little altered by the tested covariance structures. However, only 54 of 216 image-analysis taxa were direct dated-backbone tips and 162 required within-genus grafting. Independent phylogenetic-signal analyses found no FDR-supported non-circular primary trait in the direct-backbone subset, and no auxiliary involucre/spine coefficient survived deterministic PGLS FDR correction. The six randomized-tree associations are therefore historical sensitivities rather than definitive evidence of resolved phylogenetically independent effects.

## Spatial and taxonomic robustness did not overturn the headline results

Diagnostic refits of the lowest-WAIC SPDE specification closely reproduced the frozen model support, with diagnostic-refit minus frozen WAIC ranging from -2.56 to +2.50. Residual Moran's I ranged from -0.0096 to 0.0030 across the nine primary endpoints, with no permutation P < 0.05. Removing one broad geographic region at a time retained taxon trait rankings strongly, with minimum endpoint-specific Spearman rho = 0.856–0.972. The weakest sensitivities followed removal of Europe.

WCVP review identified eight synonym candidates whose accepted-name targets were already separate active source units. Collapsing all eight simultaneously reduced the balanced atlas from 216 to 211 operational units and the exhaustive primary cohort from 259 to 251, but all eight BH-supported primary component rows remained supported, no coefficient changed sign and the maximum absolute standardized-beta change was 0.000385. The minimum below-unit visible-variance fraction increased from 0.589 to 0.597.

The computational spatial, niche-null and taxonomic-robustness diagnostics therefore did not overturn the central conclusions. Independent detector accuracy and continuous-trait biological validity remain the unresolved external scientific gates.

## Legacy lability analysis is provenance only

The former RMS of unpooled absolute taxon-specific slopes was strongly confounded with slope sample size and standard error. The previously reported negative variation–association relation and median-split quadrants were withdrawn. Precision-aware reanalysis detected no common cross-module coupling. This correction remains part of the audit trail but is not a submission-facing biological headline.
