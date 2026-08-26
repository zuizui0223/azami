# Methods

## Study design and cohort separation

We used public biodiversity photographs to quantify continuous capitulum traits across *Cirsium* and related thistles while retaining repeated observations within source-assigned taxa. Taxon names attached to the frozen public observations were retained as operational taxonomic units rather than treated as a resolved genus-wide species taxonomy. This preserved exact observation provenance while allowing authority-backed lumping to be tested explicitly as sensitivity.

Two executed data streams answered complementary questions. The balanced image-comparison atlas contained 3,725 observations, 3,725 photographs, 6,626 detected capitula and 216 source-assigned taxa. Each atlas observation contributed one photograph, while photographs could contain multiple heads. This layer was used for trait extraction summaries, nested visible-variance decomposition, taxon-level trait PCA and among-taxon historical analyses.

The exhaustive stream ran detection and continuous measurement before trait-based thinning. It contained 406,582 observations with detected capitula from 286 taxa. Coordinate filtering retained 392,989 observations from 271 taxa; a positional-accuracy threshold of ≤10 km retained 297,293 observations from 259 taxa; and one observation per taxon × 0.25° cell yielded the exhaustive spatially thinned primary cohort of 46,276 observations from 259 taxa. These cohorts were not pooled under one FDR family.

## Trait scope, capitulum detection and continuous image measurements

Trait coverage was module-complete for the capitulum structures targeted by the study but explicitly tiered by what a photograph can support. For every endpoint, a frozen measurement contract records (i) the biological construct, (ii) the numerical image variable and formula, (iii) its unit and bounds, (iv) the required view, resolution and focus, (v) an endpoint-specific usable versus unassessable rule, (vi) its inferential tier and botanical-validation status and (vii) interpretations that the image variable cannot support. Fine structures such as phyllary projections, armature, indumentum, glandularity and visible exudate signals were therefore assigned continuous image variables rather than whole-head states. When the required structure or resolution was absent, the endpoint remained missing; it was never coded as absence. Absolute dimensions without a scale reference and whole-plant characters remained outside the contract.

Capitula were localized with a frozen single-class YOLO11n detector initialized from `yolo11n.pt` and trained for 40 epochs at 640-pixel input size on 270 pseudo-labelled source images. The production confidence threshold was 0.25. Detector output defined regions of interest only; downstream trait inference did not use the legacy species-level upward/nodding classifier. Each retained detection preserved source observation and photo identifiers, bounding-box coordinates and contextual crops.

The contract contains 27 continuous endpoints across orientation, visible corolla colour and display, capitulum outline, involucral architecture, armature, indumentum-related surface texture and visible exudate signal. **Continuous** means that each assessable detected head returns a numerical measurement on a defined image-derived scale rather than membership in a biological category. The inferential tiers were frozen before the extended run.

- **Primary:** signed head-axis angle relative to EXIF-oriented image vertical; CIELAB lightness and chroma; sine and cosine components of circular hue; outline aspect ratio, circularity, solidity and width-profile coefficient of variation.
- **Candidate:** visible-floret fraction; involucre length/width and apical/basal taper ratios; radial projection roughness, 95th percentile, maximum and spread fraction; projection-peak density and bilateral projection asymmetry.
- **Descriptive composition:** four mutually exclusive visible floral-pixel fractions. Because these fractions are compositional, they are not fitted as four independent univariate outcomes.
- **Validation-only:** involucre-surface edge density, normalized local-binary-pattern entropy, high-frequency energy and specular-pixel fraction.

Measurements were deterministic image-processing functions. Horizontal mirroring served as a technical replicate, and failed or unstable measurements remained missing rather than being coded as biological absence. In the 6,626-head atlas, 5,777 heads were usable for colour, 5,324 for outline and 4,585 for orientation.

The continuous design retained information that would be lost by forcing images into categories such as upright/nodding or globose/cylindrical. Categories may be used descriptively for examples or maps, but inferential analyses use the numerical endpoints. Hue sine and cosine are mathematical components of one circular response and are interpreted jointly. Image vertical is a reproducible image reference rather than a direct inclinometer measurement of gravity.

A leakage-free detector-audit packet was selected independently of detector development. All proposal/training photo and observation identifiers were excluded before sampling, yielding 1,000 source images across 323 species and 85 ten-degree spatial blocks, with 250 images assigned to a second annotator. Because adjudicated human boxes are not yet complete, detector precision and recall are not used as biological evidence.

## Exhaustive high-resolution architecture, armature and surface cohort

The earlier 3,725-photograph atlas served as a development and variance-partitioning resource but no longer defines the high-resolution inferential sample. Extended measurement starts from all retained YOLO head metadata, intersects the frozen strict spatial cohort and applies only predeclared image gates. Heads must have a minimum detected-box dimension of 150 pixels and an available source-image URL. One head per observation is retained deterministically by the largest minimum box dimension, with detector confidence and stable identifiers used only as tie-breakers. Trait values and categorical states are forbidden from cohort selection. This produces a complete attrition audit from all queued photographs through detector-positive observations, strict spatial eligibility, high-resolution eligibility and endpoint-specific usability.

The recovered involucre region is summarized continuously by length/width ratio, apical and basal taper, radial-projection roughness, 95th-percentile and maximum projection, spread fraction, local peak density and bilateral asymmetry. The former `spine_relative_length_max_proxy` is not a second trait because it duplicated the maximum-projection formula; `involucre_projection_max` is the sole canonical endpoint. These quantities describe two-dimensional contour geometry. They do not by themselves identify individual phyllaries, distinguish spine length from bract spreading or recurvature, or measure stiffness.

For crops with a minimum head dimension of 300 pixels and stricter focus criteria, surface edge density, texture entropy, high-frequency energy and specular-pixel fraction are also calculated. They remain validation-only image signals: edge density is not called hair or gland density, and specular fraction is not called mucilage presence, stickiness or secretion amount. Promotion requires a double-annotated continuous reference set containing involucre axes and polygons, visible phyllary base/tip keypoints, free-tip length, local-surface-relative phyllary angle, tip curvature and annotated hair, gland and exudate areas. Required validation includes annotator agreement, error and rank preservation, calibration slope/intercept, view- and resolution-specific performance, horizontal-mirror repeatability and taxon-held-out evaluation.

## Visible variation within and among taxa

For each assessable primary endpoint in the balanced atlas, total visible image sums of squares were decomposed exactly into differences among assigned taxon means, among photographs/observations within assigned taxa and among detected heads within photographs. Because the atlas retained one photograph per observation, photograph and observation were the same nested level. Assigned taxa were resampled as clusters 2,000 times for uncertainty intervals.

Two sensitivities reduced unequal replication and within-photograph multiplicity. First, one deterministic head was retained per photograph. Second, 10 photographs per eligible taxon were repeatedly sampled without replacement across 500 balanced replicates. Taxon medians for all nine primary endpoints were standardized and analysed by PCA to describe among-taxon multivariate trait architecture. These quantities are visible image variance and do not represent genetic variance components.

## Within-taxon environmental associations

Four CHELSA v2.1 predictors represented annual mean temperature (BIO1), temperature seasonality (BIO4), annual precipitation (BIO12) and precipitation seasonality (BIO15), using 1981–2010 normals. Environmental extraction is performed once for every observation in the frozen strict spatial cohort without phenotype-based filtering. A registry-driven analysis then evaluates every primary and candidate endpoint with sufficient complete-case support. For each endpoint–predictor pair, taxa represented by fewer than two complete observations are removed before enforcing the predeclared overall observation and taxon minima. Linear outcomes and predictors are demeaned within taxon and standardized after demeaning; ordinary least-squares coefficients are fitted without an intercept with taxon-clustered standard errors. Circular hue contributes one joint sine/cosine test per predictor using within-taxon predictor permutations. BH correction is applied separately to the primary and candidate families.

As a spatial robustness analysis, each primary endpoint was fitted separately with Gaussian SPDE-INLA under four predeclared predictor groups: climate, climate + topography, climate + soil and full environment. Predictors were centred within taxon and standardized; models included a taxon iid random intercept and a Matérn spatial field. Within each endpoint, the four predictor groups used the same complete-case cohort. Posterior two-sided tail-area values were BH-adjusted globally and within model groups, while sign consistency and 95% credible intervals were tracked across groups.

Candidate high-resolution endpoints use the same strict spatial cohort and CHELSA table but retain endpoint-specific resolution, sharpness, segmentation and repeatability missingness. Coverage is reported for every endpoint–predictor pair, including those that do not meet the inferential minimum. Validation-only surface endpoints are excluded from ecological coefficient screening until botanical calibration is passed.

## Among-taxon environmental sorting

For taxa represented by at least five observations, medians of all nine primary traits and available climate, topographic and soil predictors were calculated. Complete taxa were projected into the first three principal components of standardized environmental space, and lower versus upper trait quartiles were compared using centroid distance and Gaussian Bhattacharyya overlap.

To test whether apparent separation could arise from quartile splitting alone, each trait was permuted among the same taxa 10,000 times while environmental PCA coordinates, environmental availability and the trait distribution were held fixed. Benjamini–Hochberg correction was applied across the nine traits separately for centroid-distance and overlap tests. The procedure was repeated for the complete subset of taxa directly represented in the dated backbone.

## Among-taxon climate associations and phylogenetic sensitivity

Taxon-level medians from the balanced atlas were related to the same four standardized CHELSA predictors. Ordinary least-squares models provided uncorrected references. Historical sensitivity was evaluated with Brownian PGLS and maximum-likelihood Pagel-λ PGLS using the dated `GBOTB.extended.LCVP` backbone.

Only 54 of 216 atlas taxa are direct dated-backbone tips. Remaining taxa were therefore placed within genus using alternative `V.PhyloMaker2` scenarios rather than treating one grafted tree as known history. Deterministic scenarios S1 and S3 were compared, and scenario S2 was repeated across 50 randomized grafting trees for the nine primary endpoints. A taxon-level association was treated as robust across tree-placement uncertainty only when direction was stable, the 2.5–97.5% coefficient range excluded zero and nominal support was widespread; we additionally report the fraction of randomized-tree fits passing the within-tree FDR screen. Extended architecture, armature and surface endpoints are not promoted into historical models until the corresponding botanical-validation gate is passed.

## Taxonomic and spatial robustness

The union of 259 source names was audited against the World Checklist of Vascular Plants (WCVP). All eight WCVP synonym candidates that would merge active source units were collapsed simultaneously, after which the nested variance decomposition and all 36 primary within-taxon climate models were recomputed. This was a sensitivity analysis, not a silent replacement of source-platform identities.

Because the compact SPDE workflow had not archived observation-level fitted values, the lowest-WAIC frozen specification for each primary endpoint was refitted diagnostically using the same centering, archived predictor selection, global mesh and PC-Matérn priors solely to export fitted values and residuals. Global residual Moran's I used eight nearest neighbours on three-dimensional unit-sphere coordinates with a fixed neighbour graph and 999 permutations. Taxon trait ranks were also compared after removing each Natural Earth broad region in turn. These diagnostics did not replace accepted coefficient tables.

## Multiplicity, terminology and interpretation limits

FDR families were kept separate by analysis tier, cohort and model family. Primary registry endpoints, candidate registry endpoints, grouped SPDE-INLA posterior screens, taxon-level niche-null tests and tree-based models were not pooled under one multiplicity correction. Descriptive colour composition and validation-only surface signals are not silently added to an inferential FDR family.

Trait–environment coefficients are cross-sectional spatial associations and are not interpreted as phenotypic plasticity, local adaptation, temporal climate response, pollinator-mediated selection, evolutionary rate or adaptive radiation. Circular hue components are never counted as independent flower-colour states.

The former raw absolute-slope lability analysis and median-split quadrants are retained only as statistical provenance/QA after their precision confounding was identified. They are not part of the current submission-facing biological inference.

Independent detector accuracy and continuous-trait biological validity remain genuinely external human/reference-data gates. Production overlays and horizontal-mirror repeatability are technical checks and are not treated as accuracy validation.

