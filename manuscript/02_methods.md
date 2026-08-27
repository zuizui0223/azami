# Methods

## Study design and spatial cohort

We used public biodiversity photographs to quantify continuous capitulum phenotypes across source-assigned taxa of *Cirsium* and allied Cardueae. Source names were retained as operational taxonomic units because the photograph records do not resolve a genus-wide taxonomy. The analysis unit was a coordinate-bearing observation rather than one species-level database row.

The exhaustive source stream contained 406,582 detector-positive observations from 286 taxa. Coordinate filtering retained 392,989 observations from 271 taxa, and a positional-accuracy threshold of 10 km retained 297,293 observations from 259 taxa. Before examining trait values, we retained one observation per taxon by 0.25-degree cell, yielding the frozen spatial cohort of 46,276 observations from 259 taxa. Dates ranged from 1985 to 2026. Missing or unassessable image measurements were never encoded as biological absence or zero.

## From photographs to a 27-endpoint continuous-trait contract

A single-class YOLO11n detector localized visible capitula at a confidence threshold of 0.25. Detection defined image regions of interest only; downstream traits were calculated with deterministic image functions. The category-free contract registered 27 continuous endpoints across orientation, visible corolla colour, gross outline, floral display, involucral architecture, armature-like projection geometry and validation-only surface-image modules.

Orientation was the signed head-axis angle relative to EXIF-oriented image vertical: 0 degrees upward, 90 degrees horizontal and 180 degrees downward. Image vertical is reproducible but is not a gravity measurement. Colour was represented in visible CIELAB space by lightness, chroma and the sine and cosine of circular hue. Outline and contour variables were dimensionless ratios or normalized geometry. Fine involucral, armature-like and surface endpoints were retained as image phenotypes and were not relabelled as botanical spine length, stiffness, hair density, glands or secretion.

The final artifact-backed run measured 22 of 27 endpoints. Eighteen of 19 inferential endpoints were available; `visible_floret_fraction` and four descriptive colour-composition fractions remained explicit `unexecuted_no_measurement` rows. Hue sine and cosine formed one joint circular inferential unit, producing 26 units from the 27 registered endpoints.

## Trait geography and information retained below taxon means

For every endpoint we recorded observation and taxon counts, latitude and longitude limits, and occupancy in 5-degree geographic cells. To quantify information compressed by a one-row-per-taxon database, we decomposed total observation-level sums of squares into variation below and between source-assigned taxon means. Because abundant taxa could dominate the raw partition, we also ran 500 deterministic equal-replication resamples with two observations per eligible taxon. These quantities describe visible image-phenotype variation; they are not heritability, genetic variance or plasticity estimates.

## Six predeclared environmental blocks

Nine predictors formed six biological hypothesis blocks. Thermal conditions comprised annual mean temperature and temperature seasonality (CHELSA BIO1 and BIO4). Hydric conditions comprised annual precipitation and precipitation seasonality (BIO12 and BIO15). Radiative and atmospheric exposure comprised mean shortwave radiation and vapour-pressure deficit. Mechanical exposure was represented by mean near-surface wind. Growing-season water input was represented by climatic growing-season precipitation, and resource/productivity by modelled potential net primary productivity. Growing-season precipitation was not taxon-specific flowering-season rainfall, and potential NPP was not measured local resource supply.

All predictors had to retain at least 98% coverage in the phenotype-blind spatial cohort. The nine-variable environmental table was constructed before endpoint-specific missingness filtering; no unavailable variable could be substituted after observing trait associations.

## Within- and among-taxon environmental associations

Within-taxon models used taxon-demeaned standardized traits and predictors, no intercept, and taxon-clustered standard errors. Eligibility required at least 100 observations, ten taxa and two observations per retained taxon. Joint hue used the magnitude and direction of its standardized sine and cosine coefficients under 9,999 permutations restricted within taxa.

Among-taxon models paired taxon trait medians with taxon environmental medians. The primary scope required at least five trait observations per taxon and at least 20 taxa; a separately corrected sensitivity required two observations. Linear endpoints used standardized univariate slopes and 9,999 taxon-label permutations. Joint hue used paired permutations of the sine and cosine slopes.

Benjamini-Hochberg correction formed one global family across all successful unit by predictor rows separately for within taxa, among-taxon minimum-five and among-taxon minimum-two scopes. Historical endpoint tiers did not split multiplicity or remove rows. Each primary minimum-five row was classified as supported at both scales, within taxa only, among taxa only, neither, or not comparable.

## Sampling-composition sensitivity

Every globally FDR-supported row entered a separately frozen retrospective sensitivity; all atlas rows remained reported regardless of entry. Sampling perturbations were selected without using trait outcomes. Within-taxon rows were refitted with equal total weight per taxon. Both scales omitted each of the ten most observed taxa separately and the two most observed taxa jointly, omitted each of six broad geographic regions separately, and retained only observations explicitly classified as native by the frozen WCVP/TDWG join. Linear rows had to retain coefficient sign; circular rows had to remain within 90 degrees of the original coefficient vector. We reported direction and effect-magnitude stability only and calculated no post-selection P or q values.

## Broad-space and historical-placement sensitivity

Globally FDR-supported rows were refitted with a second-order spherical-coordinate basis. Passage required a spatial permutation P below 0.05, preserved sign for a linear endpoint and residual eight-nearest-neighbour Moran P of at least 0.05. This was a broad geographic sensitivity, not removal of every spatial process.

Among-taxon rows passing the broad-space gate were then fitted by univariate Pagel-lambda PGLS on 52 audited dated-backbone placement trees: two deterministic scenarios and 50 randomized within-genus placements. A final candidate had to retain direction and P below 0.05 on every tree. Because only 54 historical-atlas taxa were direct backbone tips, this procedure is a placement sensitivity rather than a resolved species-tree correction.

## Secondary whole-capitulum synthesis

After endpoint-level inference was complete, a secondary complete-18 analysis asked whether the measured traits reduced to one capitulum syndrome. Standardized taxon medians were analysed by principal components, and within- versus among-taxon association matrices were compared by measurement module. This synthesis did not select endpoints or alter the endpoint-level multiplicity families.

## Claim boundary

The study is a retrospective exploratory analysis of visible image phenotypes and their present spatial associations. It does not identify plasticity, local adaptation, selection, causal environmental mechanisms, botanical identity of validation-only proxies, independent evolutionary origins or convergence. The two final rows are termed adaptive-pattern candidates under current sequential controls, not demonstrations of adaptation.
