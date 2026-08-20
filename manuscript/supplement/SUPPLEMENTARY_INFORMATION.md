# Supporting Information

## Submission scope
This supplement follows the current continuous within-taxon image-phenomics manuscript. The withdrawn raw lability/quadrant result is retained only as statistical provenance/QA.

## Methods sections

### S1. Frozen cohort hierarchy and analysis permissions
The submission uses separate executed cohorts for separate inferential scales. The balanced image-comparison atlas is the head/photo/taxon layer for visible-variance decomposition and taxon-level summaries. The exhaustive post-detection stream is filtered by coordinate availability, positional accuracy and taxon × 0.25° thinning before the primary within-taxon climate analysis. Grouped SPDE-INLA uses endpoint-specific complete cases. The high-resolution involucre subset is separate because assessability depends on crop resolution. Cohorts are never merged under a single FDR family.

### S2. Operational taxonomic units and provenance
Public-platform names are retained as source-assigned operational units so that observation identity remains exactly reproducible. This is not a claim that the analysis has resolved genus-wide species limits. WCVP names are used for an authority-backed audit, and all eight synonym candidates that would merge active source units are collapsed simultaneously in a sensitivity analysis.

### S3. Capitulum localization and continuous trait extraction
A frozen single-class YOLO11n detector localizes visible capitula only. Downstream inference uses deterministic numerical image measurements rather than the legacy upward/nodding classifier. Primary endpoints are orientation relative to EXIF-oriented image vertical; CIELAB lightness and chroma; circular hue sine/cosine; and four two-dimensional outline statistics. Failed or unstable measurements remain missing rather than being coded as biological absence.

### S4. Independent detector and measurement validation
The leakage-free detector audit contains 1,000 source images sampled only after excluding every photo_id and obs_id used in detector development; 250 images are assigned for independent double annotation. Human adjudication is still pending, so detector precision/recall/F1 are not reportable. Independent biological/reference measurements for orientation, colour and outline are likewise a remaining external gate. Production overlays and horizontal-mirror repeatability are technical checks, not independent accuracy validation.

### S5. Nested visible-variance decomposition
For each primary endpoint, visible image sums of squares are decomposed into among assigned-taxon means, among photographs within assigned taxa, and among heads within photographs. Taxa are cluster-resampled 2,000 times for uncertainty. Two sensitivities reduce pseudoreplication and unequal sampling: one deterministic head per photo, and 500 balanced repeats of 10 photos per eligible taxon.

### S6. Primary within-taxon climate coefficients
Four CHELSA v2.1 1981–2010 normals (BIO1, BIO4, BIO12 and BIO15) are analysed after within-taxon demeaning. Outcomes and predictors are standardized after demeaning, and no-intercept OLS coefficients use taxon-clustered standard errors. BH correction is applied to the 36 endpoint-component tests. Hue sine and cosine are computational components of one circular response and must be interpreted jointly.

### S7. Grouped SPDE-INLA models
Each primary endpoint is fitted separately under four predeclared predictor groups: climate, climate + topography, climate + soil and full environment. For a given endpoint the groups share the same complete-case cohort. Models include a taxon iid random intercept and Matérn spatial field. Posterior two-sided tail-area summaries are BH-adjusted; direction and 95% credible intervals are tracked across groups. Coefficients are not pooled numerically with the primary OLS family.

### S8. High-resolution involucral architecture
Heads with minimum detected dimension ≥150 pixels enter the high-resolution route, followed by sharpness, segmentation and mirror-repeatability QC. Three final inferential proxies are retained: projection roughness, outward spread fraction and maximum spine-like projection relative to head radius. They are two-dimensional contour proxies and cannot distinguish bract length, spreading and recurvature as botanical mechanisms.

### S9. Among-taxon environmental sorting
For taxa with adequate observations, taxon medians are projected into the first three axes of standardized environmental PCA space. Lower and upper trait quartiles are compared using centroid distance and Gaussian Bhattacharyya overlap. Trait labels are permuted 10,000 times among the same taxa while environmental coordinates and availability remain fixed. BH correction is applied separately by metric and scope.

### S10. Historical placement sensitivity
Taxon-level medians are related to four CHELSA predictors using OLS and Brownian/Pagel-lambda PGLS on a dated synthetic backbone. Only 54 of 216 atlas taxa are direct dated-backbone tips; the remainder require within-genus grafting. Deterministic scenarios and 50 randomized scenario-2 trees quantify placement sensitivity. This layer is not a resolved Cirsium phylogeny and is not used for ancestral-state inference.

### S11. Spatial and taxonomic robustness
For each primary endpoint, the lowest-WAIC frozen SPDE specification is refitted diagnostically only to recover fitted values/residuals. Residual Moran I uses an eight-nearest-neighbour graph on unit-sphere coordinates with 999 permutations. Taxon trait rankings are also recomputed after omitting each broad geographic region. Separately, all eight WCVP synonym conflicts are collapsed simultaneously and the headline variance and 36-model climate families are recomputed.

### S12. Withdrawn lability result and current claim boundary
The former species environmental-responsiveness score was the RMS magnitude of unpooled absolute taxon-specific slopes and was strongly confounded with slope sample size and standard error. The negative variation–responsiveness relation and median-split quadrants are therefore withdrawn. A noise-adjusted association-energy summary and hierarchical variance meta-regression detect no common cross-module coupling. These analyses remain statistical QA/provenance and are not part of the biological headline. All current trait–environment results are cross-sectional spatial associations, not demonstrated plasticity, adaptation, pollinator selection or evolutionary rate.

## Tables and figures
The summary tables currently committed in `tables/` are frozen to the current claim registry. High-volume coefficient/model tables and Figures S1–S5 remain release-generation products until the submission branch passes the cleanup regression checks and the external measurement gates close.
