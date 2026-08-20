# Supplementary Information

This Supplement is paired with the Chapter 1 submission manuscript **Global citizen-science images reveal continuous within-taxon variation and environmental structure in thistle capitulum traits**. It follows the same frozen cohorts, endpoint definitions, multiplicity families and interpretation limits as `manuscript/final_claims.json` and `manuscript/COHORT_FLOW_AND_ANALYSIS_LEDGER.md`.

## Supplementary Methods

### S1. Cohort construction and non-interchangeability
The balanced image-comparison atlas, exhaustive detector-positive stream, exhaustive spatially thinned primary cohort, grouped SPDE-INLA complete-case cohorts, high-resolution involucre subset and historical-sensitivity taxon summaries are separate executed datasets. Table S1 records their sizes and permitted analyses. FDR counts are never transferred across cohorts.

### S2. Continuous image-trait endpoints
Nine primary endpoints span image-referenced orientation, visible corolla colour and gross two-dimensional outline. Orientation is measured relative to EXIF-oriented image vertical; visible colour is camera-recorded CIELAB lightness/chroma plus circular hue sine/cosine; outline endpoints are aspect ratio, circularity, solidity and width-profile CV. Unassessable measurements remain missing rather than biological absence.

### S3. Nested visible-variance decomposition
For each primary endpoint in the balanced atlas, total image sums of squares are partitioned exactly into among assigned-taxon means, among photographs within assigned taxa, and among heads within photographs. Species/taxa are resampled as clusters for uncertainty. One-head-per-photo and equal-10-photo sensitivities reduce within-photograph multiplicity and unequal taxon replication.

### S4. Exhaustive within-taxon climate coefficients
The primary 46,276-observation cohort uses one observation per source-assigned taxon × 0.25-degree cell after the ≤10 km positional-accuracy filter. Four CHELSA v2.1 predictors are demeaned within taxon and standardized. Nine endpoint components × four predictors yield 36 models with taxon-clustered standard errors and one BH family.

### S5. Grouped SPDE-INLA
Each primary endpoint is analysed separately under four predeclared predictor groups: climate, climate + topography, climate + soil, and full environment. Groups for a given endpoint share the same complete-case cohort. Models include a source-assigned taxon iid random intercept and a Matérn spatial field. Full model summaries, fixed effects, stability summaries, hyperparameters and predictor selections are provided as Tables S5–S6d.

### S6. High-resolution involucral architecture
Detected heads meeting the minimum resolution requirement are subjected to sharpness, segmentation and mirror-repeatability QC. The final inferential auxiliary endpoints are projection roughness, outward spread fraction and maximum spine-like projection relative to head radius. These are two-dimensional contour proxies and are not direct botanical measurements of bract angle, recurvature, spine length or defence.

### S7. Among-taxon environmental niche permutation
Taxa complete for all nine primary endpoints and environmental variables are projected into standardized environmental PCA space. Low- and high-trait quartile groups are compared by centroid distance and Gaussian Bhattacharyya overlap. Trait labels are permuted 10,000 times while environmental positions and the observed trait distribution are fixed. Tests are repeated for direct dated-backbone taxa.

### S8. Historical placement sensitivity
Taxon medians are analysed with OLS and Pagel-lambda PGLS on the dated GBOTB.extended.LCVP backbone. Only 54 of 216 atlas taxa are direct backbone tips; the remainder require within-genus grafting. Fifty randomized scenario-2 grafting trees are therefore used as placement sensitivity, not as a resolved species phylogeny.

### S9. Residual spatial and broad-region robustness
For each primary endpoint, the lowest-WAIC frozen SPDE specification was refitted diagnostically with the same predictor set, centering, mesh and priors only to recover observation-level fitted values and residuals. Residual Moran's I uses eight nearest neighbours on three-dimensional unit-sphere coordinates and 999 permutations. Taxon rankings are compared after omitting each broad Natural-Earth region in turn.

### S10. WCVP taxonomic sensitivity
The union of frozen source names was checked against WCVP through the GBIF checklist interface. Eight synonym conflicts would merge active operational units. All eight are collapsed simultaneously in a sensitivity rerun of nested visible variance and the full 36-model primary climate family.

### S11. Withdrawn lability analysis and precision audit
The former environmental-association score was the RMS magnitude of unpooled absolute taxon-specific slopes. It was strongly associated with slope sample size and standard error. The former negative variation–association relationship and median-split quadrants are therefore withdrawn. A precision-aware reanalysis using noise-adjusted association energy and hierarchical variance meta-regression found no common cross-module coupling. The full archived species-specific coefficient table is supplied only for statistical provenance.

### S12. Independent validation gates
The leakage-free detector audit comprises 1,000 images, including 250 double-labelled images, with zero observation/photo overlap with detector development. Independent human adjudication remains incomplete. Orientation, colour and outline also require genuinely independent reference measurements. Internal overlays and mirror repeatability are technical checks, not biological accuracy validation.

## Supplementary Tables

- **Table S1.** Canonical cohorts and analysis roles.
- **Table S2.** Trait definitions, measurement references and assessability limits.
- **Table S3.** Full nested visible-variance decomposition and replication sensitivities.
- **Table S4.** Complete 36-row exhaustive within-taxon climate coefficient family.
- **Table S5.** Grouped SPDE model-group summary and WAIC.
- **Table S6.** Full grouped SPDE fixed-effect table; S6a lists globally BH-supported rows, S6b effect stability, S6c hyperparameters and S6d predictor selection.
- **Table S7.** Among-taxon niche-permutation summary.
- **Table S8.** Randomized-tree PGLS rows retained in all 50 trees.
- **Table S9.** Residual Moran and leave-one-region-out robustness.
- **Table S10.** WCVP synonym conflicts used in simultaneous-collapse sensitivity.
- **Table S11.** Withdrawn-lability precision-confounding audit; archived taxon-specific coefficients are supplied as a separate machine-readable table.
- **Table S12.** Submission-completion gates.

## Supplementary Figures

- **Figure S1.** Nested visible-variance sensitivities.
- **Figure S2.** Complete exhaustive primary endpoint–climate coefficient map.
- **Figure S3.** SPDE model-group ΔWAIC by endpoint.
- **Figure S4.** Residual Moran's I for the nine primary endpoints.
- **Figure S5.** Leave-one-broad-region-out taxon-rank stability.

## Interpretation boundary

All coefficients are cross-sectional spatial associations. Visible image variance is not genetic variance. Camera-recorded colour is not calibrated reflectance, image vertical is not gravity, and two-dimensional outline/involucre variables remain viewpoint-dependent. Historical analyses are sensitivity tests under incomplete backbone coverage. No result in this Supplement establishes plasticity, local adaptation, selection, pollinator causation, evolutionary rate or adaptive radiation.
