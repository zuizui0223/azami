# Chapter 1 cohort flow and analysis ledger

This ledger replaces ambiguous uses of *strict*, *expanded*, *balanced* and *primary*. Each label below refers to one frozen table and one analysis purpose. Throughout the submission-facing analysis, source-platform taxon assignments are treated as operational taxonomic units; WCVP lumping is evaluated as sensitivity rather than silently rewriting observation identities.

## Exhaustive all-photo inference stream

The exhaustive stream began with all eligible public photographs before trait-based thinning. Cohorts were then derived using coordinate quality, positional accuracy, spatial cells and stable hashes only.

| Canonical cohort name | Frozen file | Observations | Taxa | Selection | Analyses permitted |
|---|---|---:|---:|---|---|
| Exhaustive detected | `all_detected_observations.csv` | 406,582 | 286 | At least one detector-positive capitulum | Detection and availability accounting only |
| Exhaustive coordinate-usable | `coordinate_usable_observations.csv` | 392,989 | 271 | Public usable coordinates | Geographic coverage summaries |
| Exhaustive ≤10 km | `strict_10km_observations.csv` | 297,293 | 259 | Positional accuracy ≤10 km | Pre-thinning sensitivity only |
| **Exhaustive spatially thinned primary** | `strict_spatial_thinned_observations.csv` | **46,276** | **259** | One observation per taxon × 0.25° cell after ≤10 km restriction | Primary within-taxon climate coefficients; source for spatial robustness and focal observation-level handoffs |
| Exhaustive between-taxon balance | `between_species_balanced_40_observations.csv` | 3,723 | 259 | Stable-hash maximum 40 observations per taxon | Exhaustive-stream between-taxon sensitivity only |

The exhaustive merge contained 777,766 photographs from 460,036 source observations, 637,745 detector-positive photographs, 406,582 observations with detected heads and 1,255,791 detected heads. These counts describe the execution stream, not the smaller image-comparison atlas.

## Balanced image-comparison atlas

A separately executed balanced image layer contained 6,626 detected capitula from 3,725 observations and 216 source-assigned operational taxa. It supports nested visible-variance partitioning, taxon-level trait PCA, among-taxon summaries and historical-placement sensitivity.

Its earlier observation-level climate screens are retained as sensitivity analyses:

| Sensitivity analysis | Endpoint rows | Typical usable rows per endpoint | BH-supported main endpoint–predictor rows |
|---|---:|---:|---:|
| Balanced atlas, all coordinate-eligible observations | 36 | 2,834–3,327 | 2 |
| Balanced atlas, positional accuracy ≤10 km | 36 | 2,029–2,388 | 0 |

These balanced-atlas screens must not be called the 46,276-observation primary cohort.

## Primary exhaustive within-taxon coefficients

The 46,276-observation spatially thinned cohort produced 36 component-wise models: nine endpoints × four CHELSA predictors. Eight component rows passed BH correction. Four were ordinary linear endpoints—orientation versus BIO1, chroma versus BIO12, and aspect ratio versus BIO4 and BIO12. Four additional rows were hue sine/cosine components. Because hue components are mathematical parts of one circular response, they are not counted as four independent colour conclusions.

Submission-facing wording is:

> In the exhaustive spatially thinned primary cohort, eight endpoint-component rows passed BH correction. Four were non-circular linear associations; four were hue components requiring joint circular interpretation.

## Grouped SPDE-INLA complete-case cohorts

The grouped spatial models use endpoint-specific complete-case cohorts of 31,666–34,472 observations and 139–141 taxa. For each endpoint, four predeclared predictor groups use the same complete-case cohort:

- climate;
- climate + topography;
- climate + soil;
- full environment.

These models provide spatially explicit robustness and broader predictor-domain sensitivity. Their coefficients are not numerically pooled with the 36-model exhaustive primary family. The strongest stable SPDE support is concentrated in orientation and visible colour; no gross-outline endpoint passes the global SPDE BH screen.

## High-resolution involucre cohort

The high-resolution layer selected 1,819 heads from 1,595 photographs and 214 taxa by detected-head resolution. After sharpness, segmentation and flip-repeatability QC, 1,443 heads from 1,292 observations and 210 taxa were usable.

The final manuscript inferential auxiliary endpoints are:

- involucre projection roughness;
- involucre spread fraction;
- maximum spine-like projection relative to head radius.

The ≤10 km auxiliary environmental family contains 904 observations from 165 taxa and 12 tests (three endpoints × four CHELSA predictors). All three final endpoints are positively associated with BIO4 temperature seasonality after auxiliary-family BH correction.

`spine_peak_count_proxy` remains in integrated derived tables for provenance/downstream exploration but is not one of the three final manuscript auxiliary inferential endpoints.

## Among-taxon environmental sorting

The primary niche-permutation analysis uses 148 taxa complete for all nine primary traits and environmental variables. The direct-backbone sensitivity uses 49 complete taxa. Trait labels are permuted among the same taxa 10,000 times while environmental PCA coordinates and availability are fixed. BH correction is applied separately across nine traits for centroid-distance and overlap metrics within each scope.

This analysis is an among-taxon environmental-sorting test and is not a within-taxon response analysis.

## Historical sensitivity cohort

Historical sensitivity uses the balanced-atlas taxon summaries. Only 54 of 216 atlas taxa are direct dated-backbone tips; remaining taxa require within-genus grafting. Deterministic scenarios and 50 randomized scenario-2 grafting trees are used to evaluate placement sensitivity.

This layer is **historical sensitivity only**. It is not a resolved phylogenetic correction and must not be used for definitive ancestral-state or adaptation claims.

## Precision/lability audit cohort — provenance only

The original 102-taxon raw absolute-slope lability analysis was withdrawn after the score was shown to depend strongly on slope sample size and standard error.

A 101-taxon complete precision-aware cohort was then used to verify that no common cross-module coupling remained after uncertainty was incorporated. It requires:

- all seven linear endpoints;
- all four predeclared CHELSA predictors for every endpoint;
- at least 10 observations for every taxon × endpoint × predictor slope;
- equal weighting of orientation, colour and shape modules;
- explicit use of every archived slope standard error.

This is now **statistical QA/provenance**, not the submission-facing biological headline. Circular hue remains in colour and PCA analyses but is excluded from this precision audit because archived joint hue vectors do not contain component standard errors.

## Independent detector audit cohort

The leakage-free detector audit contains 1,000 source images across 323 species and 85 ten-degree spatial blocks, with 250 images assigned to a second annotator. Photo and observation overlap with detector-development data is zero. Human adjudication remains incomplete, so no independent precision/recall/F1 is currently claimed.

## Naming rules

- **Image-comparison atlas** = 6,626 heads, 3,725 observations, 216 source-assigned taxa.
- **Exhaustive detected layer** = 406,582 observations, 286 taxa.
- **Exhaustive spatially thinned primary cohort** = 46,276 observations, 259 taxa.
- **Grouped SPDE complete-case cohort** = endpoint-specific 31,666–34,472 observations, 139–141 taxa.
- **High-resolution involucre cohort** = 1,443 usable heads, 1,292 observations, 210 taxa; ≤10 km environmental family = 904 observations, 165 taxa.
- **Precision-aware lability cohort** = 101 fully complete taxa, provenance/QA only.
- Do not use *expanded pooled cohort* without a filename and count.
- Do not report an FDR count without the cohort, endpoint family and number of tests.
- Do not call source-assigned operational units a resolved genus-wide species taxonomy.
