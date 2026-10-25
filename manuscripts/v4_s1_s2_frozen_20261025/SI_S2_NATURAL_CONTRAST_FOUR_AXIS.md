# Supporting Information — S2 native–introduced four-climate comparison

**Analysis status:** real source-data extension conducted on 25 October 2026 after the 24 October S1–S2 four-variable/Holm scientific freeze. The original primary statistical tests are unchanged. This supplement documents the **observational** S2 population, sampling and effect estimands, not adaptive evolution.

## SI-S2.1. Why these two Cirsium species?

The pre-existing V4 candidate pair was *C. arvense* and *C. vulgare*. An all-source-name (350 nominal thistle taxa) WCVP status audit was then run against the original source observations. Seven source names had **both** wild native and explicitly introduced regional records; only these two species had at least 100 usable original angle observations, at least 15 spatially separate 5°/sine-latitude cells, and at least 10 photographers **in both** ranges. The all-name survey was conducted **after initial candidate choice**, so this is a coverage audit, not an outcome-blind preregistered selection of independent invasion treatments.

| Taxon | WCVP-wild native records | WCVP-wild introduced records | Original angle-assessable native | Original angle-assessable introduced |
|---|---:|---:|---:|---:|
| *Cirsium arvense* | 78,646 | 70,964 | 49,072 | 43,195 |
| *Cirsium vulgare* | 62,368 | 85,644 | 33,637 | 47,833 |

These two species account for 156,608/156,750 (99.91%) explicitly introduced wild records in the source. However, the introduced records are strongly North America-weighted (69,138/70,964 and 71,667/85,644, respectively). WCVP status refers to a **regional occurrence category**; it does not identify a specimen's ancestral genotype, arrival date, founding population or independently evolving introduction event. Native S2 photographs partly overlap S1 photographs.

## SI-S2.2. Definitions and competing estimands

The four original environmental axes are CHELSA BIO18 (**mean monthly precipitation in the warmest quarter**, mm/month), CHELSA BIO05 (maximum temperature in the warmest month, °C), source RSDS (shortwave radiation) and WorldClim elevation (m). These are climatological/terrain surfaces at observed coordinates, not contemporaneous weather or independent experimental exposures. Every comparison keeps the four-axis-complete same-site information separate from the photo-angle-qualified subset.

**Observed-record marginal comparison:** compare native versus introduced medians separately for each environment from the available source records (not the same denominator for every predictor when pixels are missing).

**Spatial/photographer-balanced climate comparison:** only observations complete for all four predictors are eligible. First keep one original observation per photographer × 5°/sine-latitude approximately equal-area cell, then average retained observations within each occupied cell. Each occupied cell contributes equally to the species-range mean. Standardise the introduced-minus-native cell mean by the **native cell-mean SD** for each exposure. Resample the observed native and introduced cell means separately 399 times to display descriptive percentile ranges; these cells are **not independent evolutionary replicates**. Cell-varying climate availability is **not** modelled.

**Descriptive four-dimensional multivariate contrast:** compute four-axis native-SD-standardised cell coordinates. Estimate their native-reference covariance and a 5%-trace ridge **in dimensionless units** before calculating a native-reference Mahalanobis centroid distance and an empirical native-95th-distance threshold. These are descriptive measures of observed climate positions only; this is **not** an accessible-environment-corrected claim of niche expansion, stability or unfilling. Multivariate values may be reported only from a passing unit-scaling-invariant source audit; first-run heterogeneous-unit ridge distances are invalid for comparison.

**Photographic head-angle comparison:** select source photographs whose four environmental measurements fall inside a *rectangular intersection of the two ranges' marginal 5th–95th percentile intervals*; this is **not** a full four-dimensional nearest-neighbour match, counterfactual treatment assignment or controlled experiment. Deduplicate by photographer × spatial cell within native/introduced range. For each species, regress EXIF/image-vertical angle on range indicator, the four standardized climates, BIO18 × range and available photograph month (sine/cosine) and year. Weight spatial cells equally within native and introduced groups; use cell-clustered pointwise 95% intervals. Neither a null interaction nor an interval crossing zero is a test of equivalence.

## SI-S2.3. Actual source-real marginal and cell-weighted results

Source execution: [Actions #40164655107](https://github.com/zuizui0223/azami/actions/runs/40164655107), artifact 12587614564, ZIP SHA-256 \`024c8c14e5f0ec4a4a767ed419dbf3f2ef30caaf8156953630c4e04914547934\`. Original frozen S1–S2 member SHAs remained unchanged. Source-script scientific outputs and JSON source IDs are under \`analysis/v4\`.

The four-environment-complete records were 78,060 native/68,397 introduced for *C. arvense*, and 61,800 native/70,653 introduced for *C. vulgare*. Equal-cell climate comparisons used 198/77 and 129/91 distinct occupied native/introduced cells respectively.

| Variable | *C. arvense*: raw median difference | *C. arvense*: cell-equal mean difference in native SD [399-cell range] | *C. vulgare*: raw median difference | *C. vulgare*: cell-equal mean difference in native SD [399-cell range] |
|---|---:|---|---:|---|
| BIO18 | +48.2 mm/month | −0.172 [−0.350, −0.004] | +49.7 mm/month | +0.976 [+0.526, +1.544] |
| BIO05 | +4.5 °C | +0.182 [−0.014, +0.388] | +3.7 °C | +0.495 [+0.194, +0.819] |
| RSDS | +3.740 source units | +0.338 [+0.131, +0.549] | +3.688 source units | +0.578 [+0.334, +0.781] |
| Elevation | +115 m | +0.074 [−0.143, +0.302] | +74 m | +0.256 [−0.022, +0.574] |

All four raw source-location median contrasts are positive. However, the BIO18 **cell-equal mean** contrast reverses in *C. arvense* while remaining strongly positive in *C. vulgare*. **Both the summary statistic and sampling weights change** between the raw median and cell-equal mean estimands; the discrepancy alone cannot be attributed solely to observer clustering. Neither species's historical introduction was an independent rainfall treatment.

## SI-S2.4. Camera-angle contrasts conditional on common climate envelope

| Comparison | *C. arvense* | *C. vulgare* |
|---|---:|---:|
| Eligible native / introduced source images | 438 / 6,069 | 791 / 8,220 |
| Occupied native / introduced cells | 30 / 31 | 25 / 35 |
| Native → introduced image-angle medians (°) | 67.41 → 74.54 | 64.38 → 64.89 |
| Adjusted introduced-range indicator coefficient (°), 95% pointwise CI | −1.03 [−10.38, +8.32] | −2.52 [−8.37, +3.32] |
| BIO18 × introduction interaction (°/pooled predictor SD), 95% pointwise CI | −1.58 [−8.96, +5.80] | −0.02 [−7.30, +7.27] |

The combined S2 128-test **original** Holm family contains 16 camera-angle × four-environment tests (four taxon/range strata × four environmental variables); none passed Holm. The common-envelope model above is a separate post-freeze descriptive sensitivity, not a replacement or multiple-testing-free confirmatory result. The small native photo subgroup and large uncertainty in *C. arvense* make **phenotypic stasis** particularly indefensible. No angle was independently calibrated to gravity, flowering stage and view geometry.

## SI-S2.5. Why this is natural-experiment-like rather than causal

Introduced populations provide two historical within-species comparisons of occupied environmental space, allowing one to ask whether the among-taxon S1 climate–image association appears in a within-species range contrast. It is informative because taxon identity is shared, not because there was random assignment to dry/wet climate. Multiple plausible confounders are uncontrolled: geographic source ancestry, regional available climates, time since establishment, propagule pressure, sampling geography, observer habit, site selection, flowering stage, illumination and camera roll.

Ecological niche shift frameworks explicitly separate **expansion, stability and unfilling relative to available environments** (Guisan et al. 2014, doi:10.1016/j.tree.2014.02.009). Introductions as natural experiments can reveal differentiation, but observational contrasts can spuriously identify adaptation under nonrandom geographic sampling (Colautti & Lau 2015, doi:10.1111/mec.13162). We use these as conceptual boundaries, not as externally validated mechanisms in this dataset.

**GEB manuscript conclusion:** Native taxa show a supported BIO18–photographic-angle association; the two source-rich introduced species occupy different climatic regions, whose apparent wetness contrasts depend on weighting and taxon; and the four-axis-conditioned image-angle comparisons do not demonstrate a common introduced-range effect. This is **scale-dependent observational correspondence**, not proof that thistle capitula could not or did not evolve after introduction.

## Pending validation

- Replace any pre-correction raw-unit 4D ridge outputs with source-reexecuted native-SD-standardised distances and source-certified artifact digest.
- Environmental **availability/background** distributions and direction-agnostic niche novelty/stability/unfilling estimands.
- Independent standardized photographic/gravity measurements and equivalence-margin design.
- Photographic licensing/archival deposition; this appendix does not grant source redistribution.
