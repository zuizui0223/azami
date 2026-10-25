# S2 natural experiment-like contrasts: **source-real four-axis reanalysis complete** (25 October 2026)

## Current adopted interpretation

**Both species have higher introduced-range observed-record BIO18 medians, but they do NOT share an equivalent wetward climatic shift after spatial-cell balancing.** In *C. arvense*, the cell-weighted mean BIO18 contrast reverses direction; in *C. vulgare* it remains strongly positive. Neither has a resolved introduction-associated camera-angle effect or BIO18 × introduction interaction in the restricted four-climate support. This is a novel **observation-weight versus spatial-scale dependence**, not proof of evolution or absent phenotypic change.

Scientific parent: original 24 October frozen V4 S1/S2 four-exposure/BIO18 Holm results are **not overwritten**.

**Source-real successful execution:** [Actions run 40164655107](https://github.com/zuizui0223/azami/actions/runs/40164655107), artifact **12587614564**, ZIP SHA-256 `024c8c14e5f0ec4a4a767ed419dbf3f2ef30caaf8156953630c4e04914547934`. The computed 4D source JSON is `V4_S2_FOUR_AXIS_JOINT_NICHE_REAL_SOURCE.json` with SHA-256 `b63662536bf26309622d3fed4a0602891658b1afcd8c00770443b516bb7bbe6e`. The original V4 source-aggregate member SHA-256 values `bb01d3b5da03148d84b2521c9c92cf14f7f6177435e145118f6ff6df8846733f` and `3738292e0a64da183dfe219c854c16cda199a7f9769f4928da8ab07657d60b70` remain unchanged in the new success artifact.

## Why these two species?

The original **350 species-name** screen found **seven** *Cirsium* taxa with wild WCVP-assigned native **and** explicit introduced regions. Only *C. arvense* and *C. vulgare* met the fixed stronger per-range support gates in both domains (≥100 orientation-eligible observations, ≥15 independent 5°/equal-sine-latitude geographic cells, ≥10 observers). The two candidates were initially specified prior to the complete all-350 check; the later screen confirms there were no other qualifying taxa. This is *not* a random draw of introductions.

*C. arvense* original source totals: **78,646 native / 70,964 introduced**; *C. vulgare*: **62,368 / 85,644**. Introduced observations are dominated by North America (**69,138 of 70,964** for *C. arvense*; **71,667 of 85,644** for *C. vulgare*), limiting geographic independence.

## Same four physical environmental variables, two valid units of sampling

The exact science-frozen environmental family is BIO18 (warmest-quarter mean monthly precipitation), BIO05 (warmest-month maximum temperature), RSDS (source-scaled downward shortwave radiation) and source-pixel elevation. All four are extracted on exactly matched original WCVP-wild observations.

One estimator computes **observed-record medians** (photograph/observer effort strongly nonuniform). The independent source-row reanalysis first deduplicates photographer × broad spatial cell and **equally weights occupied cells**, comparing four-coordinate cell means standardized by the native-range cell SD of each predictor. 399 within-range cell resamplings yield descriptive 2.5–97.5% intervals. The standardized 4D centroid distance uses a native covariance matrix with 5%-trace ridge regularization. The share of introduced cell means outside the native 95th-percentile Mahalanobis ellipsoid diagnoses source climate novelty **without climate-availability correction**, not a formal niche expansion/unfilling estimate.

| Species | Four-exposure complete native/introduced records | Equal-weight native/introduced cells |
|---|---:|---:|
| *C. arvense* | 78,060 / 68,397 | 198 / 77 |
| *C. vulgare* | 61,800 / 70,653 | 129 / 91 |

## Genuine four-axis site results

| Explanatory axis | *C. arvense*, source median introduced−native | *C. arvense*, cell-balanced native-SD difference [bootstrap range] | *C. vulgare*, source median introduced−native | *C. vulgare*, cell-balanced native-SD difference [bootstrap range] |
|---|---:|---|---:|---|
| BIO18 | **+48.2 mm/month** | **−0.172 [−0.350, −0.004]** | **+49.7 mm/month** | **+0.976 [+0.526, +1.544]** |
| BIO05 | +4.5 °C | +0.182 [−0.014, +0.388] | +3.7 °C | +0.495 [+0.194, +0.819] |
| RSDS | +3.740 source units | +0.338 [+0.131, +0.549] | +3.688 | +0.578 [+0.334, +0.781] |
| Elevation | +115 m | +0.074 [−0.143, +0.302] | +74 m | +0.256 [−0.022, +0.574] |

The native-reference 4D Mahalanobis centroid distances are **0.182** (*arvense*) and **0.913** (*vulgare*). Native-ellipsoid out-of-support introduced-cell fractions are **0/77 (0%)** and **17/91 (18.7%)**, respectively, under this empirical covariance threshold only. **Do not interpret the latter as 18.7% ecological niche expansion**: available climates, spatial region confounding and establishment history are not included.

## Photo angle within four-axis overlapping climate conditions

For each species, the overlap is a four-dimensional intersection of the 5th–95th percentile **marginal bounds** for all four predictors. Source photo observations are filtered and deduplicated by observer × cell. Each fitted within-species WLS includes a range indicator, BIO18 × introduction interaction, BIO18/BIO05/RSDS/elevation main terms, month sine/cosine and year, equal range×cell weights, and cell-clustered uncertainty. This is a common-climate sensitivity, not exact pair matching or a causal estimand.

| Metric | *C. arvense* | *C. vulgare* |
|---|---:|---:|
| Valid source angles, native/introduced in overlap | **438 / 6,069** | **791 / 8,220** |
| Distinct spatial cells, native/introduced | 30 / 31 | 25 / 35 |
| Observed overlap-subset image angle medians (°), native / introduced | 67.41 / 74.54 | 64.38 / 64.89 |
| **Adjusted introduction status coefficient (°)** | **−1.03** [−10.38, +8.32] | **−2.52** [−8.37, +3.32] |
| **BIO18 × introduction (°/pooled climate SD)** | **−1.58** [−8.96, +5.80] | **−0.02** [−7.30, +7.27] |

**Both adjusted confidence intervals cross zero, but neither is a test of shape equivalence.** The *C. arvense* raw median angle difference within four-axis overlap is not trivial, and its native eligible cohort is much smaller than its introduced cohort; this rules out concluding “no shape change” or evolutionary stasis from a non-significant contrast. The older source V4 128-test Holm family independently retained 0/16 head-angle effects across four range strata; all original positive colour results remain untouched.

## GEB ecological contribution and interpretation limits

The S1 between-taxon BIO18/photo-angle association, S2 native/introduced source-domain climatic occupancy, and S2 within-species/common-climate angle contrasts are **three different estimands**. Their disagreement identifies scale- and sampling-dependent **apparent phenotype–climate coupling**, not a demonstrated mechanism.

Biological introductions are useful **natural experiment-like opportunities**, but not randomized. Source cell coverage is sparse relative to the number of photographs, climatic availability differs by region, the native and introduced sampled genotypes are not genealogically paired, and original photographs do not validate gravitational nodding or anthocyanin physiology. An ellipsoid threshold, an axis-aligned common support rectangle and 399-cell bootstrap do **not** resolve selection or true niche evolution. For background-corrected niche expansion/stability/unfilling, see Guisan et al. (2014, doi:10.1016/j.tree.2014.02.009); for invasion comparisons and nonrandom geography, see Colautti & Lau (2015, doi:10.1111/mec.13162).

**Freeze policy:** This 25 October source-real comparison is a distinct, reproducible **V4 S2 extension**; do not change any 24 October V4 S1/S2 Holm estimate or reinterpret a null as a causal nonresponse.
