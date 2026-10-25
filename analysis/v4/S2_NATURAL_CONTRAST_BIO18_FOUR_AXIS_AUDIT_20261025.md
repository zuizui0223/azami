# S2 natural experiment-like contrasts: **source-real four-axis reanalysis complete** (25 October 2026)

## Current adopted interpretation

**Both species have higher introduced-range observed-record BIO18 medians, but they do NOT share an equivalent wetward climatic shift after spatial-cell balancing.** In *C. arvense*, the cell-weighted mean BIO18 contrast reverses direction; in *C. vulgare* it remains strongly positive. Neither has a resolved introduction-associated camera-angle effect or BIO18 × introduction interaction in the restricted four-climate support. This is a novel **observation-weight versus spatial-scale dependence**, not proof of evolution or absent phenotypic change.

Scientific parent: original 24 October frozen V4 S1/S2 four-exposure/BIO18 Holm results are **not overwritten**.

**Corrected source-real execution:** [Actions #40166532404](https://github.com/zuizui0223/azami/actions/runs/40166532404) **SUCCESS**, artifact **12588158569**, ZIP SHA-256 `1fcb9ff110f2c2d6bba6afc58b45d9a38f86cf919f989f354f82054a80cec121`, corrected JSON `V4_S2_FOUR_AXIS_JOINT_NICHE_REAL_SOURCE.json` SHA-256 `800f2f5dc37d7156324a4756f372306c3f8b28cb00c88d683e86cd01c2c64063`. The earlier [Actions #40164655107](https://github.com/zuizui0223/azami/actions/runs/40164655107) was scientifically valid for marginal climate and conditional camera-angle estimates but computed its 4D ridge on raw heterogeneous physical units; its distance and outside-threshold fraction are **superseded**. Recomputing the original V4 source aggregates changed JSON bytes at floating-point roundoff level only (max absolute numeric deviation `9.93e-11` across compared members; no structural/nominal data differences). The **original frozen** October 24 input/member digest identities remain `bb01d3b5...` and `3738292e...`; do not falsely call rerun JSON byte-identical.

## Why these two species?

The original **350 species-name** screen found **seven** *Cirsium* taxa with wild WCVP-assigned native **and** explicit introduced regions. Only *C. arvense* and *C. vulgare* met the fixed stronger per-range support gates in both domains (≥100 orientation-eligible observations, ≥15 independent 5°/equal-sine-latitude geographic cells, ≥10 observers). The two candidates were initially specified prior to the complete all-350 check; the later screen confirms there were no other qualifying taxa. This is *not* a random draw of introductions.

*C. arvense* original source totals: **78,646 native / 70,964 introduced**; *C. vulgare*: **62,368 / 85,644**. Introduced observations are dominated by North America (**69,138 of 70,964** for *C. arvense*; **71,667 of 85,644** for *C. vulgare*), limiting geographic independence.

## Same four physical environmental variables, two valid units of sampling

The exact science-frozen environmental family is BIO18 (warmest-quarter mean monthly precipitation), BIO05 (warmest-month maximum temperature), RSDS (source-scaled downward shortwave radiation) and source-pixel elevation. All four are extracted on exactly matched original WCVP-wild observations.

One estimator computes **observed-record medians** (photograph/observer effort strongly nonuniform). The independent source-row reanalysis first deduplicates photographer × broad spatial cell and **equally weights occupied cells**, comparing four-coordinate cell means standardized by the native-range cell SD of each predictor. 399 within-range cell resamplings yield descriptive 2.5–97.5% intervals. The standardized 4D centroid distance first scales every exposure by the **native-cell SD**, then applies native-reference covariance with a 5%-trace ridge in the resulting dimensionless coordinates. The share of introduced cell means outside the native 95th-percentile Mahalanobis ellipsoid diagnoses source climate novelty **without climate-availability correction**, not a formal niche expansion/unfilling estimate.

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

After scaling **each of BIO18, BIO05, RSDS and elevation by the native-range cell SD before covariance regularization**, the native-reference 4D Mahalanobis centroid distances are **0.448** (*C. arvense*) and **1.791** (*C. vulgare*). Outside-native-95th-distance fractions are **1/77 (1.3%)** and **23/91 (25.3%)**, respectively. The native-SD-standardized source code passes an explicit unit-rescaling invariance regression test. These are conditional descriptive climatic **occurrence** measures, not environmental-availability-corrected niche expansion: **25.3% is not a percentage of new adaptive ecological niche**. The first-run heterogeneous-unit distances of 0.182 and 0.913 must not be reused.

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
