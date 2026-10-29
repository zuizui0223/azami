# Supporting Information: original V4 full trait-pair and configuration context

This document retains the original 24 October 2026 S1–S2 source science while
the GEB main article foregrounds independently reported climate–trait results.
All 28 exploratory visible-trait associations, including post hoc cross-domain
pairs, are retained here rather than used as additional main-text discoveries.
No new phenotype estimation, adjustment, correlation test or multiple-testing
family was performed in the editorial refocus.

## Methods for the original descriptive comparisons

**The descriptive covariance comparison is a distinct analytical question applied to existing photographed observations, not an independent source cohort.** The principal analysis uses the identical **156,683 native-source photographs of 174 taxa** with all eight visible constructs; the separate S2 species–range strata supply supplementary descriptive RV medians. We retain all original measurements from three image domains: presentation (one), colour (three) and outline (four). For biological interpretation, we focus on **post hoc, exploratory cross-domain pairs**—angle×chroma, angle×aspect and chroma×aspect—rather than treating correlated numerical descriptions of the same colour or outline as independent functions.

The 28 unordered construct pairs were estimated on the same photograph-complete native-taxon cohort. Within-taxon covariation was calculated from taxon-centred phenotypes with equal total taxon weight, whereas between-taxon covariation was based on equally weighted taxon means of the matched observations. For scalar pairs we used the squared correlation as an RV-style association magnitude; hue comparisons used the full two-dimensional circular embedding. Environmental adjustment removed the fitted fixed effects of the frozen four predictors, without removing the species mean itself. This produces conditional descriptive residual associations, **not** estimates of genetic covariance or the fraction of integration caused by environment.

Median pairwise RV was compared within and across the registered domains. For an outcome-blind descriptive reference, the eight trait labels were reassigned exhaustively to groups of sizes 1, 3 and 4 (**280** allocations), and the original within-module minus between-module contrast was ranked among the allocations. These label ranks condition on one measured covariance matrix, whose components partly share pixels and mathematical definitions; they are not phylogenetically independent replications or hypothesis-test P values for heritable modules (Pigliucci, 2003; Diggle, 2014). Raw RV estimates also have sample-size- and dimension-dependent baselines, which the exploratory equal-sampling sensitivity attempts to diagnose.

The full configuration analysis is detailed in [Supporting Information](SI_V4_FULL_COVARIATION_AND_CONFIGURATION_CONTEXT.md) and combines image angle, chroma and outline aspect ratio. One representation dichotomized each construct at its own taxon-level median, enabling enumeration of the eight visual combinations within taxa. A separate higher-order analysis contrasted observed triads with pair-constrained maximum-entropy expectations. A change in binarization rule was explicitly examined, because the direction of such parity or three-high contrasts is not invariant.

Two apparent involucre projection/spreading proxies **were included in the original native-taxon 80-test multiple-comparison family, but were interpreted only in Supporting Information**, because they were available for only **1,168 observations and 39 taxa each**, and chroma could be compared with these proxies on **1,007 common photographs and 38 taxa**. The subset-specific four-predictor VIF was **6.14**; those data therefore remain descriptive Supporting Information rather than a basis for primary four-variable mechanistic inference. No direct bract spine length, glandular secretion or defense measurement was available.

## Original data and complete tables

- Original source archive: GitHub Actions artifact **12547876704**, ZIP SHA-256
  `6e081dc081bb8ee2311a05ac5b7fcddd3f6309be2fa909e87d432789f4e0e205`.
- Common S1 cohort: **156,683 photographs, 174 taxa**, eight image constructs,
  **28** distinct pairs, with and without BIO18/BIO05/RSDS/elevation adjustment.
- The machine-generated complete paired table is
  `si/Table_S2_all_28_RV_pairs.csv`; the full comparison heatmap is
  `si/Figure_S1_S1_all_28_RV.png` / SVG, rebuilt using
  `build_full_supplement.py` after source verification.
- More intensive full-image level replications, taxon-balanced checks,
  random fixed-size trait-label allocations and the thresholded triplets
  are retained in the original frozen artifact / V4 science record.

## Why the original within-domain comparisons are supplemental

Colour lightness and chroma are different numerical descriptors of a photographed
floral colour, not automatically independent biological adaptations. In the
four-climate-adjusted data their signed r is **−0.547588** among taxa and
**−0.291968** within taxa; this can accompany a light-to-pigmented colour
gradient, as well as image illumination differences.

Outline circularity and solidity are computed from the **same segmentation**.
Their adjusted r is **+0.862672** among taxa and **+0.856496** within taxa.
Their large correlation is not independent evidence that two separate
pollinator-facing head features co-evolved.

The original four-variable-adjusted descriptive median RV values are:

| Scale | Within registered image domain | Across image domains |
|---|---:|---:|
| Within taxa | 0.085245 | 0.001419 |
| Among taxa | 0.197455 | 0.015453 |

The among-taxon RV exceeds the within-taxon RV for **24/28** original adjusted
pairs (21/28 unadjusted). The reference of all **280** allocations of eight
constructs into groups of size 1/3/4 placed the registered grouping at
**1/280** within taxa and **3/280** among taxa. These are descriptive ranks
of a covariance matrix whose variables overlap algebraically, not evidence
of separate genetic developmental modules or an adjusted phenotypic
selection test.

## Exploratory original cross-domain correlations (environment-adjusted)

These three **post hoc** descriptively selected pairs are calculated on the
same **156,683 photos / 174 native taxa**, controlling for the frozen four
environmental predictors. The **sign of the correlation** does not come from
RV (which is nonnegative): we report the original signed *r* values separately.

| Photographic construct pair | Between-taxon signed r | Within-taxon signed r |
| --- | ---: | ---: |
| Image angle × chroma | **−0.225128** | **−0.028067** |
| Image angle × 2D aspect | **−0.246946** | **−0.069387** |
| Chroma × 2D aspect | **+0.019940** | **−0.018247** |

These weak effects are not independent biological tests or evidence for a
functional presentation–colour/outline trade-off, geographic environmental
convergence, pollinator filtering, or genetically separable modules.
An image-based angle and outline aspect can share camera viewpoint artifacts.

## Original three-construct observations

Image angle × chroma × 2D outline aspect were dichotomised and enumerated;
**123/135** eligible taxa occupied all eight measured high/low combinations
after environmental adjustment (**124/135** without adjustment).
Pair-constrained max-entropy analyses returned a parity contrast **+0.03640**
[taxon-bootstrap 0.01456, 0.05762]. The within-taxon high/high/high contrast
under an alternative threshold was **−0.003595**
[−0.007336, −0.000049]; the between-taxon contrast was **+0.011494**
[−0.004680, +0.027401]. As the contrast depends on binarisation and the
probability of finding eight states depends on the sampling opportunity,
these are **not** stable evidence of an adaptive three-way interaction.

## Native and introduced descriptive photo covariation

Within- versus cross-domain adjusted RV medians were:
*C. arvense* introduced **0.0487/0.00532**, native **0.0218/0.0119**;
*C. vulgare* introduced **0.1180/0.00785**, native **0.1102/0.0125**.
These data are from observational strata with different source compositions,
not proven changes to a genetic module after introduction.

## Relationship to GEB V4 main manuscript (2026-10-29)

The GEB Research Article's principal environmental evidence is the
BIO18 × image-referenced capitulum-angle association across native taxa and
three colour × radiation/elevation associations within introduced strata.
**None** of the three selected angle–chroma–aspect cross-domain correlations
is counted as a new statistically supported ecological or evolutionary
finding. The full 28-pair matrix, within-/among-taxon signed-r subset and
old exploratory image (`figure2_cross_functional_trait_pairs.png`) remain
here as supplementary descriptive records. The historical plot filename
does not represent the current main-paper Figure 2, which is reserved for
a not-yet-built source-grounded geographic coverage map.

Original photographed trait covariation was **chosen for biological
interpretation after inspecting the results**; it was not independently
cross-validated against new photograph or taxon samples, was not tested
as genetic modularity, and cannot establish independent selection or
pollinator/predator filtering.


### Alignment with local climate-matched capitulum resemblance

The current main Results 3.3 concerns H3c: whether the joint photographic phenotype is closer at climate-matched cells of geographically separated taxa. It uses its own three-trait eligibility cohort and within-taxon cell-vector permutations. The descriptive associations retained here are not additional evidence for that null-model result, and no original estimates or test families were changed when the main Introduction, Methods and Discussion were aligned with H3c.
