# Scale-dependent climatic associations and partial phenotypic modularity of thistle capitula in public photographs

**Manuscript status:** full first scientific draft, 25 October 2026; **not submitted**.  
**Data lock:** 24 October 2026, V4 S1–S2 final freeze, commit 7e25e409d39a7e8bb8a1ce1043d6bbfd41b7d159.  
**Scope:** S1 native-taxon analysis and S2 native/introduced within-species contrasts only; T2/T3 evolutionary and cross-tribe analyses are outside the inferential scope.  
**Target format:** journal-neutral macroecology/ecology article. Author list, ethics/data deposition details and journal-specific formatting remain to be confirmed.

## Abstract

How floral traits are coordinated across species and environments remains difficult to assess at broad geographic scales, partly because species-level trait summaries discard repeated, geographically explicit observations. We used public photographs of thistle (*Cirsium*) capitula to quantify eight visible traits spanning capitulum presentation, colour and two-dimensional outline, and compared climatic associations and phenotypic covariation across taxonomic and geographic scales. Our frozen V4 study begins with 313,350 WCVP-native source observations (S1) and 297,622 wild native/introduced observations of *C. arvense* and *C. vulgare* (S2). The primary ecological model included mean monthly precipitation of the warmest quarter (BIO18), maximum temperature of the warmest month (BIO05), incident shortwave radiation (RSDS) and elevation. Among 80 S1 and 128 S2 biological endpoint–environment tests, respectively, Holm correction retained four associations: a positive between-taxon association of image-referenced capitulum angle with BIO18 in S1 (standardized β = 0.04545, 95% pointwise CI 0.02235–0.06855, adjusted P = 0.00921), and three colour–environment associations confined to introduced-range S2 strata. Median BIO18 at wild introduced-range source locations exceeded corresponding native-range medians by 48.2 mm month⁻¹ for *C. arvense* and 49.7 mm month⁻¹ for *C. vulgare*, but neither species showed a Holm-supported camera-angle association in the corresponding S2 fits. On the exact common S1 image set (156,683 photographs, 174 taxa), between-taxon RV exceeded within-taxon RV in 24 of 28 trait pairs following environmental adjustment, while median within-module RV exceeded between-module RV at both taxonomic scales. All eight median-thresholded combinations of image angle, chroma and outline aspect ratio occurred in 123 of 135 eligible taxa; however, higher-order contrasts depended on the threshold definition. These results demonstrate a heterogeneous, partly modular organisation of **photographically observed** capitulum traits and distinct associations among taxa versus within introduced/native range strata. They do not identify genetic modules, gravitropic orientation, adaptive rainfall protection or evolutionary responses to introduction.

**Keywords:** *Cirsium*; citizen science; phenotypic integration; floral traits; bioclimatic gradients; intraspecific variation; introduced ranges; image phenomics

## 1. Introduction

Functional trait databases have made large-scale comparisons of plant form and ecology possible, yet a species-level representation cannot by itself preserve the joint distribution of traits within populations, photographs or geographic environments (Kattge et al., 2020). Intraspecific variation may contribute substantially to observed plant-trait distributions (Siefert et al., 2015), but its interpretation requires care: an observed trait range combines biological heterogeneity with differences in observation, sampling and measurement. Floral displays offer a particularly interesting case, because appearance is constructed jointly by colour, spatial presentation and organ geometry rather than by a single scalar trait.

Phenotypic integration describes covariation among traits, whereas phenotypic modularity concerns subsets of traits that are more tightly associated internally than with other subsets (Pigliucci, 2003; Diggle, 2014). A flower can exhibit detectable association among constituent measurements without behaving as one invariant, indivisible phenotypic syndrome. Conversely, differences between within- and among-species covariance need not reveal genetic modularity: they may arise from development, allometry, environmental structure, species delimitation or photographic measurement. Resolving even the **observable** covariance pattern requires repeated phenotypes at a shared observational scale rather than descriptions of the modal condition of each species.

Citizen-science photographs make repeated observations across geographic settings accessible and have been used in automated trait estimation (Schiller et al., 2021). Such pictures are not standardized experiments. Camera position changes apparent orientation and outline; illumination, JPEG processing and white balance change measured colour (Laitly et al., 2021); and photographic effort is uneven among species and locations. Their value therefore lies in carefully defined image-derived estimands, independent data-quality gates and explicitly limited inference, not in treating every image measurement as a direct estimate of anatomical or genetically determined morphology.

Thistles provide numerous capitula differing visibly in presentation, colour and outline. The study asks whether these three domains act as a rigid visible package, or whether distinct combinations and partially separable covariance structure can be recovered in photographs. It also asks whether climate–trait associations differ between a comparison among native taxa and comparisons conducted inside the same species across native and introduced ranges. These are complementary descriptive scales, not independent causal replications. In particular, introduced occurrence environments can differ from native environments because of sampling and colonization history without a corresponding evolutionary or plastic phenotypic change.

We address three linked questions using a **post hoc frozen V4 analysis**, rather than claiming preregistration. First, which of the photographic traits have supported associations with precipitation, temperature, radiation and elevation at the between-taxon and within-range scales? Second, is covariation more concentrated within predefined visible domains—presentation, colour and outline—than across them? Third, do multiple combinations of these domains coexist within the same taxa, inconsistent with a uniquely fixed visible configuration? We report observed associations and covariance, while leaving selection, function, phylogenetic origin and the molecular autonomy of possible trait modules unresolved.

## 2. Materials and methods

### 2.1. Source populations, taxonomy and non-independence

The S1 population comprised 313,350 source observations assigned to native ranges with World Checklist of Vascular Plants (WCVP) status. The independent *analysis question* S2 comprised 297,622 wild source observations drawn from native and explicitly introduced ranges of *C. arvense* and *C. vulgare*. These population totals represent different source definitions and must not be presented as numbers of heads measured for every phenotype. S1 and the native part of S2 can contain the same underlying photographs; they therefore do **not** provide independent confirmation when an estimate has a similar sign.

WCVP range assignments distinguish wild native from explicitly introduced status, without silently reclassifying other statuses as introduced (Govaerts et al., 2021). We retained source photograph/observation provenance and quality eligibility internally. Public reporting uses aggregates; original observer identifiers, precise source coordinates and photographs not licensed for redistribution are not re-hosted in the analytical results. Taxa, photos and detected heads are different levels of the sampling hierarchy.

### 2.2. Visible phenotypes and common-cohort definition

Eight image-derived constructs were analysed: (1) capitulum angle relative to the image or EXIF-corrected camera vertical; (2) colour lightness; (3) chroma; (4) circular hue, represented jointly by its sine and cosine coordinates; (5) capitulum outline aspect ratio; (6) circularity; (7) solidity; and (8) width-profile coefficient of variation. The constructs form three registered observational domains: **presentation** (one), **colour** (three) and **outline** (four). Hue has two computational coordinates but is **one biological construct** for joint testing and the 28 construct-pair comparison.

Image-referenced angle is not a directly calibrated angle to gravitational vertical. Colour estimates are not measures of anthocyanin concentration, and two-dimensional outline measures are sensitive to perspective, occlusion and flowering stage. The principal integration comparison used the exact set of 156,683 photographs across 174 taxa that had all eight constructs and the needed environmental values. This is not interchangeable with the larger angle-specific sample of 168,907 photographs/observations across 178 taxa. The limited involucre (bract) proxies were analysed separately.

### 2.3. Environmental measurements and frozen model

Predictors were sampled at source-pixel resolution from CHELSA v2.1 climatological layers and WorldClim v2.1 elevation. The **single frozen primary model** contained BIO18 (mean monthly precipitation during the warmest quarter), BIO05 (maximum temperature of the warmest month), RSDS (downwelling shortwave radiation) and elevation. BIO18 is not a measure of rainfall during each actual flowering event, of individual storms, or of protection from precipitation. The original source environmental matrix had a maximum variance inflation factor of **4.9559** for this four-variable specification, not proof of causal separability.

The model was chosen after earlier V4 exploration and should not be described as prospectively preregistered or as an independently selected confirmatory model. Alternative four-predictor BIO12 models, VIF-pruned five-/seven-predictor sets and the 11 scenario explorations are preserved as **sensitivity/background evidence**, not pooled with the final Holm family or substituted into main-text conclusions.

### 2.4. Climate–trait association models and multiplicity

S1 separated among-taxon climate means from within-taxon environmental deviations using equal-taxon weighting and taxon-cluster robust uncertainty. Outcome-specific quality and availability determined the exact observations in each fit. For S2, separate estimates were obtained for each of the four combinations of species (*C. arvense*, *C. vulgare*) and status (native, introduced). Within S2 strata the workflow balanced observations by photographer × 5° spatial cell, included identifiable photographed month/year terms and used spatial-cell cluster-robust uncertainty.

All four environmental predictors were included together in the frozen principal association models. For hue, the sine and cosine components were combined in a joint two-degree-of-freedom test, preventing counting colour hue twice. Holm step-down family-wise-error adjustment was performed **separately** over the 80 biological tests in S1 (**10 constructs, including the two SI-only involucre image proxies, × four exposures × two taxonomic scales**) and the 128 tests in S2 (**eight constructs × four exposures × two species-by-range strata pairs**). Reported 95% coefficient intervals are **pointwise and unadjusted**; these must not be misrepresented as familywise simultaneous confidence intervals. Non-significance is not a formal equivalence test.

### 2.5. Within-/between-taxon phenotypic association structure

The 28 unordered construct pairs were estimated on the same photograph-complete S1 cohort. Within-taxon covariation was calculated from taxon-centred phenotypes with equal total taxon weight, whereas between-taxon covariation was based on equally weighted taxon means of the matched observations. For scalar pairs we used the squared correlation as an RV-style association magnitude; hue comparisons used the full two-dimensional circular embedding. Environmental adjustment removed the fitted fixed effects of the frozen four predictors, without removing the species mean itself. This produces conditional descriptive residual associations, **not** estimates of genetic covariance or the fraction of integration caused by environment.

Median pairwise RV was compared within and across the registered domains. For an outcome-blind descriptive reference, the eight trait labels were reassigned exhaustively to groups of sizes 1, 3 and 4 (**280** allocations), and the original within-module minus between-module contrast was ranked among the allocations. These label ranks condition on one measured covariance matrix, whose components partly share pixels and mathematical definitions; they are not phylogenetically independent replications or hypothesis-test P values for heritable modules. Raw RV estimates also have sample-size- and dimension-dependent baselines, which the exploratory equal-sampling sensitivity attempts to diagnose.

### 2.6. Three-way image combinations and supplementary involucral proxies

The three-domain configuration example combined image angle, chroma and outline aspect ratio. One representation dichotomized each construct at its own taxon-level median, enabling enumeration of the eight visual combinations within taxa. A separate higher-order analysis contrasted observed triads with pair-constrained maximum-entropy expectations. A change in binarization rule was explicitly examined, because the direction of such parity or three-high contrasts is not invariant.

Two apparent involucre projection/spreading proxies **were included in the S1 80-test multiple-comparison family, but were interpreted only in Supporting Information**, because they were available for only **1,168 observations and 39 taxa each**, and chroma could be compared with these proxies on **1,007 common photographs and 38 taxa**. The subset-specific four-predictor VIF was **6.14**; those data therefore remain descriptive Supporting Information rather than a basis for primary four-variable mechanistic inference. No direct bract spine length, glandular secretion or defense measurement was available.

### 2.7. Provenance and reproducibility

The fixed analysis is anchored to commit **7e25e409d39a7e8bb8a1ce1043d6bbfd41b7d159** and to the original-data execution at commit **d2c45b00018a0c04e7c69a29857207dc1ac06d0c** (GitHub Actions run **40062126006**; aggregate artifact **12547876704**). The aggregate artifact ZIP SHA-256 is **6e081dc081bb8ee2311a05ac5b7fcddd3f6309be2fa909e87d432789f4e0e205**. The independent frozen-result audit passed in run **40063813578** (artifact **12548870574**). The canonical manuscript facts are transcribed from **analysis/v4/FINAL_V4_BIO18_HOLM_LEGO_FREEZE_20261024_JA.md**, not reconstructed from older V3 or superseded exploratory V4 results.

## 3. Results

### 3.1. Four endpoint–environment associations remain after Holm correction

The four-variable BIO18 model retained **one of 80 S1** and **three of 128 S2** biological associations after Holm correction (Table 1). In S1, the among-taxon standardized coefficient of camera-image capitulum angle on BIO18 was **+0.045447** (95% pointwise CI **+0.022349 to +0.068545**; raw P = **0.00011508**; Holm P = **0.00920635**) using **168,907** angle-eligible records. Thus, wetter warm-quarter climates were associated with larger photographed presentation angles across the analysed native taxa after conditioning on BIO05, RSDS and elevation. This observation alone does not resolve whether the corresponding heads droop more relative to gravity or are exposed to less rain.

In introduced *C. arvense*, colour lightness was negatively associated with elevation (standardized β = **−0.173732**, CI **−0.250914 to −0.096550**; n = **26,171**; Holm P = **0.00131233**). In introduced *C. vulgare*, RSDS showed a positive association with chroma (β = **+0.134504**, CI **+0.066158 to +0.202850**; n = **27,200**; Holm P = **0.01456813**) and a negative association with lightness (β = **−0.142801**, CI **−0.216718 to −0.068884**; the same n; Holm P = **0.01925414**). Chroma and lightness from the same colour-processing pipeline are not independent replications of a radiation mechanism. No common Holm-supported S2 photographic-angle association was established across the two introduced species.

**Table 1. Frozen Holm-supported associations (pointwise coefficient intervals).**

| Stratum and scale | Image-derived trait | Exposure | n | Standardized β | 95% pointwise CI | Raw P | Holm P |
| --- | --- | --- | ---: | ---: | --- | ---: | ---: |
| S1, among native taxa | image-referenced capitulum angle | BIO18 | 168,907 | +0.045447 | +0.022349, +0.068545 | 0.00011508 | 0.00920635 |
| S2, introduced *C. arvense* | lightness | elevation | 26,171 | −0.173732 | −0.250914, −0.096550 | 0.00001025 | 0.00131233 |
| S2, introduced *C. vulgare* | chroma | RSDS | 27,200 | +0.134504 | +0.066158, +0.202850 | 0.00011471 | 0.01456813 |
| S2, introduced *C. vulgare* | lightness | RSDS | 27,200 | −0.142801 | −0.216718, −0.068884 | 0.00015281 | 0.01925414 |

### 3.2. Association geometry differs between and within taxa

The frozen S1 complete-image cohort comprised **156,683 photographs of 174 taxa**. For the **28** construct pairs, between-taxon RV exceeded within-taxon RV in **21/28** pairs before adjustment and **24/28** pairs after adjusting for BIO18, BIO05, RSDS and elevation. These raw pair-direction counts should not be interpreted as 24 statistically independent positive tests.

After four-predictor adjustment, median within-taxon RV was **0.085245** for pairs within the same registered trait domain versus **0.001419** for cross-domain pairs. At the among-taxon scale, corresponding medians were **0.197455** and **0.015453** (Table 2). Among all 280 allocations of the measured traits into one-, three- and four-trait domains, the registered within-minus-between domain contrast was matched or exceeded by **1/280** allocations at the within-taxon scale and **3/280** at the among-taxon scale. This conditional rank demonstrates a distinctive observed grouping of the **measured** axes, but part of the contrast can be produced by common image colour coordinates and derived outline geometry.

**Table 2. Four-variable-adjusted descriptive median RV by pair grouping.**

| Taxonomic scale | Pairs within a registered domain | Pairs across domains |
| --- | ---: | ---: |
| Within taxa | 0.085245 | 0.001419 |
| Among taxa | 0.197455 | 0.015453 |

### 3.3. Multiple three-domain configurations coexist, but higher-order effects are definition-sensitive

Using within-taxon median thresholds for angle, chroma and aspect ratio, **123 of 135 eligible taxa** exhibited **all eight** possible high/low combinations after the four-variable adjustment (124/135 before adjustment). This rules out describing the observed *image configurations* as one universally invariant combination, but does not measure the population frequency or functional viability of all combinations independently of photograph number.

The pair-constrained maximum-entropy reference gave a three-way parity contrast of **+0.03640** (taxon-resampling 2.5–97.5% range **+0.01456 to +0.05762**). When a different threshold and contrast defined the "all three high" configuration, its within-taxon observed frequency was **0.087729** versus a pair-only expectation of **0.091324**, a difference of **−0.003595** (resampling range **−0.007336 to −0.000049**). At the between-taxon scale the analogous difference was **+0.011494**, with a range spanning zero (**−0.004680 to +0.027401**). The sign therefore cannot be presented as a threshold-invariant three-way functional synergy.

### 3.4. Two species occupy wetter warm quarters in introduced ranges, without a common image-angle slope

In the original WCVP-qualified wild source-occurrence populations, the *C. arvense* BIO18 median was **222.2 mm month⁻¹** in its native range (78,646 locations) versus **270.4 mm month⁻¹** in its introduced range (70,964), a **+48.2 mm month⁻¹** difference. For *C. vulgare*, the corresponding medians were **218.0 mm month⁻¹** in the native range (62,368 locations) and **267.7 mm month⁻¹** in the introduced range (85,644), a **+49.7 mm month⁻¹** difference. These are **observed site-climate distribution medians**, not difference-of-slopes estimates, matched-camera image phenotypes, or inferred evolution of precipitation preference. Both introduced populations include relatively dry as well as wet sites, and the raw observation distribution is not a colonization experiment.

Despite these wet-quarter occupancy differences, the frozen four-variable S2 family did not retain a common Holm-supported image-referenced angle–environment association for the two species. This distinguishes **among-taxon environmental sorting** in S1 from **within-species geographic source-range comparisons** in S2, without demonstrating that the observed angle distributions are biologically equivalent across ranges.

After conditioning on the four environmental exposures, median within-domain versus between-domain RV in the four S2 range strata were **0.0487 versus 0.00532** for introduced *C. arvense*, **0.0218 versus 0.0119** for native *C. arvense*, **0.1180 versus 0.00785** for introduced *C. vulgare*, and **0.1102 versus 0.0125** for native *C. vulgare*. These are within-stratum photographic covariances. Differences among them cannot be assigned to introduction-driven evolution without sampling, ancestry and measurement controls, and S2 lacks the corresponding original involucre-proxy columns.

In the limited S1 shared bract–colour cohort, four-predictor-adjusted between-taxon RV for chroma versus an involucre projection proxy was **0.153365**, compared with within-taxon RV of **0.001796**. Because the cohort-specific VIF exceeded five and the measure is a photographic proxy, this is a supporting observational result, not evidence for a defense module.

## 4. Discussion

### 4.1. Partial phenotypic organisation, not an indivisible floral syndrome

The most general result is **heterogeneity of association structure**. Pairs of measurements within colour or outline show larger descriptive RV values than most pairs crossing the predefined presentation, colour and outline domains, while many different three-domain configurations occur in photographs assigned to the same taxon. Together, these observations are compatible with partially modular **visible phenotypes**, rather than a single fixed capitulum configuration. This interpretation does not require—and the evidence cannot establish—freely interchangeable genetic "LEGO" components. A single photographic head can yield correlated colour estimates by construction; multiple shape indices similarly arise from the same segmented silhouette. Shared pixels are an important measurement-level alternative explanation for the high within-domain correlations.

The large 24/28 between-over-within count is also descriptive. Between-taxon taxon means may have much less within-group noise and a different range than individual images, and RV baselines depend on dimensionality and sample size. The fixed-size 280-label reference tests how unusual one *labelled grouping* looks in this observed matrix, not how frequently the same integration architecture evolved independently. The three-way analysis reinforces this boundary: combinatorial richness is visible, but a stable higher-order excess or deficit cannot be claimed when changing the dichotomisation changes the statistic’s direction.

### 4.2. Warm-quarter rainfall corresponds to photographic presentation among taxa

The S1 BIO18 coefficient is the strongest frozen ecological association involving capitulum presentation and survives Holm correction over all 80 S1 biological tests. Its cross-taxon scale is essential. The exposure is a long-term warm-quarter rainfall climatology, not rain deposited on the capitulum while photographed. The measured angle is relative to the camera frame, not gravitational vertical. Therefore, this result proposes a geographic pattern in **image presentation** worth testing against rainfall interception, head architecture, observer viewpoint and phylogenetic structure, but is not itself a demonstration of rainfall-driven selection for nodding.

Earlier V4 exploratory models using BIO12 or additional highly collinear climate predictors sometimes retain a positive coefficient and sometimes lose corrected support. These are useful sensitivity findings, yet they must not be combined with the frozen four-variable BIO18 test as if separate independent experiments. Geographic sorting of taxa, spatial autocorrelation, unequal photo availability and shared ancestry remain alternative interpretations. A field study measuring head orientation to gravity and wetting exposure, with geographically replicated populations and independently recorded photographic roll, would resolve substantially more than an additional reparameterisation of the same images.

### 4.3. Colour associations within introduced ranges are not automatically evolutionary shifts

Three corrected S2 associations involve colour rather than a generalizable direction-angle shift: elevation and image lightness in introduced *C. arvense*, plus radiation and chroma/lightness in introduced *C. vulgare*. They indicate that the environmental association structure of photographic colour can be detectable inside specific introduced-range strata. They do **not** show that introduction caused a colour change, that anthocyanin production was altered or that a common environmentally responsive genotype occurred across the two species.

The contrast with the S1 result is nevertheless informative: introduced-range observation locations have wetter warm-quarter median climates in both focal species, yet a significant among-taxon camera-angle–BIO18 pattern need not manifest as a common, detectable, within-species photo-angle slope in the two introduced cases. A null S2 test does not prove zero slope or equivalence; its confidence intervals, environmental support and measurement precision still matter. Since the S1 native-photo archive overlaps with native S2, simple agreement in native estimates cannot be counted as an external replication.

### 4.4. What public photographs measure—and cannot measure

The gain from public photos is the ability to retain many observations and to compare the joint distributions of continuous visible traits across environmental settings. This complements traits summarized in reference databases, but does not replace standardized field phenotyping. Even correct image segmentation yields outline measures dependent on view angle and head maturity; JPEG colour is not a pigment assay; and an EXIF-corrected frame does not define gravity. The present results should accordingly be read as **image-derived phenotypic associations**.

Future mechanistic tests can be staged without rewriting the present source freeze: (i) independent orientation landmarks including a gravity marker and camera-roll record; (ii) spectral reflectance and chemical pigment assays; (iii) three-dimensional or standardized multi-view capitulum geometry, including true bract/spine structures; (iv) phylogenetic and cytotype-aware comparisons; and (v) pollinator, herbivory, rainfall and reproductive-fitness observations. Such data would permit testing whether the apparent modules are developmentally or genetically separable, and whether any configuration actually improves ecological function. These are next-study questions rather than claims about V4 S1–S2.

### 4.5. Statistical scope and generalisability

The estimates derive from a fixed, post hoc selected four-predictor model. Holm control applies to the declared 80/128 test families **within that model**, not to every environmental specification explored before the model was frozen. Pointwise intervals do not provide simultaneous coverage. Multiple photographs, observers, spatial cells and taxa still have complicated dependencies; clustered standard errors and balancing reduce but do not eliminate these concerns. A future independent sample and photo-to-field validation could change both the estimated coefficients and apparent module contrast.

The biological scope is *Cirsium* source taxa in S1 and exactly **two** introduced species in S2; extending the conclusion to all flowering plants or claiming a repeatable adaptation to introduction is premature. Strongly photographed taxa and regions, taxonomic revisions and camera behaviour could all affect which parts of morphospace enter the study. We explicitly exclude putative molecular-module reuse and evolutionary gain/loss reconstructions from the evidence for this paper.

## 5. Conclusion

Across geographically distributed *Cirsium* photographs, the registered domains of capitulum presentation, visible colour and outline show nonuniform, partly modular **phenotypic association geometry**. Four environment–image-trait associations survive the frozen S1/S2 Holm families: warm-quarter rainfall with among-taxon photographic presentation angle, and elevation/radiation with colour metrics in particular introduced-range strata. The two introduced-range source distributions show wetter warm-quarter medians than their native counterparts without a corresponding universally supported angle slope. Multiple visible trait configurations can co-occur within taxa, but their higher-order interaction is definition-sensitive.

The key ecological contribution is thus a scale-aware account of **visible floral organisation and environmental correspondence**, not proof of functional or genetic modularity. It supplies quantitative hypotheses for independent tests of rainfall exposure, pigmentation, development, ancestry and reproductive performance without converting exploratory image associations into adaptive mechanisms.

## Figure plan (draft captions; panels not yet verified for submission)

**Figure 1. V4 observation and inferential-scale design.** S1 source population (313,350 native observations), endpoint-specific angle set (168,907/178 taxa), all-eight-trait complete set (156,683/174 taxa), S2 source population (297,622 wild native/introduced observations of two species), taxon-between/taxon-within models, native/introduced strata and the no-independent-replication overlap warning. A global sampling map can be added only from license/provenance-safe aggregate location output; do not invent locations.

**Figure 2. Frozen Holm-supported environmental coefficients.** Forest plot of the four Table 1 standardized β estimates with pointwise 95% intervals; symbols split S1 among-taxon and S2 introduced strata, with adjusted P beside each estimate. Include the full association table in Supporting Information to prevent selection-of-only-successes presentation.

**Figure 3. Phenotypic integration across scales.** Exact common S1 cohort: 28 pair RV heatmap at within- and between-taxon scales before and after four-variable adjustment; same-domain and cross-domain links distinguished. Main panel can display Table 2 medians with counts, while the full heatmap should come directly from the frozen aggregate artifact, not from invented pairwise values.

**Figure 4. Three-way visible combinations and threshold sensitivity.** Counts of the eight observed configurations for angle/chroma/aspect by taxon and side-by-side pair-only reference contrasts under the two binarizations. Mark differences as phenotype-only, with no genetic or fitness interpretation.

## References (publicly verifiable starting bibliography)

- Diggle, P. K. (2014). Modularity and intra-floral integration in metameric organisms: plants are more than the sum of their parts. *Philosophical Transactions of the Royal Society B* 369:20130253. https://doi.org/10.1098/rstb.2013.0253
- Govaerts, R., Nic Lughadha, E., Black, N., Turner, R. & Paton, A. (2021). The World Checklist of Vascular Plants, a continuously updated resource for exploring global plant diversity. *Scientific Data* 8:215. https://doi.org/10.1038/s41597-021-00997-6
- Kattge, J. et al. (2020). TRY plant trait database—enhanced coverage and open access. *Global Change Biology* 26:119–188. https://doi.org/10.1111/gcb.14904
- Laitly, A., Callaghan, C. T., Delhey, K. & Cornwell, W. K. (2021). Is color data from citizen science photographs reliable for biodiversity research? *Ecology and Evolution* 11:4071–4083. https://doi.org/10.1002/ece3.7307
- Pigliucci, M. (2003). Phenotypic integration: studying the ecology and evolution of complex phenotypes. *Ecology Letters*. https://doi.org/10.1046/j.1461-0248.2003.00428.x
- Schiller, C., Schmidtlein, S., Boonman, C., Moreno-Martínez, A. & Kattenborn, T. (2021). Deep learning and citizen science enable automated plant trait predictions from photographs. *Scientific Reports* 11:16395. https://doi.org/10.1038/s41598-021-95616-0
- Siefert, A. et al. (2015). A global meta-analysis of the relative extent of intraspecific trait variation in plant communities. *Ecology Letters* 18:1406–1419. https://doi.org/10.1111/ele.12508
- CHELSA climatologies, version 2.1: https://www.chelsa-climate.org/datasets/chelsa_climatologies

## Editorial/source-integrity notes (remove before submission)

1. All four estimates, 80/128 Holm families, module medians and higher-order contrasts come from the **final 24 Oct 2026 freeze**, not an amalgam of older summaries.
2. The S2 site-climate medians (222.2→270.4 and 218.0→267.7 mm month⁻¹) are taken independently from `S1_S2_extreme_LEGO_results.json` in the same verified Actions artifact; they are source occurrence medians and **not** the image-qualified n used in S2 fitted regressions. "Standardized β" uses the source workflow's reported standardization; do not convert these numbers to image degrees or mm without checking exposure/response scaling and the retained coefficients.
3. Some methodological provenance (exact per-endpoint eligibility, original photo processing/weights, climate pixel resolutions, complete sensitivity coefficients, full RV pair heatmaps, uncertainty for integration, and residual geographic diagnostics) needs supplement-level transcription from the already archived pipeline/artifact. These are **reporting gaps**, not permission to invent results.
4. The 280 label allocations and three-way taxon-resampling ranges are **conditional descriptive diagnostics**; journal reviewers should not be told they demonstrate genetic independence or high-order biological interaction.
5. The present first draft is not a submission claim. Author approval, image licensing/attribution, deposition/access terms, exact software versions, figure QA and bibliography verification remain necessary.
