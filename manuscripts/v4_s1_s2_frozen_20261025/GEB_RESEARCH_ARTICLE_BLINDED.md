# Scale-dependent climatic associations and partial phenotypic integration in globally photographed thistles

**Running title:** Thistle traits across scales

## Abstract

**Aim:** Determine how climatic associations of photographic capitulum presentation, colour and outline vary between native taxa and within species-by-range strata, and test whether the visible floral phenotype behaves as a single fixed configuration or comprises partly separable trait domains.

**Location:** Worldwide georeferenced photographs of native *Cirsium* taxa and wild native/introduced occurrences of *C. arvense* and *C. vulgare*.

**Time period:** CHELSA 1981–2010 climatological normals; photograph dates are heterogeneous and require a final source-date audit.

**Major taxa studied:** *Cirsium* (Asteraceae).

**Methods:** On frozen source cohorts of 313,350 native observations (S1) and 297,622 wild native/introduced observations in two species (S2), we quantified eight image-derived traits. Climate–trait models jointly included warm-quarter precipitation (BIO18), maximum warm-month temperature (BIO05), radiation and elevation; 80 S1 and 128 S2 tests were separately Holm-adjusted. We compared 28 matched visible-trait pairs at within- and among-taxon scales. S2 represents stratified phenotype–environment analysis, **not a comparison of climatic niche occupancy**.

**Results:** Four associations survived Holm correction: a positive among-taxon image-frame angle–BIO18 coefficient (standardized β=+0.0454, Holm P=0.0092) and three introduced-range colour associations with elevation or radiation. In the common S1 set (156,683 photographs/174 taxa), 24/28 environmentally adjusted trait pairs had larger between- than within-taxon RV, with larger same-domain than cross-domain medians at both scales. All eight binary configurations of presentation, chroma and outline were observed in 123/135 eligible taxa; their higher-order contrasts were threshold-sensitive. Across S2 strata, none of 16 angle–environment terms passed Holm.

**Main conclusions:** Globally photographed thistle capitula show scale-dependent environmental associations and partly distinct visible phenotypic domains. These observational patterns establish neither functional or genetic modules nor adaptive rainfall-related evolution; climatic niche expansion or native–introduced rainfall-niche differences were not tested.

**Keywords:** biogeography, capitulum, citizen science, climate, floral colour, introduced ranges, phenotypic integration, trait covariance

## 1. Introduction

Functional trait databases have made large-scale comparisons of plant form and ecology possible, yet a species-level representation cannot by itself preserve the joint distribution of traits within populations, photographs or geographic environments (Kattge et al., 2020). Intraspecific variation may contribute substantially to observed plant-trait distributions (Siefert et al., 2015), but its interpretation requires care: an observed trait range combines biological heterogeneity with differences in observation, sampling and measurement. Floral displays offer a particularly interesting case, because appearance is constructed jointly by colour, spatial presentation and organ geometry rather than by a single scalar trait.

Phenotypic integration describes covariation among traits, whereas phenotypic modularity concerns subsets of traits that are more tightly associated internally than with other subsets (Pigliucci, 2003; Diggle, 2014). A flower can exhibit detectable association among constituent measurements without behaving as one invariant, indivisible phenotypic syndrome. Conversely, differences between within- and among-species covariance need not reveal genetic modularity: they may arise from development, allometry, environmental structure, species delimitation or photographic measurement. Resolving even the **observable** covariance pattern requires repeated phenotypes at a shared observational scale rather than descriptions of the modal condition of each species.

Citizen-science photographs make repeated observations across geographic settings accessible and have been used in automated trait estimation (Schiller et al., 2021). Such pictures are not standardized experiments. Camera position changes apparent orientation and outline; illumination, JPEG processing and white balance change measured colour (Laitly et al., 2021); and photographic effort is uneven among species and locations. Their value therefore lies in carefully defined image-derived estimands, independent data-quality gates and explicitly limited inference, not in treating every image measurement as a direct estimate of anatomical or genetically determined morphology.

Thistles provide numerous capitula differing visibly in presentation, colour and outline. We ask whether these domains form one rigid visible package or show distinct combinations and partly separable covariance structures in photographs. We also compare **the environmental associations of photographed traits** at two observational levels: among native taxa (S1), and within the native/introduced strata of two focal species (S2). The S2 comparison is not an experiment, and this manuscript does **not** compare the climatic niches or environmental occupancy of native and introduced ranges.

We address three linked questions using a **post hoc frozen V4 analysis**, rather than claiming preregistration. First, which of the photographic traits have supported associations with precipitation, temperature, radiation and elevation at the between-taxon and within-range scales? Second, is covariation more concentrated within predefined visible domains—presentation, colour and outline—than across them? Third, do multiple combinations of these domains coexist within the same taxa, inconsistent with a uniquely fixed visible configuration? We report observed associations and covariance, while leaving selection, function, phylogenetic origin and the molecular autonomy of possible trait modules unresolved.

## 2. Materials and methods

### 2.1. Source populations, taxonomy and non-independence

The S1 source population comprised **313,350** WCVP-assigned native observations. S2 comprised **297,622** wild observations of *Cirsium arvense* and *C. vulgare*, classified separately by WCVP native or explicitly introduced regional status. S2 uses those four **species × range-status strata** to estimate photographed phenotype–environment associations, not to infer climatic-niche shifts or changes caused by introductions. Both focal species were initially chosen for extensive, measureable records across native and introduced ranges. A subsequent screen of 350 original species names found seven species with both wild range categories; only the original two focal species met the same per-range measurement and sampling coverage gates (≥100 image-angle-eligible records, ≥15 broad spatial cells and ≥10 observers). The complete screen **followed** their initial selection and is a coverage validation, not fully prospective random selection. Together the two species account for **156,608 of 156,750** explicitly introduced wild source observations (99.91%); introduced records are strongly geographically concentrated, particularly in North America. They therefore represent a narrow observational case set rather than independent replicated invasions.

These population totals represent different source definitions and must not be presented as numbers of heads measured for every phenotype. S1 and the native part of S2 can contain the same underlying photographs; they therefore do **not** provide independent confirmation when an estimate has a similar sign.

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

### 2.5. S2 species-by-range phenotype strata (no climatic-niche comparison)

To retain original source provenance without treating introduction as a climatic treatment, each focal species was analysed in **separate native and introduced photographic strata** under the same frozen four-variable BIO18/BIO05/RSDS/elevation model. The model conditions the image-derived trait on the local environmental predictors **within that stratum**. It does not estimate native-versus-introduced differences in climate occupancy, multivariate climatic-niche expansion, introduced-minus-native niche centroid shifts, or range-status × BIO18 effects.

Range labels are coarse WCVP regional statuses, not records of founding population ancestry, introduction timing or genetic changes. Spatial cells and photographic observers remain unevenly sampled. Comparing S1 among-taxon coefficients to S2 within-stratum coefficients highlights possible **differences in observational scale**; it does not constitute a randomized natural experiment, a paired genotype design or a formal test that the coefficients differ between scales.

### 2.6. Within-/between-taxon phenotypic association structure

The 28 unordered construct pairs were estimated on the same photograph-complete S1 cohort. Within-taxon covariation was calculated from taxon-centred phenotypes with equal total taxon weight, whereas between-taxon covariation was based on equally weighted taxon means of the matched observations. For scalar pairs we used the squared correlation as an RV-style association magnitude; hue comparisons used the full two-dimensional circular embedding. Environmental adjustment removed the fitted fixed effects of the frozen four predictors, without removing the species mean itself. This produces conditional descriptive residual associations, **not** estimates of genetic covariance or the fraction of integration caused by environment.

Median pairwise RV was compared within and across the registered domains. For an outcome-blind descriptive reference, the eight trait labels were reassigned exhaustively to groups of sizes 1, 3 and 4 (**280** allocations), and the original within-module minus between-module contrast was ranked among the allocations. These label ranks condition on one measured covariance matrix, whose components partly share pixels and mathematical definitions; they are not phylogenetically independent replications or hypothesis-test P values for heritable modules. Raw RV estimates also have sample-size- and dimension-dependent baselines, which the exploratory equal-sampling sensitivity attempts to diagnose.

### 2.7. Three-way image combinations and supplementary involucral proxies

The three-domain configuration example combined image angle, chroma and outline aspect ratio. One representation dichotomized each construct at its own taxon-level median, enabling enumeration of the eight visual combinations within taxa. A separate higher-order analysis contrasted observed triads with pair-constrained maximum-entropy expectations. A change in binarization rule was explicitly examined, because the direction of such parity or three-high contrasts is not invariant.

Two apparent involucre projection/spreading proxies **were included in the S1 80-test multiple-comparison family, but were interpreted only in Supporting Information**, because they were available for only **1,168 observations and 39 taxa each**, and chroma could be compared with these proxies on **1,007 common photographs and 38 taxa**. The subset-specific four-predictor VIF was **6.14**; those data therefore remain descriptive Supporting Information rather than a basis for primary four-variable mechanistic inference. No direct bract spine length, glandular secretion or defense measurement was available.

### 2.8. Provenance and reproducibility

The frozen observational statistics were obtained from provenance-verified source observations, a fixed source-pixel climatic extraction and executable model/test-family specifications, with full corrected test results and RV matrices in the separate statistical supplement. The analysis contains source-level identifiers and third-party photographs that cannot all be openly redistributed. For submission, reviewer access to releasable source tables and code must be made identity-concealed and durable; neither access nor an archived V4 DOI is represented as complete in this draft.

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

### 3.4. S2 phenotype–environment patterns in four range-status strata

The original frozen S2 family tests eight observable constructs against the four environmental predictors within each of the **four** species × range-status strata. Among these **128** Holm-adjusted tests, three colour coefficients in introduced strata survived (Section 3.1), whereas **none of the 16** original head-image-angle × exposure terms survived Holm adjustment. This indicates that the frozen analysis detects some photographic colour–environment structure in the focal strata, without a common supported head-angle association. Failure to reject those angle tests is **not** evidence of equal or unchanged head angles.

The frozen S2 descriptive **four-variable-adjusted RV** medians also differed among image-trait domains. Same-domain versus cross-domain RV was **0.0487/0.00532** in introduced *C. arvense*, **0.0218/0.0119** in native *C. arvense*, **0.1180/0.00785** in introduced *C. vulgare*, and **0.1102/0.0125** in native *C. vulgare*. These are correlations among photographed variables, and differ in sample composition across strata. They are not independent replicated tests of genetic modules or of a phenotypic change caused by introduction.

**Out of scope by design:** no BIO18 occurrence-site median comparison, no native-versus-introduced climatic-niche analysis, no spatial-cell climate-centroid shift, and no post-freeze four-dimensional common-climate range-effect regression are presented as results in this article.

## 4. Discussion

### 4.1. Partial phenotypic organisation, not an indivisible floral syndrome

The most general result is **heterogeneity of association structure**. Pairs of measurements within colour or outline show larger descriptive RV values than most pairs crossing the predefined presentation, colour and outline domains, while many different three-domain configurations occur in photographs assigned to the same taxon. Together, these observations are compatible with partially modular **visible phenotypes**, rather than a single fixed capitulum configuration. This interpretation does not require—and the evidence cannot establish—freely interchangeable genetic "LEGO" components. A single photographic head can yield correlated colour estimates by construction; multiple shape indices similarly arise from the same segmented silhouette. Shared pixels are an important measurement-level alternative explanation for the high within-domain correlations.

The large 24/28 between-over-within count is also descriptive. Between-taxon taxon means may have much less within-group noise and a different range than individual images, and RV baselines depend on dimensionality and sample size. The fixed-size 280-label reference tests how unusual one *labelled grouping* looks in this observed matrix, not how frequently the same integration architecture evolved independently. The three-way analysis reinforces this boundary: combinatorial richness is visible, but a stable higher-order excess or deficit cannot be claimed when changing the dichotomisation changes the statistic’s direction.

### 4.2. Warm-quarter rainfall corresponds to photographic presentation among taxa

The S1 BIO18 coefficient is the strongest frozen ecological association involving capitulum presentation and survives Holm correction over all 80 S1 biological tests. Its cross-taxon scale is essential. The exposure is a long-term warm-quarter rainfall climatology, not rain deposited on the capitulum while photographed. The measured angle is relative to the camera frame, not gravitational vertical. Therefore, this result proposes a geographic pattern in **image presentation** worth testing against rainfall interception, head architecture, observer viewpoint and phylogenetic structure, but is not itself a demonstration of rainfall-driven selection for nodding.

Earlier V4 exploratory models using BIO12 or additional highly collinear climate predictors sometimes retain a positive coefficient and sometimes lose corrected support. These are useful sensitivity findings, yet they must not be combined with the frozen four-variable BIO18 test as if separate independent experiments. Geographic sorting of taxa, spatial autocorrelation, unequal photo availability and shared ancestry remain alternative interpretations. A field study measuring head orientation to gravity and wetting exposure, with geographically replicated populations and independently recorded photographic roll, would resolve substantially more than an additional reparameterisation of the same images.

### 4.3. Colour associations within introduced ranges are not automatically evolutionary shifts

Three corrected S2 associations involve colour rather than a generalizable direction-angle shift: elevation and image lightness in introduced *C. arvense*, plus radiation and chroma/lightness in introduced *C. vulgare*. They indicate that the environmental association structure of photographic colour can be detectable inside specific introduced-range strata. They do **not** show that introduction caused a colour change, that anthocyanin production was altered or that a common environmentally responsive genotype occurred across the two species.

The contrast with S1 is useful as a **scale comparison of modelled phenotype–environment relationships**, not as a climatic-niche comparison: the positive S1 among-taxon BIO18 × image-angle association does not automatically supply a supported S2 within-stratum angle effect. Importantly, the lack of S2 Holm-supported angle slopes neither establishes equivalence between range statuses nor excludes effects under different sampling and measurement designs. Native observations may overlap between S1 and S2, so parallel findings are not independent replication.

### 4.4. What S2 range stratification adds—and does not test

S1 considers the structure of repeated photographs among native *Cirsium* taxa. S2 asks a different, limited question: among two abundantly photographed thistles, which **photographic traits** show associations with the **same four environmental predictors** within native and introduced strata? The observed differences in support across analyses can motivate hypotheses about taxonomic scale, environmental association structure and sampling dependence. They cannot be interpreted as an evolutionary response to introduction or as evidence of climatic-niche conservatism.

Only the two original focal species passed retrospective coverage auditing across both statuses, and introduced records are strongly concentrated in particular geographic regions. Native and introduced observations are neither experimentally randomized nor paired by ancestral genotype. They do not provide independent transitions into comparable rainfall regimes. As a result, this study **deliberately does not estimate** wetward niche shifts, expansion, stability or unfilling, and does not assume that introduced distributions offer a natural rainfall intervention.

Nor can the unsupported S2 head-angle associations demonstrate morphological stasis. The original tests estimate **within-stratum slopes**, not native-versus-introduced phenotype differences under matched climates. Equivalence margins, physical gravity-referenced angles, matched climatic backgrounds, and source-descendant genotypes would be necessary for testing stronger adaptive or evolutionary claims. The manuscript retains S2 chiefly because its significant colour-environment coefficients and descriptive photographic trait integration provide a complementary, sharply delimited **within-stratum** view of global floral trait patterns.

### 4.5. What public photographs measure—and cannot measure

The gain from public photos is the ability to retain many observations and to compare the joint distributions of continuous visible traits across environmental settings. This complements traits summarized in reference databases, but does not replace standardized field phenotyping. Even correct image segmentation yields outline measures dependent on view angle and head maturity; JPEG colour is not a pigment assay; and an EXIF-corrected frame does not define gravity. The present results should accordingly be read as **image-derived phenotypic associations**.

Future mechanistic tests can be staged without rewriting the present source freeze: (i) independent orientation landmarks including a gravity marker and camera-roll record; (ii) spectral reflectance and chemical pigment assays; (iii) three-dimensional or standardized multi-view capitulum geometry, including true bract/spine structures; (iv) phylogenetic and cytotype-aware comparisons; and (v) pollinator, herbivory, rainfall and reproductive-fitness observations. Such data would permit testing whether the apparent modules are developmentally or genetically separable, and whether any configuration actually improves ecological function. These are next-study questions rather than claims about V4 S1–S2.

### 4.6. Statistical scope and generalisability

The estimates derive from a fixed, post hoc selected four-predictor model. Holm control applies to the declared 80/128 test families **within that model**, not to every environmental specification explored before the model was frozen. Pointwise intervals do not provide simultaneous coverage. Multiple photographs, observers, spatial cells and taxa still have complicated dependencies; clustered standard errors and balancing reduce but do not eliminate these concerns. A future independent sample and photo-to-field validation could change both the estimated coefficients and apparent module contrast.

The biological scope is *Cirsium* source taxa in S1 and exactly **two** introduced species in S2; extending the conclusion to all flowering plants or claiming a repeatable adaptation to introduction is premature. Strongly photographed taxa and regions, taxonomic revisions and camera behaviour could all affect which parts of morphospace enter the study. We explicitly exclude putative molecular-module reuse and evolutionary gain/loss reconstructions from the evidence for this paper.

## 5. Conclusion

Across globally distributed *Cirsium* photographs, capitulum presentation, colour and outline show heterogeneous, partly domain-concentrated **observable phenotypic covariance**. In the frozen V4 BIO18/BIO05/RSDS/elevation models, **one S1 among-taxon** image-angle association with warm-quarter rainfall and **three S2 within-introduced-stratum** colour associations with elevation or radiation survived their separate Holm test families. The native S1 common-image cohort also exhibited larger between- than within-taxon RV for 24 of 28 environmentally adjusted trait pairs, while several combinations of visible trait domains co-occurred within taxa.

The ecological contribution is to demonstrate **scale-dependent photographic trait–environment associations** alongside partial organisation of visible floral morphology. It is not a comparison of native-versus-introduced rainfall niches, nor a demonstration of adaptive rainfall shielding, genetic module autonomy or evolutionary change following introduction.

## References

- Diggle, P. K. (2014). Modularity and intra-floral integration in metameric organisms: plants are more than the sum of their parts. *Philosophical Transactions of the Royal Society B* 369:20130253. https://doi.org/10.1098/rstb.2013.0253
- Govaerts, R., Nic Lughadha, E., Black, N., Turner, R. & Paton, A. (2021). The World Checklist of Vascular Plants, a continuously updated resource for exploring global plant diversity. *Scientific Data* 8:215. https://doi.org/10.1038/s41597-021-00997-6
- Kattge, J. et al. (2020). TRY plant trait database—enhanced coverage and open access. *Global Change Biology* 26:119–188. https://doi.org/10.1111/gcb.14904
- Laitly, A., Callaghan, C. T., Delhey, K. & Cornwell, W. K. (2021). Is color data from citizen science photographs reliable for biodiversity research? *Ecology and Evolution* 11:4071–4083. https://doi.org/10.1002/ece3.7307
- Pigliucci, M. (2003). Phenotypic integration: studying the ecology and evolution of complex phenotypes. *Ecology Letters*. https://doi.org/10.1046/j.1461-0248.2003.00428.x
- Schiller, C., Schmidtlein, S., Boonman, C., Moreno-Martínez, A. & Kattenborn, T. (2021). Deep learning and citizen science enable automated plant trait predictions from photographs. *Scientific Reports* 11:16395. https://doi.org/10.1038/s41598-021-95616-0
- Siefert, A. et al. (2015). A global meta-analysis of the relative extent of intraspecific trait variation in plant communities. *Ecology Letters* 18:1406–1419. https://doi.org/10.1111/ele.12508
- CHELSA climatologies, version 2.1: https://www.chelsa-climate.org/datasets/chelsa_climatologies


## Data and code availability

The exact primary four-variable statistical results, the 208-test multiple-comparison table and the 28-pair matched covariance matrix have source-verified internal records and analysis code. Durable long-term deposition of redistributable data/code, a genuinely anonymous reviewer-access URL and attribution of third-party public photographs are **pending** before submission. Geographic coordinates and uploader identifiers from protected source datasets will not be redistributed; original photos are subject to source licences.

## Main display pieces

**Figure 1.** Source-observation/cohort flow and source-supported geographic coverage (original source geography and image attribution required before publication).

**Figure 2.** Four Holm-supported associations with standardized coefficient and pointwise 95% interval; full 208-test family reported in Supporting Information.

**Figure 3.** All 28 RV associations for native-taxon common photographs, raw/adjusted and within-/among-taxon scales.

**Figure 4.** Four original S2 species-by-range strata: comparison of the frozen **same-domain versus cross-domain photographic RV medians**, plotted as descriptive values without confidence-interval claims. These are image-trait association summaries, **not climate-occupancy or niche-shift results**.

**Figure 5.** Three-way trait combinations and alternative threshold definitions, observational only.

**Tables 1–2.** Frozen Holm associations and descriptive same-/between-domain RV medians.

