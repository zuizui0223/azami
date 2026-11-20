# Supporting Information — V5 annual-climate *Cirsium* manuscript

**Companion to:** Climatic matching reveals local capitulum resemblance and trait-specific gradients in *Cirsium*  
**Target:** *Global Ecology and Biogeography*, Research Article  
**Status:** 19 November 2026 working draft. Source-audited descriptions of already adopted V5 analyses; not a report of new fits. Author identifiers intentionally omitted.

The source of record for the figures, counts and estimates below is the adopted annual five-climate-axis V5 analysis on repository main, frozen for this draft at e7c57c8693478df80d6fd488d5a5ab53b41408fd. See ANALYSIS_PHOTO_COUNTS_JA.md and EVIDENCE_AND_SUBMISSION_GATES_JA.md in this manuscript directory for the corresponding source-table audits. Excluded exploratory PR #174 (within-cell climate pixels) is not part of this manuscript or Supporting Information.

## S1. Observational units, selection and provenance

The underlying source includes taxonomically assigned biological observations, links to source photographs and computer-vision-derived traits. These are **not interchangeable observation units**. One source observation need not correspond to one unique JPEG, and a photograph may show multiple detectable capitula. Likewise, a geographical cell containing many photograph-derived records is only one cell-level contribution to the relevant aggregated test. We make no claim that the sum of source observation counts equals the number of unique photographs processed, nor that the S1 and S2 sources are disjoint.

Part I starts from 313,350 records of wild observations classified as native in the S1 source. Of these, 157,822 photographic measurement records across 242 taxa satisfy the joint requirements for image-referenced capitulum angle, photographed corolla chroma, projected outline and the adopted annual climatic covariates. Geographic eligibility and the subsequent eight-cell sampling further restrict the sets contributing to each spatial-grain test. At 1°, 153,843 photographic measurement records lie within the retained taxa and occupied cells, yielding 63 taxa and 7,117 taxon–cell units (3,107 distinct geographic cells). At 5°, the corresponding values are 132,998, 31 taxa and 814 taxon–cell units (294 distinct geographic cells). In each of 24 rarefaction draws, eight occupied cells are sampled per eligible taxon. A source-pool denominator of 157,822 must not be described as 157,822 independent climatic analogues, unique image files or independent permutation units.

Part II starts from 297,622 wild observations in the native and explicitly introduced WCVP-defined regions of *Cirsium arvense* and *C. vulgare*. Annual-climate PCA occupancy uses 144,178 records after climate eligibility and observer–geographic-cell thinning: *C. arvense* native 30,523 and introduced 38,406; *C. vulgare* native 26,161 and introduced 49,088. A distinct photo-derived head-trait-complete layer contains 79,421 records: *C. arvense* native 18,505, introduced 22,269; *C. vulgare* native 13,115, introduced 25,532. The head-trait analysis is thus not a reanalysis of all 144,178 climate-only records.

**Table S1. Audit of analysis denominators.** Values below represent different stages of filtering or aggregation and must not be added as independent data.

| Analysis | Initial source | Eligible measurement records | Final units for the test |
|---|---:|---:|---|
| Part I native source | 313,350 source observation records | 157,822 joint-trait-and-climate records; 242 taxa | Common eligibility pool only |
| Part I, annual five axes, 1° | Same Part I source | 153,843 records in retained cells | 63 taxa × 7,117 taxon–cells; 24 draws |
| Part I, annual five axes, 5° | Same Part I source | 132,998 records in retained cells | 31 taxa × 814 taxon–cells; 24 draws |
| Part II climate occupancy, *C. arvense* | 78,646 native / 70,964 introduced | 30,523 native / 38,406 introduced | Occupancy densities on PCA climate grid |
| Part II climate occupancy, *C. vulgare* | 62,368 native / 85,644 introduced | 26,161 native / 49,088 introduced | Occupancy densities on PCA climate grid |
| Part II photographed heads, *C. arvense* | Same two-species source | 18,505 native / 22,269 introduced | 1,366 native / 860 introduced 1° cells |
| Part II photographed heads, *C. vulgare* | Same two-species source | 13,115 native / 25,532 introduced | 1,113 native / 1,139 introduced 1° cells |

The original archive contains 1,122,854 photograph links associated with 665,115 source records across a broader set of work. Those totals are not estimates of V5's analysed or unique JPEG files and should not be used to advertise the V5 sample size.

## S2. Taxonomic status and source-population interpretation

Species names and native/introduced geographic assignments in Part II follow the original source's World Checklist of Vascular Plants (WCVP) regional status mapping. Explicitly introduced observations are kept distinct from native and ambiguous-status observations. The geographic category 'native' is not a genetically verified source population, and 'introduced' does not identify an independent colonisation event. Analyses neither reconstruct introduction history nor compare independently replicated source–descendant population pairs.

The 63-versus-31 native-taxon sample difference between 1° and 5° reflects the existing occupied-cell eligibility criteria, not a taxonomic revision. The result at 5° is therefore conditional on its smaller set of sufficiently sampled taxa. The planned same-31-taxon grain sensitivity has **not** been run and is not reported as completed.

## S3. Photographic phenotype definitions and limitations

Three measurements were retained together on an eligible source photograph: image-referenced capitulum angle, photographed corolla CIELAB chroma C* and a projected outline aspect index. Higher coded angle means a more downward-looking head **in the photographic frame**, not a measured angle against local gravity. Chroma reflects image colour intensity and is not pigment concentration or spectral reflectance under standardised illumination. Projected outline aspect is viewpoint- and rotation-dependent and is not a direct anatomy-invariant shape index.

The source analysis retains its pre-existing image-quality rules; this supplement does **not** invent unreported detector validation metrics, photography hardware calibration, colour-card correction, interobserver error statistics, or trait-specific missingness rates. Source-image detection and continuous measurement require independent validation if biological claims about physical orientation, pigmentation or intrinsic form are to be elevated beyond visible-image phenotypes.

Cell-level trait values were constructed by taking the median of eligible photograph-derived measurements within each taxon-by-geographic-cell unit. This aggregation limits the direct weight of many repeated records in a single cell but does not eliminate uneven geographic coverage, observer selection or image-processing bias. The adopted analysis is observational and cross-sectional.

## S4. Climate variables and three distinct geographic resolutions

The annual-climate analysis adopts five CHELSA v2.1 1981–2010 normals:
(1) BIO01, annual mean temperature;
(2) BIO04, temperature seasonality;
(3) BIO12, annual precipitation;
(4) BIO15, precipitation seasonality; and
(5) RSDS, downwelling shortwave radiation.

Original values were assigned at photographic coordinates from the original 30-arcsecond CHELSA raster pixels (about 1 km north–south, varying in ground dimensions with latitude). These are long-term climate normals, not weather observed on the photo date.

Part I then aggregates photographs and their assigned climate values into angular latitude–longitude cells of 1° or 5°. The corresponding north–south angular spans are approximately 111 km and 555 km; east–west ground spans decrease with latitude. Taxon-by-cell medians represent the **photographed locations in a cell**, not a spatial mean over every high-resolution raster pixel in that cell. The spatial unit is the taxon–cell, not a 1-km or 5-km plot.

Part II represents regionally available climatic environments by a 0.5°-spaced geographic sampling mesh; the original audit reports 85,414 background mesh points in TDWG regions across the two-species scope. A climate value at each mesh point is sampled from its corresponding original CHELSA pixel. Latitude-dependent area weights are used in characterising the available environment. This background mesh is not the unit of the native 1°/5° analogue test.

## S5. Part I: five-axis climatic analogue and three-trait resemblance tests

For each spatial grain, climate coordinates and the jointly observed photographic trait vectors were standardised with taxa given equal total weight. The source annual five-axis climatic-analogue rule uses a scaled Euclidean criterion equivalent to a root-mean-square standardised climate-axis distance of 0.9 (full vector threshold 1.161895003862225). In each of 24 rarefaction draws, eight eligible occupied cells per taxon were sampled without replacement. Between-taxon comparisons considered the 64 possible cell pairs within each eligible taxon pair, and within-taxon matching was separately assessed. These cells and draws are technical sampling supports, not independent evolutionary replicates.

The response is the signed excess reduction in standardised **joint** angle–chroma–outline distance for climate-matched pairs compared with its taxon-conditioned reference and permutation expectation. The null shuffles the **entire measured three-trait vector** among sampled geographic cells **within each taxon**. This preserves the observed taxon identity, its multivariate phenotype marginal distribution and intertrait correlation but destroys its empirical allocation to climate. The adopted annual test uses 499 original taxon-preserving permutations with the existing stream and 24 rarefaction draws. Original two-grain Holm correction families are retained separately for the source's within- and between-taxon tests; repeated draws are not assigned independent P values as though they were distinct original datasets.

**Table S2. Annual five-axis multi-trait analogue tests.** Positive values represent greater similarity than the permutation reference, not estimated selection coefficients.

| Comparison | Grain | Taxa | Taxon–cells | Excess distance reduction | Original Holm P |
|---|---|---:|---:|---:|---:|
| Within taxon | 1° | 63 | 7,117 | +0.084212 | 0.004 |
| Within taxon | 5° | 31 | 814 | +0.238369 | 0.004 |
| Between taxa | 1° | 63 | 7,117 | +0.030938 | 0.076 |
| Between taxa | 5° | 31 | 814 | +0.065330 | 0.012 |

**Important untested alternative.** The permutation retains each taxon's phenotypic distribution but not spatial dependence. Climatic analogues may be geographically closer than unmatched pairs, and regionally structured photographic traits can create apparent environmental matching without climatic effects. Consequently, these tests quantify *association beyond the taxa's marginal trait distributions*, **not** association beyond spatial distance, regional history or shared environmental causes. A geography-preserving or distance-stratified null, and spatially separated climate-analogue comparisons, are submission-critical controls rather than already successful results.

## S6. Directional climate associations are a different estimand

The adopted signed regressions are the original **all-cells, within-taxon-pair** comparisons, not a regression limited to the subset of climate-matched pairs used in S5. Their source-filter identifiers are sample=all_cells and scale=within_pair in analysis/v5/results/native_seasonality/v5_signed_pairs_seasonal.csv. Analyses remove mean differences among taxon pairs when relating signed geographic cell-pair differences in climate to signed differences in angle, chroma and outline. The five annual climate variables are included together. Spatial grains share underlying photographic sources, so two point estimates are not independent validation studies.

**Table S3. Source-verified signed slopes used in the manuscript.** The intervals are taxon-bootstrap, *pointwise* 95% intervals for individual coefficients. They have not been adjusted simultaneously across all climate-axis × trait outcomes.

| Climate axis → image trait | 1° standardised β [95% interval] | 5° standardised β [95% interval] | Inferential status |
|---|---|---|---|
| BIO12 precipitation → downward-looking angle | +0.1171 [0.0154, 0.2202] | +0.1829 [0.0216, 0.3393] | Exploratory geographical association |
| BIO15 precipitation seasonality → chroma C* | +0.2059 [0.0670, 0.3476] | +0.2263 [0.0464, 0.3747] | Exploratory geographical association |

The source multivariable paired design showed maximum variance inflation factors of 1.66 (1° within-pair) and 2.20 (5° within-pair), which do not suggest severe *linear predictor collinearity*. This diagnostic does not address spatial dependence or unmeasured environmental confounding. A pointwise interval excluding zero is neither proof of rain protection nor familywise confirmatory evidence over the full search space.

## S7. Separate three-climate-axis trait tests

An earlier analysis used BIO01, BIO12 and RSDS alone, without BIO04 and BIO15, and asked whether climate matching was associated with similarity in individual photo traits. It used the existing taxon-preserving 499-permutation scheme and corrected its **six individual-trait tests (two grains × three traits)** using a separate Holm family. Chroma had excess distance reductions of +0.0300 (1°; Holm P=0.024) and +0.0291 (5°; P=0.030); image angle had +0.0246 (1°; P=0.032) and +0.0034 (5°; P=1.000). These results cannot be labelled the component-wise decomposition of the **annual five-axis** joint test without running that precise decomposition, and their P values cannot be pooled with the five-axis Holm family as an independent replication.

The source audit's 67,719 (1°) and 114,448 (5°) photographic records appearing in matched three-axis comparison sets refer specifically to this earlier model, not the eligible photo count or comparison count of the five-axis tests.

## S8. Part II: background-corrected native–introduced climate occupancy

For each of *C. arvense* and *C. vulgare*, the observed and available annual climates were projected to a common **species-specific two-dimensional PCA**. Native and introduced geographical availability were represented by the area-weighted 0.5° background mesh (S4), and the occurrence density was corrected by environmental availability. Overlap was computed on the source 80 × 80 PCA grid as Schoener D = sum of the smaller of the normalised native and introduced density values per grid square (from 0 to 1).

The adopted overlap estimates are D=0.474 for *C. arvense* and D=0.442 for *C. vulgare*. The proportion of the **introduced background-corrected occupied climate-density mass** lying outside the main 95% region of native occupancy was 18.8% and 44.0%, respectively. These are not proportions of introduced photographs or the fraction of genomes with evolved climatic tolerance. Native-to-introduced background comparisons also differ in geographic opportunity; D is conditional on the WCVP regional boundaries, environmental covariates, availability weights and smoothing.

## S9. Part II: analogous/outside-core head overlap and geographic matching

Introduced photographed head-trait cells were assigned to the operational 'inside the main 95% native climatic core' and 'outside that core' groups. This boundary is relative to a bivariate PCA density estimate of current native climatic occupancy; outside-core does not mean outside species' physiological tolerance.

The main three-trait head-space statistic is overlap against a **common native photographed head reference**. The introduced groups were balanced by selecting equal numbers of **1° cells within 20° geographical blocks**, followed by region-block label permutations. For *C. arvense*, a total of 823 introduced cells with complete traits fell in an unequivocally inside-core category and 30 outside-core; seven cells were mixed. For *C. vulgare*, the counts were 931 inside, 190 outside and 18 mixed. The final geographically paired comparisons comprised 30 cell pairs in nine blocks and 146 pairs in 14 blocks respectively. Block counts must not be construed as independent introductions.

The *C. arvense* native-reference head D was 0.740 (inside-core) versus 0.735 (outside), difference +0.005, two-species Holm P=0.953. The *C. vulgare* corresponding values were 0.847 and 0.789, difference +0.057, Holm P=0.096. For *C. vulgare*, removing the rotation-dependent outline in the two-trait angle–chroma sensitivity left a difference of +0.083 but yielded raw regional permutation P=0.066. None of these outcomes is an equivalence test establishing that head phenotypes are identical across climates.

The separately published descriptive comparison of 938 inside-core versus 201 outside-core *C. vulgare* cells assigns all 1,139 cells to a climatic category and has **different handling of mixed cells**. Its sample counts cannot be substituted for the 931/190 unequivocal head-overlap strata.

## S10. Part II: chroma contrasts and environmental sensitivity

In *C. vulgare*, the introduced-minus-native mean chroma difference at the geographic 1° cell unit was +1.7414 C*. Reweighting to shared climatic PCA bin frequencies (4 × 4 fine grid cells combined into one of the 20 × 20 coarse classes) gave +2.4646 C*. A coarser-bin alternative gave +2.3486 C*. The principal shared-support comparison retained 1,052 native and 1,015 introduced **geographic cells**, with reweighted effective sample sizes of approximately 694.2 and 691.9. These are descriptive cell-level contrasts, not a fitted evolutionary effect or a formal causal regression coefficient.

Five of seven introduced regions had chroma **distribution dissimilarity (1−D)** greater than the native regional 90th percentile; none exceeded the native regional maximum. This is evidence of unusually different *distributions* under a descriptive native-region benchmark, not evidence that all five introduced regions have higher average C*. The geographically matched native-reference chroma-only D for introduced analogous versus outside-core sites was approximately 0.965 versus 0.860, difference +0.105, but its six-comparison Holm P=0.462. It must not be described as a familywise-supported outside-core-specific chroma change.

The broad *C. vulgare* outside-core climatic group (201 1° cells) was drier and more precipitation-seasonal than the inside group (938 cells): medians BIO12 697.0 versus 897.5 mm; BIO15 60.50 versus 30.05; chroma C* 32.23 versus 31.40. This three-way descriptive contrast is **not** an adjusted model of why colour differs. In direct posthoc dry/high-seasonality versus wet/low-seasonality contrasts, the dry/high-seasonality subset instead had *lower* C* (native −1.91, introduced −3.12). In adjusted outside-core region-aware models, BIO12 β=−0.045 (raw P=0.166) and BIO15 β=−0.057 (raw P=0.432) both had Holm P=1. Therefore the manuscript does not conclude that low precipitation, precipitation seasonality or post-introduction adaptation caused the signed chroma contrast.

## S11. Claim-to-file provenance and figure-production controls

**Table S4. Source mapping; exact source paths are on repository main at the frozen analysis commit.**

| Main-paper output | Source files | Unit and plotting rule |
|---|---|---|
| Fig. 1, Table 1 cohort flow | manuscripts/v5_annual_20261111/ANALYSIS_PHOTO_COUNTS_JA.md (PR branch); analysis/v5/results/three_stages/v5_photo_to_cells.csv | Flow of original source records → eligible photo-derived records → taxon-cells; distinguish data-source photograph links |
| Fig. 2, five-axis resemblance | analysis/v5/results/native_seasonality/v5_within_analogue_seasonal.csv; v5_cross_analogue_seasonal.csv | Separate within- and between-taxon comparisons, two grains; original Holm P; no invented uncertainty bars |
| Fig. 3, signed gradients | analysis/v5/results/native_seasonality/v5_signed_pairs_seasonal.csv | Filter sample=all_cells and scale=within_pair; plot pointwise bootstrap intervals only |
| Fig. 4, climate occupancy | analysis/v5/results/niche_overlap/README_JA.md and associated source-locked density outputs | Two separate species-specific PCA environments, 80 × 80 density; D and outside-core mass are distinct quantities |
| Fig. 5, head D and chroma | analysis/v5/results/novel_head_overlap/novel_head_overlap.csv; analysis/v5/results/invasion_design/climate_balanced_heads.csv; analysis/v5/results/invasion_design/RESULTS_JA.md | Match 20° blocks for inside/outside D, but common climatic class support for C*; show distinct denominators |
| Fig. S1, original 3-axis single-trait test | analysis/v5/results/three_stages/v5_trait_convergence.csv | Draw three-axis chroma/angle components separately from Fig. 2's five-axis test |

No high-resolution figure file is asserted to be completed by this document. Figure 4's full climate-density surfaces and Figure 1's real photographs need their original authorised source artifacts. Redistribution of third-party photographs must comply with individual licences and attributions; placeholders are not evidence of permission. Vector figure production should use the named locked CSV values, declare any rounding and avoid derived invented confidence intervals.

## S12. Mandatory validation and submission gates

1. **Spatial null and taxonomy/grain**: Compare climatic analogue pairs against geographically distance-matched or region-constrained controls and test both 1° and 5° on the same eligible 31 taxa. If the excess disappears, revise the headline conclusion instead of relegating it to a footnote. Preserve the originally reported analysis as the observational baseline.
2. **Photographic measurement validity**: Report blinded detector precision/recall; independent physical or repeat-photo validation for orientation, colour and projected shape; examine whether trait assessability or camera geometry varies with geography. Report measured results, not synthetic or pseudolabel accuracy as independent validation.
3. **Climatic and taxonomic provenance**: Freeze WCVP accepted names, regional status, exclusions, CHELSA raster identifier, geographic aggregation, source collection date and original run identifiers. Repeat sensitivity to source and geographic imbalance as permitted by data provenance.
4. **Data access and originality**: Deposit rights-cleared reproducible aggregate tables and analysis code, create a persistent DOI or review-accessible archive, and verify photograph reuse permissions. This V5 draft substantially overlaps the frozen V4 GEB working manuscript in original observations and Part I analyses; treat V5 as an **integrated replacement/revision** unless distinct originality and nonoverlap are independently established.
5. **Journal materials**: Confirm original figure artwork, support files, submission-ready title pages, authorship, acknowledgements, competing interests, CRediT contributions and funders. Prepare the anonymous main text and separate identifying page under the journal's double-anonymous policy.

**No result listed as 'mandatory validation' has been silently assumed successful.** Functional rain shielding, anthocyanin physiology, pollinator responses, genetic differentiation, and evolutionary adaptation are hypotheses for later direct experiments.

## S13. Reproducibility statement

The manuscript is a synthesis of **already completed and source-audited annual V5 analyses**, not a new computational reanalysis. Original scripts, aggregate output and inspection ledgers reside in https://github.com/zuizui0223/azami under analysis/v5/ and analysis/v5/results/. The main scientific source baseline for this draft is e7c57c8693478df80d6fd488d5a5ab53b41408fd. Original observation-linked photographs, raw coordinates or credentials are not redistributed by this supplement. A review-accessible archive and rights-cleared permanent data deposit remain to be produced before submission.
