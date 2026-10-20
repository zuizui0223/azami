# Native-first 27-endpoint handoff (2026-10-20)

## Decision: what is active and what remains paused

Chapter 1's frozen v3/GEB outputs remain the historical reference; they are **not** automatically native-only by replaying them. The next primary cohort is selected from the **original pre-detector source, prior to any detector or trait selection**, using a fixed native-range definition. An all-upstream native classification and inherited-measurement linkage already exist in the archived PR142 evidence. **Reuse that completed result; do not reacquire WCVP, CHELSA rasters, original images, or remeasure missing pixels in this consolidation.** No new V4 Actions or competing analysis entry point is authorized by this note.

The complete pre-consolidation branch is preserved at [archive/pr142-v4-20261020-preconsolidation](https://github.com/zuizui0223/azami/tree/archive/pr142-v4-20261020-preconsolidation). Archived scripts and results are diagnostic evidence, not part of the currently merged v3/GEB scientific claim set.

## Cohort ledger and fixed ordering

| Stage / denominator | Count | Meaning |
| --- | ---: | --- |
| Original pre-detector metadata, 2026-07-03 | 665,115 observations; 1,122,854 photo links | Starting frame, not 665,115 trait measurements |
| Native-primary after source-level classification | 313,350 observations; 534,186 photo links; 306 accepted WCVP keys | Before phenotype, detector, colour, and trait completeness |
| Native with all nine stored CHELSA values | 313,240 observations | Environment availability, **not** measurement availability |
| Native with an inherited original measurement-table row | 217,104 observations | Join to previous image measurements, not remeasurement |
| Native lacking any inherited original measurement-table row | 96,246 observations | Preserve as missing; do not silently drop from denominator |
| Native inherited measurements plus all nine environmental values | 217,102 observations | Still not a 27-trait complete-case cohort |

Audit provenance: archived [full-source master report](https://github.com/zuizui0223/azami/blob/archive/pr142-v4-20261020-preconsolidation/analysis/v4/NATIVE_MASTER665_RESULT_20261020.md), [full-source receipt](https://github.com/zuizui0223/azami/blob/archive/pr142-v4-20261020-preconsolidation/analysis/v4/NATIVE_665115_FULLSOURCE_20261020.md), original all-source replay [run 39407663676](https://github.com/zuizui0223/azami/actions/runs/39407663676).

Native-primary requires exactly resolved accepted WCVP taxon, a matching native TDWG3 region, usable un-obscured coordinates, and a non-captive source record. Missing WCVP region, uncertain identification, introduced region, invalid/ambiguous coordinates, unknown captive status and non-species ranks are separate **exclusion/unresolved** states; they must not be recoded as an indigenous absence. TDWG3-level native membership does not prove local indigenous origin. Keep the all-source audit privately; never publish individual coordinates, source URLs, record IDs or private rights ledgers.

**Selection order (enforced conceptually by the archived `native_scope.py` / `native_665115_pipeline.py`):**

1. Original metadata → exact accepted-name/taxonomic mapping → geographic/captive native-primary flag for **every source observation**.
2. Freeze that native roster and source/authority hashes.
3. Left-join inherited observation-level detector/trait values by observation ID, verify taxon identity and one-to-one uniqueness; do not require that a detector succeeded to enter the native denominator.
4. Join already-cached environmental values and keep explicit nodata; use only environment rows available for the specific model.
5. Select eligible **per-endpoint** response rows for analysis (not a 27-column complete-case), carrying QC, missing reason, taxon, spatial/observer sampling attributes and image-level provenance.

## All 27 endpoints: registration versus measurement

The frozen inventory is [`v2_full27_endpoint_inventory.csv`](../analysis_outputs/v2_full27_environment_atlas_2026-08-27/v2_full27_endpoint_inventory.csv). Counts in its `n_observations_measured` column belong to the historical **46,276-observation V2 frame**. They are **not** counts from the 313,350-observation native master.

| Status and interpretation | Registered endpoints | Historical measured observations per endpoint | Native-first decision |
| --- | ---: | ---: | --- |
| Core photographic axes: image-vertical head angle (1); floral L*, C*, hue sine/cosine (4); 2D whole-head silhouette metrics (4) | 9 | 37,251; 40,503–40,510; 38,511 respectively | Retain separately with trait-specific QC; image proxies, not calibrated physical traits |
| Involucre body (3) and bract *contour projection* (6) | 9 | 2,326 | Exploratory/validation candidate only; not actual phyllary opening angle, spine count, or defence |
| Surface texture (3) and specularity (1) | 4 | 909 | Technical/validation-only; specular fraction is **not** sticky secretion |
| Visible-floret fraction (1) and white, red/magenta, purple, yellow floral pixel-composition fractions (4) | 5 | **0; never executed in V2** | Keep registered with NA/unknown, not biological absence; defer any new image segmentation |

Thus **22/27 endpoint definitions have at least some inherited measurements; five have none**. This is not 22/27 endpoints measured for every observation. Hue sine and cosine are one circular colour variable, and four compositional fractions are not four independent biological traits.

Archived native-first conditional-environment result [run 39416125522](https://github.com/zuizui0223/azami/actions/runs/39416125522) reported *response-specific* eligible rows: angle **170,064 / 244 taxa**; lightness/chroma **186,437 / 244 taxa**; circular hue **186,388 / 244 taxa**; bract contour projection/spread **1,303 / 111 taxa**. These are **inherited photo measurements**, not 313,350 newly analysed/masked heads, not all 27 measured, and not a confirmed V3 result. The corresponding [archived result note](https://github.com/zuizui0223/azami/blob/archive/pr142-v4-20261020-preconsolidation/analysis/v4/results/20261020_selected_environment/README.md) documents the exploratory joint-model limits.

## Native-first 27-row aggregate receipt — verified, not an extrapolation

The [machine-readable 27-endpoint coverage inventory](native_first_27_endpoint_coverage_20261020.json) was reconstructed from the **actual** aggregate `measurement_receipt.json` in [Actions run 39408603488](https://github.com/zuizui0223/azami/actions/runs/39408603488), artifact `12257893990`, receipt SHA-256 `246a5ea9643617c6dc766e0b8d27a444a5496afe150d7768032740cf6a53f0a6`. It reports each endpoint's **finite**, **eligible** and **eligible with four climate predictors** counts across the same 313,350 native records. No private observation IDs or raw photographs are included.

| Native-primary trait measurements | Per-endpoint numeric | Per-endpoint computationally eligible | Meaning |
| --- | ---: | ---: | --- |
| Image vertical orientation | 170,066 | **170,066** | Image angle, not gravity |
| Corolla L* and chroma | 186,439 each | **186,439 each** | Uncalibrated JPEG colour |
| Hue sine/cosine | 186,390 each | **186,390 each** | One joint circular colour trait |
| Whole-head outline (4) | 176,048 each | **176,048 each** | Viewpoint-dependent silhouette |
| Involucre body (3) and bract contour (6) | 1,303 each | **1,303 each** | Eligible in code, but botanical landmarks **not independently validated** |
| Involucre texture (3) and specular fraction (1) | 570 each | **0** | Values exist, but *every* native eligibility flag is false (validation-only) |
| Visible floret fraction and four colour composition fractions | 0 | **0** | Five endpoints never executed |

So **27 registered → 22 with finite values → 18 computationally eligible in the native-only roster**. Even among the 18, nine involucre candidates remain unsuitable for strong botanical-function claims. The four surface outputs are **not** eligible just because 570 numeric values exist. The five never-executed outputs have **no measurement**, not an observed zero state.

The full native roster has **313,244** observations with the four previously selected climate predictors, versus **313,240** with all nine; these are predictor-availability denominators, independent of endpoint measurability. Per-endpoint usable counts in the receipt happen to remain unchanged on the four-predictor complete-case gate. A different selected-predictor set or model design can change them and must report its own denominators.

## Next scientific analysis: bounded by existing data

- **The 27-endpoint/native-taxon ledger is now verified.** The [machine-readable inventory](native_first_27_endpoint_coverage_20261020.json) includes native-only finite/eligible/environment-complete counts and accepted-taxon coverage for each endpoint. **Not yet available** from these public aggregates: per-endpoint photographer, region, taxon-by-region and QC failure-reason counts. These must be computed only within an authorized environment from the pre-existing private native table, without another acquisition or image measurement.
- **One response-specific analysis rather than an all-traits complete-case analysis.** Use the nine existing core photographic endpoints first, preserving the circular joint hue unit. Define response-specific samples, taxon/observer balance, within- versus among-taxon associations, spatial holdouts and multiplicity before inspecting results. Treat shape and orientation as view-dependent visible phenotypes.
- **Separate supporting audits.** Involucre contour/texture/specularity can describe *image signals* only until expert landmark, scale, lighting and surface validations exist. Keep them out of unqualified ecological-function claims.
- **Only if needed later:** measure the five never-executed endpoints from independently checked floral masks and original pixels; do not reconstruct pixel-composition fractions or pixelwise a*/b* distributions from legacy channel medians. Added measurements require their own provenance and independent validation.
- **Do not overwrite the v3 manuscript.** The prior mixed-range V3 baseline stays frozen; a native-only result is a new, labelled comparison only after exact sample counts, uncertainty, dominant-taxon sensitivity, observer effects and geographic transfer have been checked.

## Operational gate

The native-primary **source selection is already completed** in archived PR142 and need not be reacquired. The native-first **27-endpoint coverage ledger and independently supported ecological results are not yet complete in current main/PR142**. The current [V3 numerical replay](CURRENT_ANALYSIS.md) validates old frozen inputs and 16 reference outputs; it **does not** implement a native-only reanalysis by itself. Do not equate green staging/integrity CI with execution of this handoff. The next executable step is reuse of the archived, checksum-verified 313,350-observation native-only prepared table to generate the 27-row eligibility ledger, without downloading or remeasuring images.


## Existing native-only joint environmental analysis — reviewed, not newly fitted

This is a **retrospective review of already-computed native-only exploratory models**, not a new model, a replacement for the frozen GEB v3 analyses, or independent validation. The source is [existing Actions run 39416125522](https://github.com/zuizui0223/azami/actions/runs/39416125522), artifact \`12261727387\`, \`environment_response.json\` SHA-256 \`707e3cf0041d127e085bfa47ef83a6ff89a1ad84740cd4ec6c09a893bbe1c517\`. Its \`reuse_inventory\` supplies accepted taxon counts for all 27 endpoints and matches every native coverage count in the earlier \`measurement_receipt.json\`.

The **environment-only** \`usdm::vifstep\` screen used 313,240 native rows with all nine environmental candidates, removed \`chelsa_bio01\` and retained **eight** variables: \`bio04\`, \`bio12\`, \`bio15\`, radiation, VPD, wind, growing-season precipitation and NPP. Selection is not causal identification. The fitted model analysed **within-taxon deviations and among-taxon means simultaneously** using approximate response-specific random-intercept REML with taxon-clustered joint tests. This is not a spatial model, random-slope reaction norm, or phylogenetically corrected model. The actual fitted row counts differ slightly from the source-eligible counts because model predictors are required.

| Pre-existing joint response | Fitted native observations | Accepted taxa | Holm-supported joint module × environment tests (between / within) |
| --- | ---: | ---: | --- |
| Floral L*, chroma, circular hue | 186,437 | 244 | **5 / 5** of eight gradients at each scale |
| Image-vertical orientation | 170,064 | 244 | **1 / 0** |
| Two bract contour projection proxies | 1,303 | 111 | **0 / 0** |
| Four head-outline metrics | **Not fitted in this result** | 243 with eligible measurements | Not evaluated, **not a null result** |

All three fitted modules comprise **48 joint tests** (3 modules × 2 scales × 8 variables), with Holm correction over the single declared family. For colour, the supported between-taxon predictors were annual precipitation, radiation, VPD, growing-season precipitation and NPP; supported within-taxon predictors were temperature seasonality, radiation, VPD, wind and growing-season precipitation. A four-coordinate joint colour test does not prove that each coordinate is separately supported, nor does it identify a pigment or adaptive mechanism.

For orientation the only Holm-supported test was **between-taxon annual precipitation**: image-referenced head angle estimate **+2.395°** per **537.774 mm** (one native-environment SD), *pointwise* 95% interval **+1.033° to +3.757°**, Holm-adjusted \`P = 0.02398\`. The within-taxon coefficient was **−0.330°** (pointwise interval **−1.609° to +0.948°**), with **no Holm support**. Thus, the detected among-taxon rain association cannot be described as an observed within-species plastic response. Four head-outline endpoints are already available in the native ledger, but this specific conditional model did **not** fit them; no inferential shape claim follows from its absence.

Separate exploratory 665,115-source master-frame diagnostics ([archived master report](https://github.com/zuizui0223/azami/blob/archive/pr142-v4-20261020-preconsolidation/analysis/v4/NATIVE_MASTER665_RESULT_20261020.md)) report native marginal chroma–radiation **β = −0.3252** (177 taxa) and orientation–precipitation **β = +0.2288** (178 taxa). These are **different marginal estimands and cohorts**, not direct coefficients from the jointly adjusted 244-taxon model. The separate master-frame geographic held-out skill was **negative** for its unregularized joint-environment-and-space model. Consequently **do not claim reliable regional transfer, universal causal environmental controls, or an established defence/pollination mechanism** based on this existing analysis.

### Next tightly bounded execution gate

Use the already verified 313,350-row encrypted native measurement input from [run 39408603488](https://github.com/zuizui0223/azami/actions/runs/39408603488) **without reclassification, new images or raster downloads**. Complete the still-missing observer/region/taxon-by-region/endpoint-failure coverage aggregate privately and publish *only* non-identifying aggregate diagnostics. The archive currently permits a full count/taxon ledger; it does **not** permit reconstructing those stratified summaries from aggregate counts alone. Any subsequent shape model or revised ecological claim must be clearly named, fit with a common predeclared estimator, and validated separately before promotion to current paper evidence.

## Cross-repository phenotype-space / functional-hypervolume handoff — existing evidence only

**Purpose:** Recover the older central question — *how do combinations of Cirsium capitulum phenotypes occupy multidimensional space, and how do this occupancy, within/among-taxon covariance and environmental regimes relate?* — without reopening upstream acquisition, adding a fourth parallel workflow or changing the frozen GEB v3 claims. This is an evidence/estimand map, **not a newly computed hypervolume**.

### Existing actual analyses (do not conflate)

| Owner | Existing output | Evidence / interpretation boundary |
| --- | --- | --- |
| **Azami v2/v3** | [PCA taxon morphospace loadings](../analysis_outputs/v2_full27_environment_atlas_2026-08-27/secondary/pca/geb_v2_taxon_morphospace_loadings.csv) and [scores](../analysis_outputs/v2_full27_environment_atlas_2026-08-27/secondary/pca/geb_v2_taxon_morphospace_scores.csv); [v3 scale-dependent construct covariance](../analysis/v3/README.md) | There **is** multivariate trait-space ordination, module structure, within/among comparisons and environmental anchors; **not** a sample-balanced n-dimensional native-only hypervolume estimate or measured ecological functions. V3's complete-18 subset is **1,734 observations / 42 taxa**, unlike the new native-only 313,350 source population. |
| **Azami PR142 archived** | [functional_pilot.py](https://github.com/zuizui0223/azami/blob/archive/pr142-v4-20261020-preconsolidation/analysis/v4/functional_pilot.py), [native_statistics.py](https://github.com/zuizui0223/azami/blob/archive/pr142-v4-20261020-preconsolidation/analysis/v4/native_statistics.py), [native-master result](https://github.com/zuizui0223/azami/blob/archive/pr142-v4-20261020-preconsolidation/analysis/v4/NATIVE_MASTER665_RESULT_20261020.md) | Multivariate within/between explained variance, block Shapley allocations, RV covariance and geographically held-out diagnostics; **not** a hypervolume or a successful geographically transferable functional model. Existing source-first native/measurement evidence should be reused; V4 Actions remain paused. |
| **EAzami functional disparity** | [seven-taxon direct four-axis PCA](https://github.com/zuizui0223/EAzami/blob/main/data/evidence/fdt4_direct_functional_disparity_pilot_v1.json) and [seven-taxon PGLS](https://github.com/zuizui0223/EAzami/blob/main/data/evidence/fdt4_direct_multitrait_7taxon_pgls_v1.json) | Four axes: measured head width, direct phyllary protrusion, binary orientation, binary flower colour. PCA variance PC1 **0.5193**, PC2 **0.2611**, PC3 **0.1980**; neither head width versus orientation (P~0.893) nor protrusion versus width (P~0.771) was supported in the tiny PGLS. **Source- and lineage-clustered pilot**, not independent functional hypervolume or trait-independent evolutionary modules. |
| **EAzami environmental replication** | [exploratory 8-taxon/63-test environment atlas](https://github.com/zuizui0223/EAzami/blob/main/data/evidence/chapter2_exploratory_reframe_integrated_result_v1.json), [colour–radiation sister-lineage screen](https://github.com/zuizui0223/EAzami/blob/main/data/evidence/chapter2_colour_rsds_focal_concordance_result_v1.json) | **0/63** BH-supported focal East-Asian traits × environment associations; of 63 global-Azami-versus-EAzami directional comparisons, **26 concordant / 37 discordant** (different taxa/estimands; **not** demonstrated opposite adaptation). White/coloured colour–radiation sign is aligned in one of two focal sister comparisons, reversed in the other; small photographic sample, not historical solar driver. |
| **EAzami current Chapter 2 V10** | [recurrent assembly / origin-maintenance separation](https://github.com/zuizui0223/EAzami/blob/main/docs/chapter2/V10_FUNCTIONAL_ASSEMBLY_ORIGIN_DECOUPLING_SPINE.md) | Repeated discrete orientation **4–6** minimum changes, phyllary posture **3**, stickiness **5**; **0/3** trait pairs have robust shared transition localization. Orientation's present BIO15-high/BIO1-low event correspondence cannot identify the historical environmental origin; one calendarized event matches present regime in **99/376** chronology×location cases and **0/324** robust broad event-climate classes. **Historical state-change minima ≠ actual functional-space evolutionary rates.** |
| **aza3** | [same-plant published morphology reanalysis](https://github.com/zuizui0223/aza3/blob/main/docs/CIRSIUM_HEAD_CORE_VS_BRACT_TIP_SPREAD_EMPIRICAL_20261020.md), [functional/genomic LEGO design](https://github.com/zuizui0223/aza3/blob/main/docs/ANTO_HIGHER_ORDER_LEGO_MEASUREMENT_PLAN_V2.md) | In an **existing external source, not the user's own field sample**, 114 non-hybrid plants in three species show flower-corolla length vs involucre length **r=+0.670** but vs relative tip width **r=−0.053** after taxon×sex adjustment. Descriptive phenotypic decoupling only; genetics/fitness/biotic function require **same-individual new validation**. Genome/cytotype/pigment/stickiness histories and functional effects are owned by aza3, not interpolated from Azami image proxies. |

**Hypervolume status:** Across the inspected current and archived code/results, PCA/morphospace, pairwise dispersion, covariance/RV, environmental models and the functional-disparity programme **exist**. A finished native-only multivariate **kernel hypervolume with matched observation/taxon/photographer sampling, overlap and uncertainty** has **not been verified**. Do not call a PCA scatterplot, convex hull, or species-median distance an estimated functional hypervolume. The methodological distinction between ordination, convex hull and density-derived hypervolume is established in [Blonder et al. (2014)](https://doi.org/10.1111/geb.12146) and [Blonder et al. (2018)](https://doi.org/10.1111/2041-210X.12865).

### One bounded next-stage estimand — conditional on private data access

The **native-first** 313,350-source selection is fixed. The only qualified broad image modules presently are **presentation (one image-vertical angle), visible floral colour (L*, C*, one *joint* circular hue), and whole-head silhouette (four perspective-dependent coordinates)**. All four image-colour quantities and four outline quantities remain distinct registered measurements, **not** independent ecological functions. Involucre contour (1,303 eligible observations), surface (0 eligible) and five never-executed metrics are retained as missing/secondary layers and **not** mandatory for this shared space.

First determine the **actual same-observation intersection** for these three core modules from the existing private native measurement table. None of the three marginal totals (170,066 angle; 186,439 lightness/chroma; 176,048 outline) identifies this intersection. The public receipt alone is insufficient to calculate a multivariate hypervolume. Preserve per-measurement QC, the accepted WCVP concept, observer/region/cell and sampling provenance; do not invent values for unmatched photos or unexecuted fractions.

Predefine an outcome-blind scaling and **balanced** reduced phenotype space (e.g., one presentation coordinate, two floral-colour coordinates including coherent circular handling, one or two whole-head shape coordinates; module contributions fixed before environment inspection). Three scientific responses are different and must be retained separately:

1. **Location** — does the centroid of the *visible phenotype space* shift along annual rainfall (BIO12), radiation (RSDS) or another already declared limited environmental axis?
2. **Breadth / occupancy** — do equal-taxon, equal-photo/observer-effort, equal-sample-size corrected within-taxon and among-taxon dispersion / Rao-type pairwise distance / low-dimensional density-derived hypervolume change across environmental regimes? Report confidence/bootstrapped uncertainty and effective coverage. An observed raw increase in volume can otherwise be a simple effect of sampling depth or camera resolution.
3. **Combination restrictions and turnover** — are some image phenotype combinations absent or disproportionately common beyond what the separate marginal distributions and taxon/clade/spatial-photo sampling imply? Separate overlap, shifts and packing from a single scalar volume. Use dimension-wise shuffles or other nulls that preserve the correct taxon / observer / region strata rather than implicitly assuming independent measurements or adaptation.

For environmental comparisons, compare broad taxon turnover and strict **within-taxon same-environment-gradient** results separately. Predeclare a **small fixed set** of nonredundant exposures and one inferential family rather than searching another 27×9 grid. Spatial block / leave-one-region-out skill, taxon-exclusion sensitivity, ascertainment/QC models and photographer/photo clustering determine whether any cross-region interpretation can be promoted. Earlier exploratory native-only geographic transfer was weak/negative for unregularized fits; there is no established universal ecological predictor of functional breadth.

The main evidence ladder remains:

**Azami**: measured *visible phenotype-space geometry* and climate correspondence → **EAzami**: whether the corresponding source-backed taxon concepts underwent independent trait-history shifts, with origin-vs-maintenance explicitly separated → **aza3**: whether linked individuals/populations show different fitness or antagonism/pollination effects and shared vs distinct genetic resources. Only the last step, with independent functional calibration, can justify upgrading a *candidate phenotypic hypervolume* into a validated **functional** hypervolume.

**Stop rule:** This note does not authorize an additional Actions workflow, source download, expanded model family, new claims in the frozen GEB v3 manuscript, or a positive hypervolume result before the native-only multivariate input and QC are actually available.


### Implemented, read-only computational feasibility gate (not yet run on private native data)

[**\`reproducibility/native_space_gate.py\`**](native_space_gate.py) is a single deterministic, **manual-only** extension of the existing analysis handoff. It requires the byte-pinned previously prepared **313,350-observation native-only input**; never downloads images/raster/source, queries a taxonomy API or decrypts a source. It rejects wrong input hashes, duplicate observation IDs, missing native taxon identity, and invalid endpoint-specific QC. No new GitHub Actions workflow is introduced.

It first reports the **actual matched same-observation overlap** of all three modules, using:
- **core3 gate**: image-vertical orientation + photographic chroma + whole-head silhouette aspect ratio;
- **full-core9 audit**: the entire nine available primary image endpoint set, with circular hue sine/cosine retained *jointly*, not treated as two biological traits.

If enough taxon × 0.5-degree spatial cells pass, it generates **exploratory** taxon-balanced 20-cell rarefaction diagnostics: Rao-type RMS pairwise dispersion and a *parametric Gaussian 95% covariance-ellipsoid volume proxy*. Climate reporting is limited to **existing pre-frozen axes BIO12 and RSDS**, and descriptive (no post-hoc significance tests); climatic taxon medians use *all* native records, not a phenotype-selected source subset. Public output is aggregate JSON only, never taxon/photographer/location/observation rows. Cell discretization is not equal-area, so this diagnostic is not a final estimator of a native-only ecological hypervolume.

The source remains an encrypted GitHub artifact inaccessible for local computation without the authorized private key. **Synthetic test success is not an empirical hypervolume result**. Before any ecological promotion, verify photographer/region-dependent missingness, independent colour/angle/outline validity and spatial transfer; a functionally calibrated hypervolume additionally requires experimental phenotype→functional-payoff links owned by \`aza3\`. The current frozen GEB v3 conclusions remain untouched.

### Covariation comes before hypervolume (2026-10-20)

**The first estimand is cross-trait covariance**, not a single volume scalar. The existing native-first manual `native_space_gate.py` now explicitly partitions the **matched three-axis core** (image orientation / flower chroma / whole-head 2D aspect) into the within-accepted-taxon covariance **W**, the among-taxon covariance **B**, and the exact equal-taxon sampling identity **T = W + B**. It reports **three named pairwise correlations at both scales**, their differences, and the share of total variance attributable to taxon means. Comparisons use one observation per taxon × 0.5-degree cell and equally weighted accepted taxa. This helps distinguish *within-taxon reorganization* from *species-assemblage turnover*, without claiming either is genetic covariance.

**Existing actual native-only result (not rerun):** [original observer-balanced exploration](https://github.com/zuizui0223/azami/blob/archive/pr142-v4-20261020-preconsolidation/analysis/v4/OBSERVER_BIAS_RESULT_20261020_JA.md) found, for **164,050 matched observations / 177 accepted taxa**, an orientation × four-dimensional colour RV of **0.006793 within**, **0.037122 among** (observer × spatial-cell balanced). Its matched-n taxon-bootstrap difference was **+0.024622**, 95% interval **[−0.020787, +0.083385]** across 400 draws. This interval includes zero; an unqualified universal stronger-between interpretation is unwarranted. The full V3 manuscript instead reports stronger among integration in the **distinct 1,734-observation / 42-taxon selected cohort**. These are **different estimands, sampling frames and technical measurements**, not direct replication.

The manual gate now additionally compares **among-taxon core3 pairwise correlations** between coarse low/high **predeclared BIO12 (annual rain) and RSDS (shortwave radiation)** groups, reporting *no value* unless **at least 10 accepted taxa per half** can be compared. These descriptive differences have **no p values**, independent phenotypic calibration, phylogenetic control or regional extrapolation. The presence of a taxon-level correlation does not establish a common within-species shift, and two ecological regimes may contain different clades. Photographic orientation/chroma/shape are proxies, not verified seed-fitness functions.

The existing Gaussian ellipsoid proxy then summarizes dispersion **conditional on covariance** and fixed-scale balanced sampling; a larger ellipsoid can result from greater variance, lower absolute pairwise correlation, or both. It is neither a kernel-estimated ecological hypervolume nor proof that occupied morphological combinations are evolutionarily or functionally accessible. Claims of functional differentiation and genetically decoupled LEGO assembly remain independently gated by `EAzami` evolutionary histories and `aza3` same-individual functional/genomic evidence. This update **does not execute** the protected source archive or change frozen V3 results.

