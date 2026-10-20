# Native-first 27-endpoint handoff (2026-10-20)

## Completed aggregate-only Gate A boundary (2026-10-20)

The [machine-readable Gate A bound receipt](native_gate_a_bounds_20261020.json) uses **only** the verified original `measurement_receipt.json`: 217,104 native observations have inherited observation-level measurements, with 170,066 eligible orientation, 186,439 chroma and 176,048 outline-aspect observations from the **same** inherited table. Sharp Fréchet inclusion–exclusion gives orientation × chroma **139,401–170,066**, orientation × outline **129,010–170,066**, chroma × outline **145,383–176,048**, and **matched all-three 98,345–170,066**. These are source-backed **bounds**, **not matched-row results**; joint nine-output lower bound remains zero.

A read-only [manual Gate A audit](native_gate_a.py) now verifies the exact historical prepared-table SHA and frozen source counts, performs **27 endpoint-specific QC checks**, audits the core three modules and their exact observation-level intersections, reports accepted-taxon and spatial-cell counts, and optionally joins private original observer IDs or a verified previously recovered five-field head→observation file. No new photo/raster collection, model fitting, publication of IDs/coordinates or public Actions workflow is involved. The current protected private input is **not present in the public repository**, so the actual intersection, independent spatial/observer support, and recovered-five native counts remain **not executed and not verified**.

## Actual native nine-endpoint observation overlap recovered from earlier aggregate audit

An **existing**, provenance-pinned [native module-overlap aggregate](native_gate_a_existing_overlap_20261020.json) (Actions [39428065666](https://github.com/zuizui0223/azami/actions/runs/39428065666), artifact \`12266812186\`, member \`native_integration_audit.json\`, SHA256 \`39d5613b611cce30be989682b3e47cdb3153aa880a0c2744cd40bb0d72b9edd4\`) contains exact observation/accepted-taxon intersections for source-matched, image-measured *modules*, irrespective of how the separately fitted ecological models behave. Its nine marginal endpoint eligible counts exactly match the previously independent 27-row receipt. These are **existing results reused**, not a new raw-data replay.

| Matched original-image measurement modules | Original native observations | Accepted taxa |
| --- | ---: | ---: |
| Presentation | 170,066 | 244 |
| Joint four-output colour (L*, C*, hue sine/cosine) | 186,390 | 244 |
| Joint four-output 2D head outline | 176,048 | 243 |
| Presentation × joint colour | 164,167 | 243 |
| Presentation × head outline | 159,967 | 242 |
| Colour × head outline | 170,227 | 243 |
| **Presentation × colour × outline (all 9 recorded outputs)** | **157,779** | **242** |

The separately fitted complete-case joint *model* may lose further taxa/observations through modelling gates and is **not** interchangeable with this measurement-only audit. Because colour C* is eligible for 186,439 observations but complete colour4 for 186,390, while the four outline outputs have an identical 176,048-observation eligible census, the simpler orientation+C*+aspect \`core3\` has a valid **157,779 to 157,828** observation range. Its **exact** identity intersection has not been reconstituted from protected row-level values and must not be invented.

The saved native measurement frame has no original \`user_id\`/\`photographer_id\` or direct regional identifiers; although source coordinates are recorded privately, original human effort and population replication cannot be inferred from image counts. Do not declare Gate A entirely passed before exact taxon×cell×observer support, season/image-quality missingness and the historical five-field head-level recovery check are verified. The manual [Gate A script](native_gate_a.py) now validates this prior receipt and tests any later private row-level replay against its exact nine-endpoint count and the bounded core3 overlap. No original identities are published and no previously excluded exploratory model is promoted into the paper.

## Canonical design (historical audit reconciled, design only)

The [single refined study design](NATIVE_CAPITULUM_SPACE_REFINED_DESIGN_20261020.md) is authoritative for future native-first work. The [historical PR/commit index](HISTORY_AUDIT_INDEX_20261020.json) preserves scoped audit metadata. This document retains the source/endpoint evidence receipts and archived analyses **as provenance, not competing execution instructions**. Follow Gate A matched-observation/QC checks before new measurements and Gate B covariance before Gate C density/KDE/Wasserstein. All Gate D environmental and phylogenetic numerical inference is owned by [PR #144](https://github.com/zuizui0223/azami/pull/144) and its [source contract](GATE_D_PR144_HANDOFF_20261021.json); this handoff cannot supersede it. No new original images, source extraction, Actions workflow, or frozen GEB V3 claim is authorized in this PR.

## Decision: what is active and what remains paused

Chapter 1's frozen v3/GEB outputs remain the historical reference; they are **not** automatically native-only by replaying them. The next primary cohort is selected from the **original pre-detector source, prior to any detector or trait selection**, using a fixed native-range definition. An all-upstream native classification and inherited-measurement linkage already exist in the archived PR142 evidence. **Reuse that completed result; do not reacquire WCVP, CHELSA rasters, original images, or remeasure missing pixels in this consolidation.** No new V4 Actions or competing analysis entry point is authorized by this note.

The complete pre-consolidation branch is preserved at [archive/pr142-v4-20261020-preconsolidation](https://github.com/zuizui0223/azami/tree/archive/pr142-v4-20261020-preconsolidation). Archived scripts and results are diagnostic evidence, not part of the currently merged v3/GEB scientific claim set.

## Cohort ledger and fixed ordering

| Stage / denominator | Count | Meaning |
| --- | ---: | --- |
| Original pre-detector metadata, 2026-07-03 | 665,115 observations; 1,122,854 photo links | Starting frame, not 665,115 trait measurements |
| Native-primary after source-level classification | 313,350 observations; 534,186 photo links; 306 accepted WCVP keys | Before phenotype, detector, colour, and trait completeness |
| Legacy raw numeric all-nine (superseded) | 313,240 observations | GSP invalid-code contamination, not valid climate completeness |
| Native with an inherited original measurement-table row | 217,104 observations | Join to previous image measurements, not remeasurement |
| Native lacking any inherited original measurement-table row | 96,246 observations | Preserve as missing; do not silently drop from denominator |
| Historic measurement rows plus raw all-nine numerics | 217,102 observations | Superseded numeric-only overlap; NOT corrected nine-climate analysis cohort |

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

The full native roster has **313,244** observations with the four previously selected climate predictors, versus **312,001** corrected-complete nine-variable observations (**313,240 was old raw numeric-only**); these are predictor-availability denominators, independent of endpoint measurability. Per-endpoint usable counts in the receipt happen to remain unchanged on the four-predictor complete-case gate. A different selected-predictor set or model design can change them and must report its own denominators.

## Next scientific analysis: bounded by existing data

- **The 27-endpoint/native-taxon ledger is now verified.** The [machine-readable inventory](native_first_27_endpoint_coverage_20261020.json) includes native-only finite/eligible/environment-complete counts and accepted-taxon coverage for each endpoint. **Not yet available** from these public aggregates: per-endpoint photographer, region, taxon-by-region and QC failure-reason counts. These must be computed only within an authorized environment from the pre-existing private native table, without another acquisition or image measurement.
- **One response-specific analysis rather than an all-traits complete-case analysis.** Use the nine existing core photographic endpoints first, preserving the circular joint hue unit. Define response-specific samples, taxon/observer balance, within- versus among-taxon associations, spatial holdouts and multiplicity before inspecting results. Treat shape and orientation as view-dependent visible phenotypes.
- **Separate supporting audits.** Involucre contour/texture/specularity can describe *image signals* only until expert landmark, scale, lighting and surface validations exist. Keep them out of unqualified ecological-function claims.
- **Audit earlier head-level measurements before contemplating new pixels:** PR #92's historical native spatial receipt reports that visible-floret fraction and four exclusive-colour fractions were **already recovered for 23,881 observations per endpoint in a separate, older 27,066-native-observation spatial cohort**, although the frozen current 27-endpoint atlas has zero eligible observations for them. Their actual native-observation intersection, head-level provenance and technical eligibility remain **unverified**. Do not infer these fractions from Lab medians or add images unless a later, separately validated acquisition gate justifies it.
- **Do not overwrite the v3 manuscript.** The prior mixed-range V3 baseline stays frozen; a native-only result is a new, labelled comparison only after exact sample counts, uncertainty, dominant-taxon sensitivity, observer effects and geographic transfer have been checked.

## Operational gate

The native-primary **source selection is already completed** in archived PR142 and need not be reacquired. The **27-endpoint aggregate coverage ledger is completed and checked in PR142**, but the **all-nine-primary photo endpoint matched intersection is now recovered from an existing aggregate receipt (157,779 native observations / 242 taxa)**. The representative core3 triple has a refined mathematical bound of 157,779–157,828 observations, while independent photographer/region missingness, actual recovered five-head fields, native KDE/Wasserstein and ecological validation **remain incomplete**. The current [V3 numerical replay](CURRENT_ANALYSIS.md) validates old frozen inputs and 16 reference outputs; it **does not** implement a native-only reanalysis by itself. Do not equate green staging/integrity CI with execution of this handoff. The next executable step is reuse of the archived, checksum-verified 313,350-observation native-only prepared table to generate the 27-row eligibility ledger, without downloading or remeasuring images.


## Retired exploratory environment-model notes

Previously fitted exploratory models remain in the Git archive, but they are not adopted as independent ecological confirmation or as a module-selection rule. Future ecology must pass source integrity, trait-specific eligibility and a frozen sampling/estimator gate.

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

### Superseded interim multivariate methods

The interim ellipsoid feasibility description is now superseded by the [canonical Gate A–E design](NATIVE_CAPITULUM_SPACE_REFINED_DESIGN_20261020.md). The original analyses remain recoverable; density overlap and Wasserstein are not interchangeable with Gaussian ellipsoid volume.

### Implemented, read-only computational feasibility gate (not yet run on private native data)

[**\`reproducibility/native_space_gate.py\`**](native_space_gate.py) is a single deterministic, **manual-only** extension of the existing analysis handoff. It requires the byte-pinned previously prepared **313,350-observation native-only input**; never downloads images/raster/source, queries a taxonomy API or decrypts a source. It rejects wrong input hashes, duplicate observation IDs, missing native taxon identity, and invalid endpoint-specific QC. No new GitHub Actions workflow is introduced.

It first reports the **actual matched same-observation overlap** of all three modules, using:
- **core3 gate**: image-vertical orientation + photographic chroma + whole-head silhouette aspect ratio;
- **full-core9 audit**: the entire nine available primary image endpoint set, with circular hue sine/cosine retained *jointly*, not treated as two biological traits.

If enough taxon × 0.5-degree spatial cells pass, it generates **exploratory** taxon-balanced 20-cell rarefaction diagnostics: Rao-type RMS pairwise dispersion and a *parametric Gaussian 95% covariance-ellipsoid volume proxy*. Climate reporting is limited to **existing pre-frozen axes BIO12 and RSDS**, and descriptive (no post-hoc significance tests); climatic taxon medians use *all* native records, not a phenotype-selected source subset. Public output is aggregate JSON only, never taxon/photographer/location/observation rows. Cell discretization is not equal-area, so this diagnostic is not a final estimator of a native-only ecological hypervolume.

The source remains private/encrypted; an authorized, explicitly source-verified native observation-level table is required before any manual feasibility calculation. **Synthetic test success is not an empirical hypervolume result**. Before any ecological promotion, verify photographer/region-dependent missingness, independent colour/angle/outline validity and spatial transfer; a functionally calibrated hypervolume additionally requires experimental phenotype→functional-payoff links owned by \`aza3\`. The current frozen GEB v3 conclusions remain untouched.

### Covariation comes before hypervolume (2026-10-20)

**Already measured pairwise construct examples in the frozen V3 *pairwise* cohort** ([canonical result](../analysis/v3/SCALE_DEPENDENT_INTEGRATION_RESULT_20260910.md); each pair has its own qualifying cohort):

| Construct pairing | Within-taxon RV | Among-taxon RV | Difference (among − within) |
| --- | ---: | ---: | ---: |
| Floral chroma × floral hue | 0.0626 | 0.3577 | +0.2951 |
| Floral chroma × visible projection prominence | 0.00002 | 0.2491 | +0.2491 |
| Floral lightness × floral chroma | 0.0981 | 0.2127 | +0.1146 |
| Projection prominence × projection pattern | 0.3038 | 0.1544 | −0.1494 |

These are *RV association magnitudes* (not signed Pearson r, adaptation effects, or percentages). Hue is circular and partially constructed from the same colour pixels, so its association with chroma cannot be interpreted as independence of ecological function. Projection variables arise from related outlines/segmentation; technical covariance may contribute. The first three strong **among** connections and the fourth stronger **within** connection illustrate why a single hypervolume number would hide scale-dependent shape. This pairwise table is a *different cohort/estimator* from the full V3 1,734/42 exact-common-cohort contrast. In the exact-common-cohort result, **33/36** relations are stronger among taxa and the bootstrap median RV contrast is **+0.065928 [0.038197, 0.119048]**, but it must **not** be transferred unchanged to the 313,350-native source frame.

**The first estimand is cross-trait covariance**, not a single volume scalar. The existing native-first manual `native_space_gate.py` now explicitly partitions the **matched three-axis core** (image orientation / flower chroma / whole-head 2D aspect) into the within-accepted-taxon covariance **W**, the among-taxon covariance **B**, and the exact equal-taxon sampling identity **T = W + B**. It reports **three named pairwise correlations at both scales**, their differences, and the share of total variance attributable to taxon means. Comparisons use one observation per taxon × 0.5-degree cell and equally weighted accepted taxa. This helps distinguish *within-taxon reorganization* from *species-assemblage turnover*, without claiming either is genetic covariance.

**Existing actual native-only result (not rerun):** [original observer-balanced exploration](https://github.com/zuizui0223/azami/blob/archive/pr142-v4-20261020-preconsolidation/analysis/v4/OBSERVER_BIAS_RESULT_20261020_JA.md) found, for **164,050 matched observations / 177 accepted taxa**, an orientation × four-dimensional colour RV of **0.006793 within**, **0.037122 among** (observer × spatial-cell balanced). Its matched-n taxon-bootstrap difference was **+0.024622**, 95% interval **[−0.020787, +0.083385]** across 400 draws. This interval includes zero; an unqualified universal stronger-between interpretation is unwarranted. The full V3 manuscript instead reports stronger among integration in the **distinct 1,734-observation / 42-taxon selected cohort**. These are **different estimands, sampling frames and technical measurements**, not direct replication.

The manual gate now additionally compares **among-taxon core3 pairwise correlations** between coarse low/high **predeclared BIO12 (annual rain) and RSDS (shortwave radiation)** groups, reporting *no value* unless **at least 10 accepted taxa per half** can be compared. These descriptive differences have **no p values**, independent phenotypic calibration, phylogenetic control or regional extrapolation. The presence of a taxon-level correlation does not establish a common within-species shift, and two ecological regimes may contain different clades. Photographic orientation/chroma/shape are proxies, not verified seed-fitness functions.

The existing Gaussian ellipsoid proxy then summarizes dispersion **conditional on covariance** and fixed-scale balanced sampling; a larger ellipsoid can result from greater variance, lower absolute pairwise correlation, or both. It is neither a kernel-estimated ecological hypervolume nor proof that occupied morphological combinations are evolutionarily or functionally accessible. Claims of functional differentiation and genetically decoupled LEGO assembly remain independently gated by `EAzami` evolutionary histories and `aza3` same-individual functional/genomic evidence. This update **does not execute** the protected source archive or change frozen V3 results.

