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
