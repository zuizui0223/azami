# Chapter 1 v4: function-led exploratory modules

**Exploratory only. Do NOT replace frozen v3 / GEB results.** No original photos are re-downloaded or reanalysed.

## Three biological questions, not 22 arbitrary endpoints

| Module | Existing image observables | Hypothesized function | Not proven |
|---|---|---|---|
| Presentation | Head principal-axis angle vs **image** vertical | Rain/wind exposure, display to pollinators | Gravity-referenced head or peduncle angle, protection |
| Visible flower colour | L*, C*, hue sin/cos; approximate Lab a*, b* reconstructed from C* and hue | Pollinator signal, light/thermal environment | Anthocyanin, calibrated reflectance, UV response |
| Involucre | Length/width, basal/apical taper, projection p95 & peak density | Protection and physical architecture | Actual spine number/length, gland/mucilage, herbivory defence |

Do not promote projection silhouette roughness into botanical spine length. Approximate Lab a*/b* recombine observation medians: this is not the median of calibrated pixel-level Lab, and low-chroma hue requires its own reliability gate. The next measurement version should report *actual visible colour distributions* (white/cream, pink, violet etc. as probability + continuous Lab), separating flower from bracts. Mucilage and indumentum need macro reference data, not shiny pixel proxies.

## Implemented and executed pilot

Executable: analysis/v4/functional_pilot.py

Inputs match the exact current v3 numerical input archives:

- frozen traits Actions artifact 9612943217; ZIP SHA256 101e996b638996a0c5ae79d358bf51293c3585f0e84c4a961b91dcbedf96211e
- nine-predictor environment artifact 9633419268; ZIP SHA256 d7c0c466f55b67695d06ae46c21a6452dbe6cfd92a52db8042caa200429e97f4
- paths inside ZIP: universe/continuous_trait_universe_observation_long.csv and process_environment/strict_spatial_chelsa_process.csv

Run: python analysis/v4/functional_pilot.py --traits PATH_TO_TRAITS_CSV --environment PATH_TO_ENV_CSV --out-dir work/v4-functional

Pilot calculations:

1. Influence of the most observed taxa: observation-weighted and equal-total-weight-per-taxon demeaned chroma/radiation coefficient; remove Cirsium vulgare, C. arvense, top two, top five.
2. Common complete three-module cohort, at least five observations per source-assigned taxon. **Not the old complete-18 cohort.**
3. For each module among taxa, a full joint *multivariate* model with five environmental blocks: radiation (1 predictor), thermal (2), water demand (4: rainfall/seasonality/VPD/GSP), wind (1), productivity (1). Report in-sample R², order-independent Shapley R² allocation, unique drop-one-block ΔR². Repeat with latitude and longitude in the baseline. Correlated predictors share explained variation; allocated shares are *not independent process effects*. Block size and 42-taxon overfitting matter.
4. **Pooled within-taxon multivariate environment contrast on the same complete cohort:** demean response and environment within each source taxon; use equal total weight per taxon; report Shapley and unique conditional R² in the same five exploratory environmental blocks. These cross-sectional slopes are not plasticity, and their in-sample R² are not directly comparable as effects with species-median R².
5. Within-taxon vs among-taxon three-module RV. Because raw RV baseline depends on n, compare against 300 repeats of one taxon-centred observation per taxon. This is sensitivity, not a bootstrap p-value.

### Observed 2026-10-20 numerical pilot

- Broad environment cohort: 46,276 observations, 259 source taxa.
- C. vulgare: **13,207** observations; C. arvense: **11,795**. Together **54.03%**.
- New three-module complete cohort: **1,777 observations, 42 taxa**. The different 1,777 versus v3's 1,734 arises from different endpoint eligibility.
- **Species-median v3 chroma × radiation:** baseline standardized β **−0.345372** on 143 taxa, **−0.342726** excluding C. vulgare, **−0.352802** excluding C. arvense, **−0.350163** excluding both, **−0.336055** excluding the top five. Thus the observed broad taxon-median negative association is **not** carried by those dominant species alone; this does not resolve native-range or spatial confounding.
- Within-taxon marginal radiation×chroma standardized β: **+0.0010** observation-weighted; **−0.0169** under equal taxon weight. Remove both dominant taxa: **−0.0058** observation-weighted. All are very small exploratory associations.
- Among-taxon in-sample multivariate R²: **0.262** presentation, **0.253** colour, **0.252** involucre. Water-demand block has the greatest Shapley allocation across all three modules, including the linear lat/lon sensitivity. This is **not evidence that water has the largest independent causal effect**; unequal predictor count and ecological/spatial covariance are unresolved.
- Taxon-balanced within-taxon **in-sample** five-block multivariate R²: presentation **0.0391**, colour **0.0175**, involucre **0.0326**. The strongest allocated block is water demand for presentation/involucre, but **wind** for colour; the small R² and collinearity preclude a mechanistic ranking. Both models use the exact new 1,777/42 complete cohort and different within/between estimands.
- Presentation × colour RV within / among: **0.030 / 0.211**; matched-n median difference **+0.172**, positive in **99.3%** of resamples.
- Presentation × involucre: **0.002 / 0.043**; matched-n difference **−0.001**, positive **47.0%**. No stable direction.
- Colour × involucre: **0.010 / 0.107**; matched-n difference **+0.037**, positive **85.3%**. Uncertain.
- The three pairs are **not** a universal law of integration. The result is more biologically specific than the raw v3 33/36 descriptive count.

This pilot uses **source-assigned names**, not accepted-name/lineage proofs. Lat/lon linear controls are NOT a spatial model, and neither native-only nor phylogenetic results are part of this new pilot.

## v4 function-led hypotheses (independent data required to test mechanisms)

**H1 Rain × presentation:** A measured gravity/flower-peduncle inclination predicts wetting/flower protection after accounting for opening stage, height and position. Test paired photography + inclinometer plus rain exposure.

**H2 Colour:** Image C* declines toward high solar radiation at the species level because taxon or regional sorting changes visible colour, not necessarily because sunlight suppresses anthocyanin within flowers. Test colour-chart/RAW photography, spectra (visible+UV), pigment amounts, native-only stratification and temperature × radiation.

**H3 Involucre:** Anatomically measured free bract tips, curvature, gland/sticky material and spine lengths can replace image contour proxies only after multi-view/macro landmark calibration. Link to herbivore damage, wetting and mechanical exposure with new field data.

**H4 Environment importance:** Report shared Shapley, unique conditional ΔR² and spatial-block held-out predictive ΔR² together. Do not rank independent *causes* from correlated cross-sectional predictors. Next fit regularized equal-complexity block models and uncertainty-aware within/between hierarchical models; label non-identifiable cases.

**H5 Scale covariance:** Under equal taxa/sample-size comparisons and taxon-bootstrap uncertainty, identify which *functional module pairs*, not all endpoint pairs, maintain different integration across within/among scales. Include native-only/WCVP resolution and phylogenetic placements **as sensitivity analyses**, not as the primary identification strategy.

## Next gates, not yet completed

1. Independently annotate a stratified physical-measurement gold set: actual head/stem orientation, colour standards and floral mask, anatomical phyllaries and glands; benchmark against image metrics.
2. Distinguish annual climate from flowering-season exposures and temperature/radiation/VPD interactions; score collinearity and partial identifiability.
3. Run spatial-blocked predictions and regional omission with taxon-level uncertainty; compare the leading environmental blocks using cross-validation.
4. Repeat leverage and integration on *native-only* and accepted-taxonomy strata; handle observation, photo and photographer pseudoreplication.
5. Treat 52-tree historical placement as a phylogenetic **sensitivity**, not 52 independent replications or a substitute for reticulate phylogenomics.

### Evidence rule

New models were designed **after** seeing the v3 results. Numerical CI validates code and data linkage; neither CI nor a low regression P establishes measurement accuracy or adaptation. v3 frozen claims and manuscript remain unchanged.