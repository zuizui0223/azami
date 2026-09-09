# Chapter 1 v3 — shortest reuse route

The active v3 route is deliberately small. It repairs upstream provenance and aggregation, then reuses the frozen v2 ecological machinery rather than rebuilding a second analysis system.

## Current path

`frozen source provenance -> frozen strict-spatial cohort -> native-range restriction -> historical measurement reuse -> phenotype-blind VIFstep<=10 predictor filter -> frozen v2 within/among estimators -> existing broad-spatial sensitivity -> existing 52-tree historical-placement sensitivity`

### Cohort

- Frozen v2 strict-spatial universe: **46,276 observations / 259 source-assigned taxa**.
- Current native-only reanalysis: **27,066 observations / 241 taxa with environment data**.
- Native membership is deterministically regenerated from frozen v2 auxiliary tables because the historical native-status Git LFS object is no longer byte-recoverable.

### Traits

- All **27 registered endpoints** are retained.
- Hue sine/cosine remain one joint circular inferential unit, giving **26 inferential units**.
- Five fields already measured historically at head level but omitted from the old observation aggregation are restored; no new image measurement is performed.
- Closed colour fractions remain descriptive and are not narrated as independent primary biological discoveries.

### Within and among scales

This distinction is inherited from v2, not introduced by v3.

- **Within taxon:** taxon-demeaned standardized marginal slope, taxon-clustered uncertainty; circular traits use the registered joint test and within-taxon predictor permutations.
- **Among taxon:** taxon-level trait medians and environment medians, standardized marginal slope, predictor permutation across taxa; min5 primary and min2 sensitivity scopes are retained.
- BH-FDR families remain separate for within, among-min5 and among-min2.

### Environment

The frozen v2 nine-predictor atlas is retained as the historical benchmark. A user-requested comparative lane applies phenotype-blind VIFstep<=10 before rerunning the same v2 marginal models.

For the native cohort, VIFstep removes **BIO1 only** at both scales. The retained eight predictors are:

- BIO4
- BIO12
- BIO15
- shortwave radiation (`rsds`)
- VPD
- surface wind
- growing-season precipitation (GSP)
- NPP

Final maximum VIF is **1.961** for within-taxon demeaned environment and **6.977** for among-taxon median environment. VIF is used only to filter the family of marginal tests; it does not convert the analysis to a simultaneous multivariable regression.

## Current shortest-route result

VIFstep<=10 native reanalysis:

- within FDR: **37**
- among-min5 FDR: **27**
- among-min2 FDR: **16**
- within broad-spatial passes: **6**
- among broad-spatial passes: **4**
- among pairs passing all 52 historical-placement trees: **4/4**

The two primary among-taxon survivors are:

1. `corolla_lab_chroma x chelsa_rsds_mean` — negative association.
2. `corolla_lab_chroma x chelsa_npp` — positive association.

Two additional closed-composition descriptive rows also survive:

- `corolla_yellow_pixel_fraction x chelsa_gsp`
- `corolla_yellow_pixel_fraction x chelsa_npp`

`orientation_image_vertical_angle x chelsa_bio12` remains in the retained predictor family but does **not** pass the existing broad-spatial sensitivity because residual Moran structure remains detectable. VIF filtering does not rescue it. `orientation x BIO1` is not tested in this lane because BIO1 is the phenotype-blind VIFstep exclusion.

## Receipts

- `native_primary_spatial_historical_receipt_20260909.json` — native-only nine-predictor benchmark.
- `native_primary_vifstep10_receipt_20260909.json` — current VIFstep<=10 comparative result and exact workflow/artifact provenance.

## What v3 no longer requires

The following are preserved only on `archive/pr92-preclean-20260909` and do not block the current route:

- Wave B/C image remeasurement;
- independent assessability research by taxon/region/environment;
- alternative uniform-chroma qualification;
- tail/repair/full-family/final calibration branches;
- new 36-slot calibration expansions;
- hypervolume/breadth synthesis;
- a new 319,244-observation primary ecological cohort.

## Boundary

This remains a pattern-first observational analysis. Broad-spatial sensitivity is not a full SPDE process model, the 52 trees are historical-placement sensitivity rather than a resolved *Cirsium* species tree, and no surviving association establishes causation, adaptation, selection or mechanism.

The principle is: **reuse v2 wherever it is still valid, patch only the parts that v2 could not support cleanly.**
