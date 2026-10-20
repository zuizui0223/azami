# v4 native-first FULL upstream analysis — 2026-10-20

## Scope and executed work

The user requested native-range restriction **before** selecting images or estimating phenotypes. Classification began from all **460,036 observations / 288 source taxon names** in the 777,766-photo upstream queue, not from v3's 46,276-observation subset. The final replay reclassified this complete roster, joined the archived measurements only afterwards, and reproduced the results below from clean CSV inputs.

- WCVP source snapshot: artifact **12241384633**, run **39373263671**, archive SHA256 `4924a2f3064ecd407f370d3fe77df1167ca4aecd9045477cef740b3ba6530fd7`.
- WCVP source JSON SHA256: `d3bf06b5237a371094ed4c52322d270fb653b689056886bbd62145105ec3d046`.
- All 288 source names queried; 272 uniquely resolved, 265 accepted keys before geographic restrictions.
- Native means an existing WCVP native distribution entry in the observation's uniquely assigned TDWG level-3 polygon. This is a **regional native-range filter**, not proof of indigenous ancestry for each photographed individual.
- Introduced, doubtful, unlisted, unresolved and unmapped records are never rescued as native. The WCVP-specific blank-establishment convention is applied only to an existing distribution row, never to an absent row.
- Coordinate usability, obscuration and captive/cultivated flags are separate gates. The frozen upstream captive flags are all false; this does not independently validate that every photograph is wild.

## Cohort accounting

| Step | Observations | Taxa |
|---|---:|---:|
| ALL upstream observations | 460,036 | 288 source labels |
| In a documented native region, before coordinate/flag restrictions | 253,609 | — |
| Eligible native-only primary roster | **244,570** | **256 source labels / 249 accepted keys** |
| Native observations with archived detector-positive measurement rows | **217,104** | 247 accepted keys |
| Same-record presentation + colour comparison, min 5 records per taxon | **164,052** | **177 accepted keys** |
| Presentation environment model | **169,950** | **178** |
| Colour environment model | **186,264** | **177** |

Range-status counts before eligibility: native 253,609; introduced 137,907; unresolved taxon 35,923; unmapped TDWG 25,947; unlisted 6,198; invalid coordinates 452. These sum to 460,036. The eligibility exclusions are recorded separately rather than conflating invalid coordinates with introduced ranges.

The old 46,276-record native-status table was used only for validation: all **46,276 classifications agreed**. It was NOT the source population for the new analysis.

## Full climate coverage is now implemented, not merely planned

All nine original CHELSA v2.1 1981–2010 predictors were freshly acquired from the source rasters. Run **39374865573** completed all nine raster jobs. The 316 occupied 5-degree cells are only **storage groups**: original raster pixels, resolution, nodata and scale/offset are retained. There is no 5-degree averaging or nearest-observation imputation.

- All native targets: **244,570**.
- Complete nine-predictor records: **244,568 (99.9992%)**.
- Remaining original raster nodata: **2 observations, wind only**; not filled from neighbours.
- Exact agreement at every finite old/native overlap: **27,066/27,066 per variable**, except wind **27,064/27,064** because the same two values are nodata.
- No fallback to the old sparse environment table was used.

## Models and comparisons

Source synonyms are collapsed to accepted WCVP keys. Within each eligible model cohort, every taxon receives equal total weight and its occupied equal-area longitude/sine-latitude cells receive equal total weight. Trait and environmental aggregation use **the same observations**. Spatial terms are an eight-column spherical basis, not raw latitude/longitude only. Phylogenetic placement remains a separate sensitivity, not the primary model.

Presentation uses the archived image-vertical angle. Colour is a multivariate L*, C*, hue-sine, hue-cosine response; hue is jointly scaled as one circular component in environmental models. This does not recover calibrated visible/UV reflectance or anthocyanin concentration.

Environmental blocks: radiation (1 column), thermal (2), water-demand (4: annual rain, rain seasonality, VPD and growing-season rain), wind (1), productivity (1). Block sizes and correlations differ. Shapley shares distribute shared **in-sample** explanatory information; unique drop-block increments condition on other environment columns and space. Neither statistic ranks independent causes.

### Joint environmental organization

| Quantity | Presentation | Visible colour |
|---|---:|---:|
| Among-taxon environment-only trace R² | **0.1391** | **0.2345** |
| Among-taxon environment + space R² | 0.2339 | 0.3370 |
| Unique environment given space | 0.0753 | 0.1221 |
| Shared environment/space component | 0.0638 | 0.1124 |
| Within-taxon environment-only trace R² | **0.00438** | **0.00486** |
| Among-taxon held-out predictive R² | **0.0198** | **0.1381** |
| Within-taxon held-out predictive R² | **−0.00243** | **−0.00178** |

Among-taxon validation holds out complete 10-degree **taxon-centroid** blocks; taxon ranges may still overlap folds. Within-taxon validation holds out actual 10-degree observation-location blocks, uses training-only taxon means and scores only taxa adequately represented in training. It tests new locations of known taxa, not new taxa. Scored records: presentation 164,679 (5,271 unseen/sparse-taxon test rows excluded), colour 180,765 (5,499 excluded). Neither validation uses a geographic buffer. Ridge alpha=0.1 is fixed, not chosen by test performance. Negative CV R² is retained.

Water-demand receives the largest in-sample Shapley share in both modules (0.0754 presentation; 0.0690 colour). That does **not** translate into one universal strongest driver. Among-taxon predictive drop-block differences favour wind for presentation (+0.0622 relative to omitting wind) and radiation for colour (+0.0527). The overall predictive performance of presentation is nevertheless weak. Environmental organization is much more evident among taxa than within taxa in this dataset; this is not proof of plasticity or adaptation.

### Recheck of the two historical ecological patterns

Exploratory standardized species-level coefficients, using native-only same-record, equal-area-weighted taxon means; HC3 intervals. These are not re-labelled frozen v3 results and do not use exactly v3's median estimand.

| Association | Marginal beta | Joint environment beta | Environment + space beta (95% interval) |
|---|---:|---:|---:|
| Chroma × radiation | **−0.3361** | −0.4757 | **−0.2108 (−0.8040, +0.3824)** |
| Presentation × annual precipitation | **+0.2708** | +0.2692 | **−0.1449 (−0.4663, +0.1765)** |

Removing the two most-observed taxa gives marginal coefficients −0.3390 and +0.2702, respectively. Removing the top five gives −0.3435 and +0.2709. Thus the marginal native-range patterns are not caused solely by *C. arvense* and *C. vulgare*. Joint spatial/environment models do not isolate a stable independent effect of either historical predictor. Report shared spatial/environment information rather than interpreting failure to separate it as proof of no association.

### Within/among presentation–colour covariance

Same native cohort: **164,052 observations / 177 accepted taxa**.

- Raw component-scaled RV: within **0.0062396**, among **0.0367894**.
- Adding taxon resampling to matched-n resampling: median difference **+0.02569**, 95% interval **−0.02152 to +0.08927**, positive fraction **0.867** (1,000 replicates).
- Joint-hue rotationally invariant sensitivity: median difference **+0.01830**, 95% interval **−0.03012 to +0.08051**, positive fraction **0.787**.
- Minimum per-taxon counts 5/10/20/50 are all retained as sensitivities, not selected by their outcome. All nested 95% intervals include zero.

These are bootstrap-plus-matched-n resampling intervals, not permutation P values. The native-wide data do not support promoting v3's small-cohort result into a general claim of stronger among-taxon integration. A tendency remains, with substantial taxon uncertainty.

## Involucre and remeasurement boundary

Only **1,303 native observations** currently have archived involucre descriptors; the common three-module comparison has **949 observations / 37 taxa**. It is a legacy-selected bridge, not the native-wide answer. A native-only high-resolution remeasurement queue of **110,515 observations** is saved. The full native-range analysis here uses archived image measurements; **raw photographs remeasured in this numerical run = 0**. Do not conflate concurrent original-photo pilot work with full-queue completion.

## Reproduction and verification

Canonical numerical command:

```bash
python -m analysis.v4.run_native_analysis \
  --roster data/source_observations.csv.gz \
  --measurements data/exhaustive_continuous_observation_level.csv \
  --evidence data/native_sources.json \
  --environment data/native_full_environment.csv.gz \
  --out-dir work/native-replay --replicates 1000
```

A clean local replay **PASS** reproduced classification denominators, all eight integration sensitivities and both multivariate environment outputs. New local climate/statistics tests: **8 passed** before repository integration. GitHub source/raster jobs validate actual remote source acquisition; general CI validates code. Do not describe those jobs as a full raw-photograph remeasurement.

Frozen numerical inputs:

- source observation roster gzip SHA256 `a05709396bf5bba642b561f2781c6806fbe2d74fac8876b4a5398e047e7d51e3`;
- all upstream observation-level measurement CSV SHA256 `b49a093c849928c456d6747ec51772b5a6b78debb2ae5eb5dbe2a91b08f66dd7`;
- native source JSON SHA256 `d3bf06b5237a371094ed4c52322d270fb653b689056886bbd62145105ec3d046`;
- complete native climate table gzip SHA256 `51304f84fc2d1674c3b3eb6f64fc1bd7bfead02baea45b07863a2e8b89bf77c1`.

This v4 exploration leaves main, v3 results and GEB claims unchanged. All new analyses are retrospective, not a new confirmatory study.
