# V4 native-first 665,115-observation MASTER execution — 2026-10-20

## Status and boundaries

The complete historical pre-detector metadata snapshot was classified **before** joining any phenotype, detector output or climate value. The full classification and analyses below ran locally on recovered checksum-verified inputs. GitHub source-completion run **39375870888** passed and supplied the all-350-name authority snapshot. CI tests numerical/source contracts; it does not establish botanical accuracy.

This is a retrospective exploratory lane. Main, the frozen v3 results and the GEB manuscript are unchanged. **No new pixel/anatomical measurements are claimed by this execution.** The ecological calculations are explicitly labelled a diagnostic reuse of inherited upstream measurements, not a new v4 measurement result.

## Relationship to the concurrent 460,036-observation result

`NATIVE_FULL_RESULT_20261020.md`, added separately while this master-frame run was executing, uses the older 460,036-observation / 288-source-label frame and reports 244,570 native observations. It is preserved unchanged. **This file instead starts at the earlier 665,115-observation / all-350-species pre-detector metadata snapshot.** Do not combine the denominators or transplant environmental centroids between them. The larger native master has 313,350 eligible observations; it is the frame used here to satisfy the all-upstream-first request. The inherited phenotype coverage overlaps strongly, but native exposure medians and estimators differ. Consequently the coefficients and geographical-holdout scores are separate exploratory estimates, not contradictory reruns of one frozen model.

## 1. Full-source native-first cohort

Source artifact **8066010557**, snapshot **2026-07-03**, was recovered from its byte-identical owner archive. It is not the current entire iNaturalist database.

| Stage | Observations | Meaning |
|---|---:|---|
| Entire pre-detector source | 665,115 | 1,122,854 photo rows, all ranks retained in audit |
| Source species rank | 580,066 | 350 source species names queried; no genus/subspecies coercion |
| Native primary | **313,350** | Noncaptive, usable non-obscured coordinates, WCVP TDWG3 native |
| All nine environmental variables available | **313,240** | **99.9649%** of native primary, no imputation |
| Legacy upstream measurements joined | 217,104 | No selection on the frozen 46,276-observation v3 subset |
| Legacy measurements plus all nine environmental variables | 217,102 | Diagnostic input only |

The native primary includes **313 source species labels / 306 accepted WCVP keys** and retains **534,186 photo rows**. Native inclusion never uses photo licensing, colour, detector success or trait completeness. Licensing remains a separate acquisition/redistribution constraint: a retained photo row is not permission to reuse its bytes.

The earlier authority roster covered 288 source names; another 62 exact source species were queried. All 350 species names are now represented in the source snapshot, including unresolved names as unresolved. Exact-name records must converge on one accepted WCVP key. Analyses collapse synonymous labels via that key.

WCVP's native status is at **TDWG level 3**, not proof that a particular local population is indigenous. An existing WCVP distribution row with blank or NATIVE establishment means is native under that source-specific convention. A missing row is unlisted, never native. Boundaries, ambiguous locations, unresolved names, captive observations and private/obscured coordinates are not promoted into the primary cohort.

Every source observation has an audit status. Ineligible rows remain in the audit; private or unusable coordinates are redacted. Non-species ranks account for **85,049 observations**, not hidden losses.

The two most observed native taxa are C. arvense (78,646) and C. vulgare (62,368): **45.0021%** jointly. They are not automatically removed; equal-taxon weights and explicit omissions are separate diagnostics.

## 2. Environment extraction on the native master frame

The nine CHELSA V2.1 1981–2010 raster caches from run **39374865573** were downloaded and checksum-verified. Source resolution is **0.0083333333 degrees**. Five-degree cells are storage partitions only; no five-degree environmental averaging occurs. Stored scale/offset values are applied explicitly. Missing unit labels are retained as missing rather than invented.

**313,240 / 313,350 observations have all nine variables.** Of the 110 incomplete rows, 106 lie outside the earlier cache footprint and four additional rows have wind nodata within it. They remain listed in native_environment_gaps.csv. The unfilled footprint and uncertainty in source coordinate accuracy remain limitations.

## 3. Native-only ecological diagnostics using inherited measurements

The exact complete observation-level member from the historical 406,582-observation measurement archive was recovered and SHA-256 verified. Matching native observation IDs and source taxon labels were checked before fitting. This does **not** substitute for the new anatomical/pixel measurement workflow.

### Marginal among-taxon anchors

Each accepted taxon requires at least five usable phenotype observations. Phenotype medians are related to environmental medians from the full eligible native frame. The intervals below resample taxa 2,000 times; they do not simultaneously control space or phylogeny.

| Association | Measured observations | Accepted taxa | Standardized beta | Taxon-bootstrap 95% interval | Beta excluding C. vulgare and C. arvense |
|---|---:|---:|---:|---:|---:|
| Image chroma × radiation | 186,313 | 177 | **-0.32516** | -0.46213 to -0.17784 | **-0.33484** |
| Image presentation angle × annual precipitation | 169,950 | 178 | **+0.22881** | +0.10482 to +0.34943 | **+0.22901** |

The native-range marginal patterns are not driven solely by those two widespread taxa. These are cross-sectional geographical associations, not causal light/rain responses, pigment concentrations or validated gravity angles.

### Joint environment, space and regional transfer

The same nine predictors are grouped into radiation, thermal conditions, water-related conditions, wind and productivity. Outcomes are either presentation angle or a joint colour representation (L*, C*, hue sine/cosine). Each hue coordinate has half a response weight so the circular hue unit does not count as two complete traits. Within-taxon models demean outcomes and predictors and assign equal total weight to each accepted taxon. Among-taxon models use medians. Eight spherical-coordinate terms describe broad space; phylogeny is not forced into this primary exploratory model.

| Module | Among-taxon environmental R² | Unique environment after broad space | Shared environment–space R² | Within-taxon environmental R² | Geographical holdout skill, environment-only |
|---|---:|---:|---:|---:|---:|
| Presentation | 0.1196 | 0.0899 | 0.0298 | **0.0026** | **-0.1416** |
| Colour | 0.2315 | 0.1307 | 0.1009 | **0.0032** | **-0.0138** |

Shapley allocation conditional on broad space assigns the largest in-sample share to wind for presentation and the water-related block for colour. These are not established causal rankings: blocks differ in size, share environmental information, and the model transfers poorly to held-out geographical clusters.

Six spherical-geography clusters of taxon native centroids define the holdout folds. Scaling is trained within each fold; skill is relative to predicting the held-out taxa with the training-taxon mean. **Negative skill means the tested model did not improve on that baseline.** The joint environment-plus-space model also has negative skill (-0.4506 for presentation; -0.2644 for colour). Adding covariates is not by itself methodological progress. These unregularized linear diagnostics motivate, but do not replace, a separately evaluated predictive model.

### Within–among presentation–colour covariance

The common native cohort contains **164,050 observations / 177 accepted taxa**. Among-taxon RV is **0.03398**. Drawing one centred observation per taxon to match the row counts, among-minus-within RV has median **+0.02099**, conditional 95% range **-0.01698 to +0.03251**, and is positive in **91.0%** of 1,000 draws.

Taxa and the among-taxon medians remain fixed in this diagnostic. This range is **not a taxon-bootstrap confidence interval or a P value**. Because it crosses zero, the larger native-only cohort does not support an unqualified claim that this module pair is always more strongly integrated among taxa. It is not equivalent to the earlier selected 42-taxon result.

## 4. Reproduction and evidence identities

Master metadata ZIP SHA-256:
`2f216652602a0b3c42b64780d9b31a3ecb2e7fc1300b28b303e8c369b31884ae`

All-350 authority artifact **12242753035** (run **39375870888**) SHA-256:
`9590aaab804bf1d77601b2352fc5ab8085b61b1226c5da10c52609fd89366370`

native_sources_350.json SHA-256:
`077709f7a94718ac4d19b5062a840c43cda46c9e893ebcbb60fc8b7aba016952`

Inherited observation-level measurement CSV SHA-256:
`b49a093c849928c456d6747ec51772b5a6b78debb2ae5eb5dbe2a91b08f66dd7`

The WCVP dataset version remains 2026-06-10T05:03:50.116+00:00 (GBIF dataset f382f0ce-323a-4091-bb9f-add557f3a9a2). The complete execution receipts include climate raster/cache hashes, scale/offsets, exclusion counts and analysis definitions.

Commands from repository root:

```bash
python -m analysis.v4.native_full --metadata-zip META.zip --native-sources native_sources_350.json --climate-root CLIMATE_ROOT --out-dir work/v4-native
python -m analysis.v4.native_exploration --native-dir work/v4-native --legacy exhaustive_continuous_observation_level.csv --out-dir work/v4-native
python -m pytest -q tests
```

Local full suite after implementation: **159 passed, 3 subtests passed**. Ten new tests cover native-first eligibility, conflicting metadata, unknown taxa, privacy, one-to-one linkage, raster scale/offset and multivariate decomposition invariants. Source-acquisition CI and mathematical tests are not independent biological validation; full master classification and the numbers here were executed locally.

## 5. Explicit remaining work

- New floral-region/peduncle/phyllary measurements at scale remain incomplete. Existing original-photo/annotation examples are not the full native cohort.
- Involucre biological validation and full-native coverage remain incomplete; no new native-wide involucre result is inferred from missing data.
- Fill the 106-observation climate cache gap and retain true nodata separately.
- Evaluate regularization, biological seasonal exposures and measurement uncertainty with held-out regions. Do not choose complexity merely to retain significance.
- Maintain native-only primary scope; assess coordinate uncertainty, native-region resolution, widespread taxa and phylogenetic placements as labelled sensitivities.
