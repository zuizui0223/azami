# GEB v2 multilevel capitulum-space final results — 2026-08-27

## Provenance

- PR: `zuizui0223/azami#72`
- successful workflow run: `33035785120`
- PR head used by the run: `227c0e7b8c338894806785b8545c7c77c8724de1`
- PR synthetic merge SHA recorded inside the artifact: `523c4ca0a1e6988a2b62f3d2f53df8b7442dff72`
- artifact ID: `9632715852`
- artifact digest: `sha256:51e7a26b5bd09e030b67b9342586699abaaf46e630f45b6bb4ee7bfc9152ced6`
- source GEB-v2 artifact: `9612943217`
- source GEB-v2 ZIP SHA-256: `101e996b638996a0c5ae79d358bf51293c3585f0e84c4a961b91dcbedf96211e`
- permutations per inferential test: 10,000
- taxon bootstrap replicates for multivariate geometry: 1,000

The source images were not remeasured. The analysis consumes the frozen final GEB-v2 continuous-trait universe and its strict spatial CHELSA table.

## 1. Symmetric within- and among-taxon environmental comparison

The 18 measured endpoints form 17 inferential units because circular hue is treated jointly.

Across the same four frozen CHELSA predictors, the main >=5 among-taxon scope contains 68 endpoint × predictor models:

- 4/68 among-taxon rows pass FDR;
- all four are in the established primary tier;
- no candidate endpoint passes among-taxon FDR;
- cross-scale classification: 3 `both_scales`, 8 `within_only`, 1 `among_only`, 56 `neither`.

Orientation is scale specific: its established within-taxon association is with BIO1, whereas its among-taxon association is with BIO12. All three final candidate atlas rows are `within_only`: involucre length/width–BIO15, apical taper–BIO12 and bract projection p95–BIO15.

## 2. Expanded among-taxon environmental sorting

The common-taxon >=5 scope retained 44 taxa complete for every measured endpoint. No linear endpoint passed either FDR sorting metric; joint circular hue was supported. The >=2 sensitivity retained 78 taxa; orientation was the only linear endpoint passing both centroid and overlap FDR, and joint hue was supported. No candidate linear endpoint passed the expanded sorting test at either threshold.

## 3. Multilevel organization of the whole capitulum phenotype

All 18 endpoints were complete for 1,874 observations from 124 taxa. The main >=5 scope retained 1,734 observations from 42 taxa; the >=2 sensitivity retained 1,825 observations from 75 taxa.

Main >=5 scope:

- within-taxon registered-module contrast: `0.164502`, bootstrap 95% CI `0.130693–0.179475`;
- among-taxon registered-module contrast: `0.088475`, CI `0.024942–0.126171`;
- within-vs-among 17-unit association-matrix Spearman: `0.366299`, CI `0.118457–0.399817`.

The >=2 sensitivity gives the same qualitative result:

- within contrast `0.157688`;
- among contrast `0.083662`;
- cross-scale matrix Spearman `0.377272`.

The registered modules are therefore more internally organized than expected under label permutation at both scales, while the full covariance geometry is only partly shared. This is evidence for **partial phenotypic / measurement-module organization**, not validated functional or genetic modularity.

## 4. Process-environment extraction

The complete-18 environmental cohort retained all 1,874 observations and 124 taxa. CHELSA coverage was 1.000 for:

- mean shortwave radiation;
- mean vapour-pressure deficit;
- mean near-surface wind;
- growing-season precipitation;
- potential NPP.

All variables and block definitions were frozen before process-extension outcomes were inspected.

## 5. Stand-alone whole-capitulum environment blocks

Six predeclared blocks were tested separately at within- and among-taxon scales. No block passed BH correction in either replication scope. The main within-taxon thermal block was nominally supported (`R2 = 0.008510`, `P = 0.01510`) but did not pass the six-block family (`Q = 0.09059`). These stand-alone tests are not promoted as new evidence-atlas discoveries.

## 6. Are BIO1/BIO4/BIO12/BIO15 sufficient?

The predeclared nested test used the frozen four predictors as the reduced model and added radiation, VPD, wind, growing-season precipitation and potential NPP.

### Within taxa

The omnibus extension was unsupported:

- >=5: `delta_R2 = 0.013350`, `partial_R2 = 0.013535`, `P = 0.24468`;
- >=2: `delta_R2 = 0.016314`, `partial_R2 = 0.016512`, `P = 0.13079`.

No within-taxon process block passed its BH family. The four-variable core is therefore adequate for the current within-taxon whole-capitulum estimand.

### Among taxa

The omnibus extension was supported at both thresholds:

- >=5: core `R2 = 0.079680`, full `R2 = 0.277528`, `delta_R2 = 0.197849`, `partial_R2 = 0.214978`, `P = 0.000700`;
- >=2: core `R2 = 0.069350`, full `R2 = 0.156691`, `delta_R2 = 0.087340`, `partial_R2 = 0.093849`, `P = 0.03020`.

Growing-season precipitation was the only stable block-specific increment:

- >=5: `partial_R2 = 0.078652`, `Q = 0.003200`;
- >=2: `partial_R2 = 0.031017`, `Q = 0.031997`.

Radiation + VPD (`partial_R2 = 0.083285`, `Q = 0.03253`) and potential NPP (`partial_R2 = 0.059760`, `Q = 0.014999`) were supported only in the main >=5 scope. Wind was unsupported in both scopes.

## Final environmental conclusion

> **The frozen four CHELSA variables are adequate for the current within-taxon 18D analysis, but they do not exhaust among-taxon environmental structure. Growing-season precipitation is the most replication-stable added environmental representation.**

This does not show that growing-season precipitation caused the phenotype. It may represent phenological water input, correlated climatic geography or historical turnover. The result is a scale-specific observational redundancy diagnosis.

## Chapter-1 synthesis

The final Chapter-1 pattern is not a single climatic capitulum syndrome. It is:

1. extensive visible variation below taxon means;
2. a high-dimensional taxon morphospace;
3. detectable internal organization of registered trait modules;
4. only partial correspondence between within- and among-taxon phenotype geometry;
5. endpoint-level climate associations that differ by trait and scale;
6. a richer process-environment representation required among taxa but not within taxa.

A submission-safe synthesis is:

> **The thistle capitulum is a multidimensional, partially organized phenotype whose environmental alignment changes across biological scales.**

## Azami → EAzami handoff

The artifact exports 62 provenance-gated observational targets:

- 6 structure targets;
- 24 environment-block R2 targets;
- 12 descriptive cross-scale coefficient-geometry targets;
- 20 incremental environment targets.

EAzami must retain all of them as `unscored_observational_target` until a future model emits statistically matched covariance and nested-environment estimands with a scoring rule frozen before model comparison.

## Claim boundary

These analyses do not establish functional or genetic modularity, plasticity, local adaptation, opposing selection, pollinator-mediated selection, defensive efficacy, direct rain/UV/thermal/aerodynamic protection, or a resolved causal environmental mechanism. Taxon-label permutation is not phylogenetic correction.
