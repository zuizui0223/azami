# Core-four environmental sufficiency hypothesis — Chapter 1

Status: design frozen before process-extension outcomes; final artifact-backed outcome appended, 2026-08-27

## Question

The established GEB-v2 environmental core uses four CHELSA predictors:

- BIO1: annual mean temperature;
- BIO4: temperature seasonality;
- BIO12: annual precipitation;
- BIO15: precipitation seasonality.

The process extension adds mean shortwave radiation, mean vapour-pressure deficit, mean near-surface wind, growing-season precipitation and potential NPP. The central question is not merely whether any of those added variables is associated with the 18-endpoint capitulum phenotype. It is:

> **Do the process-extension variables contain spatial information about the observed 18D capitulum phenotype that is not already represented by BIO1/BIO4/BIO12/BIO15?**

This is the analysis that determines whether the four-variable environmental core is adequate for the Chapter-1 spatial phenotype estimand.

## Frozen estimand

The response is the same complete-18 observation cohort used by the whole-capitulum environmental analysis. Environment extraction is restricted to that predeclared complete-case cohort only to reduce remote raster I/O; observations are selected by endpoint measurement completeness, not by trait magnitude or environmental outcome.

Two biological scales are evaluated separately:

1. within taxon: taxon-centred traits and environmental predictors, with inverse taxon sample-size weights so each retained taxon has equal total weight;
2. among taxon: taxon medians.

The main replication scope requires at least five complete-18 observations per taxon. A >=2-observation sensitivity is retained exactly as in the other capitulum-space analyses.

## Frozen nested tests

The reduced model always contains the four core predictors. Full models add predeclared extension predictors.

### Omnibus test

`all_process_extension_beyond_core4`

Adds simultaneously:

- mean shortwave radiation;
- mean VPD;
- mean near-surface wind;
- growing-season precipitation;
- potential NPP.

This is one predeclared omnibus test per replication scope and biological scale. It is not included in the block-specific BH family.

### Block-specific tests

Four nested tests are evaluated separately:

- radiative + atmospheric drying beyond core4;
- mechanical exposure beyond core4;
- growing-season water input beyond core4;
- climatic productivity beyond core4.

BH correction is applied across these four block-specific tests separately within each replication scope and biological scale.

## Statistics

For every nested comparison report:

- `R2_core4`;
- `R2_core4_plus_extension`;
- `delta_R2 = R2_full - R2_core`;
- `partial_R2 = delta_R2 / (1 - R2_core)`.

Inference uses a Freedman–Lane reduced-model residual permutation:

- within taxon: the complete 18-response residual vector for an observation may move only among observations of the same taxon;
- among taxon: residual vectors are permuted among taxa.

The whole response vector is permuted as a row, preserving covariance among the 18 measured endpoints under the null.

## Frozen interpretation rule

The environmental core will be described as **adequate for the current 18D Chapter-1 estimand** only if the process-extension omnibus test is unsupported in the main >=5 scope and no stable block-specific incremental signal survives the predeclared multiplicity/sensitivity interpretation.

If the omnibus test is supported, or a process block shows a stable supported incremental contribution, the conclusion will instead be that the four core variables are a useful low-dimensional baseline but are **not sufficient to exhaust environmental structure in the capitulum phenotype space**.

A supported incremental test does not establish that the named environmental process caused the phenotype. It only shows that the added environmental representation contains spatial information not captured by the four core predictors in this observational dataset.

## Relation to the existing evidence atlas

This nested analysis does not reopen or alter the established endpoint-level A/B/C evidence grading. The original four-variable primary family remains frozen. The incremental analysis is a higher-level multivariate diagnostic answering whether a broader process-based environmental representation improves the whole-capitulum spatial description.

## EAzami handoff

EAzami receives the incremental `partial_R2`, `delta_R2`, permutation P, block-specific Q and support state as provenance-gated **unscored observational targets**. They cannot enter a mechanism-family score until a future generator emits a statistically matched nested core-vs-extension estimand and its scoring rule is frozen before model comparison.

## Final artifact-backed outcome

The final successful execution was GitHub Actions run `33035785120` on PR-head `227c0e7b8c338894806785b8545c7c77c8724de1`. Artifact `9632715852` has digest `sha256:51e7a26b5bd09e030b67b9342586699abaaf46e630f45b6bb4ee7bfc9152ced6`. It consumes the completed GEB-v2 source artifact `9612943217` and retains the same 18 measured endpoints.

The process-environment cohort contained 1,874 complete-18 observations from 124 taxa. Coverage was 1.000 for mean shortwave radiation, mean VPD, mean near-surface wind, growing-season precipitation and potential NPP. The main >=5 scope retained 1,734 observations from 42 taxa; the >=2 sensitivity retained 1,825 observations from 75 taxa.

### Stand-alone environmental blocks

When each of the six predeclared blocks was tested separately against the whole 18D response, no block passed BH correction at either biological scale or replication threshold. In the main scope, the within-taxon thermal block was nominally supported (`R2 = 0.008510`, permutation `P = 0.01510`) but failed the six-block BH family (`Q = 0.09059`). Stand-alone block tests therefore do not provide a promoted whole-capitulum result.

### Increment beyond BIO1/BIO4/BIO12/BIO15

The nested comparison produced a clear scale contrast.

At the within-taxon scale, adding all five process variables to the frozen four-variable core was unsupported in both scopes:

- >=5: `delta_R2 = 0.013350`, `partial_R2 = 0.013535`, `P = 0.24468`;
- >=2: `delta_R2 = 0.016314`, `partial_R2 = 0.016512`, `P = 0.13079`.

No within-taxon block-specific extension passed the predeclared BH rule. Thus the four-variable core is an adequate low-dimensional representation for the present **within-taxon whole-capitulum** estimand, even though individual endpoints retain trait-specific associations.

At the among-taxon scale, the five process variables added substantial information beyond the core in both scopes:

- >=5: `R2_core4 = 0.079680`, `R2_full = 0.277528`, `delta_R2 = 0.197849`, `partial_R2 = 0.214978`, `P = 0.000700`;
- >=2: `R2_core4 = 0.069350`, `R2_full = 0.156691`, `delta_R2 = 0.087340`, `partial_R2 = 0.093849`, `P = 0.03020`.

Growing-season water input was the only block-specific increment supported at both thresholds:

- >=5: `delta_R2 = 0.072385`, `partial_R2 = 0.078652`, `Q = 0.003200`;
- >=2: `delta_R2 = 0.028866`, `partial_R2 = 0.031017`, `Q = 0.031997`.

Radiative/atmospheric drying (`partial_R2 = 0.083285`, `Q = 0.03253`) and climatic productivity (`partial_R2 = 0.059760`, `Q = 0.014999`) were supported in the main >=5 scope but not after the >=2 sensitivity (`Q = 0.12012` and `0.07279`, respectively). Mechanical exposure was unsupported in both scopes.

## Final decision on the four-variable question

The answer is scale dependent:

> **BIO1/BIO4/BIO12/BIO15 are adequate for the current within-taxon 18D capitulum analysis, but they are not sufficient to exhaust among-taxon environmental structure.**

The most replication-stable missing information is growing-season precipitation, not another undifferentiated expansion of the BIOCLIM screen. Radiation/VPD and potential productivity remain plausible but threshold-sensitive among-taxon extensions. This result strengthens the Chapter-1 argument that environmental organization rotates across biological scales.

The result remains observational. Growing-season precipitation may summarize phenological water input, correlated climatic geography or historical turnover; it is not demonstrated rain protection, water-balance selection, plasticity or adaptation. The process-extension result also does not alter the frozen endpoint-level evidence atlas.
