# Core-four environmental sufficiency hypothesis — Chapter 1

Status: frozen before process-extension outcomes, 2026-08-27

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
