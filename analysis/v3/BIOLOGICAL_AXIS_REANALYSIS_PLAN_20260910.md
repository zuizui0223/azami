# Chapter 1 v3 biological-axis reanalysis plan — 2026-09-10

## Purpose

The frozen v2 analysis remains the scientific baseline. This v3 extension does not replace its endpoint-level results or change the interpretation of the two leading v2 candidates. It asks whether the same 22 measured endpoints can be represented by a smaller set of biologically interpretable trait constructs before the environment association analysis is repeated.

The aim is to reduce measurement-level redundancy without forcing biologically distinct dimensions into a single global PCA.

## Starting universe

Use the 22 endpoints that were already measured in frozen v2. Do not introduce the five later-restored display / colour-composition endpoints into this analysis family.

The 22 endpoints collapse to 11 named constructs in `biological_trait_axis_contract_20260910.csv`. Ten constructs are testable ecological units; `surface_specularity` is retained descriptively because the current among-taxon min5 cohort has no taxon-median variation.

## Aggregation principles

1. Biological meaning precedes covariance. Endpoints are grouped only when they measure facets of the same biological construct.
2. A composite is used only where component direction can be stated biologically in advance.
3. Joint multivariate responses are used when related endpoints are not adequately one-dimensional.
4. CIELAB lightness, chroma and hue remain distinct because they are biologically different colour dimensions. Hue sine and cosine are one circular trait.
5. v2 endpoint-level candidate selection is frozen. Any new association created by the smaller v3 multiplicity family is exploratory rather than a retroactive v2 discovery.

## Trait-only structure audit motivating the aggregation

Using taxon medians with the frozen min5 rule and no environmental outcomes:

- capitulum circularity and solidity are strongly aligned (r about 0.73), while width-profile CV varies in the opposite direction; these three support a compactness / regularity composite. Aspect ratio remains a separate elongation dimension.
- projection roughness, p95, maximum and spread fraction form a coherent projection-prominence cluster (for example roughness-p95 r about 0.80 and p95-maximum r about 0.78).
- projection peak density is not part of that magnitude cluster and asymmetry contains a distinct spatial-pattern component, so these are retained jointly as `projection_pattern`.
- involucre length/width, apical taper and basal taper show weak mutual correlations; they are therefore retained as one three-dimensional biological construct rather than forced into a scalar score.
- surface edge density and high-frequency energy are related, but LBP entropy is largely distinct; surface texture is therefore a joint three-dimensional construct rather than a single score.
- the taxon-median specularity endpoint is constant under the min5 among-taxon rule and cannot support an among-taxon regression.

These decisions are based on endpoint semantics and trait-only covariance. Environmental coefficients are not used to choose construct membership or sign.

## Reanalysis

### Among taxa

For each scalar construct:

- require at least five original observations per taxon for every member endpoint;
- calculate endpoint taxon medians;
- standardize member endpoints across eligible taxa;
- calculate the frozen signed equal-weight composite where applicable;
- fit the same marginal standardized association against each of the nine frozen CHELSA predictors;
- use 9,999 taxon-label permutations for p values.

For circular or multivariate constructs, use the Euclidean norm of the standardized component coefficient vector as the joint test statistic and permute taxon labels 9,999 times.

Apply one Benjamini-Hochberg family across all successful biological-construct x predictor tests at the among-taxon min5 scale.

### Within taxa

Project the same construct definitions to observation-level endpoint values. Scalar components are standardized with fixed construct weights and then taxon-centred. Single scalar constructs use the same taxon-clustered OLS logic as v2. Joint constructs use within-taxon predictor permutations while preserving component vectors.

Apply one BH family across all successful biological-construct x predictor tests at the within-taxon scale.

### Cross-scale interpretation

Compare sign / joint coefficient geometry across the two scales. Do not call a mismatch opposing selection or plasticity. The allowed language is scale-dependent association geometry.

## Expected role in the manuscript

The v2 endpoint atlas remains the primary evidence map. The v3 biological-axis analysis is a higher-level synthesis designed to answer a different question: whether ecologically interpretable phenotype constructs recover the same leading environmental patterns while reducing repeated measurement facets.

A successful v3 result should strengthen, not redefine, the v2 conclusion. In particular:

- `floral_chroma` retains direct comparability with the v2 chroma-radiation result;
- `presentation_angle` retains direct comparability with the v2 orientation-precipitation result;
- new v3-only FDR hits are reported as exploratory consequences of the reduced construct family, not promoted above the frozen v2 candidates.
