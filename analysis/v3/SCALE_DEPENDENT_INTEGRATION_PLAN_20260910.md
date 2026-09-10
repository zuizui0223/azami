# Chapter 1 v3 scale-dependent construct integration

## Question

The reviewer-highlighted whole-capitulum perspective is not simply whether any one within-taxon coefficient is significant. The biological question is whether the capitulum has the same internal trait organization within taxa as it has among taxa.

The frozen v2 complete-18 synthesis already suggested only partial alignment between within- and among-taxon association geometries. v3 sharpens that test by replacing the 18 raw measurement endpoints with nine biological constructs while preserving the same endpoint scope.

## Fixed construct scope

The analysis uses the 18 non-surface endpoints from the original complete-18 synthesis, represented as nine constructs:

1. presentation angle;
2. floral lightness;
3. floral chroma;
4. floral hue;
5. head elongation;
6. head compactness;
7. involucre form;
8. projection prominence;
9. projection pattern.

Surface texture and surface specularity are excluded from this integration test because they were not part of the original complete-18 whole-capitulum synthesis and remain validation-only image proxies.

## Pairwise scale comparison

For every pair of constructs, the exact same merged observation scope is used for both scales. A taxon must contribute at least five paired observations, and at least 20 taxa must be available.

- **Within-taxon:** each construct component is taxon-centred and observations are weighted so every retained taxon contributes equal total weight.
- **Among-taxon:** taxon medians are calculated from the exact same paired observation scope.
- **Integration strength:** RV coefficient between the two construct blocks. This accommodates scalar and multivariate constructs without forcing hue or other multidimensional constructs onto arbitrary PCs.

The resulting nine-by-nine within and among integration matrices are compared by Spearman correlation of their upper triangles. A 9,999-permutation QAP-style construct-label test asks whether the named relational organization is more similar across scales than expected under relabelling.

## Environmental signature comparison

The already-computed construct-by-environment atlases are also compared descriptively. For each construct, its nine effect magnitudes are normalized within scale and the relative environmental signature is compared between within and among scales. Because the environmental predictors are correlated, no predictor-label permutation is used as confirmatory inference.

## Interpretation boundary

A difference between within- and among-taxon integration is not direct evidence for plasticity, selection, genetic modularity or developmental modularity. It shows that the organization of visible image phenotypes is scale dependent. The robustness chain for chroma-radiation and presentation-angle-precipitation remains a separate test of specific environmental patterns.
