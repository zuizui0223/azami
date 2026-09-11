# Chapter 1 v3 — scale-dependent biological-construct integration

## Why this analysis exists

The manuscript's earlier `Secondary whole-capitulum synthesis` asked whether the measured capitulum traits reduce to one syndrome. The complete-18 endpoint analysis found only partial alignment between within- and among-taxon association geometries (Spearman rho = 0.3663). The v3 analysis asks the same question at a more biologically interpretable level: does that partial alignment persist after the 18 original non-surface endpoints are reorganized into nine named biological constructs?

This analysis therefore develops the whole-capitulum perspective rather than treating within-taxon significance as a separate headline.

## Exact scope

The original complete-18 endpoint scope is retained and represented as nine constructs:

1. presentation angle;
2. floral lightness;
3. floral chroma;
4. floral hue;
5. head elongation;
6. head compactness;
7. involucre form;
8. projection prominence;
9. projection pattern.

Surface texture and surface specularity are excluded because they were not part of the original complete-18 whole-capitulum synthesis and remain validation-only image proxies.

For each of the 36 construct pairs, the exact same merged observations and taxa are used at both scales. Taxa require at least five paired observations, and every pair has at least 20 taxa. Within-taxon values are taxon-centred with equal total weight per taxon; among-taxon values are taxon medians calculated from that identical pairwise observation scope. Integration strength is the RV coefficient so multivariate constructs are not forced onto arbitrary scalar PCs.

## Main result: partial integration is reproduced after biological aggregation

All 36 construct pairs were evaluable.

The upper triangles of the within- and among-taxon construct-integration matrices were positively but only moderately aligned:

- Spearman rho = **0.3747748**;
- QAP-style one-sided P = **0.0131** from 9,999 construct-label permutations;
- QAP null mean = **0.001995**;
- QAP null 95% interval = **-0.3179 to +0.3295**.

This is strikingly close to the frozen complete-18 endpoint/unit result (rho = **0.3663**). Thus the earlier conclusion of partial, scale-dependent capitulum organization is not removed by replacing raw image endpoints with biologically named constructs.

The result should not be described as genetic, developmental or functional modularity. It is a scale-dependent organization of visible image-phenotype constructs.

## Integration strength differs strongly between scales

The median pairwise RV coefficient was:

- within taxa: **0.002352**;
- among taxa: **0.037031**.

Thirty of the 36 construct relations were stronger among taxa than within taxa. Because matrix edges share nodes, this count and the median contrast are descriptive rather than independent pairwise tests.

The largest among-taxon strengthening occurred in colour-related couplings:

| Construct relation | Within RV | Among RV | Among - within |
|---|---:|---:|---:|
| floral chroma — floral hue | 0.0626 | 0.3577 | +0.2951 |
| floral lightness — floral hue | 0.1089 | 0.3964 | +0.2875 |
| floral chroma — projection prominence | 0.00002 | 0.2491 | +0.2491 |
| floral lightness — floral chroma | 0.0981 | 0.2127 | +0.1146 |
| floral lightness — projection pattern | 0.00090 | 0.1081 | +0.1073 |

The clearest within-taxon strengthening was instead the armature/projection relation:

- projection prominence — projection pattern: within RV = **0.3038**, among RV = **0.1544**, delta = **-0.1494**.

A second within-strengthened relation was involucre form — projection prominence (within RV = 0.0396; among RV = 0.0176).

This suggests a useful biological interpretation for future testing: taxon-level differentiation integrates visible colour dimensions much more strongly, whereas repeated observations within taxa show particularly strong local co-variation among projection/armature dimensions. The present data do not identify the developmental or selective processes producing those differences.

## Relationship to the current manuscript synthesis

The pairwise result is retained as the first construct-level whole-capitulum check. The manuscript's stronger primary synthesis then forces all nine constructs and all 36 relations onto one exact complete-18 common cohort, adds taxon-bootstrap uncertainty and module-cohesion tests, and summarizes environmental organization with the six predeclared biological predictor groups. Those common-cohort results are reported in `CONSTRUCT_SCALE_UPGRADE_RESULT_20260910.md`.

Against that broader whole-capitulum reorganization, the two frozen-v2 among-taxon anchors — lower floral chroma under higher shortwave radiation and larger image-referenced presentation angle under higher annual precipitation — remain unusually robust, surviving construct aggregation and the full sampling -> spatial/residual -> 52-tree historical sensitivity chain.

## Reproducibility

- successful GitHub Actions run: `34421253852`;
- workflow job: `102696956553`;
- branch head: `858f7e681d901d1fd6099101611b63089934bde9`;
- output artifact: `10131007603`;
- artifact digest: `sha256:a483e116c46211023df95bef9388aa89deaef3441d85701fc93792867884514a`;
- QAP permutations: `9999`;
- source trait artifact: `9612943217`;
- source process-environment artifact: `9633419268`.

## Claim boundary

These results concern the organization of visible image-phenotype constructs. They do not establish plasticity, genetic or developmental modularity, causal environmental response, local adaptation, selection or convergence.
