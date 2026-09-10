# Chapter 1 v3 construct-scale upgrade — 2026-09-10

## Purpose

This layer strengthens the reviewer-highlighted whole-capitulum synthesis without changing the frozen v2 endpoint atlas, multiplicity families, or the two manuscript-level ecological headline candidates.

It asks whether the partial within-vs-among organization remains when the v3 biological constructs are forced back onto the exact complete-18 repeated-observation cohort used by the frozen v2 whole-capitulum synthesis, and whether biologically defined modules remain internally cohesive at both scales.

## Exact complete-18 common cohort

The 18 non-surface endpoints yield exactly the frozen v2 repeated-observation cohort:

- 1,874 complete observations across 124 taxa before the minimum-five requirement;
- **1,734 observations across 42 taxa** after requiring at least five complete observations per taxon.

All nine biological constructs and all 36 construct pairs use this same 1,734-observation / 42-taxon cohort at both scales.

## Scale-dependent integration survives the common-cohort restriction

The common-cohort within- and among-taxon construct integration matrices are positively but incompletely aligned:

- upper-triangle Spearman rho = **0.439125**;
- QAP one-sided P = **0.0041** from 9,999 construct-label permutations.

For comparison:

- frozen v2 endpoint/unit complete-18 synthesis: rho = **0.3663**;
- pairwise-cohort v3 construct synthesis: rho = **0.374775**, QAP P = **0.0131**;
- common-cohort v3 construct synthesis: rho = **0.439125**, QAP P = **0.0041**.

Thus the conclusion of **partial, scale-dependent whole-capitulum organization is not created by pairwise cohort differences**. It persists, and is slightly stronger, when every construct pair is forced onto the exact same complete-18 cohort.

### Taxon-bootstrap uncertainty

A 1,000-replicate taxon bootstrap gives:

- median rho = **0.31210**;
- 95% bootstrap interval = **0.01708 to 0.53517**;
- fraction of bootstrap replicates with rho > 0 = **0.979**.

The alignment is therefore uncertain in magnitude but consistently positive under taxon resampling. This supports partial conservation rather than either complete invariance or complete reorganization.

## Integration is biologically modular at both scales

Constructs were grouped a priori into presentation, colour, head form, and involucre/armature modules. Module labels are biological definitions rather than environmental-result labels.

### Within taxa

- mean within-module RV = **0.09616**;
- mean between-module RV = **0.00387**;
- difference = **+0.09228**;
- label-permutation P = **0.0013**.

### Among taxa

- mean within-module RV = **0.16238**;
- mean between-module RV = **0.07297**;
- difference = **+0.08941**;
- label-permutation P = **0.0365**.

Therefore the capitulum is neither one undifferentiated syndrome nor a collection of unrelated image measurements. **Biologically related dimensions are more internally integrated than unrelated dimensions at both scales, while the detailed relation matrix is only partially conserved across scales.**

This construct-level result recovers and sharpens the frozen v2 whole-capitulum interpretation.

## The strength of integration increases mainly among taxa

Across the 36 construct relations:

- median within-taxon RV = **0.00224**;
- median among-taxon RV = **0.04321**;
- **33/36** relations are stronger among taxa than within taxa.

This is a descriptive scale contrast, not a claim that evolutionary processes necessarily increase integration. Photography, taxonomic sorting, history, ecological sorting and biological differentiation can all contribute to the among-taxon structure.

## Six predeclared environmental blocks

The nine v2 predictors were returned to the six biological blocks already used in the manuscript:

1. thermal: BIO1 + BIO4;
2. hydric: BIO12 + BIO15;
3. radiative/atmospheric: shortwave radiation + VPD;
4. mechanical: wind;
5. growing-season water: GSP;
6. resource/productivity: NPP.

Block scores use size-balanced RMS effect magnitude and are normalized within each construct. This is descriptive because correlations remain within and among blocks.

The block-level within-vs-among environmental signatures remain only weakly rank-aligned:

- flattened normalized signature Spearman rho = **0.10349**;
- cosine similarity = **0.85263**;
- only **2/9** constructs have the same strongest biological block at the two scales.

Thus coarsening the correlated predictors into predeclared biological blocks reduces the apparent disagreement compared with the raw nine-predictor comparison, but the dominant environmental organization is still usually scale-specific.

## Phenotypic integration is not simply environmental-profile similarity

An exploratory QAP compared each scale's construct-integration matrix with pairwise similarity of six-block environmental profiles.

- among taxa: rho = **0.04299**, QAP P = **0.8475**;
- within taxa: rho = **0.23449**, QAP P = **0.1652**.

There is therefore no evidence here that constructs are integrated merely because they share similar broad environmental association profiles. This negative result helps separate the two manuscript ideas: **phenotypic organization** and **environmental sorting/covariation** are related questions, not interchangeable statistics.

## Relationship to the frozen v2 positive results

None of these analyses changes the frozen endpoint-level positive results. In particular, the two frozen-v2 final among-taxon candidates remain:

- lower floral chroma with higher shortwave radiation;
- larger image-referenced presentation angle with higher annual precipitation.

The separate construct-level robustness chain already showed that both survive biological reaggregation, sampling composition, broad/residual spatial checks and all 52 historical placement trees.

The added result here is conceptual rather than a new headline association: **the same multidimensional capitulum has a partially conserved but scale-dependent internal organization, and the two frozen ecological anchors remain robust against that broader reorganization.**

## Reproducibility

- GitHub Actions run: `34432661967`
- job: `102731234004`
- artifact: `10135139012`
- artifact digest: `sha256:65438082976af8aa0135edc964a59f531e8763fde25f92739ae31c5056429f7e`
- bootstrap replicates: `1000`
- QAP/module permutations: `9999`
- frozen trait artifact: `9612943217`
- frozen environment artifact: `9633419268`

Claim boundary: this is a secondary biological-construct synthesis. It is not evidence of functional, developmental or genetic modularity, not proof of plasticity, and not a replacement for the frozen v2 endpoint-level inference.
