# Chapter 1 v3 RV estimator-validity sensitivity — 2026-09-15

## Why this check was added

The frozen common-cohort scale contrast uses the ordinary RV coefficient on two representations with different row counts: **1,734 within-taxon observations** versus **42 taxon medians**. Ordinary RV has a positive finite-sample baseline that can depend on row count and block dimensionality. During manuscript audit this raised a submission-critical question: does the headline statement that visible-phenotype integration is stronger overall among taxa survive when the scale contrast is not allowed to benefit from this raw-RV baseline difference?

This is explicitly a **post-hoc estimator-validity sensitivity**. The raw result was already known before this check was designed. It is not preregistration and it does not reopen ecological discovery.

The exact frozen trait/environment inputs and exact **1,734-observation / 42-taxon / nine-construct / 36-relation** common cohort were reused.

## Frozen raw result

The original common-cohort estimate remains unchanged:

- median RV within taxa = **0.0022379**;
- median RV among taxon medians = **0.0432123**;
- **33/36** raw construct relations are stronger among taxa.

These values remain descriptive outputs of the frozen estimator. The new analysis asks how much of their magnitude survives estimator-focused sensitivities.

## 1. Equal-n within resampling — primary estimator-validity gate

Within-taxon construct features were first centred using the full common cohort. One centred observation was then sampled per taxon, so every within-scale replicate had **42 rows**, matching the **42 rows** in the among-taxon median calculation. The among-taxon matrix was held fixed and this sampling was repeated **1,000 times**.

Results:

- matched within-taxon median RV: median **0.023897**, 95% interval **0.012861–0.041329**;
- among-minus-within median RV: median **+0.019316**, 95% interval **+0.001884 to +0.030352**;
- the median-RV difference was positive in **98.3%** of replicates;
- relations stronger among taxa: median **23/36**, 95% interval **18–27**;
- a majority (>18/36) of relations was stronger among taxa in **97.3%** of replicates.

The frozen decision rule therefore **passes**. The qualitative manuscript statement that integration is **stronger overall among taxa** is retained after matching the nominal row count used by the two RV calculations.

The raw **33/36** count is *not* promoted as a bias-robust headline. Under equal-n resampling the typical count is 23/36, not 33/36.

## 2. Permutation-null centring — diagnostic chance-baseline check

For each construct pair, chance RV was estimated separately at each scale with **499 permutations**. Among taxa, the right construct was permuted across taxon medians. Within taxa, the right construct was permuted within taxon after centring, preserving taxon sizes and equal-total-taxon weighting. The mean null RV was subtracted from observed RV.

Results:

- median null-centred RV within taxa = **−0.001068**;
- median null-centred RV among taxa = **+0.006950**;
- median among-minus-within contrast = **+0.008019**;
- **19/36** relations are stronger among taxa after null centring.

This diagnostic reaches the same overall direction but shows why the raw 33/36 relation count should not be treated as estimator-invariant evidence. Negative null-centred RV values simply indicate an observed RV below the estimated permutation baseline; they are not negative biological integration.

No new significance test is defined from this diagnostic.

## 3. Coordinate-standardization sensitivity

All construct coordinates were standardized by their across-taxon-median standard deviation before recomputing RV. This primarily changes the relative weighting of coordinates inside multivariate constructs.

Results:

- median RV within taxa = **0.002245**;
- median RV among taxa = **0.049800**;
- **34/36** relations stronger among taxa;
- cross-scale matrix alignment: **rho = 0.40849**, QAP **P = 0.0067**;
- module cohesion remains supported within taxa (**P = 0.0010**) and among taxa (**P = 0.0369**).

Thus the partial cross-scale geometry and module result are not artifacts of leaving the joint construct coordinates on their original numerical scales.

## Manuscript consequence

The estimator audit changes **emphasis**, not the central conclusion.

Supported wording:

> Visible-phenotype integration is stronger overall among taxon medians than within taxa on the common cohort. This direction persists when the within-scale RV calculation is restricted to the same nominal row count as the among-taxon calculation, although the raw 33/36 relation count and raw RV difference are partly sensitive to the finite-sample baseline of RV.

For the Main paper:

- retain the raw matrices and raw RV values as the frozen descriptive analysis;
- retain **stronger overall among taxa** as the conceptual result;
- do **not** present raw `33/36` or the raw taxon-bootstrap `+0.0659 / 100%` result as if they alone establish estimator-robust strength;
- show or report the equal-n sensitivity: median difference **+0.0193**, 95% **+0.0019 to +0.0304**, **98.3% positive**, and typical relation count **23/36**;
- retain the matrix-alignment and module-cohesion claims, which also survive the coordinate-standardization sensitivity.

## Reproducibility

- workflow: `Chapter 1 RV estimator validity`;
- run: **34933875866**;
- job: **104267550836**;
- artifact: **10382387052**;
- artifact digest: `sha256:345866a7e333f78677ad3797811e2cbe82d3e09ee061f7642df7f4c4d5ec008e`;
- frozen trait artifact: `9612943217`;
- frozen environment artifact: `9633419268`;
- equal-n replicates: `1000`;
- null permutations per pair per scale: `499`;
- coordinate-standardized QAP/module permutations: `9999`.

The compact frozen summary is `reproducibility/current_reference/estimator_validity/rv_estimator_validity_summary.json`.

## Claim boundary

This sensitivity supports only a descriptive scale contrast in **visible image phenotypes**. It does not establish that evolution increases integration, genetic or developmental modularity, plasticity, adaptation, or a causal mechanism.
