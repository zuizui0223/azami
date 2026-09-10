# Chapter 1 v3 colour-space interpretation defense — 2026-09-10

## Purpose

This post-hoc defense asks what the frozen v2 negative floral-chroma × shortwave-radiation association means in CIELAB colour space. It does **not** alter the frozen v2 endpoint family, q-values, or headline association.

The key distinction is that lower C* means lower chroma, not necessarily darker colour. Darkening would require an accompanying decrease in L*.

## Frozen joint-colour cohort

The same minimum-five taxon rule was applied simultaneously to:

- `corolla_lab_lightness` (L*),
- `corolla_lab_chroma` (C*),
- `corolla_hue_sin`,
- `corolla_hue_cos`.

The resulting common cohort contains **143 taxa**.

## Radiation displacement in colour space

Against `chelsa_rsds_mean`:

- **L***: standardized beta = **+0.086009**, permutation P = **0.2973**;
- **C***: standardized beta = **−0.345372**, permutation P = **0.0001**;
- hue sin/cos joint magnitude = **0.500371**, permutation P = **0.0001**;
- joint L* + C* + hue displacement magnitude = **0.614045**, permutation P = **0.0001**.

The hue component vector is:

- sin beta = **+0.252489**;
- cos beta = **+0.431996**.

## Bootstrap direction check

Across **2,000 taxon bootstrap replicates**:

- L* beta median = **+0.08577**, 95% interval **−0.07894 to +0.25296**;
- C* beta median = **−0.34755**, 95% interval **−0.48543 to −0.19830**;
- P(C* beta < 0) = **1.000**;
- P(L* beta < 0) = **0.1545**;
- P(L* < 0 and C* < 0) = **0.1545**.

## Interpretation

The frozen v2 result is therefore best described as a **radiation-associated decrease in floral chroma embedded in a broader colour-space shift**. The present data do **not** support calling that shift darker: L* does not decrease detectably and its point estimate is slightly positive.

Accordingly:

- supported: **higher radiation → lower photographed floral chroma**;
- supported: radiation is also associated with **hue reorganization**;
- not supported: **higher radiation → darker flowers**;
- not demonstrated: **higher radiation → greater anthocyanin concentration**.

Anthocyanin remains a biologically plausible mechanistic hypothesis for pink–purple thistle coloration, but CIELAB displacement alone cannot distinguish pigment concentration, pigment composition, vacuolar pH/co-pigmentation, structural optical effects, or image-colour calibration effects. Biochemical or calibrated spectral validation is required before translating the chroma association into an anthocyanin-quantity claim.

## Claim boundary

This is an interpretation defense, not a new confirmatory endpoint family. The manuscript should retain the frozen wording `lower floral chroma under higher radiation` and avoid replacing it with `darker`, `more pigmented`, or `more anthocyanic` without external calibration.

Reproduction workflow: `Chapter 1 v3 colour-space defense`, run `34438065808`, artifact `10136908709`.
