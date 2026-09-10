# Chapter 1 v3 L*a*b* colour interpretation defense — 2026-09-10

## Purpose

This post-hoc defense does **not** replace the frozen-v2 `floral chroma × shortwave radiation` result or its multiplicity family. It asks what that colour shift means in a more intuitive CIELAB decomposition.

The exact 143-taxon colour cohort used by the earlier colour-space defense was retained. `a*` and `b*` were reconstructed from frozen taxon-level `C*` and hue direction:

- `a* = C* cos(h)` — red/magenta direction when positive;
- `b* = C* sin(h)` — yellow direction when positive, blue direction when negative.

These reconstructed axes are interpretive derivatives, not new measured endpoints and not biochemical anthocyanin concentrations.

## Radiation slopes

| axis | standardized beta | permutation P | interpretation |
|---|---:|---:|---|
| L* | +0.0860 | 0.2973 | no supported darkening; if anything the point estimate is slightly lighter |
| C* | **-0.3454** | **0.0001** | strong desaturation, reproducing the frozen-v2 anchor |
| a* | -0.0980 | 0.2430 | weak/uncertain shift away from red-magenta |
| b* | **+0.2481** | **0.0026** | supported shift toward the yellow/warmer side and away from the blue/purple side |

The joint `L*–a*–b*` displacement has magnitude **0.2803** with permutation **P = 0.0301**.

## Taxon bootstrap (2,000 replicates)

- L*: median +0.0858, 95% interval **-0.0789 to +0.2530**; negative in 15.45% of replicates.
- C*: median -0.3475, 95% interval **-0.4854 to -0.1983**; negative in 100% of replicates.
- a*: median -0.1012, 95% interval **-0.2692 to +0.0671**; negative in 88.5% but interval crosses zero.
- b*: median +0.2477, 95% interval **+0.0938 to +0.3927**; positive in 99.85% of replicates.

## Interpretation

The frozen radiation association is therefore **not a simple darkening axis**. Higher-radiation taxa are consistently lower in chroma and are shifted in hue mainly through the `b*` direction, with no supported decrease in L* and no robust `a*` increase.

A defensible manuscript interpretation is:

> Higher shortwave radiation is associated with a less chromatic and warmer-shifted floral colour state, rather than demonstrated darkening.

The result does **not** show that anthocyanin concentration increases with radiation. Anthocyanin amount, anthocyanin composition, vacuolar pH, co-pigmentation and other optical effects remain candidate mechanisms requiring biochemical calibration.

## Reproducibility

- workflow run: `34444840554`
- job: `102767294595`
- artifact: `10139260925`
- artifact digest: `sha256:cc08d7e9e88d56c32741b75880d32c43dffb975a3d75e52d323f224a6979f7f4`
- permutations: `9999`
- taxon bootstrap replicates: `2000`
