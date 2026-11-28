# Part 2 standard R Gaussian hypervolume pilot — audited results (2026-11-28)

## Provenance and completion

- Full run [46106255971](https://github.com/zuizui0223/azami/actions/runs/46106255971), commit 65dd6e11ab87a16c4fe13f5afa32f65e835545d1, success; artifact 15250546191, actual ZIP SHA256 4f2d47748273aeee02d10cabdfc8ec75b1f934a50103dc4732caac9d7cac0fac, matches GitHub digest.
- Priority run [46109779446](https://github.com/zuizui0223/azami/actions/runs/46109779446), commit 3eb112f479258255564d63839b010a7bdd11b9dd, success; artifact 15248972769, ZIP SHA256 0966054331a64d90c85ba92ce7dc3fa6b226c7939e2a1638f5be9405beddaddd, matches GitHub digest.
- R 4.4.3, hypervolume 3.1.7. Gaussian bandwidth 0.5 native SD, 95% probability content, common scale within each species, frozen cell-median cohorts. All 108 planned comparisons are present. All four priority primary fits reproduce the corresponding full-run volumes, intersection, ratios, Jaccard and Sørensen exactly. Recalculated maximum absolute set partition errors < 8.1e-12, Jaccard error < 5.0e-15. This is an independent audit of aggregate outputs, not an independent fit of confidential source observations.

## Primary comparisons

| Species | Space and comparison | Reference n | Comparison n | Volume ratio | Jaccard | Intersection/comparison volume |
|---|---|---:|---:|---:|---:|---:|
| Cirsium arvense | Raw observed climate, native vs all introduced, 5D | 1366 | 860 | .911608 | .467730 | .668251 |
| Cirsium vulgare | Raw observed climate, native vs all introduced, 5D | 1113 | 1139 | 2.158066 | .279814 | .319947 |
| Cirsium arvense | Trait, native vs introduced novel climate, 3D | 1366 | 30 | .334951 | .327428 | .983080 |
| Cirsium vulgare | Trait, native vs introduced novel climate, 3D | 1113 | 190 | .866290 | .608861 | .815296 |
| Cirsium arvense | Trait, native vs all introduced, 3D | 1366 | 860 | .644067 | .610379 | .967521 |
| Cirsium vulgare | Trait, native vs all introduced, 3D | 1113 | 1139 | .934968 | .759719 | .893484 |

The climate package volumes are **background uncorrected**. Their projection is not the same as the prior PCA-env density cloud. The previous background-corrected grid climate volume ratio .614247 / 1.636041 remains a separate coarse-grid descriptive result.

## Numerical and sample-size sensitivity

Higher-point-density repeats of primary 3D native-vs-novel trait volume ratios .336022 / .866924 and shared fractions .984524 / .809955 closely agree. At equal cell counts (eight fixed native subsamples, 30 and 190 each, respectively), ratios .472213–.878074 / .956945–1.260951; shared fractions .723140–.882457 / .635042–.756945. Shared volume remains a majority, but the magnitude depends strongly on reference coverage. These ranges are not confidence intervals, and there is no geographic matching of native to introduced in these subsamples. Under the seven full-sample settings (including 2D sensitivity), shared fractions .959373–.986419 / .788177–.862468. The full-sample comparison does not establish robust volumetric contraction or equivalence.

Raw climate equal-n volume ratios .886737–1.107483 / 2.009244–2.198940: only C. vulgare consistently shows expansion under these settings. Geographically matched analogous-to-novel trait ratios .701232 / 1.250834 for the core95 definition, and 1.005659 / .897167 for the density-surplus definition. Do not infer different radiation modes or significant expansion/turnover. No new inferential tests were run; previous block-randomization tests of grid-based geometry remain non-significant after correction.

## Interpretation and figure adoption

The standard package reproduces the descriptive coexistence of broad novel-climate trait-space overlap and raw climate expansion in C. vulgare, with reduced shared fractions after equal-n subsampling. Inference remains limited by only two focal taxa, geographic and climate opportunity differences, image measurements, and native coverage. The previous independently calibrated near-neighbor support proportions 96.2% / 92.4% are fractions of cell-draws in novel climate, not hypervolume intersection fractions. The previous positive regional rank association for C. vulgare (rho=.149545, Holm4=.018) and nonlinear continuous chroma comparisons remain separate analyses.

Use native and introduced-novel modeled set-operation point clouds (not raw observations) to show true 3D Gaussian hypervolume geometry. The poster retains background-corrected 2D PCA-env climate contours for clarity and labels the coarse 5D corrected ratio as a separate approximation; the optional figure displays both package estimates, with raw 5D climate projected into two dimensions. The main claim is broad overlap of observed trait spaces, not proof of adaptive use of ancestral alleles. Preserve all negative results. No change to main branch or underlying raw data. 
