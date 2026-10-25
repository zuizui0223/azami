# Fixed spatial diagnostic

Defined before spatial outputs are inspected. The all-observation model remains
the reference. Native-range filtering does not remove spatial dependence.

Use a single global spherical equal-area cell grid: longitude intervals of
5 degrees from -180 and sin(latitude) intervals of 0.05 from -1. Cells have equal
spherical area except numerical boundary handling; their shapes differ with
latitude. This is a broad sampling diagnostic, not an estimated correlation
range, a climate-grid resolution, or proof of independent observations.

For every endpoint separately:

1. Reconstruct fixed-effect predictions on its exact main-model rows. Reconstruct
   the fitted species BLUP from the saved variance ratio and group residual sums.
   Assess conditional residuals, not raw traits.
2. Average residuals within occupied cells. Calculate Moran I on a symmetric
   eight-nearest-cell-centroid graph using spherical coordinates; 999 fixed-seed
   residual-label permutations. These are aggregate descriptive diagnostics,
   not a cell-level error model. A non-significant result is not proof of no
   spatial dependence. Report species/observations/cells and the residual scale.
3. Choose one eligible observation per species x cell, uniformly using seed
   20261020 after sorting by observation identity. Refit the same environment
   model; freeze full-population environment means/scales. No trait-dependent
   choice, search across cell sizes, or repeatedly thinning until significance.
4. Compare all-data and balanced-sample slopes and intervals. Sampling intervals
   remain species-cluster intervals; thinning does not independently establish
   spatially robust inference.

This comparison evaluates sensitivity to local sampling concentration. It does
not fully equalize species abundances across their geographic ranges or account
for photographer effects (observer identities are absent from the prepared
table). Do not use the diagnostic to declare that no spatial model is needed
solely because signs or significance agree. Remaining spatial dependence limits
the final claim and must be reported alongside the all-data estimates.
