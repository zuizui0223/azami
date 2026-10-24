# V4 phylogenetic sensitivity implementation plan

**Goal:** Test the existing V4 environmental coefficients against inherited
historical-placement trees without replacing observations by species means.

**Architecture:** Reuse authenticated native inputs and fixed five-predictor
selection. Fit Gaussian observation-level REML with covariance
`sigma_e^2 I + Z (sigma_phy^2 A + sigma_species^2 I) Z'`.
Species sums are exact sufficient statistics, not a species-mean analysis.
Compare phylogenetic and independent-species fits on identical tree-covered
observations. Report missing taxa and all 52 placement scenarios.

**Constraints:** Actions-only computation; no new environmental selection;
no fabricated placements; no edits to frozen V3. These are placement
sensitivities, not molecular species-tree confirmation. Conditional
model-based intervals are not the earlier species-cluster intervals.
No joint colour or evolutionary covariance claim from univariate fits.

- [x] Verify pinned Zenodo tree archive and tip coverage.
- [x] Check sufficient-statistic REML against dense observation covariance.
- [x] Fit nine main coordinates across the frozen placement scenarios.
- [x] Archive coverage, convergence, variance boundaries and coefficients.
- [x] Interpret rainfall-orientation first; report other components without
      selecting a favourable tree or changing predictor sets.

## Execution ledger

- New analysis is retrospective sensitivity to the already reported V4
  coefficients, not an independently preregistered confirmation.
- The frozen tree resource contains grafted placements. The historical
  V3 report explicitly says most taxa are grafted; no claim of resolved
  Cirsium species history is permitted.
- Authority identifiers are not tree labels. Exact source names map each
  authority group to one tip; absent or ambiguous matches are excluded
  and listed, and many-authority-to-one-tip matches stop execution.
- First Actions attempt stopped on a duplicate Biopython version constraint.
  Reused the repository's pinned dependency without changing analysis.
- Second attempt passed likelihood tests but exposed pandas underscore-column
  attribute access; changed to explicit bracket access.
- Run 39435174409 executes the corrected model. Integrity run 39435183498
  passed 160 tests and 3 subtests.
- No species means replace individual responses. Grouped crossproducts give
  the exact Gaussian observation likelihood; a dense-covariance test covers
  positive and zero phylogenetic/independent variance components.
- Secondary coordinate intervals are pointwise, not whole-colour tests.
  A phylogenetic covariance analysis between traits is not implemented here.
- Run 39435174409 completed successfully: 468 phylogenetic fits plus nine
  matched-cohort independent-species fits, no recorded failures.
  Tree coverage is 194 taxa. Orientation retains 169912/170066 observations.
  The between-species rainfall coefficient is positive on all 52 placements
  (0.4722 to 0.5078 degrees/100 mm); every corresponding conditional interval
  excludes zero. Every within-species rainfall interval includes zero.
