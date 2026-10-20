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

- [ ] Verify pinned Zenodo tree archive and tip coverage.
- [ ] Check sufficient-statistic REML against dense observation covariance.
- [ ] Fit nine main coordinates across the frozen placement scenarios.
- [ ] Archive coverage, convergence, variance boundaries and coefficients.
- [ ] Interpret rainfall-orientation first; report other components without
      selecting a favourable tree or changing predictor sets.
