# Native multivariate covariance decomposition

Retrospective extension requested after the environmental phylogenetic check.
Reuse the native common nine-coordinate cohort, at least five complete
observations per authority species, five fixed environmental predictors and
the original coordinate scaling (shared scale for the two hue coordinates).
Use exact, previously audited tree matching; do not graft additional taxa.

For raw and environmentally adjusted observations, estimate
`Cov(y_i,y_j) = A_species(i),species(j) P
               + 1[same species] U + 1[same observation] R`.
P is phylogenetic-placement covariance, U independent species covariance,
R pooled within-species covariance. Neither P nor U identifies genetic
covariance, natural selection, speciation or a molecular evolutionary history.

The environmental adjustment subtracts the previously estimated fixed
within/between effects, not species BLUPs. Estimates and simulations condition
on those coefficients. Raw and adjusted use exactly the same observations.

Estimation: multivariate quadratic method of moments with positive-semidefinite
constraints, NOT REML/MCMC. Species means and within-species crossproducts are
sufficient summaries here, not an omission of within-species variation.
For H = I - 11'/g, C1=HAH, C2=H, D=H diag(1/n_species) H:
S_l = M' C_l M - trace(C_l D) R; G_lk=trace(C_l C_k).
Fit P,U by minimizing the convex quadratic
0.5 sum_lk G_lk trace(V_l V_k) - sum_l trace(S_l V_l).
R is pooled within-species residual covariance (N-g degrees of freedom).
Keep unconstrained estimates and eigenvalues as diagnostics; never silently
clip an unconstrained matrix and call it an unbiased estimate. Constrained
estimates are biased at boundaries. Gram condition >1e8 stops identification.

All 52 historical placements are evaluated. For each placement and version,
200 parametric bootstrap replicates simulate species means using fitted
P,U,R and exact observation counts, plus Wishart within-species sums with
N-g degrees of freedom. This is conditional Gaussian-model uncertainty,
not robust species-cluster inference or a posterior. The 2.5/97.5 percentile
intervals are approximate, pointwise and may have boundary undercoverage.
Report across-tree ranges, not pooled confidence intervals.

Report 28 construct pairs using the established RV definition and signed
correlations for scalar pairs, plus variance fractions and full matrices.
Undefined correlation is retained when a component has negligible variance.
Within-species R pools observations and differs from the earlier equal-species
descriptive covariance. A changed estimate must not be attributed solely
to phylogeny. Shared pixel/outline construction remains a limitation.

Validation: exact expected-moment recovery with unequal observation counts,
PSD solver convergence, and nonidentifiability of proportional kernels.
All computations and tests run on Actions. No changes to V3, upstream
measurement, environmental selection or the prior REML estimates.
