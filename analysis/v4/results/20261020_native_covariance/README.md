# Joint native phenotype covariance decomposition

Numerical Actions run:
https://github.com/zuizui0223/azami/actions/runs/39437849165

All 104 decompositions completed (52 placements × raw/adjusted), with 200
conditional Gaussian parametric bootstrap replicates each. Three focused
recovery/identifiability tests passed. The cohort is 157561 observations from
170 tree-covered taxa, drawn from the existing common cohort of 157654/174.

The estimator is **PSD-constrained quadratic moments, not REML**. Full
definitions, assumptions and boundaries are in
`analysis/v4/COVARIANCE_CONTRACT.md`. Environmental coefficients remain fixed
at their existing estimates. This is a conditional sensitivity analysis,
not a new confirmatory test of evolutionary modules.

## Principal result

Phylogenetic correlation directions are unstable across placements and
bootstrap samples. Every decomposition has at least one phylogenetic
covariance eigenvalue on the PSD boundary. Extreme correlations near ±1
must not be described as strong evolutionary integration.

After environmental adjustment, the phylogenetic variance fraction is
0.01–1.79% for orientation and 0.01–6.93% for chroma across placements;
lightness ranges 0.01–11.77%. These are point-estimate ranges, not confidence
intervals or evidence that phylogeny has no effect.

The independent-species component retains several associations:

- Lightness–chroma r = −0.606 to −0.466; conditional intervals exclude zero
  on 51/52 placements.
- Orientation–chroma r = −0.518 to −0.203; intervals exclude zero on 51/52.
- Circularity–solidity r = 0.853 to 0.930; positive intervals on 52/52.

Within-species pooled residual correlations are approximately −0.294 for
lightness–chroma, −0.048 for orientation–chroma and +0.853 for
circularity–solidity. Shared measurement geometry contributes to the last
association. These are not the previous equal-species descriptive estimates.

Most variance is assigned to the within-species residual component under
this model. That component combines biological variation with photographic
and measurement variation; it is not a direct estimate of plasticity.

Thus the data support observed and model-conditional phenotype associations,
but do not establish a stable evolutionary covariance structure from these
grafted trees. Preserve uncertainty rather than presenting the bright
phylogenetic RV heatmap as biological evidence.
