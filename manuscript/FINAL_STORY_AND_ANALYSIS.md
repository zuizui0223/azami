# Chapter 1 final story and analysis boundary

## Frozen result

The strict within-species analysis contains 46,276 spatially thinned observations
from 259 taxa. The final species-level lability analysis requires at least 10
observations, six eligible traits and three climate predictors per trait; 102 taxa
pass this completeness rule.

The central result is that two quantities often called “lability” are not the same:

- **species within-variation**: robust, trait-standardized dispersion of the
  image-derived phenotype within a species;
- **species environmental responsiveness**: the root-mean-square magnitude of
  species-specific standardized trait–climate slopes.

Across the 102 complete taxa, these axes are moderately negatively associated
(Spearman rho = -0.333). High visible variation therefore does not imply strong
tracking of the four measured CHELSA gradients. The four-quadrant plot contains
22 high-variation/high-response, 29 low-variation/high-response, 29
high-variation/low-response and 22 low-variation/low-response taxa.

## Module contrast

Median species-level indices differ among modules:

| Module | Median within-variation | Median environmental responsiveness | Species |
|---|---:|---:|---:|
| Colour | 0.553 | 0.175 | 104 |
| Orientation | 0.630 | 0.129 | 101 |
| Shape | 0.665 | 0.153 | 103 |

Colour has the strongest median climate tracking but the lowest median total
variation. Shape has the largest visible variation without the largest climate
responsiveness. This contrast is more informative than ranking traits on a single
“labile–conserved” axis.

## Global within-species coefficients

Four linear endpoint–predictor combinations pass Benjamini–Hochberg FDR at 0.05,
but all standardized effects are small:

- orientation angle × BIO01: beta = 0.0171, q = 0.0363;
- corolla chroma × BIO12: beta = 0.0393, q = 0.0363;
- shape aspect ratio × BIO04: beta = 0.0137, q = 0.0207;
- shape aspect ratio × BIO12: beta = -0.0108, q = 0.0363.

The result must therefore be described as weak but detectable pooled climatic
structure, not as strong within-species adaptation. Circular hue is retained as a
joint sine/cosine vector; its magnitude is descriptive unless a joint uncertainty
test is explicitly added.

## Sensitivity

Changing the minimum observations per species preserves rankings:

- n >= 5: variation rho = 0.939; responsiveness rho = 0.999 relative to primary;
- n >= 20: variation rho = 0.914; responsiveness rho = 0.987 relative to primary.

The two-axis conclusion is therefore not driven by the chosen minimum-n threshold.
Mesh, correlation-cutoff and prior sensitivity belong to the spatial/Bayesian
between-species analysis and must not be presented as sensitivity of these
frequentist species-specific slopes.

## Novelty

The strongest novelty is the combination of:

1. continuous capitulum colour, orientation and outline measured from public
   photographs across a large taxonomic and geographic sample;
2. explicit analysis of measurement assessability rather than silently discarding
   failed images;
3. separation of among-species climatic structure from within-species association;
4. decomposition of species lability into total variation and climate tracking;
5. module-level comparison showing that colour, orientation and shape occupy
   different positions on the variation–responsiveness plane;
6. circular treatment of hue and cautious historical sensitivity under uncertain
   tree placement.

## Recommended headline

> **Visible phenotypic variation and climate tracking are distinct dimensions of
> floral-trait lability across thistles.**

A strong result-led title is:

> **Public photographs reveal decoupled variation and climate tracking in thistle
> capitulum traits**

A safer macroecological title is:

> **Scale-dependent climatic structure of image-derived capitulum traits in
> Cirsium**

## Interpretation boundary

The indices describe image-derived variation and observational climate association.
They do not demonstrate plasticity, local adaptation, pollinator selection,
evolutionary rate, rain protection or a causal environmental response. Public-image
colour is not spectrally calibrated, image vertical is a proxy for gravity, outline
traits are viewpoint dependent and species-specific slopes can contain geographic
or unmeasured environmental structure.

## Frozen outputs

- `species_lability_axes.csv`: final species-level axes and quadrants;
- `species_trait_within_variation.csv`: trait-level robust dispersion;
- `species_trait_environmental_responsiveness.csv`: trait-level slope magnitude;
- `species_module_lability.csv` and `module_lability_summary.csv`;
- `table_s1_all_global_coefficients.csv`;
- `table_s2_species_specific_coefficients.csv`;
- `sensitivity_minimum_sample_size.csv`;
- Figure S1 effect-size heatmap;
- species-level four-quadrant figure;
- `analysis_metadata.json`.
