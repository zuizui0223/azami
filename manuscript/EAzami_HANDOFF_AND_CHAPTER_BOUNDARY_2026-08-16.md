# Chapter 1 boundary and EAzami handoff — 2026-08-16

## Decision

Chapter 1 remains the **global macro-scale continuous image-phenomics layer**. Its methodological contribution is to convert public photographs into repeated numerical capitulum-trait observations and retain within-taxon distributions that would be lost by coarse categories or one taxon mean. It already includes primary orientation/colour/outline traits and a high-resolution auxiliary involucre/spine-like layer. It does not need a resolved *Cirsium* nuclear history or causal adaptation tests before submission.

```text
Azami / Chapter 1
Global public-image continuous phenomics
        ↓ trait distributions + macro-scale hypotheses
EAzami
East Asian nuclear tree + explicit trait-transition history
        ↓ replicated focal transitions
Focal systems
Population genomics + mechanism + ecological fitness tests
```

## What Chapter 1 already contains

### Primary continuous trait modules

1. orientation angle relative to EXIF-oriented image vertical;
2. visible corolla colour: Lab lightness/chroma and circular hue;
3. capitulum outline: aspect ratio, circularity, solidity and width-profile variation.

These are numerical image-derived endpoints, not categorical classifier outputs. The Chapter 1 inference retains repeated head- and observation-level values rather than treating each taxon as one fixed state.

### High-resolution auxiliary involucre/spine-like layer

Chapter 1 has already measured high-resolution image-geometry variables for:

- `involucre_projection_roughness`;
- `involucre_spread_fraction`;
- `spine_peak_count_proxy`;
- `spine_relative_length_max_proxy`.

The **three final manuscript inferential auxiliary endpoints** are:

1. involucre projection roughness;
2. involucre spread fraction;
3. maximum relative spine-like projection.

`spine_peak_count_proxy` remains in the integrated derived tables and is preserved for provenance/downstream exploration, but it is not one of the three final auxiliary inferential endpoints reported in the submission-facing main results.

The three final auxiliary endpoints have already entered the high-resolution within-taxon environmental analysis. Their main manuscript role is explicitly auxiliary: they appear in Methods, Results and Discussion because they extend the continuous-phenomics framework to finer capitulum architecture, but they are not promoted to the primary nine-endpoint confirmatory family. Auxiliary historical/PGLS mapping remains sensitivity/downstream material.

They must therefore **not** be described as wholly new or previously unanalysed traits.

## What remains unresolved about the involucre/spine traits

The current variables are image-geometry proxies. They are not yet equivalent to botanically validated measurements such as:

- actual phyllary spreading/recurvature angle;
- actual spine length;
- spine orientation or stiffness;
- homologous phyllary rank/position;
- direct florivore or pollinator accessibility.

The downstream question is therefore not “add spines for the first time,” but:

> Do the existing global involucre/spine proxy patterns correspond to reproducible botanical states, and how did those states evolve on the resolved nuclear tree?

Visible involucre stickiness/glandularity remains a genuinely separate trait requiring its own assessability and provenance work.

## Chapter 1 scientific boundary

### Keep in Chapter 1

- global public-image phenomics;
- continuous head-level orientation/colour/outline measurements;
- nested head → photograph → taxon variance structure;
- within-taxon climatic associations;
- grouped SPDE-INLA spatial models;
- among-taxon environmental niche permutation tests;
- existing high-resolution involucre/spine-like proxy analysis;
- spatial, broad-region and taxonomic robustness diagnostics;
- historical-tree analysis as sensitivity only;
- explicit measurement and sampling limitations.

### Hand off to EAzami

- accepted East Asian/Japanese nuclear species-tree ensemble;
- explicit transition-history and ancestral-state reconstruction;
- repeated white-loss / candidate-regain inference;
- Japan-38 and other focal lineage crosswalks;
- mapping of existing continuous orientation/colour/outline distributions onto the resolved nuclear history;
- phylogenetic mapping of the existing involucre/spine-like proxies;
- tests of whether image proxies correspond to explicit botanical phyllary/spine states;
- correlated or decoupled evolution among colour, orientation, outline and defensive architecture;
- reticulation, ploidy and cytonuclear uncertainty.

### Hand off to focal empirical studies

- pollinator preference and effective pollination;
- rain/UV/thermal protection experiments for orientation;
- florivore/seed-predator tests for phyllary/spine architecture;
- sticky-involucre defence tests;
- population ancestry and gene flow;
- floral RNA, pigment chemistry and candidate-locus tests;
- genotype → expression → pigment → phenotype → interaction → fitness chains.

## Cross-repository trait bridge

The Azami→EAzami bridge must preserve more than one taxon mean **when the downstream question concerns within-taxon structure**. Two bridge layers are therefore distinguished.

### A. Species/taxon summary bridge

Useful for broad topology mapping and coverage audits. For each endpoint retain:

- source and accepted taxon names;
- trait module and endpoint;
- continuous/circular/proxy state type;
- median estimate plus MAD or other uncertainty/range summary;
- number of observations and usable images;
- assessability/QC status;
- direct nuclear-tip match status;
- evidence provenance and claim boundary.

The current `azami_ch1_eazami_trait_handoff_v1` exporter is this **summary bridge**. It must not be described as preserving the full within-taxon distribution.

### B. Observation/distribution bridge

Required whenever EAzami tests population structure, multimodality, environmental gradients or within-taxon transition heterogeneity. It should preserve observation-level or distribution-resolving values from the exhaustive Chapter 1 source with:

- taxon and observation identifiers;
- continuous endpoint values;
- coordinate/QC provenance;
- assessability status;
- spatial block;
- source artifact checksum and frozen cohort identity.

PR57 already adds a Japan-38 exhaustive trait export from the frozen observation-level Chapter 1 layer. Broad image records are reduced only to species binomials for coverage and are **not** silently assigned to varieties or subspecies. This is the correct pattern for downstream focal exports.

For the auxiliary involucre/spine layer, every bridge must explicitly label values as **image-geometry proxies** so later botanical validation does not silently reinterpret them as direct phyllary measurements.

## Downstream analysis sequence

Once the accepted EAzami nuclear topology ensemble exists:

1. flower-colour transition-history analysis;
2. orientation ancestral/transition analysis;
3. continuous outline history;
4. map existing involucre spread/roughness/spine-like proxies onto the same tree;
5. quantify whether proxy patterns are phylogenetically structured and whether they covary with other capitulum modules;
6. validate or replace proxy endpoints with explicit botanical phyllary/spine measurements in targeted material;
7. add genuinely new traits such as stickiness only after assessability/provenance are established;
8. promote high-information transitions to population and ecological experiments.

## Stop rules

- No definitive ancestral-state reconstruction from the current grafted Chapter 1 tree.
- No adaptation or adaptive-radiation claim from macro correlation alone.
- No pollinator/herbivore causation without interaction and fitness data.
- No molecular regain claim from Chapter 1 images.
- Do not call the already-analysed involucre/spine-like proxy layer “new” or “unanalysed.”
- Do not equate image proxies with botanical phyllary/spine characters without validation.
- Do not replace continuous within-taxon distributions with a single species mean when the downstream question explicitly concerns population or within-taxon structure.

## Current Chapter 1 submission status

The repository-computable taxonomic, residual-spatial, broad-region and environmental-niche-null gates are complete. The two remaining scientific submission blockers are external:

1. adjudicated human boxes for the independent detector audit;
2. independent reference measurements for orientation, colour and outline.

Chapter 1 should remain focused on those credibility gates while EAzami proceeds in parallel with the nuclear tree and provenance-preserving trait bridge.
