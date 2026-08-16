# azami — global continuous capitulum phenomics in *Cirsium*

This repository contains a multi-chapter project on the ecology and evolution of *Cirsium* floral architecture. Chapter 1 is the submission-focused component: a global analysis that converts public biodiversity photographs into repeated **continuous head-level trait measurements** rather than one categorical state or one mean value per taxon.

## Chapter 1 in one sentence

Global public photographs recover continuous within-taxon distributions of capitulum orientation, visible colour and outline that are largely hidden by taxon means, and part of that variation is structured along environmental gradients in a trait- and scale-specific manner.

The analysis is observational and non-causal. Spatial climatic association is not proof of temporal response, local adaptation, phenotypic plasticity, pollinator selection, evolutionary rate or adaptive radiation.

## Methodological contribution

The core contribution is not simply large image volume. The production workflow is:

```text
public biodiversity photograph
        ↓
YOLO localization of visible capitula
        ↓
deterministic continuous measurement per assessable head
        ↓
orientation angle + colour coordinates + outline geometry
        ↓
repeated within-taxon trait distributions
        ↓
nested variance + within-taxon environment association + among-taxon sorting
```

The inferential traits are **not categorical classifier outputs**. Each assessable head contributes numerical values such as an orientation angle, CIELAB lightness/chroma and circular hue components, or outline statistics. Failed or unassessable measurements remain missing rather than being coded as biological absence.

Across the nine primary endpoints, 0.589–0.931 of visible image variance occurs below source-assigned taxon means. One-head-per-photo and equal-replication sensitivities retain the same qualitative conclusion. These are visible image-variance components, not genetic variance estimates.

## Why *Cirsium*

*Cirsium* is useful because the same reproductive structure expresses conspicuous diversity in orientation, visible corolla colour, gross outline and involucral architecture. The lineage also includes recent regional radiations and has a complex history involving hybridization, incomplete lineage sorting, polyploidy and phenotypic convergence. That combination makes it a strong system for asking whether visible phenotype is environmentally structured at the level of repeated observations, among taxa, or both, without assuming that one species mean or one resolved bifurcating tree is sufficient.

Chapter 1 does **not** itself demonstrate adaptive radiation. Rapid regional diversification is evolutionary context and motivation; explicit transition histories and resolved nuclear phylogeny are downstream EAzami tasks.

## Role in the wider research program

Chapter 1 is the **global macro-scale hypothesis-generation layer**. The next evolutionary-resolution stage is maintained in `zuizui0223/EAzami`:

```text
azami / Chapter 1
Global public-image continuous phenomics
        ↓ trait distributions and macro-scale hypotheses
EAzami
East Asian nuclear phylogeny + explicit transition histories
        ↓ replicated focal transitions
population / mechanism studies
Ancestry + expression + pigment + biological interaction + fitness
```

Therefore Chapter 1 does **not** need to absorb definitive ancestral-state reconstruction, transition counts, pollinator causation or molecular regain mechanisms before submission. Those are downstream tests. The explicit boundary and handoff are recorded in `manuscript/EAzami_HANDOFF_AND_CHAPTER_BOUNDARY_2026-08-16.md`.

## Canonical executed cohorts

| Cohort | Observations | Taxa | Purpose |
|---|---:|---:|---|
| Balanced image-comparison atlas | 3,725 | 216 | 6,626-head nested variance, PCA and among-taxon summaries |
| Exhaustive detector-positive layer | 406,582 | 286 | Post-detection source layer |
| Exhaustive spatially thinned primary | 46,276 | 259 | Primary within-taxon climate coefficients |
| Grouped SPDE-INLA complete cases | 31,666–34,472 per endpoint | 139–141 | Spatially explicit within-taxon models |
| High-resolution involucre subset | 1,292 observations / 1,443 usable heads | 210 | Auxiliary involucre/spine-like contour proxies |
| Species-level historical sensitivity | 213–214 primary taxa | 216 atlas taxa; 54 direct tips | PGLS across alternative grafting trees |
| Independent detector audit | 1,000 source images | 323 | Human validation packet; annotation pending |

The full flow from 777,766 photographs to derived tables is recorded in `manuscript/COHORT_FLOW_AND_ANALYSIS_LEDGER.md`.

## Primary and auxiliary traits

### Primary continuous endpoints

- orientation angle relative to EXIF-oriented image vertical;
- corolla Lab lightness and chroma;
- circular hue sine/cosine;
- capitulum aspect ratio, circularity, solidity and width-profile variation.

### Auxiliary high-resolution involucre layer

The final manuscript retains three inferential auxiliary contour proxies:

- involucre projection roughness;
- involucre spread fraction;
- maximum relative spine-like projection length.

`spine_peak_count_proxy` also exists in the integrated derived table and is preserved for downstream provenance, but it is **not one of the three final manuscript auxiliary inferential endpoints**.

All involucre/spine variables are image-geometry proxies, not direct botanical measurements of phyllary angle, spine length, orientation, stiffness or defensive function.

## Main retained results

- Below-taxon visible variance: 0.589–0.931 across nine primary endpoints.
- Among-photograph within-taxon component: 0.440–0.691; among-head within-photo component: 0.143–0.379.
- Primary within-taxon climate family: 8/36 BH-supported component rows, including orientation × BIO1, chroma × BIO12, aspect ratio × BIO4/BIO12 and four circular-hue component rows.
- Grouped SPDE-INLA most consistently retains temperature associations with orientation and visible colour, precipitation with chroma and soil-pH associations with visible colour; gross-outline effects do not survive the global SPDE screen.
- High-resolution auxiliary layer: roughness, outward spread and maximum spine-like projection all increase with temperature seasonality in the ≤10 km cohort.
- Permutation-supported among-taxon environmental sorting is concentrated mainly in orientation and visible colour.
- Residual Moran's I is near zero for all nine primary endpoints; leave-one-broad-region-out taxon rankings remain strong.
- Simultaneous collapse of all eight WCVP synonym conflicts changes neither primary coefficient signs nor FDR decisions.
- Historical/PGLS results remain sensitivity analyses because only 54/216 atlas taxa are direct dated-backbone tips.

The withdrawn raw absolute-slope lability relationship and median-split quadrants are retained only as statistical provenance/QA and are **not** part of the current headline story.

## Canonical manuscript path

Start here:

1. `manuscript/SUBMISSION_MANUSCRIPT.md` — current submission-facing story and section order;
2. `manuscript/COHORT_FLOW_AND_ANALYSIS_LEDGER.md` — immutable cohort names, counts and analysis permissions;
3. `manuscript/final_claims.json` — machine-readable frozen claims;
4. `manuscript/results/nested_visible_variance_summary.csv` — image-hierarchy decomposition;
5. `manuscript/results/NICHE_PERMUTATION_RESULTS_2026-08-07.md` — niche-null test;
6. `manuscript/results/SPATIAL_ROBUSTNESS_RESULTS_2026-08-07.md` — residual and region-omission diagnostics;
7. `manuscript/results/WCVP_TAXONOMIC_SENSITIVITY_2026-08-07.md` — synonym-collapse sensitivity;
8. `analysis/ch1/pipeline.json` — canonical executable stages;
9. `manuscript/EXTERNAL_COMPLETION_GATES.md` — genuinely external remaining validation;
10. `manuscript/EAzami_HANDOFF_AND_CHAPTER_BOUNDARY_2026-08-16.md` — macro → phylogeny → mechanism handoff.

## Current submission gates

The computational spatial, niche-null and taxonomic-robustness gates are complete for the frozen operational-unit analysis. The two remaining **scientific** submission blockers are genuinely external:

1. adjudicated human boxes for the independent detector audit;
2. genuinely independent reference measurements for orientation, colour and outline.

Nomenclatural notes for the 24 high-priority WCVP rows and authorship/administrative metadata also require human confirmation. Durable archive/DOI release is the final step after the measurement gates close.

## Reproducibility and interpretation limits

- Every FDR count must name its cohort, endpoint family and number of tests.
- Multiple heads from one photograph remain a nested level or are reduced by one-head-per-photo sensitivity.
- Circular hue components are interpreted jointly and never as separate named flower colours.
- Public photographs are not colour calibrated; image vertical is a reproducible image reference, not gravity; outline measures remain viewpoint dependent.
- Production overlays and mirror repeatability are technical checks, not independent biological accuracy validation.
- The grafted mega-tree is a historical sensitivity analysis, not a resolved *Cirsium* species tree.
