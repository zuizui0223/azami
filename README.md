# azami — submission repository for global continuous thistle image phenomics

This branch is organized around the Chapter 1 submission manuscript:

**Global citizen-science images reveal continuous within-taxon variation and environmental structure in thistle capitulum traits**

## Scientific scope
Public biodiversity photographs are converted into repeated continuous head-level phenotype measurements after YOLO localization. The primary endpoints are orientation, visible corolla colour and two-dimensional capitulum outline; a high-resolution layer adds three involucral contour proxies.

The current headline is deliberately scale-explicit:

> public photographs → detected capitula → continuous head-level measurements → within-taxon trait distributions → within- and among-taxon environmental structure

Across nine primary endpoints, 0.589–0.931 of visible image variance occurs below source-assigned taxon means. Orientation and visible colour show the clearest environmental structure after spatial modelling; most gross-outline effects are weaker or model-dependent. These are observational image-derived associations, not genetic variance, plasticity, local adaptation, pollinator causation or evolutionary rate.

## Azami → EAzami boundary

Azami is the **global observational discovery layer**. It ends after freezing:

`continuous visible capitulum variation -> within/among-taxon environmental structure`

Azami does **not** own inference about:

- pollinator- or antagonist-mediated selection;
- defensive function of spine/phyllary or sticky traits;
- East Asian/Japanese phylogenetic or biogeographic history;
- repeated loss/regain of floral colour;
- repeated/parallel erect↔nodding transitions;
- trait-specific evolutionary rates or modular evolvability;
- rapid/adaptive radiation.

Those questions belong to `zuizui0223/EAzami`, the downstream **mechanism + evolutionary-history zoom**. EAzami combines the frozen Azami pattern bundle with quantitative ecological literature, East Asian evolutionary history and ancestry-resolved focal experiments. Its current doctoral architecture is:

`Azami global phenotype landscape`

`-> mechanism meta-analysis`

`-> East Asian/Japanese rapid-radiation history`

`-> repeated-state / parallel-evolution tests`

`-> focal trait -> interaction/protection -> reproductive fitness`

Do not move EAzami mechanism-reduction, East Asian repeated-state analysis or adaptive-radiation inference back into Chapter 1. The original observation→mechanism handoff is documented in:

`manuscript/EAzAMI_MECHANISM_REDUCTION_BOUNDARY_2026-08-20.md`.

## Start here
1. `manuscript/SUBMISSION_MANUSCRIPT.md` — manuscript-facing story and section order.
2. `manuscript/COHORT_FLOW_AND_ANALYSIS_LEDGER.md` — frozen cohort names, counts and permitted analyses.
3. `manuscript/final_claims.json` — machine-readable numerical claim registry.
4. `analysis/ch1/pipeline.json` — active submission analysis map.
5. `analysis/ch1/run_submission.py` — structural/claim checks and status summary.
6. `manuscript/supplement/` — Supporting Information source and tables.
7. `manuscript/EXTERNAL_COMPLETION_GATES.md` — remaining external validation.
8. `manuscript/EAzAMI_MECHANISM_REDUCTION_BOUNDARY_2026-08-20.md` — explicit observation → mechanism handoff.

## Frozen implementation
Exact numbered files under `ch1_global/v2/` remain in place because executed workflows and artifact provenance cite those paths. Do not rename or rewrite them merely to simplify the tree. The active author/reviewer interface is `analysis/ch1/`.

## Current cohorts
- balanced image-comparison atlas: 3,725 observations / 6,626 heads / 216 source-assigned taxa;
- exhaustive detector-positive layer: 406,582 observations / 286 taxa;
- exhaustive spatially thinned primary: 46,276 observations / 259 taxa;
- grouped SPDE-INLA complete cases: 31,666–34,472 observations and 139–141 taxa per endpoint;
- high-resolution involucre layer: 1,443 usable heads / 1,292 observations / 210 taxa;
- independent detector audit: 1,000 source images / 323 taxa / 250 double-labelled images.

## Remaining scientific submission gates
Only two scientific blockers remain external to repository computation:

1. adjudicated human boxes for the independent detector audit;
2. genuinely independent reference measurements for orientation, colour and outline.

Taxonomic synonym-collapse sensitivity, environmental-niche permutation tests, residual Moran diagnostics and broad-region omission are complete for the frozen operational-unit analysis.

## Withdrawn result
The former raw absolute-slope lability relation (`rho = -0.333`) and median-split quadrants were withdrawn after precision-confounding audit. They remain provenance/QA only and must not re-enter the biological headline.
