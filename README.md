# azami — global continuous thistle image phenomics

> **Submission status: HOLD (2026-08-26).** The global pattern analyses and major computational sensitivities are complete for the frozen primary trait set. An expanded continuous-trait run is in progress for explicitly defined involucral architecture, armature and surface-image metrics. Independent measurement validation remains open.

This repository is organized around the Chapter 1 submission manuscript:

**Global citizen-science images reveal continuous within-taxon variation and environmental structure in thistle capitulum traits**

The canonical manuscript entry point is [`manuscript/SUBMISSION_MANUSCRIPT.md`](manuscript/SUBMISSION_MANUSCRIPT.md).

## Current paper in one sentence

Public biodiversity photographs are used to move broad-scale thistle trait ecology beyond species means and coarse morphological categories by converting visible capitula into repeated continuous phenotype measurements, mapping their global spatial diversity and testing how different trait components align with environmental gradients.

## Scientific position

The Chapter 1 target is **pattern first, mechanism later**.

The paper asks:

1. what the global visible capitulum morphospace looks like when traits are quantified continuously rather than assigned to coarse states;
2. how much phenotype varies within source-assigned taxa as well as among taxa;
3. which environmental gradients align with orientation, colour, outline and fine capitulum architecture;
4. how those associations differ among within-taxon, spatial and among-taxon analytical scales.

The paper does **not** attempt to prove the causal mechanism behind each spatial association. Plasticity, adaptation, pollinator selection, antagonist defence and evolutionary-rate hypotheses are downstream tests motivated by the observed patterns.

## Analysis route

`public photographs -> detected capitula -> explicitly defined continuous traits -> repeated within-taxon distributions -> global phenotypic geography -> trait- and scale-specific environmental structure`

YOLO11n is used only to localize visible capitula. Deterministic image-processing functions then return numerical measurements for each assessable head. Unassessable traits remain missing rather than being converted into categorical absence.

### Frozen primary modules

- **Orientation** — signed head-axis angle relative to image vertical;
- **Visible corolla colour** — CIELAB lightness, chroma and circular hue;
- **Capitulum outline** — aspect ratio, circularity, solidity and width-profile CV.

### Expanded continuous modules under current reanalysis

The category-free continuous-trait contract additionally defines continuous measurements for:

- visible floral display;
- involucre length/width ratio;
- apical and basal involucre taper;
- radial bract projection roughness, upper-tail projection and maximum projection;
- outward-projection fraction;
- projection-peak density;
- bilateral projection asymmetry;
- high-resolution surface edge/texture/specular signals.

Surface signals remain validation-only and cannot currently be called hair density, gland density or mucilage.

## Current empirical pattern

### Continuous diversity below taxon means

Across the nine frozen primary endpoints, **0.589–0.931** of visible image variance occurs below source-assigned taxon means.

- among photographs within taxa: **0.440–0.691**;
- among heads within photographs: **0.143–0.379**;
- one-head-per-photo sensitivity: **0.582–0.899** below taxon means;
- equal-10-photo-per-taxon sensitivity: median **0.528–0.879**.

These are visible image-phenotype components, not genetic variance components. The result nevertheless shows how much phenotypic information is discarded by representing each taxon with one mean.

### Trait-specific environmental structure

The frozen 46,276-observation global within-taxon analysis identifies environmental associations involving:

- **orientation**;
- **visible colour** including chroma and circular hue;
- **outline aspect ratio**.

Grouped spatial models concentrate the strongest support in **orientation and visible colour**. Among-taxon permutation analyses likewise show the clearest environmental sorting for orientation and colour. Selected outline associations are strongest in the within-taxon or historical-sensitivity layers rather than uniformly across models.

### Current strongest and sensitive associations

Orientation–BIO1 and aspect ratio–BIO4 have the strongest robustness across seasonal, dominant-taxon and native-range sensitivities.

Chroma–BIO12 and aspect ratio–BIO12 remain same-signed global associations but lose multiplicity support under the restricted native-only analysis. They are therefore treated as **supported but range-sensitive**, not erased from the global pattern.

The original high-resolution involucral screen detected coherent seasonality-associated contour signals. Resolution/sharpness adjustment weakened them and two metrics reversed direction in the lowest-resolution stratum, so they are treated as **exploratory continuous signals** rather than confirmatory botanical-function results.

## Evidence tiers

Sensitivity analyses grade evidence strength rather than acting as one universal veto. The policy is fixed in [`analysis/ch1/evidence_tiering_policy.json`](analysis/ch1/evidence_tiering_policy.json).

- **A — robust:** global support plus consistency across relevant sensitivities;
- **B — supported but sensitive:** global support, but weakened in a restricted domain while retaining direction;
- **C — exploratory continuous signal:** coherent image-derived association requiring stronger image-quality or botanical validation;
- **D — unresolved/artifact-prone:** important sign reversal or measurement-integrity problem;
- **validation-only:** numerical image signals not yet linked securely to a botanical quantity.

## Why *Cirsium*

A single thistle capitulum combines orientation, colour, gross outline and involucral architecture on one reproductive structure. This makes it possible to compare several phenotype modules across the same global observational network rather than treating each as a separate categorical character.

The genus also contains substantial ecological and evolutionary diversity, but Chapter 1 uses that history as context rather than as evidence that the measured patterns are adaptive.

## Claim boundary

Azami Chapter 1 does **not** claim demonstrated:

- phenotypic plasticity;
- local adaptation or selection;
- pollinator or antagonist causation;
- rain-protection or defensive adaptation;
- evolutionary rate differences or adaptive radiation;
- a resolved *Cirsium* species tree;
- direct botanical spine/phyllary measurements from unvalidated image proxies;
- hair, gland or mucilage identity from generic surface-image metrics.

Image vertical is a reproducible image reference rather than a direct inclinometer measurement of gravity.

## Azami -> EAzami boundary

Azami is the global observational discovery layer. It ends after freezing:

`continuous visible capitulum diversity -> spatial pattern -> within/among-taxon environmental structure`

Mechanism experiments, East Asian/Japanese evolutionary history, repeated-state analyses and adaptive-radiation inference belong to [`zuizui0223/EAzami`](https://github.com/zuizui0223/EAzami).

## Start here

1. [`manuscript/SUBMISSION_MANUSCRIPT.md`](manuscript/SUBMISSION_MANUSCRIPT.md) — canonical paper story and submission boundary.
2. [`manuscript/00_title_abstract.md`](manuscript/00_title_abstract.md) through [`manuscript/07_data_code_availability.md`](manuscript/07_data_code_availability.md) — manuscript text and bibliography.
3. [`analysis/ch1/evidence_tiering_policy.json`](analysis/ch1/evidence_tiering_policy.json) — interpretation of global results versus sensitivity analyses.
4. [`ch1_global/v2/ontology/ch1_continuous_trait_contract.csv`](ch1_global/v2/ontology/ch1_continuous_trait_contract.csv) — category-free endpoint definitions.
5. [`manuscript/COHORT_FLOW_AND_ANALYSIS_LEDGER.md`](manuscript/COHORT_FLOW_AND_ANALYSIS_LEDGER.md) — frozen cohorts and permitted analyses.
6. [`manuscript/final_claims.json`](manuscript/final_claims.json) — numerical claim registry; narrative interpretation is governed by the evidence-tier policy.
7. [`analysis/ch1/pipeline.json`](analysis/ch1/pipeline.json) — active executable submission map.
8. [`analysis_outputs/README.md`](analysis_outputs/README.md) — auditable analysis snapshot.

## Current cohorts

- balanced image-comparison atlas: 3,725 observations/photos, 6,626 heads, 216 taxa;
- exhaustive detector-positive layer: 406,582 observations, 286 taxa;
- exhaustive spatially thinned primary: 46,276 observations, 259 taxa;
- grouped SPDE-INLA complete cases: 31,666–34,472 observations, 139–141 taxa per endpoint;
- original high-resolution layer: 1,443 heads, 1,292 observations, 210 taxa;
- original high-resolution <=10 km family: 904 observations, 165 taxa;
- independent detector audit: 1,000 source images, 323 taxa, 250 double-labelled images.

## Remaining scientific gates

The package is not yet submission-ready. Important open checks are:

1. developmental-stage labels or an anthesis-restricted sensitivity;
2. paired flower/background visible-colour negative controls;
3. repeat-photo trait remeasurement and variance decomposition;
4. adjudicated human boxes for the independent detector audit;
5. independent reference measurements for orientation, colour and outline;
6. botanical reference validation before expanded architecture/armature/surface image proxies are renamed as direct botanical traits.

These gates constrain interpretation and evidence strength. They do not change the paper's primary purpose: **quantifying global thistle phenotype diversity and its spatial environmental structure from public photographs.**
