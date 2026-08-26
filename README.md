# azami — global continuous thistle image phenomics

> **Submission status: HOLD (2026-08-26).** Seasonal, dominant-taxon,
> hemisphere-aware and native-range controls are complete. The native-only
> audit retained two of four non-circular primary rows. The former involucre
> headline failed its locked image-quality rule and was withdrawn. Image
> remeasurement and independent validation remain open.

This repository is organized around the Chapter 1 submission manuscript:

**Global citizen-science images reveal continuous within-taxon variation and environmental structure in thistle capitulum traits**

The canonical manuscript entry point is
[`manuscript/SUBMISSION_MANUSCRIPT.md`](manuscript/SUBMISSION_MANUSCRIPT.md).
Machine-readable claim status is frozen in
[`manuscript/final_claims.json`](manuscript/final_claims.json).

## Current paper in one sentence

Public biodiversity photographs provide repeated continuous head-level image
phenotypes, but after predeclared bias controls the present biological headline
is limited to two small observational climate rows—orientation–BIO1 and aspect
ratio–BIO4—both still awaiting developmental-stage control and independent
measurement validation.

## Scientific scope

YOLO11n is used only to localize visible capitula. Deterministic image-processing
functions then return continuous numerical measurements for each assessable
detected head.

Primary modules:

- **Orientation** — signed head-axis angle relative to image vertical;
- **Visible corolla colour** — CIELAB lightness, chroma and circular hue components;
- **Capitulum outline** — aspect ratio, circularity, solidity and width-profile CV.

A high-resolution layer also measures three exploratory two-dimensional
involucral contour proxies. Those proxies are retained as withdrawn sensitivity
provenance, not as current supported biological results.

The declared route is:

`public photographs -> detected capitula -> continuous head-level measurements -> repeated within-taxon distributions -> multiscale environmental analysis`

## Current results and reviewer controls

### 1. Visible variation below taxon means

Across the nine primary endpoints, **0.589–0.931** of visible image variance
occurs below source-assigned taxon means.

- among photographs within taxa: **0.440–0.691**;
- among heads within photographs: **0.143–0.379**;
- one-head-per-photo sensitivity: **0.582–0.899** below taxon means;
- equal-10-photo-per-taxon sensitivity: median **0.528–0.879**.

These are visible image-phenotype variance components. They do not separate
biological variation from imaging conditions or measurement error. An
outcome-blind preflight identified 20,073 observations with multiple public
photographs; trait remeasurement is still pending.

### 2. Primary within-taxon climate rows were narrowed

The frozen 46,276-observation analysis contained eight BH-supported rows among
36 endpoint-component × CHELSA tests. Four were non-circular primary rows.

All four non-circular rows retained their signs and BH support after both
raw-calendar and hemisphere-aware cyclic-season definitions. They also retained
their signs across the predeclared dominant-taxon omissions, although coefficient
magnitudes varied. *Cirsium vulgare* and *C. arvense* together account for 54.0%
of the strict cohort.

Under the separately locked WCVP/TDWG native-only rule, only two rows retained
the frozen sign and BH support:

- orientation–BIO1: native-only beta **0.02265**, BH q **0.0311**;
- aspect ratio–BIO4: native-only beta **0.02842**, BH q **1.25 × 10⁻17**.

The following rows were withdrawn from the current biological headline:

- chroma–BIO12: native-only BH q **0.0754**;
- aspect ratio–BIO12: native-only BH q **0.133**.

The two retained effects remain small, observational and developmentally
unlabelled. Native-only restriction does not establish adaptation, niche response
or introduction history.

### 3. The involucre headline was withdrawn

Resolution and sharpness were added as covariates and the BIO4 rows were checked
within fixed head-resolution strata. All three adjusted rows failed the locked
BH rule (**q = 0.0696–0.0730**), and two reversed sign in the 150–199 px stratum.
The former headline is therefore withdrawn. `involucre_spread_fraction` remained
positive in every resolution stratum and is preserved only as a lead for a
resolution-designed follow-up.

### 4. Spatial and historical layers are sensitivities

Residual Moran diagnostics, broad-region omission, grouped SPDE-INLA models,
taxonomic synonym collapse and trait-label permutation tests remain supporting
analyses. The phylogenetic layer is confined to Supporting Information because
only 54 of 216 atlas taxa are direct dated-backbone tips and retained Pagel
lambda estimates are zero. It is historical-placement sensitivity, not a
resolved species-tree correction.

### 5. The lability result remains withdrawn

The former raw absolute-slope lability relation (`rho = -0.333`) and median-split
quadrants were withdrawn after a precision-confounding audit. Precision-aware
reanalysis found no common cross-module coupling. The former result remains
statistical provenance/QA only.

## Claim boundary

Azami Chapter 1 does **not** claim demonstrated:

- phenotypic plasticity;
- local adaptation or selection;
- pollinator or antagonist causation;
- rain-protection or defensive adaptation;
- evolutionary rate differences or adaptive radiation;
- a resolved *Cirsium* species tree;
- direct botanical spine/phyllary measurements from image proxies.

Image vertical is a reproducible image reference rather than a direct
inclinometer measurement of gravity.

## Azami -> EAzami boundary

Azami is the global observational discovery layer. It ends after freezing:

`continuous visible capitulum variation -> within/among-taxon environmental structure`

Mechanism experiments, East Asian/Japanese evolutionary history, repeated-state
analyses and adaptive-radiation inference belong to
[`zuizui0223/EAzami`](https://github.com/zuizui0223/EAzami). The handoff is fixed
in
[`manuscript/EAzAMI_MECHANISM_REDUCTION_BOUNDARY_2026-08-20.md`](manuscript/EAzAMI_MECHANISM_REDUCTION_BOUNDARY_2026-08-20.md).

## Start here

1. [`manuscript/SUBMISSION_MANUSCRIPT.md`](manuscript/SUBMISSION_MANUSCRIPT.md) — manuscript order, current headline and blockers.
2. [`manuscript/00_title_abstract.md`](manuscript/00_title_abstract.md) through [`manuscript/07_data_code_availability.md`](manuscript/07_data_code_availability.md) — submission text and bibliography.
3. [`manuscript/COHORT_FLOW_AND_ANALYSIS_LEDGER.md`](manuscript/COHORT_FLOW_AND_ANALYSIS_LEDGER.md) — frozen cohorts and permitted analyses.
4. [`manuscript/final_claims.json`](manuscript/final_claims.json) — machine-readable numerical claim registry.
5. [`manuscript/BIAS_CONTROL_REANALYSIS_PROTOCOL.md`](manuscript/BIAS_CONTROL_REANALYSIS_PROTOCOL.md) — locked reviewer-control sequence and decision rules.
6. [`manuscript/FIGURE_TABLE_MAP.md`](manuscript/FIGURE_TABLE_MAP.md) — authoritative main/Supplement crosswalk.
7. [`analysis/ch1/pipeline.json`](analysis/ch1/pipeline.json) and [`analysis/ch1/run_submission.py`](analysis/ch1/run_submission.py) — active executable submission interface.
8. [`analysis_outputs/README.md`](analysis_outputs/README.md) — auditable result snapshot included with this review branch.
9. [`manuscript/EXTERNAL_COMPLETION_GATES.md`](manuscript/EXTERNAL_COMPLETION_GATES.md) — completed and open validation gates.

## Current main figures

1. image -> phenotype -> ecological scale;
2. geographic sampling and analytical domain;
3. nested visible variance;
4. taxon-level trait architecture;
5. environmental effect sizes and bias-control decisions.

The phylogenetic/historical-placement layer is Supporting Information only.

## Frozen implementation and provenance

Exact numbered files under `ch1_global/v2/` remain in place because executed
workflows and artifact provenance cite those paths. The active reviewer interface
is `analysis/ch1/`; manuscript-facing controls live under `manuscript/`.

This review branch includes a bounded audit snapshot under `analysis_outputs/`,
including locked model summaries, public-record audit matrices and PNG/PDF review
figures. GitHub Actions artifacts are not treated as a durable archive. The full
data/code release still requires an immutable DOI-backed repository.

## Current cohorts

- balanced image-comparison atlas: 3,725 observations/photos, 6,626 heads, 216 taxa;
- exhaustive detector-positive layer: 406,582 observations, 286 taxa;
- exhaustive spatially thinned primary: 46,276 observations, 259 taxa;
- grouped SPDE-INLA complete cases: 31,666–34,472 observations, 139–141 taxa per endpoint;
- high-resolution layer: 1,443 heads, 1,292 observations, 210 taxa;
- high-resolution <=10 km family: 904 observations, 165 taxa;
- independent detector audit: 1,000 source images, 323 taxa, 250 double-labelled images.

## Remaining submission gates

The package is not submission-ready. Open scientific gates are:

1. developmental-stage labels or an anthesis-restricted sensitivity;
2. paired flower/background visible-colour negative controls;
3. repeat-photo trait remeasurement and variance decomposition;
4. adjudicated human boxes for the independent detector audit;
5. independent reference measurements for orientation, colour and outline.

Human/administrative finalization also includes 24 nomenclatural notes,
authorship/CRediT metadata, funding and acknowledgements, selection of a
repository `LICENSE`, and a durable DOI-backed release.
