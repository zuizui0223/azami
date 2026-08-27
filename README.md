# azami — the present phenotypic field of thistle capitula

> **Submission status: HOLD FOR INDEPENDENT VALIDATION (2026-08-27).** The accepted computational Chapter 1 pattern analysis is complete, including the expanded continuous-trait universe, matched within/among comparison, whole-capitulum organization and process-environment sufficiency tests. The remaining scientific blockers are independent detector/trait validation, developmental-stage and repeat-photo controls, paired colour controls and botanical calibration.

This repository is organized around the Chapter 1 manuscript:

**Global citizen-science images reveal multidimensional and scale-dependent environmental organization of thistle capitula**

The canonical submission entry is [`manuscript/SUBMISSION_MANUSCRIPT.md`](manuscript/SUBMISSION_MANUSCRIPT.md). The explicit hypothesis-recovery and completeness audit is [`manuscript/AZAMI_CH1_HYPOTHESIS_RECOVERY_AND_COMPLETENESS_2026-08-27.md`](manuscript/AZAMI_CH1_HYPOTHESIS_RECOVERY_AND_COMPLETENESS_2026-08-27.md).

## Chapter 1 in one sentence

> **The thistle capitulum is a multidimensional, partially organized phenotype whose environmental alignment changes across biological scales.**

## Scientific position

Azami reconstructs the **present phenotypic field** of thistle capitula from public photographs. It asks:

1. what phenotype diversity is lost when taxa are represented by coarse states or one mean;
2. how a single reproductive structure occupies a continuous multidimensional trait space;
3. how that space is organized within and among taxa;
4. which environmental representations align with the phenotype at each biological scale.

The chapter is pattern first. It does not infer plasticity, adaptation, pollinator or antagonist causation, defensive efficacy or a unique evolutionary history from cross-sectional image associations.

## Analysis route

`public photographs → detected capitula → explicitly defined continuous traits → repeated phenotype distributions → below-taxon diversity → taxon morphospace → within/among phenotype geometry → scale-specific environmental organization`

YOLO11n localizes visible capitula. Deterministic image functions then return numerical measurements for each assessable head. Unassessable traits remain missing rather than being converted to categorical absence.

## Frozen continuous-trait universe

The category-free contract defines 27 endpoints across orientation, visible colour, floral display, outline, involucre architecture, armature and validation-only surface/exudate signals.

Final execution:

- 27 registered endpoints;
- 22 measured endpoints;
- 19 inferential endpoints, of which 18 were measured;
- 9 established primary + 9 expanded candidate endpoints;
- `visible_floret_fraction` is the only unexecuted inferential endpoint;
- 46,276 spatially thinned observations from 259 source-assigned taxa;
- 374,255 analysis-eligible observation × endpoint measurements.

Surface metrics remain validation-only and cannot be called hair density, gland density, mucilage or secretion without targeted botanical validation.

## Current empirical pattern

### Taxon means conceal most visible variation

Across the nine established primary endpoints, **0.589–0.931** of visible image variance occurs below source-assigned taxon means.

- among photographs within taxa: **0.440–0.691**;
- among heads within photographs: **0.143–0.379**;
- one-head-per-photo sensitivity: **0.582–0.899**;
- equal-10-photo-per-taxon sensitivity: median **0.528–0.879**.

These are visible image-phenotype components, not genetic variance components.

### The capitulum is multidimensional rather than one syndrome

The expanded taxon morphospace contains 127 taxa complete for all 18 measured inferential endpoints.

- PC1: **18.49%**;
- PC2: **12.01%**;
- PC3: **11.78%**;
- PC1–PC3 cumulative: **42.28%**.

Finer measurement therefore reveals several partly independent phenotype dimensions rather than one hidden master axis.

### The whole phenotype is partially organized

Among 1,734 complete-18 observations from 42 taxa in the main replication scope:

- within-taxon registered-module contrast: **0.164502**;
- among-taxon registered-module contrast: **0.088475**;
- within/among 17-unit matrix Spearman: **0.366299**.

Registered modules are more internally organized than expected under label permutation, but their covariance geometry is only partly shared between scales. This is partial phenotypic / measurement-module organization, not validated functional or genetic modularity.

### Environmental structure is trait specific and scale specific

Across the same 17 inferential units and four frozen CHELSA predictors, the 68 matched rows classify as:

- 3 supported at both scales;
- 8 within taxa only;
- 1 among taxa only;
- 56 at neither scale.

Orientation aligns with annual mean temperature within taxa but annual precipitation among taxa. All three final expanded candidate rows are within-taxon only.

The strongest endpoint-level evidence remains:

- **A-grade:** orientation–BIO1 and outline aspect ratio–BIO4;
- **B-grade/range-sensitive:** chroma–BIO12 and aspect ratio–BIO12;
- **C-grade quality-robust:** involucre length/width–BIO15 and apical taper–BIO12;
- **C-grade image-sensitive:** `bract_projection_p95`–BIO15.

### Four climate variables are sufficient only at the within-taxon scale

The frozen environmental core is BIO1, BIO4, BIO12 and BIO15. Adding shortwave radiation, VPD, wind, growing-season precipitation and potential NPP gives:

- within taxa: unsupported incremental information (`partial R² = 0.013535`, `P = 0.24468` at the main threshold);
- among taxa: substantial incremental information (`partial R² = 0.214978`, `P = 0.000700`);
- growing-season precipitation is the only block-specific increment supported at both replication thresholds.

The core four are therefore an adequate low-dimensional representation for the current within-taxon 18D estimand, but they do not exhaust among-taxon environmental structure. This is an observational redundancy result, not proof of direct precipitation causation.

## Hypothesis recovery

The central hypotheses are resolved as follows:

- repeated continuous image phenotypes: computationally established, independently unvalidated;
- large below-taxon visible variation: supported for the primary endpoints;
- one universal capitulum syndrome: rejected;
- partial measurement-module organization: supported within and among taxa;
- identical within/among geometry: rejected; correspondence is partial;
- trait- and scale-specific environmental organization: supported;
- stable cross-scale coefficient-vector rotation: inconclusive;
- negative lability–responsiveness relationship: withdrawn;
- direct functional mechanisms: not tested in Azami.

Full definitions, numbers and claim boundaries are in the [hypothesis-recovery audit](manuscript/AZAMI_CH1_HYPOTHESIS_RECOVERY_AND_COMPLETENESS_2026-08-27.md).

## Evidence tiers

Sensitivity analyses grade evidence strength rather than acting as one universal veto. The policy is fixed in [`analysis/ch1/evidence_tiering_policy.json`](analysis/ch1/evidence_tiering_policy.json).

- **A — robust:** global support plus consistency across relevant sensitivities;
- **B — supported but sensitive:** global support weakened in a restricted domain while retaining direction, or calibration-limited;
- **C — exploratory continuous signal:** coherent candidate image phenotype requiring botanical validation;
- **image-sensitive C:** candidate signal with sign instability across plausible resolution strata;
- **validation-only:** numerical image signals not yet securely linked to a botanical quantity.

## Azami → EAzami: one vertical research line

Azami and EAzami are not parallel trait projects.

`Azami = present phenotypic field`

`continuous phenotype → repeated distribution → spatial hierarchy → frozen present-day geometry`

`↓`

`EAzami = generative-history discrimination`

`candidate function → interaction and fitness structure → historical models → admissible/rejected generators → next experiment`

Azami defines what currently exists and how it is organized. [`EAzami`](https://github.com/zuizui0223/EAzami) asks which declared histories and mechanism families can generate that frozen present. EAzami reproduction cannot promote Azami observations to causal results.

## Computational completion

Complete for the frozen Chapter 1 estimand:

- continuous-trait contract and final measurement workflow;
- primary nested variance and sampling sensitivities;
- expanded morphospace;
- within-taxon core-climate models;
- matched among-taxon models and environmental sorting;
- whole-capitulum within/among organization;
- process-environment extraction and core-sufficiency tests;
- seasonal, hemisphere, dominant-taxon, native-range, spatial, taxonomic and image-quality audits;
- evidence atlas and machine-readable EAzami handoff.

No additional uncontracted correlational raster screen is required to complete the current Chapter 1 question.

## Remaining scientific gates

The package is not submission-ready until the following are completed or formally dispositioned:

1. adjudicated human boxes for independent detector precision/recall;
2. independent reference measurements for orientation, colour and outline;
3. botanical calibration of expanded involucral/armature image phenotypes;
4. developmental-stage labels or an anthesis-restricted sensitivity;
5. paired flower/background visible-colour negative controls;
6. repeat-photo trait remeasurement and bounded variance decomposition.

Actual final multilevel figures and expanded Supplement exports must also be frozen after the submission narrative and validation outcomes are synchronized.

## Start here

1. [`manuscript/SUBMISSION_MANUSCRIPT.md`](manuscript/SUBMISSION_MANUSCRIPT.md) — canonical paper story and claim boundary.
2. [`manuscript/AZAMI_CH1_HYPOTHESIS_RECOVERY_AND_COMPLETENESS_2026-08-27.md`](manuscript/AZAMI_CH1_HYPOTHESIS_RECOVERY_AND_COMPLETENESS_2026-08-27.md) — hypothesis recovery and completion audit.
3. [`manuscript/GEB_V2_CROSS_SCALE_FINAL_RESULTS_2026-08-27.md`](manuscript/GEB_V2_CROSS_SCALE_FINAL_RESULTS_2026-08-27.md) — final multilevel result ledger.
4. [`manuscript/results/GEB_V2_FINAL_RESULTS_2026-08-26.md`](manuscript/results/GEB_V2_FINAL_RESULTS_2026-08-26.md) — final continuous-trait measurement/evidence-atlas ledger.
5. [`ch1_global/v2/ontology/ch1_continuous_trait_contract.csv`](ch1_global/v2/ontology/ch1_continuous_trait_contract.csv) — endpoint definitions.
6. [`analysis/ch1/evidence_tiering_policy.json`](analysis/ch1/evidence_tiering_policy.json) — evidence interpretation.
7. [`manuscript/EXTERNAL_COMPLETION_GATES.md`](manuscript/EXTERNAL_COMPLETION_GATES.md) — completed versus open validation gates.
8. [`analysis/ch1/pipeline.json`](analysis/ch1/pipeline.json) — executable submission map.

## Claim boundary

Azami Chapter 1 does **not** claim demonstrated:

- phenotypic plasticity, heritability or genetic variance;
- local adaptation or selection;
- pollinator or antagonist causation;
- rain, radiation, VPD, wind or productivity causation;
- defensive efficacy of image-derived involucral/armature geometry;
- functional or genetic modularity;
- evolutionary-rate differences or adaptive radiation;
- a resolved *Cirsium* species tree;
- floral-display structure from the unexecuted `visible_floret_fraction` endpoint.

Image vertical is a reproducible image reference rather than a direct inclinometer measurement of gravity.
