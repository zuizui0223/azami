# External completion gates for Chapter 1

Chapter 1 is on **submission hold for independent validation and final presentation**, not for additional uncontracted correlational analysis. The accepted computational estimand is complete through continuous-trait measurement, below-taxon diversity, expanded morphospace, matched within/among geometry and process-environment sufficiency. Completed diagnostics are rerun only if a frozen input, operational taxonomic unit, measurement definition or accepted model changes materially.

The current submission-facing story is a global reconstruction of **present continuous-trait geography in thistle capitula**: public photographs become repeated continuous measurements, each trait is paired with predeclared environmental gradients, and its support is classified across nested within- and among-taxon scales. Whole-capitulum organization is a secondary synthesis. The remaining scientific gates ask whether the image-derived measurements and below-taxon components are independently valid.

## 1. Independent detector audit — MANUAL ANNOTATION REMAINS

The earlier detector-development packet cannot serve as independent validation because it contributed proposals and pseudo-labelled training/validation images. The leakage-free replacement packet and evaluation infrastructure are complete. The workflow:

- excludes prior proposal/training `photo_id` and `obs_id` values;
- freezes 1,000 taxonomically and geographically spread source images;
- blinds taxonomy, coordinates and detector predictions;
- assigns 250 images to a second annotator;
- fixes the production weight hash and confidence threshold 0.25;
- generates hidden low-threshold predictions for later matching.

Precision, recall and F1 remain unreported until independent boxes are completed and disagreements adjudicated. No repository-only automated re-interpretation can substitute for this gate.

## 2. Independent primary-trait validity — REFERENCE DATA REMAIN

Horizontal-mirror stability and production overlays are technical checks, not independent accuracy validation. Genuinely external reference measurements are required for:

- orientation: human landmarks/axes and camera-roll assessment;
- colour: calibration-supported references or a justified repeat-image error design;
- outline: independent manual masks/contours and continuous-trait agreement.

A material failure requires a new analysis version and uncertainty propagation into the headline sensitivities. This is the most important methodological gate because Chapter 1 treats each assessable head as a numerical phenotype observation.

## 3. Botanical calibration of candidate geometry — REFERENCE DATA REMAIN

The expanded contract produced two quality-robust C-grade candidate rows and one image-sensitive C-grade row, but image-derived involucre and armature measurements are not yet calibrated against direct botanical measurements.

Required reference work should compare image endpoints with quantities such as direct involucre dimensions, bract-tip projection/spreading geometry and independently traced contours. It should estimate agreement and failure conditions rather than assign categorical states after inspection. Until then, candidate endpoints remain image phenotypes and cannot be named as actual spine length, phyllary recurvature, stiffness or defensive efficacy.

## 4. Developmental stage — IMAGE LABELLING REMAINS

Raw-calendar and hemisphere-aware cyclic analyses retained the established non-circular directions, but calendar timing does not distinguish bud, anthesis and post-anthesis heads. A stage-labelled or anthesis-restricted sensitivity remains submission blocking, especially for orientation and outline.

## 5. Visible-colour negative control — IMAGE REMEASUREMENT REMAINS

The paired flower-versus-background `region × BIO12` implementation is tested and fixed. It measures outer-context and nonfloral-green controls and does not infer flower specificity from a nonsignificant background slope. Production-scale context-image remeasurement remains necessary. Chroma–BIO12 is already B-grade/range-sensitive rather than A-grade, but the control remains required for broader corolla-specific colour interpretation.

## 6. Repeat photographs — COHORT FROZEN, MEASUREMENTS REMAIN

An outcome-blind preflight identified 20,073 strict-cohort observations from 236 taxa with at least two photographs: 58,748 photo rows and 38,675 photographs beyond the first. Traits must be remeasured on this cohort. Between-photo variance will still combine camera, illumination, viewpoint and possible subject differences unless same-individual status is reviewed, but it will provide a bounded observation-layer sensitivity unavailable from the one-photo stream.

## 7. Primary nested variance — COMPUTATION COMPLETE, INTERPRETATION BOUNDED

For the nine established primary endpoints, taxon → photograph → head decomposition and one-head/equal-replication sensitivities are complete. The result establishes substantial visible variance below taxon means. It does not identify genetic variance or plasticity.

Observation-level below-taxon fractions exist for the expanded trait universe, but the full photograph/head decomposition has not been extended to every candidate endpoint. This is a secondary completeness gap; the Chapter 1 headline about nested visible variation remains explicitly limited to the primary nine.

## 8. Taxonomic robustness — COMPLETED FOR OPERATIONAL UNITS

Chapter 1 uses exact source-platform taxon assignments as operational units. A WCVP audit identified eight synonym candidates, which were collapsed simultaneously in sensitivity analysis. The collapse changed the balanced atlas from 216 to 211 units and the primary cohort from 259 to 251, but all eight primary BH-supported component rows remained supported, no coefficient changed sign and below-unit variance changed little.

The scientific taxonomic-robustness gate is closed for the operational-unit analysis. Twenty-four high-priority rows still require nomenclatural notes but are no longer an untested explanation of the ecological patterns.

## 9. Spatial and region diagnostics — COMPLETED

Observation-level residual export reproduced the frozen best-WAIC SPDE specifications. Across the nine primary endpoints, global Moran's I ranged from −0.0096 to 0.0030 and no permutation P was below 0.05. Leave-one-region-out taxon-rank Spearman minima ranged from 0.856 to 0.972. This gate is closed unless cohort or model specification changes materially.

## 10. Environmental-sorting nulls — COMPLETED

The established nine-trait environmental sorting and the expanded common-taxon sorting used 10,000 trait-label permutations. The matched 17-unit × four-core-predictor analysis and 68-row scale classification are complete. No additional all-BIOCLIM fishing screen is required for the current Chapter 1 estimand.

## 11. Process-environment sufficiency — COMPLETED

Shortwave radiation, VPD, wind, growing-season precipitation and potential NPP were extracted with complete coverage for the 1,874-observation complete-18 cohort. The six stand-alone blocks and nested comparisons against BIO1/BIO4/BIO12/BIO15 are complete.

The frozen decision is scale specific: the four-variable core is adequate for the current within-taxon 18D estimand but does not exhaust among-taxon environmental structure. Growing-season precipitation is the only stable block-specific increment across both replication thresholds. This is an observational redundancy result, not direct environmental causation.

## 12. Historical sensitivity — COMPLETED AS SENSITIVITY, HISTORY UNRESOLVED

Only 54 of 216 atlas taxa are direct dated-backbone tips; most taxa require placement and Pagel-λ estimates were zero. The randomized-tree PGLS layer is complete as a placement sensitivity but cannot be promoted to resolved phylogenetic correction. A resolved species history, ploidy and hybridization information belong to the downstream generative-history programme.

## 13. Final figure and Supplement build — PRESENTATION REMAINS

The final narrative now requires six main figures. Figures 1–3 retain the workflow, geography and nested variance. Figures 4–6 must be rendered or rebuilt to show:

- the 18-endpoint taxon morphospace;
- within/among 17-unit association matrices and module contrasts;
- endpoint evidence grades, matched scale classes and core-four sufficiency.

The expanded Supplement must export Tables S1.15–S1.22 and Figures S1.10–S1.17 from frozen artifact tables, followed by visual QA, manifest regeneration and output hashes. This is presentation completion, not a new scientific analysis.

## 14. Authorship, licence and durable release — HUMAN FINALIZATION REMAINS

Author order, affiliations, CRediT roles, acknowledgements, funding and corresponding-author details require team confirmation. After validation outcomes are final, the project must select licences, create an immutable release tag, generate SHA-256 manifests and deposit a licence-safe bundle in a durable repository with a persistent DOI. GitHub Actions artifacts are not the durable archive.

## Current completion decision

### Computational pattern estimand — COMPLETE

Completed:

- 27-endpoint contract and final 18-inferential-endpoint measurement layer;
- primary nested visible variance and sampling sensitivities;
- expanded taxon morphospace;
- primary/candidate within-taxon environment models;
- matched within/among models and expanded sorting;
- whole-capitulum measurement-module organization;
- process-environment extraction and core-four sufficiency;
- seasonal, hemisphere, dominant-taxon, native-range, image-quality, spatial and taxonomic audits;
- evidence atlas, hypothesis recovery and EAzami observational handoff.

### Submission-blocking scientific gates — SIX OPEN

1. adjudicated independent detector boxes;
2. independent orientation/colour/outline reference measurements;
3. botanical calibration of candidate involucre/armature image phenotypes;
4. developmental-stage labels or anthesis restriction;
5. paired flower/background colour controls;
6. repeat-photo trait remeasurement and bounded variance decomposition.

### Presentation and administrative work — OPEN

- final Figures 4–6 and expanded Supplement rendering/QA;
- nomenclatural notes;
- authorship and metadata;
- licence and durable DOI release.

The next scientific value lies in independent validation and discriminating measurements, not in adding another unregistered set of environmental correlations.
