# Azami Chapter 1 hypothesis recovery and analysis-completeness audit — 2026-08-27

## Purpose

This ledger separates four questions that had become mixed during the successive Chapter 1 revisions:

1. Which empirical hypotheses have been supported, rejected or left unresolved?
2. Which analyses are computationally complete?
3. Which scientific conclusions still depend on independent validation or new data?
4. Where does Azami end and the EAzami generative-history programme begin?

The audit does not redefine any estimand or reopen a completed multiplicity family. Numerical results remain tied to the immutable artifact-backed executions below.

## Frozen numerical sources

### Continuous-trait measurement and GEB v2 evidence atlas

- workflow run `32975451732`;
- analysis head `f4a6fd5e01a2befd4f49174984a99e53856c2330`;
- artifact `9612943217`;
- digest `sha256:101e996b638996a0c5ae79d358bf51293c3585f0e84c4a961b91dcbedf96211e`.

### Matched within/among capitulum-space and process-environment analysis

- workflow run `33035785120`;
- analysis head `227c0e7b8c338894806785b8545c7c77c8724de1`;
- artifact `9632715852`;
- digest `sha256:51e7a26b5bd09e030b67b9342586699abaaf46e630f45b6bb4ee7bfc9152ced6`;
- 10,000 permutations per inferential test;
- 1,000 taxon bootstraps for multivariate geometry.

## Status vocabulary

- **SUPPORTED** — the predeclared analysis directly supports the hypothesis within its stated observational scope.
- **PARTIALLY SUPPORTED** — one part, scale or evidence tier is supported while another is not.
- **REJECTED / WITHDRAWN** — the preregistered or audited result does not support the former claim and it is excluded from the Chapter 1 conclusion.
- **INCONCLUSIVE** — the implemented analysis does not distinguish the proposed alternatives with adequate uncertainty control.
- **NOT TESTED IN AZAMI** — the hypothesis is mechanistic, historical or functional and lies beyond the observational Chapter 1 estimand.
- **OPERATIONALLY COMPLETE, EXTERNALLY UNVALIDATED** — the computation is complete but the image-derived measurement still lacks genuinely independent reference validation.

# 1. Foundational Chapter 1 hypotheses

| Hypothesis | Recovery status | Evidence and interpretation |
|---|---|---|
| Public images can be converted into repeated continuous capitulum phenotypes at global scale | **OPERATIONALLY COMPLETE, EXTERNALLY UNVALIDATED** | The frozen trait universe contains 46,276 observations from 259 source-assigned taxa and 374,255 analysis-eligible observation × endpoint measurements. The pipeline is reproducible, but independent detector boxes and independent trait reference measurements remain open. |
| Taxon means conceal substantial visible phenotype diversity | **SUPPORTED for the nine established primary endpoints** | The fraction of visible image variance below taxon means is 0.589–0.931; one-head-per-photo and equal-replication sensitivities retain the result. This is visible image variance, not genetic variance or plasticity. |
| The same nested taxon → photograph → head decomposition is established for all expanded candidate endpoints | **PARTIALLY SUPPORTED / INCOMPLETE** | Observation-level below-taxon fractions exist for the expanded trait universe, but the full photograph/head nested decomposition remains completed only for the nine primary endpoints. |
| Finer continuous measurement reveals a multidimensional phenotype rather than one hidden master syndrome | **SUPPORTED** | In the 127-taxon × 18-endpoint morphospace, PC1 explains 18.49% and PC1–PC3 explain 42.28%. No single capitulum syndrome captures most taxon-level variation. |
| Environmental organization is trait specific rather than one warm-, wet- or dry-climate morphology | **SUPPORTED** | Established support is distributed across orientation, visible colour and selected outline traits; candidate involucral rows differ in quality grade. Most endpoint × environment rows are unsupported, which is part of the result rather than missing reporting. |
| Environmental organization is scale specific | **SUPPORTED** | Across 68 matched endpoint × core-climate rows: 3 are supported at both scales, 8 within only, 1 among only and 56 neither. Orientation aligns with BIO1 within taxa but BIO12 among taxa. |
| Fine involucral geometry contains environmental information missed by coarse categories | **PARTIALLY SUPPORTED** | Involucre length/width–BIO15 and apical taper–BIO12 are quality-robust C-grade within-taxon signals; `bract_projection_p95`–BIO15 is image-sensitive. No candidate endpoint passes matched among-taxon FDR. |
| A universal negative relationship exists between within-taxon variation and environmental responsiveness | **REJECTED / WITHDRAWN** | The former lability/responsiveness index was sample-size and uncertainty sensitive. The rho/quadrant story is provenance only and must not re-enter the manuscript conclusion. |

# 2. Recovery of the frozen whole-capitulum hypotheses H1–H5

The original pre-result contract is preserved in `CAPITULUM_FUNCTIONAL_SPACE_HYPOTHESES_2026-08-27.md`. This outcome ledger records recovery without rewriting that preregistration.

## H1 — Measurement-module organization exists at both biological scales

**Status: SUPPORTED.**

Main complete-18, >=5 observations/taxon scope:

- within-taxon registered-module contrast: `0.164502`, bootstrap 95% CI `0.130693–0.179475`;
- among-taxon registered-module contrast: `0.088475`, CI `0.024942–0.126171`.

The >=2 sensitivity gives the same qualitative result. The supported claim is partial phenotypic / measurement-module organization. It is not evidence that the modules are functional, developmental or genetic modules.

## H2 — Within- and among-taxon organization is only partially shared

**Status: SUPPORTED.**

- within-vs-among 17-unit association-matrix Spearman: `0.366299`;
- bootstrap 95% CI: `0.118457–0.399817`;
- >=2 sensitivity: `0.377272`.

The matrices are positively related but far from identical. Within-taxon organization is therefore not simply a small-scale copy of among-taxon organization.

## H3 — Environment should be represented by predeclared process blocks rather than a large post hoc BIOCLIM screen

**Status: DESIGN RECOVERED; BIOLOGICAL RESULT IS SCALE SPECIFIC.**

The six frozen blocks were extracted with complete coverage and analysed without expanding to all 19 BIOCLIM variables. Stand-alone block tests did not survive the six-block BH family. The more informative predeclared test was the nested core-four sufficiency comparison.

- within taxa: all five process variables beyond BIO1/BIO4/BIO12/BIO15 were unsupported (`partial R2 = 0.013535`, `P = 0.24468` at >=5);
- among taxa: the same extension was supported (`partial R2 = 0.214978`, `P = 0.000700` at >=5; supported again at >=2);
- growing-season precipitation was the only block-specific increment supported at both replication thresholds.

Therefore the four-variable core is adequate for the current within-taxon 18D estimand but does not exhaust among-taxon environmental structure.

## H4 — Environmental effect geometry may rotate across biological scales

**Status: INCONCLUSIVE AS A COEFFICIENT-VECTOR HYPOTHESIS.**

Cross-scale coefficient cosines were computed, but their bootstrap intervals include zero and they are retained as descriptive geometry only. The broader proposition that environmental organization differs across scales is supported by matched endpoint support, module contrasts and the nested environment-sufficiency result. The stronger directional claim that coefficient vectors rotate in a stable, interpretable way is not recovered.

## H5 — Azami should export phenotype-space targets, not causal mechanisms, to EAzami

**Status: SUPPORTED AS A WORKFLOW AND CHAPTER-BOUNDARY RESULT.**

The final artifact exports 62 provenance-gated observational targets:

- 6 structure targets;
- 24 environment-block R2 targets;
- 12 descriptive coefficient-geometry targets;
- 20 incremental environment targets.

EAzami imports these as observational, noncausal targets. Model reproduction cannot upgrade Azami correlations to adaptation, selection, plasticity, defence or direct climatic causation.

# 3. Functional hypotheses motivated by literature

| Functional hypothesis | Azami recovery status | Correct destination |
|---|---|---|
| Orientation, involucre architecture and outline jointly regulate rain, radiation, thermal or wind exposure | **NOT TESTED IN AZAMI** | EAzami structural-sufficiency models plus orientation manipulation, microclimate and reproductive-performance experiments |
| Colour and display form an advertisement axis affecting pollinators | **NOT TESTED IN AZAMI** | Calibrated colour/contrast, pollinator observations and preference experiments |
| Advertisement also increases florivore or predispersal seed-predator exposure | **NOT TESTED IN AZAMI** | EAzami interaction trade-off models and field antagonist data |
| Involucral/armature geometry has defensive efficacy | **NOT TESTED IN AZAMI** | Botanical calibration, stiffness/penetration measurements and antagonist exclusion/manipulation |
| Different trait combinations provide many-to-one functional solutions | **NOT TESTED IN AZAMI** | Mechanistic model and fitness-equivalence tests; current PCA/covariance alone cannot establish this |
| Registered modules are true functional or genetic modules | **INCONCLUSIVE / NOT IDENTIFIED** | Repeat-photo and developmental controls, independent botanical measurements, common-garden/genetic covariance and resolved history |

The functional literature is used to annotate hypotheses, not to convert image-derived associations into demonstrated functions.

# 4. Analysis-completeness audit

## A. Computational pattern layer — COMPLETE

The following accepted computational components are complete for the frozen operational-unit analysis:

- global image and cohort assembly;
- 27-endpoint category-free contract;
- final measurement run with 22 measured endpoints;
- 18/19 inferential endpoints measured;
- primary nested visible-variance decomposition and sampling sensitivities;
- 18-endpoint taxon morphospace and module-specific PCA;
- registry-driven within-taxon core-climate models;
- matched among-taxon models and 68-row cross-scale classification;
- expanded common-taxon environmental sorting;
- whole-capitulum within/among association geometry;
- CHELSA process-variable extraction;
- stand-alone environment-block tests;
- nested core-four sufficiency tests;
- seasonal, hemisphere, dominant-taxon and native-range sensitivities;
- residual-spatial, region-omission and taxonomic sensitivities;
- candidate resolution/sharpness and fixed-resolution-stratum audit;
- evidence atlas and EAzami machine-readable handoff.

No additional correlational raster screen is required to complete the current Chapter 1 estimand.

## B. Computational components that are intentionally incomplete or partial

- `visible_floret_fraction` is the only unexecuted inferential endpoint; no floral-display conclusion is permitted.
- Four descriptive colour-composition pixel fractions are unexecuted.
- The full photograph/head nested variance decomposition has not been extended from the primary nine endpoints to every candidate endpoint.
- Repeat-photo records are identified, but traits have not been remeasured across the frozen repeat-photo cohort.
- The phylogenetic analysis remains a randomized-placement sensitivity, not resolved phylogenetic correction.
- Soil/topographic predictors exist in the broader spatial layer, but a new all-environment whole-capitulum screen is not required for the frozen Chapter 1 claim and would be exploratory unless separately contracted.

## C. Manuscript and figure integration — PARTIAL AT THE START OF THIS AUDIT

Before this audit, the abstract, conclusion, canonical submission entry and final result ledger contained the multilevel H1–H5 outcomes, but the main `03_results.md`, `04_discussion.md`, README and active figure crosswalk still emphasized the earlier endpoint-level v2 story. Therefore the numerical analysis was ahead of the submission narrative.

The completion target of the present branch is to:

1. integrate the final module-organization, matched-scale and core-sufficiency results into Results;
2. interpret the whole-capitulum and scale-specific environment result in Discussion;
3. state the Azami → EAzami line as **present phenotypic field → generative-history discrimination**;
4. update repository status from “expanded run in progress” to “computational analysis complete; external validation open”;
5. make the hypothesis recovery ledger a canonical control file.

Actual final Figure 4–6 rendering and the expanded Supplement export remain a separate presentation build after the figure specification is synchronized.

## D. Independent scientific validation — OPEN AND SUBMISSION BLOCKING

The central remaining blockers are not additional model fitting. They are:

1. adjudicated human boxes for independent detector precision/recall;
2. genuinely independent reference measurements for orientation, colour and outline;
3. independent botanical calibration of candidate involucral/armature image phenotypes;
4. developmental-stage labels or an anthesis-restricted sensitivity;
5. paired flower-versus-background colour negative controls;
6. repeat-photo trait remeasurement and bounded image/biological variance decomposition.

A material failure in independent validation requires a versioned remeasurement and propagation of uncertainty. Until those gates are closed, the computational pattern is complete but the submission remains on hold.

# 5. Overall completion decision

## What is complete

The **Azami scientific pattern estimand is computationally complete**:

> repeated continuous image phenotypes → below-taxon diversity → multidimensional taxon morphospace → partial whole-capitulum organization → matched within/among environmental structure → scale-specific environmental sufficiency.

The final Chapter 1 synthesis is supported:

> **The thistle capitulum is a multidimensional, partially organized phenotype whose environmental alignment changes across biological scales.**

## What is not complete

Azami has not established:

- the biological-versus-photographic share of all below-taxon variation;
- direct functional meaning for the registered modules;
- plasticity, local adaptation or selection;
- pollinator or antagonist causation;
- direct rain, radiation, VPD, wind or productivity effects;
- a resolved evolutionary history.

These are not missing correlational analyses. They require independent measurement, experiments and historical data.

# 6. Vertical research line to EAzami

Azami and EAzami should be read as one directed programme rather than two parallel trait projects:

`Azami: present phenotypic field`

`continuous traits → repeated distributions → spatial hierarchy → frozen phenotype/environment geometry`

`↓`

`EAzami: generative-history discrimination`

`candidate function → interaction/fitness structure → historical models → rejected/admissible generators → next experiment`

Azami answers **what currently exists and how it is organized across space and hierarchy**. EAzami asks **which declared histories and mechanism families can generate that frozen present**. EAzami outcomes do not alter the Azami observation layer; they constrain the set of compatible generative histories.

# 7. Priority order from here

1. Finish submission-text and figure-specification synchronization without changing frozen numerical results.
2. Complete independent detector and trait validation.
3. Remeasure repeat photographs and label developmental stage.
4. Execute the paired colour negative control.
5. Calibrate candidate involucral/armature endpoints botanically.
6. Freeze the submission figures, Supplement exports and durable release only after material validation outcomes are known.
7. Continue EAzami as a separate generative-history discrimination layer, not as post hoc causal interpretation of Chapter 1.
