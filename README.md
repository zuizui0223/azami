# azami — present spatial geography of continuous thistle capitulum traits

> **Submission status: scientific HOLD for independent validation.** The frozen v2 computational analysis is complete. Independent detector and trait-reference validation, developmental-stage and colour controls, repeat-photo remeasurement, licence selection and durable DOI release remain open.

This `main` branch is the clean submission-facing tree for:

**Global image phenomics reveals hidden variation and scale-specific environmental geography in thistle capitula**

## Result in one line

Public photographs recover repeated continuous phenotypes that taxon means compress; environmental alignment differs by endpoint and biological scale, and only chroma-radiation and orientation-annual-precipitation pass the current sequential among-taxon spatial and historical-placement gates while retaining direction under declared sampling perturbations.

## Why thistles

Thistles are the comparative model system, not the genus-specific endpoint of the paper. The Carduus–Cirsium group places multiple visible phenotype dimensions on one repeatedly photographed capitulum and contains geographically replicated rapid regional radiations, including separate Pleistocene radiations in Japan and North America. That combination allows the paper to ask whether repeated phenotypic diversification across distinct regional lineages is lost when morphology is reduced to categories or taxon means. The current chapter does **not** assume that every regional radiation is a demonstrated adaptive radiation and does not infer repeated origins or convergence from present spatial associations.

## Frozen v2 lane

- 27 continuous endpoints registered;
- 22 endpoints measured and five retained as explicit unexecuted rows;
- 46,276 spatially thinned observations from 259 source-assigned taxa;
- 81.4-98.6% of raw visible variation below taxon means;
- 28.7-58.5% below taxon means after equalizing to two observations per eligible taxon;
- nine environmental predictors in six predeclared blocks;
- seven both-scale, 18 within-only, three among-only, 152 neither and 54 not-comparable rows;
- 674 sampling-composition scenarios, all evaluable;
- five within-taxon broad-space passes, four also directionally stable throughout sampling perturbations;
- two among-taxon candidates passing broad-space control, all relevant sampling directions and all 52 placement trees.

The public-image sample is extensive but uneven: Europe plus North America contributed 92.26% and the two most observed taxa 54.03%. The sampling audit annotates that dependence; it does not convert the data into a probability sample.

## Main versus Supporting Information

Main contains the general inferential story in five figures and two tables: a compact analytical roadmap plus the exact evidence chain for the two candidate patterns. Supporting Information is the audit layer: complete endpoint coverage, full supported/unsupported model grids, all sampling perturbations, spatial diagnostics, all tree-level historical-placement results and the full roadmap. Main Figures 1–5 are not reused as Supporting Figures; S1–S4 are distinct coverage, sampling, spatial-diagnostic and historical-placement displays, and Tables S1–S10 are supplied as a machine-readable workbook.

## Candidate boundaries

Higher-radiation taxa had lower camera-recorded CIELAB corolla chroma. Radiation-responsive pigment regulation, including anthocyanin synthesis, is a candidate hypothesis, but chroma is not pigment concentration or UV reflectance.

Higher-annual-precipitation taxa had a larger angle on the frozen 0-up, 90-horizontal, 180-down EXIF-image scale. Rain shielding is a candidate hypothesis, but image vertical is not gravity and rain interception, pollen viability and fitness were not measured.

Both are **adaptive-pattern candidates under current controls**, not demonstrated adaptations or mechanisms.

## Start here

1. [`manuscript/SUBMISSION_MANUSCRIPT.md`](manuscript/SUBMISSION_MANUSCRIPT.md) — active paper order, model-system rationale and Main/SI boundary.
2. [`manuscript/results/GEB_V2_FULL27_FULL_ENVIRONMENT_RESULTS_2026-08-27.md`](manuscript/results/GEB_V2_FULL27_FULL_ENVIRONMENT_RESULTS_2026-08-27.md) — numerical result ledger.
3. [`manuscript/AZAMI_CH1_V2_DEFENSE_2026-08-27.md`](manuscript/AZAMI_CH1_V2_DEFENSE_2026-08-27.md) — method and ecological defence.
4. [`manuscript/final_claims.json`](manuscript/final_claims.json) — machine-readable active claims.
5. [`analysis/ch1/pipeline.json`](analysis/ch1/pipeline.json) — canonical executable stages.
6. [`analysis_outputs/v2_full27_environment_atlas_2026-08-27/`](analysis_outputs/v2_full27_environment_atlas_2026-08-27/) — frozen v2 tables.
7. [`submission/`](submission/) — blinded Word/PDF files and the hashed submission-candidate bundle.

## Reproduce the submission-facing checks

```text
python analysis/validate_final_claims.py
python analysis/validate_ch1_hypothesis_recovery.py
python analysis/validate_manuscript_citations.py
python analysis/ch1/run_submission.py check
python analysis/validate_v2_submission_artifacts.py
python analysis/build_v2_submission_figure_inputs.py
python analysis/build_v2_submission_figures.py
python analysis/build_v2_submission_tables.py
python analysis/build_v2_submission_package.py
```

The old lability/quadrant result, PR69-era nine-endpoint hierarchy and superseded v1 output directories are absent from the active tree. They remain recoverable in Git history and the immutable scientific tag `azami-ch1-v2-2026-08-27`.

## Claim ceiling

This repository does not demonstrate plasticity, heritability, selection, causal environmental mechanisms, botanical identity of validation-only image proxies, a resolved *Cirsium* species tree, independent origins, adaptive radiation or convergence.
