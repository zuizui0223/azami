# Chapter 1 submission runbook

## Canonical execution

1. Use the reviewed submission branch or a reviewed descendant of the frozen
   analysis commit.
2. Run `.github/workflows/ch1-phase1-lability.yml`.
3. Confirm both `ch1-phase1-lability` and `Ch1 submission code CI` succeed.
4. Run `.github/workflows/ch1-phylogenetic-evidence-summary.yml` only when the
   pinned phylogenetic-signal artifact or its decision rules change.
5. Run `.github/workflows/ch1-final-manuscript-claims.yml` and require success.
6. Download the validated lability, phylogenetic-evidence and claim outputs.
7. Record artifact digests, workflow run IDs and commit SHA before copying any
   manuscript table or figure.

The pinned climate artifact was generated from the frozen strict spatial cohort.
Rebuild it only when the observation cohort, CHELSA source or climate-extraction
code changes. Do not silently mix climate tables and coefficient files from
different runs.

## Frozen validation values

- strict within-species input observations: 46,276;
- strict within-species input taxa: 259;
- species with primary-complete lability axes: 102;
- primary minimum observations per trait-species combination: 10;
- minimum eligible traits per species: 6;
- minimum climate predictors per trait: 3;
- sensitivity thresholds: 5 and 20 observations;
- hue representation: joint sine/cosine endpoint;
- variation-responsiveness Spearman rho: -0.333;
- direct backbone tips: 54;
- phylogenetic-signal fits: 636 successful, 0 failed.

## Required lability outputs

- `analysis_metadata.json`;
- `species_lability_axes.csv`;
- `species_trait_within_variation.csv`;
- `species_trait_environmental_responsiveness.csv`;
- `species_module_lability.csv`;
- `module_lability_summary.csv`;
- `table_s1_all_global_coefficients.csv`;
- `table_s2_species_specific_coefficients.csv`;
- `sensitivity_minimum_sample_size.csv`;
- both PNG and PDF versions of the heatmap and four-quadrant figure.

## Required manuscript-control files

- `manuscript/FINAL_MANUSCRIPT_STRATEGY.md`;
- `manuscript/final_claims.json`;
- `manuscript/FINAL_FREEZE_CHECKLIST.md`;
- `analysis/validate_final_claims.py`;
- `ch1_global/SUBMISSION_READINESS.md`.

## Manuscript checks

Before submission, verify that:

- the text reports 102 complete taxa, not all 259 input taxa;
- every FDR count names its cohort and reporting scope;
- the strict <=10 km primary cohort is reported as zero BH-FDR-significant primary
  effects;
- the expanded pooled coefficient table is reported as four small
  BH-FDR-significant effects;
- the variation-responsiveness Spearman correlation is reported as -0.333;
- module summaries use medians and retain their sample sizes;
- hue magnitude is not given a scalar-angle p-value;
- phylogenetic signal is described as grafting-sensitive and not as a resolved
  species-history result;
- `plasticity`, `adaptation`, `selection` and `evolutionary rate` are not claimed as
  demonstrated results;
- all headline numbers match `manuscript/final_claims.json`.

## Remaining credibility gates

The core analysis is frozen. The following may still change presentation or
confidence, but should not trigger redesign of the core lability model:

1. measurement demonstration panel;
2. independent detector precision/recall audit;
3. accepted-name and synonym decision table;
4. residual spatial and broad-region sensitivity summary;
5. final captions, supplement and manuscript text;
6. permanent archive with licence-safe identifiers and checksums.

## Submission bundle

Create a durable bundle containing:

- final CSV tables and PNG/PDF figures;
- the pinned climate report and CHELSA source metadata;
- phylogenetic-evidence and molecular-database summaries;
- `analysis_metadata.json` and `final_claims.json`;
- Python package versions and all relevant R `sessionInfo()` files;
- Git commit SHA and workflow run IDs;
- SHA-256 checksum manifest;
- a file-level README mapping outputs to manuscript tables and figures;
- source-observation identifiers and licence fields without redistributing images
  beyond their licence terms.

GitHub Actions artifacts are temporary working products. Deposit the exact frozen
bundle in a durable archive such as Zenodo, OSF or Dryad before submission.
