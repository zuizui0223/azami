# Repository review for Chapter 1 submission

## Submission-critical path

The authoritative path is now:

```text
frozen strict observation cohort
  -> ch1_global/v2/75_run_exhaustive_within_species_climate_analysis.py
  -> ch1_global/v2/77_build_lability_axes_and_supplement.py
  -> ch1-final-lability-results
  -> manuscript/RUNBOOK.md
```

## Findings corrected in PR #34

- The first lability implementation plotted traits rather than species.
- Its responsiveness index used pooled global coefficients rather than
  species-specific slopes.
- Species with one measured trait could be ranked beside species with complete
  trait coverage.
- The workflow mixed a current script with coefficients downloaded from an older
  artifact.
- Hue magnitude was displayed without a clear distinction between descriptive
  joint magnitude and inferential scalar coefficients.

The final implementation rebuilds its own CHELSA table and pooled coefficients,
fits species-specific associations, imposes completeness rules and records
minimum-sample-size sensitivity.

## Keep active

- Chapter 1 v2 continuous measurement pipeline;
- strict spatial cohort and CHELSA join;
- QC-retention analyses;
- within-versus-between scale comparison;
- alternative-tree historical sensitivity;
- final two-axis lability analysis;
- submission CI and runbook.

## Move to the archive after output freeze

- superseded orientation-only analyses not cited in Methods;
- dated one-off workflows after their run IDs and checksums are recorded;
- intermediate lability figures generated before species-level completeness
  rules;
- duplicate reports and recovery scripts whose products are represented in the
  final provenance manifest.

## Keep exploratory and outside the headline claim

- Random Forest variable-importance analyses;
- involucre/spine contour proxies;
- mucilage and hair evidence queues;
- MaxEnt, hindcasting and palaeoclimate analyses;
- pollinator and adaptive-function hypotheses;
- unresolved ancestral-state claims.

## Remaining submission risks

1. Detector precision and false-negative audit must be reported.
2. Colour remains camera dependent and requires restrained interpretation.
3. Species-specific slopes may reflect geographic structure or unmeasured
   covariates rather than plasticity.
4. The grafted historical tree is a sensitivity device, not a resolved phylogeny.
5. Final species rankings must be checked for stability across minimum sample
   sizes before appearing in the main text.
