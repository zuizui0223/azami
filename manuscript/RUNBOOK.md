# Chapter 1 submission runbook

## 1. Control checks
```bash
python analysis/ch1/run_submission.py check
python analysis/ch1/run_submission.py claims
python analysis/ch1/run_submission.py summary
```

## 2. Frozen result rule
Do not rerun or alter accepted scientific inputs during repository cleanup. Exact numbered `ch1_global/v2` paths remain stable. Rerun only under the conditions declared in `analysis/ch1/pipeline.json`.

## 3. Remaining scientific gates
- adjudicated independent human boxes for the 1,000-image detector audit;
- independent reference measurements for orientation, colour and outline.

A material failure creates a new analysis version and a full downstream rerun. Internal overlays, pseudo-label metrics and mirror repeatability cannot substitute for these gates.

## 4. Completed robustness gates
Taxonomic synonym-collapse sensitivity, residual Moran/broad-region omission and 10,000-permutation niche-null tests are complete for the frozen operational-unit analysis. For the full-27 lane, the retrospective sampling-composition audit is also complete: it evaluates equal-taxon weighting, dominant-taxon omission, leave-one-region-out and native-only subsets for the preselected q-supported rows, reports direction stability without new p/q values and is validated by `analysis/validate_geb_v2_full27_sampling_composition.py`.

## 5. Supplement and release
Build the Supplement from the same frozen tables as the manuscript, include complete machine-readable coefficient/model tables, then create a final checksummed release manifest only after the external measurement gates close. Archive source identifiers/licence metadata without redistributing photographs beyond licence terms.

## Forbidden resurrection
Do not restore the raw variation–responsiveness rho=-0.333 or median-split quadrants as a biological conclusion; do not call cross-sectional associations climate tracking, plasticity, adaptation or evolutionary rate.
