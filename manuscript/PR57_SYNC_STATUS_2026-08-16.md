# PR57 sync status — 2026-08-16

PR57 now synchronizes the Chapter 1 repository-facing story with the current submission manuscript and the EAzami handoff.

## Canonical Chapter 1 story

The Chapter 1 headline is no longer the precision-corrected lability analysis. The submission-facing story is:

> Global public photographs are converted into repeated **continuous head-level phenotype measurements** rather than categorical trait labels; these retain within-taxon trait distributions that are largely hidden by taxon means and reveal trait- and scale-specific environmental structure.

The former raw absolute-slope lability result and quadrants remain statistical provenance/QA only.

## Files synchronized in PR57

- `README.md`
- `manuscript/SUBMISSION_MANUSCRIPT.md`
- `manuscript/00_title_abstract.md`
- `manuscript/01_introduction.md`
- `manuscript/02_methods.md`
- `manuscript/03_results.md`
- `manuscript/04_discussion.md`
- `manuscript/05_conclusion_and_declarations.md`
- `manuscript/COHORT_FLOW_AND_ANALYSIS_LEDGER.md`
- `manuscript/final_claims.json`
- `manuscript/EXTERNAL_COMPLETION_GATES.md`
- `analysis/ch1/pipeline.json`
- `manuscript/EAzami_HANDOFF_AND_CHAPTER_BOUNDARY_2026-08-16.md`
- `analysis/export_eazami_trait_handoff.py`
- `tests/test_export_eazami_trait_handoff.py`

## Resolved scope mismatches

1. Taxonomic, residual-spatial, broad-region and niche-null diagnostics are now marked complete; only detector human boxes and independent orientation/colour/outline reference measurements remain scientific submission blockers.
2. The high-resolution involucre layer is represented consistently: roughness, spread fraction and maximum spine-like projection are the three final manuscript auxiliary inferential endpoints; `spine_peak_count_proxy` is retained for provenance/downstream exploration only.
3. Auxiliary involucre results remain in the main manuscript as an auxiliary extension of the continuous-phenomics framework, while their resolved historical interpretation is handed to EAzami.
4. The EAzami summary bridge is explicitly distinguished from a full within-taxon distribution bridge. Median/MAD taxon summaries are suitable for topology mapping and coverage, but observation-level exports are required whenever downstream inference concerns within-taxon structure.
5. The reason for focusing on *Cirsium* is aligned with the manuscript: conspicuous multi-module capitulum disparity in a lineage containing recent regional radiations and complex reticulate/phylogenetic history, without claiming that Chapter 1 itself demonstrates adaptive radiation.

## Scientific stop rule

No frozen coefficient, cohort count or accepted result is changed by this synchronization. The changes align manuscript framing, machine-readable claims, pipeline status and cross-repository handoff semantics with already-completed analyses.
