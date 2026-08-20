# Chapter 1 manuscript workspace

This directory contains the active submission manuscript, frozen claim controls, result ledgers and Supplement. Do not copy numbers from legacy notes, superseded figures or workflow logs.

## Canonical entry point

1. `SUBMISSION_MANUSCRIPT.md` — current manuscript story, section order and figure priorities;
2. `00_title_abstract.md` through `05_conclusion_and_declarations.md` — active text sections;
3. `COHORT_FLOW_AND_ANALYSIS_LEDGER.md` — canonical cohort names, counts and permitted analyses;
4. `final_claims.json` — frozen machine-readable claims;
5. `FIGURE_TABLE_MAP.md` — current main and supplementary figure/table map;
6. `results/` — frozen result summaries used by the manuscript;
7. `supplement/` — Supplementary Information and machine-readable tables;
8. `DETECTOR_INDEPENDENT_AUDIT_PROTOCOL.md` and `FIGURE1_MEASUREMENT_AUDIT_PROTOCOL.md` — remaining independent-validation designs;
9. `EXTERNAL_COMPLETION_GATES.md` — completed versus external gates.

## Current scientific story

The submission is a global continuous within-taxon public-image phenomics analysis. Public photographs are converted into repeated numerical head-level measurements of orientation, visible colour and outline; those distributions are analysed at nested, within-taxon environmental, among-taxon and historical-sensitivity scales.

The former raw absolute-slope lability relationship and median-split quadrants were withdrawn after precision confounding was demonstrated. They remain only in the statistical QA/provenance record and are not part of the biological headline.

## Drafting rules

- Write numerical claims from `final_claims.json` and the current result ledgers.
- Name the exact cohort and multiplicity family for every FDR count.
- Keep the 6,626-head image atlas, 406,582-observation detector-positive layer and 46,276-observation spatially thinned primary cohort separate.
- Interpret hue sine/cosine jointly as a circular response.
- Treat visible image variance as image-level variation, not genetic variance.
- Do not convert cross-sectional associations into climate tracking, plasticity, adaptation, selection, pollinator causation or evolutionary rate.
- Treat the grafted historical tree as sensitivity analysis rather than a resolved Cirsium phylogeny.
- Do not report detector accuracy until the independent human audit is adjudicated.

## Reproducibility

Executable stages and rerun policy are registered in `../analysis/ch1/pipeline.json`. Submission controls are run with `../analysis/ch1/run_submission.py`. Exact executed numbered scripts under `../ch1_global/v2/` remain in place for provenance.
