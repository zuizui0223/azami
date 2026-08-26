# Canonical Chapter 1 submission entry point

Use this directory to navigate the frozen continuous within-taxon image-phenomics analysis and its predeclared reviewer bias controls. The current package is on submission hold.

## Source of truth
- `pipeline.json` — stage map, rerun policy and remaining gates.
- `manuscript/final_claims.json` — frozen numerical claims.
- `manuscript/COHORT_FLOW_AND_ANALYSIS_LEDGER.md` — cohort names/counts and permitted analyses.
- `manuscript/SUBMISSION_MANUSCRIPT.md` — current narrative and figure priorities.
- `manuscript/EXTERNAL_COMPLETION_GATES.md` — genuinely external pending validation.
- `bias_control_contract.json` — locked raw-calendar seasonal, taxon-omission, image-quality and image-remeasurement rules.
- `hemisphere_season_sensitivity_contract.json` — separately locked half-cycle and hemisphere-specific seasonal sensitivity; it preserves rather than replaces the raw-calendar family.
- `native_range_sensitivity_contract.json` — separately locked WCVP/TDWG source versions, fail-closed mapping rules and native-only retention decision.

## Commands
```bash
python analysis/ch1/run_submission.py check
python analysis/ch1/run_submission.py claims
python analysis/ch1/run_submission.py summary
```

The runner intentionally has no `lability` command. The raw lability/quadrant result was withdrawn after precision-confounding audit and is not a submission headline.

Exact numbered scripts in `ch1_global/v2/` remain where they were executed; submission cleanup must not rename or move those provenance targets.
