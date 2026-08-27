# Canonical Chapter 1 submission entry point

Use this directory to navigate the frozen continuous within-taxon image-phenomics analysis and its predeclared reviewer bias controls. The current package is on submission hold.

## Source of truth
- `pipeline.json` — stage map, rerun policy and remaining gates.
- `manuscript/final_claims.json` — frozen numerical claims.
- `manuscript/COHORT_FLOW_AND_ANALYSIS_LEDGER.md` — cohort names/counts and permitted analyses.
- `manuscript/SUBMISSION_MANUSCRIPT.md` — current narrative and figure priorities.
- `manuscript/AZAMI_CH1_SPATIAL_SCALE_CANONICAL_2026-08-27.md` — PR #71 trait-by-gradient-by-scale mainline and EAzami temporal handoff.
- `v2_full27_environment_atlas_contract.json` — canonical 27-endpoint × nine-predictor analysis and gate contract.
- `v2_full27_sampling_composition_sensitivity_contract.json` — independent retrospective contract for equal-taxon-weight, dominant-taxon, broad-region and native-only direction checks.
- `../../manuscript/AZAMI_CH1_V2_DEFENSE_2026-08-27.md` — frozen v2 method/ecology defence and candidate-mechanism claim boundaries.
- `manuscript/results/GEB_V2_FULL27_FULL_ENVIRONMENT_RESULTS_2026-08-27.md` — validated full-27 result and provenance ledger.
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
