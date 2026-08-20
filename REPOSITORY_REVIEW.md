# Repository review for Chapter 1 submission

## Authoritative path
`README.md` → `manuscript/SUBMISSION_MANUSCRIPT.md` → `manuscript/final_claims.json` + `manuscript/COHORT_FLOW_AND_ANALYSIS_LEDGER.md` → `analysis/ch1/pipeline.json` → exact frozen scripts/artifacts.

The raw lability/quadrant result is withdrawn and must not appear in the active submission path.

## Keep
- exact `ch1_global/v2` frozen implementation paths;
- current `analysis/` audit and robustness scripts referenced by `pipeline.json`;
- current manuscript sections, claims, cohort ledger, results and supplement;
- submission utilities, tests and relevant workflows;
- image-trait dictionary and licence/provenance tables required by the final release.

## Remove from active submission branch
- Japan/pollinator and future-chapter trees;
- v1 categorical/orientation-only code;
- legacy model checkpoints and categorical direction classifier;
- old planning/review documents once no active links remain;
- workflows not reachable from the current pipeline or remaining validation gates.

Git history remains the archive. Do not move/rename executed v2 numbered scripts before durable release.
