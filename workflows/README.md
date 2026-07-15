# Workflows

The executable GitHub Actions files remain in `.github/workflows/`, as required by
GitHub. This directory documents which workflows are authoritative.

## Submission-facing

- `.github/workflows/ch1-phase1-lability.yml`
  rebuilds the strict CHELSA table, reruns the global within-species models and
  generates the final species-level lability bundle.
- the Chapter 1 submission CI validates deterministic code and output contracts.

## Research-scale workflow policy

Heavy workflows are manual or path-gated. One-off dated workflows are retained
only for provenance and must not be cited as the current runbook. The current
workflow must rebuild its own upstream statistical inputs rather than silently
mixing artifacts produced by different commits.
