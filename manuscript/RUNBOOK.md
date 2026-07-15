# Chapter 1 submission runbook

## Canonical execution

1. Dispatch or update `.github/workflows/ch1-phase1-lability.yml`.
2. The workflow downloads the frozen strict spatial cohort.
3. Script 75 resamples the four CHELSA V2.1 variables and rebuilds the global
   within-species models.
4. Script 77 generates species-level variation and responsiveness indices,
   sensitivity analyses, tables and figures.
5. Confirm both the submission CI and lability workflow succeed.
6. Download `ch1-final-lability-results` and freeze it with the exact commit SHA.

## Required validation

- input observations: 46,276;
- input species: 259;
- primary quadrant analysis contains at least 20 species;
- all required CSV, PNG, PDF and JSON outputs are non-empty;
- the primary minimum sample size is 10;
- sensitivity runs use minimum sample sizes 5 and 20;
- hue remains a joint sine/cosine endpoint;
- no claim uses the words plasticity, adaptation or evolutionary rate as an
  inferred result.

## Submission bundle

The durable bundle must include:

- final CSV tables;
- PNG and PDF figures;
- `analysis_metadata.json`;
- CHELSA source metadata and URLs;
- Python and R package versions where applicable;
- `sessionInfo()` for R-derived outputs;
- Git commit SHA;
- workflow run ID;
- SHA-256 checksum manifest;
- README describing each file and its role in the manuscript.

GitHub Actions artifacts are temporary working products. The final bundle should
be deposited in a durable archive such as Zenodo when the manuscript is frozen.
