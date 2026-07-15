# Canonical Chapter 1 analysis entry point

This directory provides a stable submission-facing interface over the numbered
`ch1_global/v2` scripts that generated the frozen results.

The numbered scripts are intentionally **not moved or renamed**. Their paths are
part of the analysis provenance recorded in workflows, artifacts and manifests.
New users should start here rather than selecting scripts by number.

## Commands

Verify that all canonical scripts and manuscript controls exist:

```bash
python analysis/ch1/run_submission.py check
```

Run the frozen two-axis lability stage from the pinned CHELSA observation and
coefficient tables:

```bash
python analysis/ch1/run_submission.py lability \
  --observations strict_spatial_thinned_with_chelsa.csv \
  --coefficients exhaustive_within_species_climate_models.csv \
  --out-dir results/lability
```

Validate the frozen manuscript claims:

```bash
python analysis/ch1/run_submission.py claims
```

Build the final checksum and software manifest after the complete submission
bundle has passed the existing output validator:

```bash
python analysis/ch1/run_submission.py manifest \
  --bundle-root submission_bundle \
  --validation-report submission_bundle/validation_report.json \
  --out-dir submission_bundle/provenance \
  --r-session-info submission_bundle/provenance/R_sessionInfo.txt
```

## Source of truth

- `pipeline.json` maps scientific stages to exact provenance scripts and outputs.
- `manuscript/final_claims.json` is the numeric claim registry.
- `manuscript/FINAL_MANUSCRIPT_STRATEGY.md` is the narrative source of truth.
- `manuscript/FIGURE_TABLE_MAP.md` maps frozen outputs to manuscript items.
- `manuscript/RUNBOOK.md` defines the release and archive procedure.

## Change policy

Changes that alter accepted taxa, trait definitions, thresholds, environmental
layers or tree scenarios require a new analysis version and full rerun. Changes
to wording, captions or file organization must not alter frozen result tables in
place.
