# Consolidated native V4 evidence

Numerical run: https://github.com/zuizui0223/azami/actions/runs/39428065666
Numerical commit: `e5db09e102d6d5bf48038ee8469008bee033a69b`.
Focused verification: 28 tests passed. All numerical stages completed.

The archived JSON/CSV files are aggregate artifacts, preserved byte-for-byte.
The input is the authenticated, corrected native observation table, SHA-256
`342d78e9a8a960e8d9bdcf0bc2cdacd13bee48214a8b962419d5a71d5aef46ed`.
Raw input source identities and transformations are in
`environment_decoding_receipt.json`. Raw private records are not redistributed.

The pipeline reuses previous measured values and the verified native roster.
It does not remeasure photos, establish physical colour/angle calibration,
or reproduce frozen V3 on a new cohort without methodological changes.

The result families and interpretation limits are specified in the environment,
model, integration and space contracts one directory above `results/`.
The report renderer checks input and model hashes across stages and scopes its
interpretation to this numerical commit.

To render these archived aggregates on Actions, the report-only workflow uses
run `39428065666`. The same renderer can consume this directory directly:

```bash
python -m analysis.v4.render_native_report \
  --source analysis/v4/results/20261020_native_integration \
  --output work/native-report
```

No refit is performed by rendering. This command documents the implementation;
the project workflow executes numerical work and rendering in GitHub Actions.

## Supported interpretation

- The adjusted rainfall–orientation association is positive between species;
  the within-species interval includes zero. Fixed spatial balancing retains
  the positive between-species direction.
- Colour associations differ by component and scale. Species occupying higher
  radiation environments are brighter and less chromatic; within species,
  warmer environments are associated with brighter, less chromatic flowers.
- Equal-n RV diagnostics retain a positive median between-minus-within
  difference before and after fixed environmental adjustment. These are
  descriptive scale comparisons, not latent genetic covariance estimates.
- The current bract SI joint tests do not show clear environmental associations.
- Coarse residual-cell diagnostics do not establish fine-scale spatial
  independence. Species bootstrap/cluster uncertainty is not phylogenetic or
  photographer adjustment. Adjusted integration intervals condition on estimated
  environmental coefficients.

A mechanistic adaptation claim requires independent biological validation;
these outputs provide phenotype/environment patterns and testable hypotheses.
