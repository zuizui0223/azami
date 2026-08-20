# Supplement status

The Supplement is part of the active submission branch. Frozen summary tables are committed under `tables/`; large release-only machine-readable tables are copied from the pinned workflow artifacts into the final checksummed archive rather than hand-edited.

Generate the five Supplement figures from the committed frozen summary tables with:

```bash
python analysis/build_supplement_figures.py \
  --tables-dir manuscript/supplement/tables \
  --out-dir manuscript/supplement/figures
```

The figure builder is presentation-only: it does not refit models or alter cohorts. It writes PNG and SVG versions of Figures S1–S5. The generator has been smoke-tested against the frozen Supplement tables.

The Supplement must remain synchronized with `manuscript/final_claims.json`, `manuscript/COHORT_FLOW_AND_ANALYSIS_LEDGER.md` and `analysis/ch1/pipeline.json`.

Release-only full tables include the complete grouped-SPDE fixed effects/hyperparameters and the archived taxon-specific slope table used solely for the withdrawn-lability precision audit. The latter is provenance, not biological headline evidence.
