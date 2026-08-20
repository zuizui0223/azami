# Supplement status

The Supplement is part of the active submission branch. Frozen summary tables and the full grouped-SPDE machine-readable tables are committed under `tables/`; they are assembled from pinned workflow artifacts rather than hand-edited. `SUPPLEMENT_MACHINE_TABLE_PROVENANCE.json` records the source artifact IDs and SHA-256 checksums.

Generate the five Supplement figures from the committed frozen summary tables with:

```bash
python analysis/build_supplement_figures.py \
  --tables-dir manuscript/supplement/tables \
  --out-dir manuscript/supplement/figures
```

The figure builder is presentation-only: it does not refit models or alter cohorts. It writes PNG and SVG versions of Figures S1–S5 and is smoke-tested by the submission CI.

Rebuild the machine-readable Supplement package from the pinned SPDE and precision-audit artifacts with `analysis/assemble_supplement_machine_tables.py` via `.github/workflows/ch1-assemble-supplement-machine-tables.yml`. The committed S06 family contains the complete grouped-SPDE fixed effects, stability summaries, hyperparameters, predictor selection and run metadata.

The 628,518-byte taxon-specific slope table used solely for the withdrawn-lability precision audit remains **release-only provenance** and is uploaded in the Supplement machine-table workflow artifact rather than committed as active biological evidence.

The Supplement must remain synchronized with `manuscript/final_claims.json`, `manuscript/COHORT_FLOW_AND_ANALYSIS_LEDGER.md` and `analysis/ch1/pipeline.json`.
