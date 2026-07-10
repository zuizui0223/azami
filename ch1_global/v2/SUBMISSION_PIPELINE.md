# Chapter 1 submission pipeline

This is the canonical runbook for the continuous-trait manuscript analysis.
Older categorical CLIP holdout and categorical trait-integration paths are not
part of the submission.

## 1. Inputs

The pipeline requires four frozen input groups:

1. YOLO-detected head metadata and crop URLs;
2. reconstructed head/context crops;
3. accepted taxon and observation metadata;
4. CHELSA observation-level climate values.

Every final run must record input file SHA-256 checksums and the Git commit. Do
not use mutable filenames without a manifest.

## 2. Continuous primary traits

Run the corrected production entrypoint:

```bash
python ch1_global/v2/56_run_primary_traits_continuous_v2.py \
  --manifest INPUT/continuous_trait_manifest.csv \
  --packet-root INPUT/rebuilt_packet \
  --metadata INPUT/global_ai_merged_yolo_crop_metadata.csv \
  --out-dir work/continuous_head

python ch1_global/v2/53_aggregate_primary_trait_continuous.py \
  --measurements work/continuous_head/primary_trait_continuous_head_measurements.csv \
  --metadata INPUT/global_ai_merged_yolo_crop_metadata.csv \
  --out-dir work/continuous_aggregated
```

Use only `*_status=usable` values. QC failure is missing measurement, not trait
absence.

## 3. QC-retention and environment join

```bash
python ch1_global/v2/57_qc_bias_join_primary_continuous.py \
  --head work/continuous_aggregated/primary_trait_continuous_head_with_metadata.csv \
  --observation work/continuous_aggregated/primary_trait_continuous_observation_level.csv \
  --species work/continuous_aggregated/primary_trait_continuous_species_level.csv \
  --environment INPUT/global_ai_observation_environment_metadata.csv \
  --out-dir work/qc_environment_join
```

The QC GEE results diagnose environmentally selective assessability. They are
not trait–environment effects.

## 4. High-resolution involucre supplement

```bash
python ch1_global/v2/58_prepare_high_resolution_involucre_subset.py \
  --yolo-metadata INPUT/global_ai_merged_yolo_crop_metadata.csv \
  --out-csv work/involucre/high_resolution_yolo_metadata.csv \
  --report work/involucre/high_resolution_subset_report.json \
  --min-bbox-dimension 150

python ch1_global/v2/54_rebuild_global_continuous_trait_packet.py \
  --yolo-metadata work/involucre/high_resolution_yolo_metadata.csv \
  --out-dir work/involucre/rebuilt \
  --context-pad-ratio 1.5

python ch1_global/v2/59_measure_involucre_auxiliary_continuous.py \
  --manifest work/involucre/rebuilt/continuous_trait_manifest.csv \
  --packet-root work/involucre/rebuilt \
  --metadata work/involucre/high_resolution_yolo_metadata.csv \
  --out-dir work/involucre/measurements \
  --min-dimension 150 \
  --min-sharpness 55 \
  --min-mask-quality 0.55
```

The output variables are image-contour proxies. Do not relabel them as
categorical bract orientation or confirmed spine states.

## 5. Evidence queues

```bash
python ch1_global/v2/60_prepare_microtrait_literature_queue.py \
  --species work/qc_environment_join/primary_continuous_species_environment.csv \
  --out-csv work/evidence/hair_mucilage_evidence_queue.csv \
  --schema-json work/evidence/hair_mucilage_evidence_schema.json
```

Ploidy/hybrid evidence is prepared later by script 67. Unknown values must remain
unknown until a reviewed source, page/table/figure and quotation are present.

## 6. Integrate and run pre-tree models

```bash
python ch1_global/v2/62_join_primary_auxiliary_continuous.py \
  --primary-observation work/qc_environment_join/primary_continuous_observation_environment.csv \
  --primary-species work/qc_environment_join/primary_continuous_species_environment.csv \
  --primary-registry work/qc_environment_join/primary_continuous_endpoint_registry.csv \
  --auxiliary-head work/involucre/measurements/involucre_auxiliary_head_level.csv \
  --out-dir work/integrated

python ch1_global/v2/63_run_prephylogenetic_continuous_environment_models.py \
  --observation work/integrated/integrated_primary_auxiliary_observation.csv \
  --species work/integrated/integrated_primary_auxiliary_species.csv \
  --registry work/integrated/integrated_continuous_endpoint_registry.csv \
  --out-dir work/prephylogenetic_models
```

Interpretation hierarchy:

- within-species, coordinate-precision <=10 km: strict local-association check;
- within-species, full coordinate cohort: sensitivity only;
- between-species OLS: screening only until historical sensitivity is complete.

## 7. Historical constraints

```bash
Rscript ch1_global/v2/65_audit_and_build_historical_trees.R \
  --species work/integrated/integrated_primary_auxiliary_species.csv \
  --out-dir work/historical_trees \
  --n-random 50 \
  --seed 20260710

Rscript ch1_global/v2/66_run_historical_constraint_pgls.R \
  --species work/integrated/integrated_primary_auxiliary_species.csv \
  --registry work/integrated/integrated_continuous_endpoint_registry.csv \
  --tree-dir work/historical_trees \
  --out-dir work/historical_models \
  --min-species 30

python ch1_global/v2/67_prepare_ploidy_hybrid_evidence_queue.py \
  --species work/integrated/integrated_primary_auxiliary_species.csv \
  --out-csv work/ploidy_hybrid/ploidy_hybrid_evidence_queue.csv \
  --schema-json work/ploidy_hybrid/ploidy_hybrid_evidence_schema.json

python ch1_global/v2/68_finalize_historical_sensitivity.py \
  --tree-audit-dir work/historical_trees \
  --model-dir work/historical_models \
  --out-dir work/historical_summary
```

Pagel lambda is bounded to [0, 1]. A random-tree result is complete only when all
50 tree replicates succeed.

## 8. Assemble the submission bundle

Create one directory containing the final tables, reports, figures, tree audits
and evidence queues. Do not include detector/model weights or temporary image
packets.

Example:

```bash
mkdir -p submission_bundle
cp -r work/qc_environment_join submission_bundle/
cp -r work/integrated submission_bundle/
cp -r work/prephylogenetic_models submission_bundle/
cp -r work/historical_trees submission_bundle/
cp -r work/historical_models submission_bundle/
cp -r work/historical_summary submission_bundle/
cp -r work/ploidy_hybrid submission_bundle/
cp -r work/evidence submission_bundle/
```

## 9. Validate and write provenance

```bash
python ch1_global/v2/69_validate_submission_outputs.py \
  --bundle-root submission_bundle \
  --config ch1_global/v2/submission_config.json \
  --out-json submission_provenance/submission_validation.json

Rscript -e 'writeLines(capture.output(sessionInfo()), "submission_provenance/R_sessionInfo.txt")'

python ch1_global/v2/70_build_submission_manifest.py \
  --bundle-root submission_bundle \
  --validation-report submission_provenance/submission_validation.json \
  --config ch1_global/v2/submission_config.json \
  --out-dir submission_provenance \
  --git-sha "$(git rev-parse HEAD)" \
  --r-session-info submission_provenance/R_sessionInfo.txt
```

Archive both `submission_bundle/` and `submission_provenance/`. A result is not
considered frozen without a passing validation report and checksummed manifest.

## 10. Versioning

Any change to a primary output table, taxon decision, endpoint definition, QC
threshold or tree scenario requires:

1. a new `analysis_version` in `submission_config.json`;
2. a full rerun of validation and manifest generation;
3. a new archive version rather than overwriting the submitted bundle.
