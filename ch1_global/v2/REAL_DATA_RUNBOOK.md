# Ch.1 v2 — real-data annotation runbook

Run this after the existing global iNaturalist+YOLO crop metadata has been generated locally. Do not use the small fixture in `tests/` for biological annotation.

## 1. Build a locked manifest

```powershell
python ch1_global\v2\00_build_annotation_manifest.py `
  --input "C:\Users\zuizui\cirsium_inat\global_cirsium_flowering_annotation_then_yolo\flowering_annotation_yolo_crop_metadata.csv" `
  --out-dir "C:\Users\zuizui\cirsium_inat\ch1_v2_annotation\manifest" `
  --n-calibration 720 `
  --n-main 3000 `
  --n-segmentation 1000
```

Inspect `annotation_species_qc.csv` before annotation. The minimum per-species quota, the cap on common species, and the spatial-block width are part of the study design. If you change them, save the previous manifest and make a new dated output directory rather than overwriting it.

**Stop here** if the output contains fewer than 36 represented species, if a few species dominate the pool, or if image paths visibly include a taxon name. Resolve those before any human annotation.

## 2. Build blinded packets

```powershell
python ch1_global\v2\01_build_blinded_annotation_packets.py `
  --manifest "C:\Users\zuizui\cirsium_inat\ch1_v2_annotation\manifest\annotation_manifest.csv" `
  --out-dir "C:\Users\zuizui\cirsium_inat\ch1_v2_annotation\packets"
```

Give `annotator_1_packet.csv` and `annotator_2_packet.csv` to different annotators. The second packet contains the calibration subset only. Keep the original manifest private: it contains taxonomy, coordinates, and CV group fields.

## 3. Complete and compile the two calibration packets

Save completed copies without changing column names. Each annotator must mark a row `annotation_complete = yes` only when every assessable field has been considered.

```powershell
python ch1_global\v2\02_compile_double_annotations.py `
  --primary "C:\Users\zuizui\cirsium_inat\ch1_v2_annotation\packets\annotator_1_packet_completed.csv" `
  --secondary "C:\Users\zuizui\cirsium_inat\ch1_v2_annotation\packets\annotator_2_packet_completed.csv" `
  --out-dir "C:\Users\zuizui\cirsium_inat\ch1_v2_annotation\agreement"
```

Read these in order:

1. `annotation_validation_errors.csv` — fix invalid values or missing columns before interpreting agreement.
2. `annotation_agreement_summary.csv` — report the completed-pair count, assessable-pair count, percent agreement, and kappa for each trait. For orientation, colour, cover, and display, kappa is calculated only where both annotators regarded the trait as assessable.
3. `annotation_adjudication_queue.csv` — a blinded queue of every disagreement. A third annotator or predeclared adjudication rule supplies `adjudicated_value`.

## 4. Gate before training models

Do not train the orientation classifier or segmentation model from the raw first-pass labels. First freeze an adjudicated calibration table, then set aside the evaluation partitions using the manifest's `observation_group`, `species_cv_fold`, and `spatial_cv_fold`. A random image split is only a descriptive baseline and is never the reported generalization result.

## 5. What GitHub Actions checks

The `Ch1 annotation foundation CI` workflow is deliberately lightweight. On every relevant branch/PR change it compiles the scripts, runs the manifest/packet/adjudication integration test, creates a small smoke dataset, and saves it as an artifact for 14 days. It does not download iNaturalist images or train YOLO models.
