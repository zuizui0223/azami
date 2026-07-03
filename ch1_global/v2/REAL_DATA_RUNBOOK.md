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

## 3. Annotate in the local browser app

Install the annotation-only dependencies once:

```powershell
python -m pip install -r ch1_global\v2\requirements_annotation_app.txt
```

Run one app instance per annotator and give each instance its own output directory. Do not point two annotators at the same response directory.

```powershell
streamlit run ch1_global\v2\03_trait_annotation_app.py -- `
  --packet "C:\Users\zuizui\cirsium_inat\ch1_v2_annotation\packets\annotator_1_packet.csv" `
  --out-dir "C:\Users\zuizui\cirsium_inat\ch1_v2_annotation\responses\annotator_1" `
  --annotator "annotator_1"
```

The app shows the head crop and source image, but not taxon, coordinates, model probability, or raw metadata. It saves `annotation_responses.csv` atomically after every Save action and can be stopped and restarted without losing progress. Complete all fields using either a biological state or `not_assessable`; a row becomes complete only when every trait has been considered.

The first app must annotate all 3,000 heads. The second app is run on the 720-head calibration packet only. This app is for **human-scored trait labels and image-quality flags**; detector boxes, keypoints, and segmentation polygons are separate annotation phases.

## 4. Compile the two calibration packets

Use the two response files created by the app. Do not rename or edit their columns manually.

```powershell
python ch1_global\v2\02_compile_double_annotations.py `
  --primary "C:\Users\zuizui\cirsium_inat\ch1_v2_annotation\responses\annotator_1\annotation_responses.csv" `
  --secondary "C:\Users\zuizui\cirsium_inat\ch1_v2_annotation\responses\annotator_2\annotation_responses.csv" `
  --out-dir "C:\Users\zuizui\cirsium_inat\ch1_v2_annotation\agreement"
```

Read these in order:

1. `annotation_validation_errors.csv` — fix invalid values or missing columns before interpreting agreement.
2. `annotation_agreement_summary.csv` — report the completed-pair count, assessable-pair count, percent agreement, and kappa for each trait. For orientation, colour, cover, and display, kappa is calculated only where both annotators regarded the trait as assessable.
3. `annotation_adjudication_queue.csv` — a blinded queue of every disagreement. A third annotator or predeclared adjudication rule supplies `adjudicated_value`.

## 5. Gate before training models

Do not train the orientation classifier or segmentation model from the raw first-pass labels. First freeze an adjudicated calibration table, then set aside the evaluation partitions using the manifest's `observation_group`, `species_cv_fold`, and `spatial_cv_fold`. A random image split is only a descriptive baseline and is never the reported generalization result.

## 6. What GitHub Actions checks

The `Ch1 annotation foundation CI` workflow is deliberately lightweight. On every relevant branch/PR change it compiles the scripts, runs the manifest/packet/adjudication/app state integration tests, creates a small smoke dataset, and saves it as an artifact for 14 days. It does not download iNaturalist images or train YOLO models.
