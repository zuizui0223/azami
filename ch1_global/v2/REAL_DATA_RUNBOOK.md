# Ch.1 v2 — real-data collection and annotation runbook

This is the operational path from a worldwide iNaturalist sampling frame to blinded human trait labels. Large images, raw observation JSON, and completed annotation files remain **local research outputs**; do not commit them to Git.

The key v2 rule is: **collect metadata first, select candidate images second, run only the detector third, annotate traits fourth.** The old species-inherited upward/nodding classifier is never used to decide what enters the v2 dataset.

## 0. Install the two local environments

Metadata and queue construction:

```powershell
python -m pip install -r ch1_global\v2\requirements_data_collection.txt
```

Detector-only screening, on the machine that has the Cirsium YOLO weights:

```powershell
python -m pip install -r ch1_global\v2\requirements_yolo_screening.txt
```

## 1. Collect a pilot global metadata inventory

This uses the iNaturalist API to resolve the active `Cirsium` genus taxon, retrieve photo-bearing observations in stable ID order, and save raw response objects plus one row per photo. It retains flower annotations, licences, and obscured coordinates as metadata, but does **not** filter to flower-annotated records.

```powershell
python ch1_global\v2\04_collect_inat_cirsium_metadata.py `
  --out-dir "C:\Users\zuizui\cirsium_inat\ch1_v2_raw_pilot" `
  --max-observations 5000
```

Inspect these outputs before scaling up:

- `collection_provenance.json` — resolved taxon, endpoint, API parameters, timestamps, and restart state.
- `collection_summary.csv` — observation and photo counts.
- `photo_metadata.csv` — one photo per row, including taxonomy, licence, coordinates, coordinate usability, and iNaturalist flower annotation status.
- `observation_raw.ndjson` — raw API objects needed for later auditing.

The command is resumable: rerun the **same command with the same output directory** after an interruption. Do not use `--restart` in an existing directory.

## 2. Complete the global metadata inventory

After the pilot columns and taxonomic scope look correct, use a distinct output folder for the complete inventory. `0` removes the artificial observation limit while preserving the API delay and checkpoint.

```powershell
python ch1_global\v2\04_collect_inat_cirsium_metadata.py `
  --out-dir "C:\Users\zuizui\cirsium_inat\ch1_v2_raw_full" `
  --max-observations 0
```

Do not start image downloading until this metadata inventory has finished or reached a predeclared versioned cutoff. This protects the global species × spatial sampling frame from being determined by download order.

## 3. Build the first species × spatially balanced screening queue

The queue defaults to open, georeferenced, non-captive, species-rank records and retains only one photo per observation. It first spreads selection across species × 10-degree blocks, then fills remaining quota by the least represented species. Flower annotations are not selection criteria.

```powershell
python ch1_global\v2\05_build_image_screening_queue.py `
  --input "C:\Users\zuizui\cirsium_inat\ch1_v2_raw_full\photo_metadata.csv" `
  --out-dir "C:\Users\zuizui\cirsium_inat\ch1_v2_screen_queue_pilot" `
  --max-photos 5000 `
  --min-photos-per-species 10 `
  --max-photos-per-species 150
```

Inspect before download:

- `queue_species_qc.csv`: no one species should dominate merely because it is common on iNaturalist.
- `queue_spatial_qc.csv`: check that global blocks are represented rather than concentrated in North America, Europe, or a few Japanese sites.
- `queue_provenance.json`: save the seed and metadata SHA as the exact sampling record.

For the final detector pool, make a new dated queue with `--max-photos 25000` after the 5,000-image pilot confirms reasonable head-detection yield. Keep both queue versions; the pilot is a documented QC stage, not discarded history.

## 4. Download only the chosen image queue

Images are stored with photo-ID filenames, not species names, so later human annotation remains blind to taxonomy. Downloads resume by `queue_id`.

```powershell
python ch1_global\v2\06_download_inat_screening_queue.py `
  --queue "C:\Users\zuizui\cirsium_inat\ch1_v2_screen_queue_pilot\image_screening_queue.csv" `
  --out-dir "C:\Users\zuizui\cirsium_inat\ch1_v2_screen_images_pilot" `
  --image-size small
```

Check `download_results.csv`. Resolve systematic HTTP 429, 404, or non-image errors before increasing scale. Rerun the same command to retry only rows without a successful download.

## 5. Screen downloaded images with YOLO — detector only

Use the existing head detector weights but do not run the legacy orientation model. This produces the crop metadata that is required by the v2 annotation manifest.

```powershell
python ch1_global\v2\07_screen_downloaded_images_with_yolo.py `
  --queue "C:\Users\zuizui\cirsium_inat\ch1_v2_screen_queue_pilot\image_screening_queue.csv" `
  --downloads "C:\Users\zuizui\cirsium_inat\ch1_v2_screen_images_pilot\download_results.csv" `
  --weights "C:\Users\zuizui\mc\runs\cirsium_detect\weights\best.pt" `
  --out-dir "C:\Users\zuizui\cirsium_inat\ch1_v2_yolo_pilot"
```

Inspect `screening_results.csv`: detector-positive rate, counts of `no_detection`, failed files, and high-confidence false positives sampled manually. The only file passed forward is `yolo_crop_metadata.csv`. It contains source/crop paths and the metadata required for grouped validation, but no trait prediction.

## 6. Build a locked annotation manifest

After the detector pilot is acceptable, run the same collection/queue/download/detector sequence for the final queue. Then use the detector crop metadata to make the 3,000-head annotation pool, 720-head double-label calibration set, and 1,000-head segmentation subset.

```powershell
python ch1_global\v2\00_build_annotation_manifest.py `
  --input "C:\Users\zuizui\cirsium_inat\ch1_v2_yolo_final\yolo_crop_metadata.csv" `
  --out-dir "C:\Users\zuizui\cirsium_inat\ch1_v2_annotation\manifest" `
  --n-calibration 720 `
  --n-main 3000 `
  --n-segmentation 1000
```

Inspect `annotation_species_qc.csv` before annotation. The minimum per-species quota, common-species cap, and spatial-block width are part of the study design. Save each changed configuration to a new dated directory rather than overwriting a previous manifest.

**Stop here** if fewer than 36 species are represented, a small set of species dominates, or image paths visibly include taxon names. Resolve those issues before human annotation.

## 7. Build blinded packets and annotate locally

```powershell
python ch1_global\v2\01_build_blinded_annotation_packets.py `
  --manifest "C:\Users\zuizui\cirsium_inat\ch1_v2_annotation\manifest\annotation_manifest.csv" `
  --out-dir "C:\Users\zuizui\cirsium_inat\ch1_v2_annotation\packets"

python -m pip install -r ch1_global\v2\requirements_annotation_app.txt

streamlit run ch1_global\v2\03_trait_annotation_app.py -- `
  --packet "C:\Users\zuizui\cirsium_inat\ch1_v2_annotation\packets\annotator_1_packet.csv" `
  --out-dir "C:\Users\zuizui\cirsium_inat\ch1_v2_annotation\responses\annotator_1" `
  --annotator "annotator_1"
```

Run a separate app instance and output directory for the second annotator. The first annotator scores all 3,000 heads; the second scores only the 720-head calibration packet. The app saves atomically after every action and hides taxonomy, coordinates, model predictions, and raw metadata.

## 8. Compile double annotations and adjudicate

```powershell
python ch1_global\v2\02_compile_double_annotations.py `
  --primary "C:\Users\zuizui\cirsium_inat\ch1_v2_annotation\responses\annotator_1\annotation_responses.csv" `
  --secondary "C:\Users\zuizui\cirsium_inat\ch1_v2_annotation\responses\annotator_2\annotation_responses.csv" `
  --out-dir "C:\Users\zuizui\cirsium_inat\ch1_v2_annotation\agreement"
```

Read outputs in order:

1. `annotation_validation_errors.csv` — fix invalid values or missing columns.
2. `annotation_agreement_summary.csv` — retain completed-pair count, assessable-pair count, percent agreement, and kappa for each trait.
3. `annotation_adjudication_queue.csv` — adjudicate every disagreement under a predeclared rule or by a third blinded annotator.

## 9. Gate before model training and macroecology

Do not train orientation, colour, cover/display, shape, or segmentation models from raw first-pass labels. Freeze the adjudicated calibration table first. All reported performance must use `observation_group`, `species_cv_fold`, and `spatial_cv_fold`; a random image split is only a descriptive baseline.

## 10. What GitHub Actions checks

The `Ch1 annotation foundation CI` workflow is deliberately lightweight. It compiles scripts, runs unit/integration tests over synthetic inputs, and saves a tiny smoke artifact. It does not call iNaturalist, download images, access local photo files, or train YOLO. Real data collection happens locally so raw media and study data stay under researcher control.
