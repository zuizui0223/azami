# YOLO-positive image and head pool

`07_screen_downloaded_images_with_yolo.py` is the **operational flower/head presence gate** for Chapter 1.

For every downloaded source image it records:

- `screen_status`: `detected`, `no_detection`, `missing_image`, or `error`
- `n_detections`: number of valid detector boxes in that photograph
- `max_yolo_conf`: strongest detector confidence in that photograph

Use `09_build_yolo_positive_pool.py` immediately after screening:

```powershell
python ch1_global\v2\09_build_yolo_positive_pool.py `
  --queue "...\image_screening_queue.csv" `
  --screening "...\screening_results.csv" `
  --crops "...\yolo_crop_metadata.csv" `
  --out-dir "...\yolo_positive_pool"
```

It writes two distinct datasets:

1. `yolo_image_summary.csv` — one row per photographed observation. Its `head_present_yolo` column decides whether the photograph contains at least one visible detected head. `visible_head_count_yolo` is the number of detected heads in that photo.
2. `yolo_positive_heads.csv` — one row per head crop. This is the source table for later blind head-trait annotation. It keeps each crop's `yolo_conf` and repeats the original photograph's visible head count as `source_image_head_count`.

## Interpretation rule

`visible_head_count_yolo` is a **photo-level visible-head count**, not an estimate of every flower on a plant. It is affected by framing, occlusion, and whether the observer photographed one head or an entire individual. It can be used to choose images, flag multi-head photos, and describe image composition.

Call it a *flower count* only after checking the YOLO training labels: if training boxes included unopened buds or closed capitula, the correct term remains *visible capitulum/head count*.

For directional, colour, involucre, and shape labels, the correct unit is the individual head crop, not the source photo. Multi-head photos can therefore contribute several independent crop candidates, while their shared `obs_id` must remain available for grouped cross-validation.
