# Chapter 1 — v1 (frozen orientation-only baseline)

This directory is the **frozen v1** of the global Chapter 1 pipeline: the
single-trait (head-orientation) analysis exactly as it stood before the
multi-trait redesign in [`../README.md`](../README.md). It is preserved, not
deleted, so that:

1. the redesign (v2) has a concrete, reproducible baseline to be compared
   against and to inherit working code from, and
2. the known limitations of v1 (see [`../METHODS_AUDIT.md`](../METHODS_AUDIT.md))
   are documented against a fixed target rather than a moving one.

**v1 is not the final Chapter 1.** Its headline number (orientation
classifier Acc 0.896 / macro-F1 0.861 on a random split) is **not
peer-review-defensible as reported**, primarily because of species-level
label leakage under a random train/test split (Audit finding A1). Treat every
v1 result as provisional until it is re-validated under the v2 grouping and
error-propagation rules.

## Pipeline (as-is)

```text
iNaturalist global images  (../../ch1_shared/download_inat_images_japan.py)
        |
v1/01_annotate_and_train_yolo.py     --> models/best.pt, best.onnx   (capitulum detector, single class)
        |
v1/02_crop_heads_with_yolo.py        --> data/head_crop_metadata.csv  (one crop per detection)
        |
v1/03_train_head_direction_classifier.py
        --> models/cirsium_head_direction_yolo_crop_model.keras
        --> data/validation_predictions.csv   (random split; see Audit A1/A2)
        |
v1/04_predict_head_direction_global.py
        |
v1/05_rf_glm_gam_env_analysis.R      <-- CHELSA + topography + SoilGrids (pre-extracted)
v1/05b_rf_only_env_analysis.R        <-- RF-only resume variant
```

`03b_*_legacy.py` and `04b_*_legacy.py` are earlier CSV-driven variants of the
classifier train/predict steps, kept for provenance only. The labels they
consume trace back to `../../ch1_shared/make_training_data_from_labels.py`,
which joins **species-level** orientation descriptions from Kahaku specimen
text to iNaturalist images **by species name** — the root cause of Audit
finding A1.

## Known non-portability

Every v1 script hardcodes Windows paths (`C:\Users\zuizui\cirsium_inat\...`)
and assumes local intermediate files (the CHELSA/SoilGrids enriched CSV is
produced by an extraction step that is *not* in this repository). v1 is
therefore a record of what was run, not a turnkey pipeline. Making it
path-portable and adding the missing environmental-extraction step are v2
tasks (Audit findings B3, C2).

## What v2 inherits vs. replaces

| v1 component | v2 disposition |
|---|---|
| Single-class capitulum detector (`01`) | Extended to multi-class ROI-A/ROI-B detector (`../README.md` §2-3). |
| Head crop (`02`) | Replaced by standardized ROI-A + segmentation-derived ROI-B. |
| Orientation classifier (`03`/`04`) | Reused as one trait module, but re-validated under Group-CV/LOSO and paired with a continuous angle estimate. |
| Species-name label join (`ch1_shared`) | Replaced by image-level human labels with inter-annotator agreement (Audit A1). |
| Random train/test split | Replaced by grouped CV (photo / species / geographic block) + LOSO (Audit A2). |
| Binary `y_nodding` RF/GLM/GAM (`05`/`05b`) | Generalized to per-trait models + trait-syndrome PCA/clustering, with aggregation-level sensitivity and error-in-variables handling (Audit B1-B4, D1). |
