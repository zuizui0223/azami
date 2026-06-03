# azami — *Cirsium* head orientation: evolution, ecology & function

Code repository for a PhD project investigating the **adaptive significance of nodding (downward-facing) capitula** in *Cirsium*, integrating macroecology, phylogenetics, and functional experiments.

> **Central question**: When, how many times, and in which environments did nodding capitula evolve — and what do they actually do?

---

## Research overview

```
Macro pattern --> Phylogeny --> Niche evolution --> Function
    Ch. 1            Ch. 2           Ch. 3            Ch. 4
```

| Chapter | Question | Approach |
|---------|----------|----------|
| Ch.1-JP | Which environments and pollinator assemblages predict nodding capitula in Japanese *Cirsium*? | GBIF occurrences x GloBI pollinator SDM x env PCA -> GLM |
| Ch.1-GL | What is the global environmental niche of nodding capitula? | iNat images -> YOLO detection -> Keras classifier -> RF/GLM/GAM |
| Ch.2 | When and how many times was the nodding trait gained? | Molecular phylogeny + ancestral state reconstruction |
| Ch.3 | Does nodding co-occur with humid/cold niche shifts? | PGLS / phylogenetic logistic regression |
| Ch.4 | Does head orientation affect pollination, rain/wind protection, or seed set? | Field observation + manipulation experiments |

---

## Repository structure

```
azami/
|
+-- data/
|   +-- labels.csv                          # Manual labels (upright / nodding)
|   +-- head_crop_metadata.csv              # Metadata for YOLO-cropped heads
|   +-- training_table_yolo_crops.csv       # Keras training table
|   +-- validation_predictions.csv          # Model validation predictions
|   +-- validation_report.txt
|
+-- models/
|   +-- best.pt                             # YOLO best weights (PyTorch)
|   +-- best.onnx                           # YOLO best weights (ONNX)
|   +-- last.pt                             # YOLO last epoch weights
|   +-- cirsium_head_direction_yolo_crop_model.keras
|
+-- ch1_shared/
|   +-- scrape_kahaku_specimen_images.py    # Scrape specimen images from NSMT
|   +-- download_inat_images_japan.py       # Bulk download iNat images (Japan)
|   +-- make_training_data_from_labels.py   # Convert NSMT labels to iNat training set
|
+-- ch1_japan/
|   +-- 01_get_gbif_occurrences.R           # Fetch GBIF records and aggregate to grid
|   +-- 02_build_pollinator_sdm.R           # GloBI fetch + per-order SSDM via ENMeval
|   +-- 02b_build_pollinator_sdm_enmeval.R  # ENMeval pipeline (extended version)
|   +-- 03_env_pollinator_pca_glm.R         # Env PCA x pollinator PCA -> GLM
|   +-- 04_sensitivity_model_comparison.R   # Pollinator grouping sensitivity analysis
|   +-- 05_export_figures.R                 # Export publication figures
|
+-- ch1_global/
|   +-- 01_annotate_and_train_yolo.py       # Annotate flowering images and train YOLO
|   +-- 02_crop_heads_with_yolo.py          # Crop capitulum regions with YOLO
|   +-- 03_train_head_direction_classifier.py
|   +-- 03b_train_head_direction_classifier_legacy.py
|   +-- 04_predict_head_direction_global.py # Apply classifier to global images
|   +-- 04b_predict_head_direction_global_legacy.py
|   +-- 05_rf_glm_gam_env_analysis.R        # CHELSA env vars -> RF/GLM/GAM + figures
|   +-- 05b_rf_only_env_analysis.R
|
+-- reports/
|   +-- ch1_global_japan_analysis_report.Rmd
|
+-- README.md
```

---

## Analysis pipeline (Ch.1)

### Japan

```
scrape_kahaku_specimen_images.py -+
download_inat_images_japan.py    -+-> make_training_data_from_labels.py
                                        |
                               labels.csv / training_table_yolo_crops.csv
                                        |
01_get_gbif_occurrences.R           <-- GBIF
02_build_pollinator_sdm.R           <-- GloBI + ENMeval (Hymenoptera / Lepidoptera / Diptera)
03_env_pollinator_pca_glm.R         <-- CHELSA climate variables
04_sensitivity_model_comparison.R
05_export_figures.R
```

### Global (AI pipeline)

```
iNaturalist global images
        |
01_annotate_and_train_yolo.py       --> models/best.pt, best.onnx
        |
02_crop_heads_with_yolo.py          --> data/head_crop_metadata.csv
        |
03_train_head_direction_classifier.py
        --> models/cirsium_head_direction_yolo_crop_model.keras
        --> data/validation_predictions.csv  (Acc=0.896, Macro F1=0.861)
        |
04_predict_head_direction_global.py
        |
05_rf_glm_gam_env_analysis.R        <-- CHELSA + SoilGrids
```

---

## Model performance (Ch.1-Global)

| Metric | Value |
|--------|-------|
| Accuracy | 0.896 |
| Macro F1 | 0.861 |
| Down Precision | 0.823 |
| Down Recall | 0.761 |

181/192 upright and 51/67 nodding heads correctly classified. Class imbalance corrected with class weights.

---

## Dependencies

### R
```r
tidyverse, sf, terra, ENMeval, dismo, ggplot2, patchwork, MuMIn, lme4
```

### Python
```
ultralytics    # YOLO
tensorflow / keras
pandas, numpy, Pillow
selenium       # NSMT scraping
pyinaturalist  # iNat API
```

---

## Hypotheses

| | Hypothesis | Chapter |
|---|---|---|
| H1 | Nodding capitula evolved independently in multiple lineages | Ch.2 |
| H2 | Nodding is more likely to evolve and persist in rainy, seasonally cold environments | Ch.1, Ch.3 |
| H3 | Nodding alters pollinator landing posture and contact site, affecting pollination efficiency | Ch.4 |
| H4 | Nodding is associated with higher selfing rates and seed set in pollinator-limited environments | Ch.4 |

---
