# azami — *Cirsium* floral trait architecture: evolution, ecology & function

Code repository for a PhD project investigating the adaptive significance of
nodding (downward-facing) capitula in *Cirsium*. Head orientation remains the
focal trait, now studied as part of a flower-centred architecture that can also
include floral display, involucral barriers, floral antagonism, and abiotic
exposure.

> **Central question**: When, how many times, and in which environments did
> nodding capitula evolve — and how do head orientation and associated floral
> traits jointly affect pollination, floral antagonism, rain/wind exposure, and
> reproductive outcomes?

---

## Research overview

```text
Macro pattern --> Phylogeny --> Trait architecture --> Function
    Ch. 1            Ch. 2            Ch. 3            Ch. 4
```

| Chapter | Question | Approach |
|---------|----------|----------|
| Ch.1-JP | Which environments and pollinator assemblages predict nodding capitula in Japanese *Cirsium*? | GBIF occurrences x GloBI pollinator SDM x env PCA -> GLM |
| Ch.1-GL | What is the global environmental niche of nodding capitula? | iNat images -> YOLO detection -> Keras classifier -> RF/GLM/GAM |
| Ch.2 | When and how many times was the nodding trait gained? | Molecular phylogeny + ancestral state reconstruction |
| Ch.3 | Does head orientation belong to recurrent floral trait syndromes, and how do those syndromes relate to niche shifts? | Trait dictionary + comparative data matrix + PGLS / phylogenetic logistic regression |
| Ch.4 | Do head orientation and associated floral traits affect pollination, floral damage, rain/wind protection, mating, or seed set? | Field observation + manipulation experiments + reproductive outcomes |

---

## Floral trait architecture extension

`polliTrait` was an unused R-package template. Its intended flower-trait idea is
therefore migrated here rather than maintained as a separate repository.

The initial extension remains flower-centred:

```text
A  attraction and pollinator access
D  floral barrier or resistance to floral antagonists
X  abiotic exposure and protection
O  measured interaction and reproductive outcomes
```

Head orientation is not assigned a single role in advance. It may alter
pollinator posture, floral-antagonist access, and exposure to rain or wind. A
trait association becomes a functional claim only when paired with a declared
mechanism and an interaction or fitness outcome.

See:

- `ch3_trait_architecture/README.md`
- `data/trait_architecture/cirsium_floral_trait_dictionary.csv`
- `data/trait_architecture/cirsium_field_panel_template.csv`
- `protocols/floral_trait_architecture_protocol.md`

Leaf resource traits, leaf defence, and leaf herbivory are intentionally outside
the initial module. They require an explicit cross-organ hypothesis before being
combined with floral inference.

---

## Repository structure

```text
azami/
|
+-- data/
|   +-- labels.csv                          # Manual labels (upright / nodding)
|   +-- head_crop_metadata.csv              # Metadata for YOLO-cropped heads
|   +-- training_table_yolo_crops.csv       # Keras training table
|   +-- validation_predictions.csv          # Model validation predictions
|   +-- validation_report.txt
|   +-- trait_architecture/
|       +-- cirsium_floral_trait_dictionary.csv
|       +-- cirsium_field_panel_template.csv
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
+-- ch3_trait_architecture/
|   +-- README.md                           # Scope, boundaries, and planned comparisons
|
+-- protocols/
|   +-- floral_trait_architecture_protocol.md
|
+-- reports/
|   +-- ch1_global_japan_analysis_report.Rmd
|
+-- README.md
```

---

## Analysis pipeline (Ch.1)

### Japan

```text
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

```text
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

181/192 upright and 51/67 nodding heads correctly classified. Class imbalance
was corrected with class weights.

---

## Dependencies

### R
```r
tidyverse, sf, terra, ENMeval, dismo, ggplot2, patchwork, MuMIn, lme4
```

### Python
```text
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
| H1 | Nodding capitula evolved independently in multiple lineages. | Ch.2 |
| H2 | Nodding is more likely to evolve and persist in rainy, seasonally cold environments. | Ch.1, Ch.3 |
| H3 | Nodding alters pollinator landing posture and contact site, affecting pollination efficiency. | Ch.4 |
| H4 | Nodding is associated with higher selfing rates and seed set in pollinator-limited environments. | Ch.4 |
| H5 | Head orientation co-occurs with floral display and involucral-barrier traits in repeatable, phylogenetically structured syndromes. | Ch.3 |
| H6 | The effects of head orientation on fitness can be mediated by pollination, floral antagonism, abiotic exposure, or more than one route. | Ch.4 |
