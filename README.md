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
| Ch.1 | Across global *Cirsium*, what repeated combinations do image-derived capitulum traits (orientation, colour, involucral cover, head shape, corolla:involucre display) form, how coordinated/constrained is that trait space, and how do the traits and syndromes correspond to climate, terrain, and geography? | iNat images -> multi-class YOLO ROI detection (ROI-A/ROI-B) -> per-trait extraction -> human-label + Group CV/LOSO validation -> RF/GLM/GAM + PCA/clustering syndrome analysis + trait integration/constraint. Two tiers: a global tier and a Japan deep-dive tier (justified by Japan being a *Cirsium* diversity hotspot, a scale-consistency test, and continuity with the Ch.2/4 study system — not by data density alone). Non-independence handled **phylogeny-free** (spatial + taxonomic-rank control), because *Cirsium* polyploidy/hybridization make a reliable tree intractable. The pollinator *hypothesis* is proposed as a candidate mechanism (the hook into Ch.4), but no pollinator *variable* is analysed in Ch.1 (a stacked-SDM proxy is climate-derived and cannot separate pollinator from climate effects). See `ch1_global/README.md` + `ch1_global/METHODS_AUDIT.md`. Macroecological pattern only; no causal claim. |
| Ch.2 | When and how many times was the nodding trait gained? | Molecular phylogeny + ancestral state reconstruction |
| Ch.3 | Does head orientation belong to recurrent floral trait syndromes, and how do those syndromes relate to niche shifts? | Trait dictionary + comparative data matrix + PGLS / phylogenetic logistic regression |
| Ch.4 | Do head orientation and associated floral traits affect pollination, floral damage, rain/wind protection, mating, or seed set? | Field observation + manipulation experiments + reproductive outcomes |

`ch1_japan` (GBIF occurrences x GloBI pollinator SDM x env PCA -> GLM)
contributes the **Japan deep-dive tier** of Chapter 1 via its
occurrence/environmental machinery — justified because Japan is a *Cirsium*
diversity hotspot, gives a fine-grain scale-consistency test of the global
pattern, and is the study system continued in Ch.2/4 (not by data density
alone, and demotable to a regional sensitivity section if a reviewer prefers).
Its **pollinator SDM is deferred** out of Chapter 1 (the pollinator hypothesis
is still *named* as the hook into Ch.4, but no pollinator variable is analysed;
a stacked-SDM proxy is climate-derived and cannot separate pollinator from
climate effects). No folder move is required. See `ch1_global/README.md`
§3.4, §5-6.

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
|   +-- ch1_image_traits/
|       +-- cirsium_ch1_image_trait_dictionary.csv  # ROI-A/ROI-B trait defs (Ch.1)
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
|   +-- README.md                           # Ch.1 redesign (v2): multi-trait macroecology plan
|   +-- METHODS_AUDIT.md                     # Peer-review / robustness audit of v1 + v2 gates
|   +-- v1/                                  # Frozen orientation-only baseline (preserved, not deleted)
|   |   +-- README.md                        # What v1 is, its limits, what v2 inherits vs replaces
|   |   +-- 01_annotate_and_train_yolo.py    # Annotate flowering images and train YOLO
|   |   +-- 02_crop_heads_with_yolo.py       # Crop capitulum regions with YOLO
|   |   +-- 03_train_head_direction_classifier.py
|   |   +-- 03b_train_head_direction_classifier_legacy.py
|   |   +-- 04_predict_head_direction_global.py
|   |   +-- 04b_predict_head_direction_global_legacy.py
|   |   +-- 05_rf_glm_gam_env_analysis.R     # CHELSA env vars -> RF/GLM/GAM + figures
|   |   +-- 05b_rf_only_env_analysis.R
|   # v2 (phases 0-8 in ch1_global/README.md) extends v1 into ROI-A/ROI-B
|   # detection, per-trait extraction, Group CV/LOSO validation, and
|   # per-trait + trait-syndrome environmental analysis.
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

Current state (orientation-only; see `ch1_global/README.md` for the
multi-trait redesign and the phase-by-phase roadmap that extends this):

```text
iNaturalist global images
        |
v1/01_annotate_and_train_yolo.py    --> models/best.pt, best.onnx
        |
v1/02_crop_heads_with_yolo.py       --> data/head_crop_metadata.csv
        |
v1/03_train_head_direction_classifier.py
        --> models/cirsium_head_direction_yolo_crop_model.keras
        --> data/validation_predictions.csv  (Acc=0.896, Macro F1=0.861;
                                               random split — see ch1_global/METHODS_AUDIT.md A1/A2)
        |
v1/04_predict_head_direction_global.py
        |
v1/05_rf_glm_gam_env_analysis.R     <-- CHELSA + SoilGrids
```

Planned extension (`ch1_global/README.md` §4): multi-class detector -> ROI-A
(head+involucre+pedicel) / ROI-B (corolla-only) -> per-trait extraction
(orientation, colour, involucral cover, head shape, corolla:involucre
display) -> human-label validation -> Group CV / leave-one-species-out ->
per-trait + trait-syndrome RF/GLM/GAM.

---

## Model performance (Ch.1-Global, orientation trait)

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
| H2 | Nodding is more likely to evolve and persist in rainy, seasonally cold environments. | Ch.3 (adaptive/evolutionary claim; Ch.1 only supplies the non-causal pattern that motivates it, see H1-Ch1b below) |
| H3 | Nodding alters pollinator landing posture and contact site, affecting pollination efficiency. | Ch.4 |
| H4 | Nodding is associated with higher selfing rates and seed set in pollinator-limited environments. | Ch.4 |
| H5 | Head orientation co-occurs with floral display and involucral-barrier traits in repeatable, phylogenetically structured syndromes. | Ch.3 |
| H6 | The effects of head orientation on fitness can be mediated by pollination, floral antagonism, abiotic exposure, or more than one route. | Ch.4 |

Chapter-1-specific, pattern-level (non-adaptive) hypotheses are in
`ch1_global/README.md` §1.2 (H1-Ch1a-d): they concern trait-syndrome
structure, environmental/geographic correspondence, and leave-one-species-out
generalization of image-derived traits, without claiming adaptive function.
