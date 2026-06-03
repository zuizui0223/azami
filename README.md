# azami 窶・*Cirsium* head orientation: evolution, ecology & function

繧｢繧ｶ繝溷ｱ槭↓縺翫￠繧・*荳句髄縺埼ｭ闃ｱ縺ｮ驕ｩ蠢懃噪諢冗ｾｩ**繧偵√・繧ｯ繝ｭ逕滓・繝ｻ邉ｻ邨ｱ繝ｻ讖溯・螳滄ｨ薙ｒ邨ｱ蜷医＠縺ｦ讀懆ｨｼ縺吶ｋ蜊壼｣ｫ遐皮ｩｶ縺ｮ繧ｳ繝ｼ繝峨Μ繝昴ず繝医Μ縺ｧ縺吶・
> **荳ｭ蠢・撫縺・*・壻ｸ句髄縺埼ｭ闃ｱ縺ｯ縲√＞縺､繝ｻ菴募屓繝ｻ縺ｩ縺ｮ迺ｰ蠅・〒騾ｲ蛹悶＠縲∝ｮ滄圀縺ｫ菴輔↓蜉ｹ縺上・縺具ｼ・
---

## Research overview

```
Macro pattern 笏笏竊・Phylogeny 笏笏竊・Niche evolution 笏笏竊・Function
   Ch. 1               Ch. 2          Ch. 3              Ch. 4
```

| Chapter | 蝠上＞ | 謇区ｳ・|
|---------|------|------|
| Ch.1-JP | 譌･譛ｬ逕｣繧｢繧ｶ繝溘〒荳句髄縺埼ｭ闃ｱ縺悟､壹＞迺ｰ蠅・・騾∫ｲ芽・嶌縺ｯ・・| GBIF 蛻・ｸ・ﾃ・GloBI 騾∫ｲ芽・SDM ﾃ・迺ｰ蠅・PCA 竊・GLM |
| Ch.1-GL | 荳也阜繧ｹ繧ｱ繝ｼ繝ｫ縺ｮ迺ｰ蠅・ル繝・メ縺ｯ・・| iNat 逕ｻ蜒・竊・YOLO 讀懷・ 竊・Keras 蛻・｡槫勣 竊・RF/GLM/GAM |
| Ch.2 | 荳句髄縺阪・縺・▽繝ｻ菴募屓迯ｲ蠕励＆繧後◆・・| 蛻・ｭ千ｳｻ邨ｱ讓ｹ + 逾門・迥ｶ諷区耳螳・|
| Ch.3 | 貉ｿ貎､繝ｻ蜀ｷ豸ｼ繝九ャ繝√∈縺ｮ騾ｲ蜃ｺ縺ｨ騾｣蜍輔☆繧九°・・| PGLS / phylogenetic logistic |
| Ch.4 | 闃ｱ蜷代″縺ｯ騾∫ｲ峨・髮ｨ鬚ｨ繝ｻ邨仙ｮ溘↓蜉ｹ縺上°・・| 驥主､冶ｦｳ蟇・+ 謫堺ｽ懷ｮ滄ｨ・|

---

## Repository structure

```
azami/
笏・笏懌楳笏 data/
笏・  笏懌楳笏 labels.csv                          # 謇句虚繝ｩ繝吶Ν・井ｸ雁髄縺・荳句髄縺搾ｼ・笏・  笏懌楳笏 head_crop_metadata.csv              # YOLO 繧ｯ繝ｭ繝・・縺ｮ繝｡繧ｿ諠・ｱ
笏・  笏懌楳笏 training_table_yolo_crops.csv       # Keras 險鍋ｷｴ繝・・繝悶Ν
笏・  笏懌楳笏 validation_predictions.csv          # 繝｢繝・Ν讀懆ｨｼ縺ｮ莠域ｸｬ邨先棡
笏・  笏披楳笏 validation_report.txt
笏・笏懌楳笏 models/
笏・  笏懌楳笏 best.pt                             # YOLO 譛濶ｯ驥阪∩ (PyTorch)
笏・  笏懌楳笏 best.onnx                           # YOLO 譛濶ｯ驥阪∩ (ONNX)
笏・  笏懌楳笏 last.pt                             # YOLO 譛邨ゅお繝昴ャ繧ｯ驥阪∩
笏・  笏披楳笏 cirsium_head_direction_yolo_crop_model.keras
笏・笏懌楳笏 ch1_shared/
笏・  笏懌楳笏 scrape_kahaku_specimen_images.py    # 遘大忽繧ｦ繧ｧ繝悶°繧画ｨ呎悽逕ｻ蜒上せ繧ｯ繝ｬ繧､繝・笏・  笏懌楳笏 download_inat_images_japan.py       # iNat 譌･譛ｬ逕｣逕ｻ蜒上・荳諡ｬ繝繧ｦ繝ｳ繝ｭ繝ｼ繝・笏・  笏披楳笏 make_training_data_from_labels.py   # 遘大忽繝ｩ繝吶Ν 竊・iNat 險鍋ｷｴ繧ｻ繝・ヨ螟画鋤
笏・笏懌楳笏 ch1_japan/
笏・  笏懌楳笏 01_get_gbif_occurrences.R           # GBIF 蛻・ｸ・Ξ繧ｳ繝ｼ繝牙叙蠕励・繧ｰ繝ｪ繝・ラ髮・ｨ・笏・  笏懌楳笏 02_build_pollinator_sdm.R           # GloBI 蜿門ｾ・(Genus/Guild)
笏・  笏懌楳笏 02b_build_pollinator_sdm_enmeval.R  # ENMeval 縺ｧ逶ｮ蛻･ SSDM 讒狗ｯ・笏・  笏懌楳笏 03_env_pollinator_pca_glm.R         # 迺ｰ蠅・PCA ﾃ・騾∫ｲ芽・PCA 竊・GLM
笏・  笏懌楳笏 04_sensitivity_model_comparison.R   # 騾∫ｲ芽・げ繝ｫ繝ｼ繝斐Φ繧ｰ諢溷ｺｦ蛻・梵繝ｻ繝｢繝・Ν豈碑ｼ・笏・  笏披楳笏 05_export_figures.R                 # 謚慕ｨｿ逕ｨ蝗ｳ迚亥・蜉・笏・笏懌楳笏 ch1_global/
笏・  笏懌楳笏 01_annotate_and_train_yolo.py       # 髢玖干逕ｻ蜒上い繝弱ユ繝ｼ繧ｷ繝ｧ繝ｳ 竊・YOLO 蟄ｦ鄙・笏・  笏懌楳笏 02_crop_heads_with_yolo.py          # YOLO 縺ｧ鬆ｭ闃ｱ鬆伜沺繧偵け繝ｭ繝・・
笏・  笏懌楳笏 03_train_head_direction_classifier.py
笏・  笏懌楳笏 03b_train_head_direction_classifier_legacy.py
笏・  笏懌楳笏 04_predict_head_direction_global.py # 蜈ｨ逅・判蜒上∈縺ｮ莠域ｸｬ驕ｩ逕ｨ
笏・  笏懌楳笏 04b_predict_head_direction_global_legacy.py
笏・  笏懌楳笏 05_rf_glm_gam_env_analysis.R        # CHELSA 迺ｰ蠅・､画焚 竊・RF/GLM/GAM 隗｣譫舌・蝗ｳ迚・笏・  笏披楳笏 05b_rf_only_env_analysis.R
笏・笏懌楳笏 reports/
笏・  笏披楳笏 ch1_global_japan_analysis_report.Rmd
笏・笏披楳笏 README.md
```

---

## Analysis pipeline (Ch.1)

### Japan

```
scrape_kahaku_specimen_images.py 笏笏・download_inat_images_japan.py    笏笏､竊・make_training_data_from_labels.py
                                   竊・                          labels.csv / training_table_yolo_crops.csv
                                   竊・01_get_gbif_occurrences.R          竊絶楳笏 GBIF
02_build_pollinator_sdm.R          竊絶楳笏 GloBI + ENMeval (Hymenoptera / Lepidoptera / Diptera)
03_env_pollinator_pca_glm.R        竊絶楳笏 CHELSA 迺ｰ蠅・､画焚
04_sensitivity_model_comparison.R
05_export_figures.R
```

### Global (AI pipeline)

```
iNaturalist 蜈ｨ逅・判蜒・      竊・01_annotate_and_train_yolo.py        竊・models/best.pt, best.onnx
      竊・02_crop_heads_with_yolo.py           竊・data/head_crop_metadata.csv
      竊・03_train_head_direction_classifier.py
      竊・models/cirsium_head_direction_yolo_crop_model.keras
      竊・data/validation_predictions.csv  (Acc=0.896, Macro F1=0.861)
      竊・04_predict_head_direction_global.py
      竊・05_rf_glm_gam_env_analysis.R         竊絶楳笏 CHELSA + SoilGrids
```

---

## Model performance (Ch.1-Global)

| Metric | Value |
|--------|-------|
| Accuracy | 0.896 |
| Macro F1 | 0.861 |
| Down Precision | 0.823 |
| Down Recall | 0.761 |

荳雁髄縺・181/192縲∽ｸ句髄縺・51/67 繧呈ｭ｣蛻・｡槭ゆｸ句髄縺阪け繝ｩ繧ｹ縺ｮ荳崎ｶｳ縺ｯ class weight 縺ｧ陬懈ｭ｣縲・
---

## Dependencies

### R
```r
tidyverse, sf, terra, ENMeval, dismo, ggplot2, patchwork, MuMIn, lme4
```

### Python
```
ultralytics   # YOLO
tensorflow / keras
pandas, numpy, Pillow
selenium      # 遘大忽繧ｹ繧ｯ繝ｬ繧､繝・pyinaturalist # iNat API
```

---

## Hypotheses

| | 莉ｮ隱ｬ | 讀懆ｨｼ遶 |
|---|---|---|
| H1 | 荳句髄縺埼ｭ闃ｱ縺ｯ隍・焚邉ｻ邨ｱ縺ｧ迢ｬ遶九↓蜿榊ｾｩ騾ｲ蛹悶＠縺・| Ch.2 |
| H2 | 螟夐岑繝ｻ蟄｣遽諤ｧ縺ｮ蠑ｷ縺・・豸ｼ迺ｰ蠅・〒迯ｲ蠕励・邯ｭ謖√＆繧後ｄ縺吶＞ | Ch.1, Ch.3 |
| H3 | 繝上メ逶ｮ縺ｮ險ｪ闃ｱ蟋ｿ蜍｢繝ｻ謗･隗ｦ驛ｨ菴阪ｒ螟牙喧縺輔○騾∫ｲ牙柑邇・↓髢｢繧上ｋ | Ch.4 |
| H4 | 騾∫ｲ芽・ｸ榊ｮ牙ｮ夂腸蠅・〒閾ｪ谿也紫繝ｻ邨仙ｮ溽紫縺ｮ螟牙喧縺ｨ騾｣蜍輔☆繧・| Ch.4 |

---

## Author

ZHANG Rachel
