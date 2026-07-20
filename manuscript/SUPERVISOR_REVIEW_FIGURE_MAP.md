# Supervisor-review figure map

These figures are presentation products generated from frozen Chapter 1 outputs. They do not refit models or change any manuscript claim. The final journal submission may consolidate panels to meet figure limits.

| Figure | File | Frozen source |
|---|---|---|
| Figure 1 | `figures/supervisor_review/figure1_real_photo_yolo_pipeline.png` and `.svg` | actual executed atlas photographs, stored YOLO boxes/confidences and final v2 continuous measurements |
| Figure 2 | `figures/supervisor_review/figure2_variance_partition.svg` | `trait_lability_conservatism_screen.csv` |
| Figure 3 | `figures/supervisor_review/figure3_pca_loadings.svg` | `species_trait_pca_loadings.csv`; `species_trait_pca_variance.csv` |
| Figure 4 | `figures/supervisor_review/figure4_lability_quadrants.svg` | `species_lability_axes.csv` |
| Figure 5 | `figures/supervisor_review/figure5_module_lability.svg` | `module_lability_summary.csv` |
| Figure 6 | `figures/supervisor_review/figure6_niche_contrasts.svg` | `trait_environmental_niche_contrasts.csv` |
| Figure 7 | `figures/supervisor_review/figure7_sample_sensitivity.svg` | `sensitivity_minimum_sample_size.csv` |
| Figure 8 | `figures/supervisor_review/figure8_historical_molecular_coverage.svg` | NCBI coverage audit and direct-backbone count |

## Frozen artifact provenance

- Actual-photo Figure 1: source artifacts `8099953404` (executed atlas and stored YOLO metadata), `8225059018` (final v2 measurements) and `8066010557` (licence and attribution metadata); build workflow run `29725188014`. The repository also stores `figure1_real_photo_provenance.csv` and `figure1_real_photo_build.json` beside the rendered figure.
- Two-axis lability: workflow run `29382767192`, artifact `8330350031`, digest `sha256:f7757b412d8390858e409696f9665bb51f8f1fd765028e290d2bd2f8f74b8a2c`.
- Grouped SPDE and module/niche outputs: workflow run `29343722222`, artifact `8315492671`, digest `sha256:f97fd9cd3adb120e5640b3a8de738f1bc3a22341850f99112a4277b679e697ce`.
- Molecular coverage and phylogenetic signal: workflow run `29391037056`, artifact `8333456619`, digest `sha256:ff453278cf0e48ed6c71ca2d02d08371eb98936608529beea8cb489306e430ea`.

`analysis/build_real_image_figure1.py` reconstructs Figure 1 from actual production records and imports the final v2 measurement functions. `analysis/build_supervisor_review_figures.py` builds Figures 2–8 from frozen result tables. Figure 1 demonstrates production provenance and measurement face validity; it does not replace the pending independent detector precision/recall audit.
