# Submission figure map

The final figure set is being rebuilt around the reviewer-revised claims. Legacy median-split quadrants and raw absolute-slope RMS summaries are not submission figures.

| Proposed figure | Content | Frozen or revised source |
|---|---|---|
| Figure 1 | Actual open-licensed photographs, stored YOLO boxes, reconstructed crops and deterministic orientation/colour/outline overlays | artifacts `8099953404`, `8225059018`, `8066010557`; workflow run `29725188014` |
| Figure 2 | Visible within/among-assigned-species variance plus complete nine-endpoint PCA loadings | `trait_lability_conservatism_screen.csv`; `species_trait_pca_loadings.csv`; `species_trait_pca_variance.csv` |
| Figure 3 | Legacy absolute-slope RMS versus median sample size and median slope SE | `legacy_precision_diagnostic_species.csv` from precision-reanalysis workflow |
| Figure 4 | Equal-module visible variation versus sampling-noise-adjusted association energy, without quadrants | `revised_species_axes.csv` from precision-reanalysis workflow |
| Figure 5 | Named-cohort within-species coefficient summary: exhaustive primary versus balanced-atlas sensitivities | exhaustive climate report plus balanced prephylogenetic coefficient table |
| Figure 6 | Among-species environmental niche contrasts and historical/data-coverage limitations | niche contrast, direct-backbone and NCBI coverage outputs |

## Revised precision provenance

- Input artifact: `8330350031`, digest `sha256:f7757b412d8390858e409696f9665bb51f8f1fd765028e290d2bd2f8f74b8a2c`.
- Script: `analysis/reanalyze_lability_precision.py`.
- Workflow: `.github/workflows/ch1-reviewer-precision-reanalysis.yml`.
- Primary revised cohort: 101 taxa, seven linear endpoints, four predictors per endpoint, n ≥ 10 for every slope.

## Excluded legacy figures

- `figure_species_lability_quadrants.*` and derivatives are legacy provenance only.
- Module panels based on median raw absolute-slope RMS are excluded.
- Minimum-n rank stability of the legacy responsiveness index is excluded because stability does not remove its precision bias.

Figure 1 demonstrates production provenance and face validity. It does not replace independent detector or continuous-measurement validation.
