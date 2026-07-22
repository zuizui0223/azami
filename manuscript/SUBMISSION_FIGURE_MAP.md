# Submission figure map

The final figure set is rebuilt around the reviewer-revised claims. Legacy median-split quadrants, two-level variance bars and raw absolute-slope RMS summaries are not submission figures.

| Proposed figure | Content | Frozen or revised source |
|---|---|---|
| Figure 1 | Actual open-licensed photographs, stored YOLO boxes, reconstructed crops and deterministic orientation/colour/outline overlays | artifacts `8099953404`, `8225059018`, `8066010557`; workflow run `29725188014` |
| Figure 2 | Nested variance across assigned species, photographs and heads, with one-head/equal-replication sensitivity in the caption or supplement | `nested_visible_variance_summary.csv` from nested-variance workflow |
| Figure 3 | Complete nine-endpoint species PCA loadings and explained variance | `species_trait_pca_loadings.csv`; `species_trait_pca_variance.csv`; submission builder must assert nine endpoints including width-profile CV |
| Figure 4 | Legacy absolute-slope RMS versus median sample size and median slope SE | `legacy_precision_diagnostic_species.csv` from precision-reanalysis workflow |
| Figure 5 | Equal-module visible variation versus sampling-noise-adjusted association energy, without quadrants | `revised_species_axes.csv` from precision-reanalysis workflow |
| Figure 6 | Named-cohort within-species coefficient summary: exhaustive primary versus balanced-atlas sensitivities | exhaustive climate report plus balanced prephylogenetic coefficient table |
| Figure 7 | Among-species environmental niche contrasts and historical/data-coverage limitations | niche contrast, direct-backbone and NCBI coverage outputs |

## Nested variance provenance

- Input artifact: `8225059018`, final v2 head-level table.
- Script: `analysis/decompose_nested_visible_variance.py`.
- Workflow: `.github/workflows/ch1-nested-visible-variance.yml`.
- Atlas hierarchy: one photograph per public observation; one or more detected heads per photograph.
- Sensitivities: one deterministic head per photograph and 10 photographs per eligible species across 500 repeats.

## Revised precision provenance

- Input artifact: `8330350031`, digest `sha256:f7757b412d8390858e409696f9665bb51f8f1fd765028e290d2bd2f8f74b8a2c`.
- Script: `analysis/reanalyze_lability_precision.py`.
- Workflow: `.github/workflows/ch1-reviewer-precision-reanalysis.yml`.
- Primary revised cohort: 101 taxa, seven linear endpoints, four predictors per endpoint, n ≥ 10 for every slope.

## Excluded legacy figures

- Two-level within/among-species bars based on the earlier observation-level decomposition are superseded by the nested head–photograph–species figure.
- `figure_species_lability_quadrants.*` and derivatives are legacy provenance only.
- Module panels based on median raw absolute-slope RMS are excluded.
- Minimum-n rank stability of the legacy responsiveness index is excluded because stability does not remove its precision bias.

Figure 1 demonstrates production provenance and face validity. It does not replace independent detector or continuous-measurement validation.
