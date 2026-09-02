# azami analysis code

This repository is intentionally limited to reproducibility code for the image-phenomics analysis pipeline.

Prepublication manuscripts, Supporting Information, reviewer materials, generated figures/tables, analysis outputs, and research data are intentionally **not tracked here**.

## Contents

- `ch1_global/v2/` — image collection, screening, trait-extraction, and QC pipeline code.
- `analysis/` — statistical analysis and sensitivity-analysis scripts.
- `analysis/figures/` — figure-generation scripts. Generated figures must be written outside the tracked source tree.
- `src/azami_ch1/` — reusable Python utilities.

## Reproducibility

Inputs and generated outputs are supplied outside this repository. Use explicit input/output paths when running scripts. The repository should remain code-only before publication.

## Prepublication hygiene

Do not commit manuscripts, submission packages, result tables, figures, raw/derived data, model weights, or review correspondence. The `.gitignore` blocks common prepublication artifacts.
