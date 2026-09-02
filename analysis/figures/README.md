# Canonical Chapter 1 figure rendering

This directory contains **code only**. Frozen analysis outputs, source-provenance
assets, manuscript files and rendered figures must remain outside the
prepublication repository.

## Single source of truth

`azami_figstyle.py` owns:

- final review-canvas widths;
- absolute font sizes;
- Okabe-Ito semantic colours;
- hatch semantics for `not comparable` and direction instability;
- the canonical endpoint labels and module order;
- the 26 analysis-unit, 22 measured-endpoint, 18-endpoint PCA and 17-unit S7 contracts;
- fixed-canvas 600-dpi PNG and vector-PDF export.

Do not pass a raw submission `figsize`, independent colour palette, or
`bbox_inches="tight"` from an individual figure script.

## External inputs

`build_v2_figures.py` requires three explicit locations:

```text
--input-dir   frozen v2 analysis-output directory
--source-dir  Figure 1 provenance composite/CSV and Figure 2 map-display assets
--output-dir  external destination for rendered figures
```

The figure-source directory is expected to contain:

```text
Figure_1_real_photo_measurement_provenance.png
Figure_1_real_photo_measurement_provenance.csv
Figure_2_v2_realized_sampling_cells.csv
natural_earth_110m_land_outline.json
```

The main builder writes Main Figures 1–5 and Supporting Figures S3–S7. The
technical-audit builder writes Supporting Figures S1–S2 from its external frozen
summary JSON.

## Render

```bash
python analysis/figures/build_v2_figures.py \
  --input-dir /path/to/frozen-v2-results \
  --source-dir /path/to/figure-source-assets \
  --output-dir /path/to/rendered-figures

python analysis/figures/run_image_to_trait_perturbation_audit.py \
  --summary /path/to/image_to_trait_automated_technical_audit_summary.json \
  --output-dir /path/to/rendered-figures
```

The frozen technical-audit summary records nine perturbations but does not retain
scenario-level numerical values for all nine. Figure S2 therefore visualizes only
the detailed numerical metrics actually retained in that summary. Missing values
must never be reconstructed from prose.

## Validate

After all 12 figures have been rendered:

```bash
python analysis/figures/validate_v2_figure_style.py \
  --figure-dir /path/to/rendered-figures
```

The validator checks:

- all 12 PNG/PDF pairs are present;
- every PNG has the same 6.27-in / 600-dpi physical width;
- no figure exceeds the 7.60-in image-height ceiling;
- the figure builders do not reintroduce tight bounding-box or 300-dpi exports;
- endpoint registry counts remain 26 / 22 / 18 / 17;
- required hatch encodings remain active.

## Figure-specific presentation rules

- **Figure 1**: eight image/glyph panels; no prose-only panel and no pixel table.
  Source-photo subtitles use continuous orientation/chroma values rather than
  categorical colour/orientation labels.
- **Figure 2**: the 2° map cell is explicitly labelled `display only`; cohort flow
  is horizontal to reduce height.
- **Figure 3**: both panels use one canonical endpoint order.
- **Figure 4**: `not comparable` is white + hatch + border, never invisible white.
- **Figure 5**: candidate coefficients and gate evidence are graphical; the old
  large prose block is not reproduced.
- **Figure S4**: within/among panels are side-by-side; direction instability uses
  orange + hatch, not red/green.
- **Figure S5**: both P = 0.05 thresholds are labelled directly and retained rows
  use leader lines.
- **Figure S6**: randomized placement results are distributions rather than 50
  overplotted points; deterministic placements remain explicit.
- **Figure S7**: PCA loadings are shown as a biplot and the matrices use the same
  17-unit registry/order as the rest of the figure system.
