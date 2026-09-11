# Reproducibility

## Current Chapter 1 analysis

Use [CURRENT_ANALYSIS.md](CURRENT_ANALYSIS.md). The entry point is `python -m reproducibility.run_current_analysis`. It verifies exact frozen input identities and runs existing models without changing the source cohort, construct definitions, predictors, random seeds or test families.

The [2026-09-11 local execution receipt](current_replay_execution.json) records all seven numerical stages, input verification and a [15-file current-scope numerical comparison](current_replay_validation.json). This is an actual numerical replay, distinct from tests that only inspect file presence and distinct from a new public Zenodo release.

Numerical reproduction begins with frozen measurements, not mutable source photographs. Repeating image acquisition and measurement from scratch is a distinct upstream task: the model, acquisition and measurement provenance is retained under `legacy/ch1_global/v2/` and [actions_artifact_catalog.json](actions_artifact_catalog.json). The numerical package must not be described as a complete, permanent archive of every original photograph or as an independent physical-trait validation.

## Frozen v2 reproduction

The original [v2 runbook](../legacy/v2/ORIGINAL_REPRODUCTION.md) applies to immutable code revision `584af97b050d15701f26ce1facea212d5b648d4d`, paired with Zenodo DOI [10.5281/zenodo.22295791](https://doi.org/10.5281/zenodo.22295791). Do not use its old paths as current entry points. The [path map](legacy_path_map.json) locates historical implementations in the cleaned tree.

The v2 `public_release_manifest.json`, `material_availability.json`, `recovery_inventory.json` and related receipts describe that historical release. Their readiness statements do not certify current v3 public availability. Archived source paths are resolved using the path map rather than changing recorded scientific inputs.

## Figure layout revisions

`python -m reproducibility.render_layout_revisions` renders the presentation-only Figure 1 and Figure S5 revisions into `work/layout-revisions/`, using the existing matplotlib/numpy/pandas figure environment. It does not refit models or replace the frozen figures. Figure 1 separates the measurement rows from the construct summary and replaces a conflicting historical angle overlay with an image-vertical guide. The production CSV value remains 0.732009 degrees (displayed as 0.7). Figure S5 moves one long label inward without changing points or statistics. The output receipt records file hashes; it does not certify Word pagination. These are two layout revisions, not a complete current manuscript figure rebuild.

## Complete construct environment figure

Main Figure 4 is reproduced separately from the complete current numerical atlas:

```bash
python -m pip install Pillow==12.3.0
python -m reproducibility.render_construct_environment --font /path/to/arial.ttf
```

Run this after the current numerical replay, or pass `--axis-dir` containing its two `biological_axes_*.csv` inputs. On Windows the manuscript font is `C:/Windows/Fonts/arial.ttf`; supply your own licensed copy on another platform. The font is not redistributed. The renderer checks the complete unique 90-row family at each scale, displays all nine non-surface constructs without outcome-based filtering, and records input/font/output hashes. With the original font and Pillow runtime, its PNG is byte-identical to the manuscript Figure 4. Scalar cells show signed slopes; joint cells show non-negative vector magnitudes, which must not be read as signed or dimension-comparable effects. No models or multiple-testing corrections are rerun.

## Permanent archive

The existing Zenodo record is unchanged. A new version is needed for current inputs, code and outputs: see [ZENODO_UPDATE_AUDIT.md](ZENODO_UPDATE_AUDIT.md). An audit or local reconstruction does not mean a new DOI has been published. Claim credential-free current reproduction only after the new version is downloaded anonymously, its hashes are verified, and the numerical runner passes its reference checks.
