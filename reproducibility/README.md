# Reproducibility

## Current Chapter 1 analysis

Use [CURRENT_ANALYSIS.md](CURRENT_ANALYSIS.md). The entry point is `python -m reproducibility.run_current_analysis`. It verifies exact frozen input identities and runs existing models without changing the source cohort, construct definitions, predictors or the frozen test families. The current runner has eight numerical stages; the eighth is the explicitly post-hoc RV estimator-validity sensitivity added during submission audit.

The [2026-09-15 execution receipt](current_replay_execution.json) and [16-file numerical comparison](current_replay_validation.json) record a completed eight-stage replay with `PASS`. The current reference manifest contains 16 checksum-locked files. This closes the current replay-execution gate; it does **not** by itself publish the current v3 archive.

Numerical reproduction begins with frozen measurements, not mutable source photographs. Repeating image acquisition and measurement from scratch is a distinct upstream task: the model, acquisition and measurement provenance is retained under `legacy/ch1_global/v2/` and [actions_artifact_catalog.json](actions_artifact_catalog.json). The numerical package must not be described as a complete, permanent archive of every original photograph or as an independent physical-trait validation. The exact layer boundary is frozen in [UPSTREAM_DATA_BOUNDARY_20260926.md](UPSTREAM_DATA_BOUNDARY_20260926.md).

## Frozen v2 reproduction

The original [v2 runbook](../legacy/v2/ORIGINAL_REPRODUCTION.md) applies to immutable code revision `584af97b050d15701f26ce1facea212d5b648d4d`, paired with Zenodo DOI [10.5281/zenodo.22295791](https://doi.org/10.5281/zenodo.22295791). Do not use its old paths as current entry points. The [path map](legacy_path_map.json) locates historical implementations in the cleaned tree.

The v2 `public_release_manifest.json`, `material_availability.json`, `recovery_inventory.json` and related receipts describe that historical release. Their readiness statements do not certify current v3 public availability. Archived source paths are resolved using the path map rather than changing recorded scientific inputs.

## Current Zenodo data staging

Zenodo is used only as **durable storage for processed numerical analysis inputs** that would otherwise depend on expiring GitHub Actions artifacts. GitHub remains the authority for analysis code, runbooks, current reference outputs, figures and replay receipts.

The current data-only archive contains exactly five analysis inputs:

1. continuous image-derived measurements — exact artifact `9612943217`;
2. nine-predictor environment input — exact artifact `9633419268`;
3. broad-region lookup — exact artifact `8983877726`;
4. 52 historical-placement trees — exact artifact `8227254443`;
5. frozen `observation_native_status.csv`, checksum-recovered from taxonomy artifact `10292140117`.

The Zenodo archive intentionally does **not** contain analysis code, manuscript/SI files, figures, fitted/reference outputs, replay receipts, original photographs, or detector-training material. The archive manifest stores the GitHub repository URL and exact code commit needed to interpret the inputs.

## Current data-only bundle builder

`python -m reproducibility.build_current_release_bundle` is the offline, fail-closed packager for this numerical-input archive. It verifies the exact four Actions artifact ZIP hashes and their required members, recovers the exact native-status CSV, writes a compact input contract/README/checksum manifest, and creates a deterministic ZIP.

Staging build:

```bash
python -m reproducibility.build_current_release_bundle \
  --input-dir /path/to/verified-input-artifacts \
  --expected-head <GITHUB_CODE_COMMIT> \
  --out /path/to/azami_ch1_v3_analysis_inputs.zip
```

Final build requires only the approved release-metadata contract in addition to the verified inputs:

```bash
python -m reproducibility.build_current_release_bundle \
  --input-dir /path/to/verified-input-artifacts \
  --release-metadata /path/to/zenodo_release_metadata.json \
  --expected-head <FINAL_GITHUB_COMMIT_SHA> \
  --final \
  --out /path/to/azami_ch1_v3_analysis_inputs.zip
```

For numerical replay, download the Zenodo data package and use the separately versioned GitHub code at the pinned commit with [CURRENT_ANALYSIS.md](CURRENT_ANALYSIS.md). This separation is deliberate: **Zenodo preserves bytes; GitHub preserves executable analysis history.**

## Current scale-dependent integration figure

`python -m reproducibility.render_scale_integration` renders the current construct-level cross-scale result directly from checksum-verified current-reference files. It performs no fitting or resampling itself; the displayed equal-n sensitivity is read from the frozen estimator-validity summary produced by `analysis.v3.run_rv_estimator_validity`.

```bash
python -m reproducibility.render_scale_integration \
  --out-dir work/scale-integration-figure
```

Panels (a) and (b) show the raw within- and among-taxon 9 × 9 integration matrices for the exact 1,734-observation / 42-taxon common cohort. Panel (c) shows all 36 raw relation-wise contrasts and explicitly separates the frozen raw count (`33/36`) from the equal-n sensitivity median (`23/36`). Panel (d) contrasts the frozen raw taxon-bootstrap strength difference with the estimator-validity calculation in which one centred observation per taxon is drawn so both scales use 42 rows.

The estimator-validity audit is post-hoc and was motivated after the raw result was known. Its equal-n result retains the qualitative scale conclusion: among-minus-within median RV has median `+0.019316`, 95% interval `+0.001884` to `+0.030352`, and is positive in `98.3%` of 1,000 replicates. The median number of relations stronger among taxa is `23/36`, and a majority of relations is stronger among taxa in `97.3%` of replicates. Pairwise permutation-null centring reduces the relation count to `19/36` but retains a larger median null-centred RV among taxa; coordinate standardization retains `34/36`, matrix alignment `rho = 0.40849` with QAP `P = 0.0067`, and module cohesion at both scales.

Therefore the current supported claim is that visible-phenotype integration is **stronger overall among taxa**, not that the raw `33/36` count is estimator-invariant. The figure and its provenance preserve both the frozen raw analysis and the estimator-validity boundary. None of these results establishes that evolution increases integration or demonstrates genetic, developmental, functional or causal modularity.

The `Render current scale integration figure` workflow runs the focused renderer test, builds PNG/PDF outputs and uploads the figure plus provenance as a GitHub Actions artifact for visual QA. The generic stem `Figure_v3_scale_integration` is intentional until final manuscript figure numbering is frozen.

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

The existing Zenodo record is unchanged. A new version is needed for the current **numerical analysis inputs only**: see [ZENODO_UPDATE_AUDIT.md](ZENODO_UPDATE_AUDIT.md). Code and outputs remain on GitHub. Claim credential-free current reproduction only after the new data version is downloaded without owner credentials, its hashes are verified, and the GitHub numerical runner passes against those archived inputs.
