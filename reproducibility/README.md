# Reproducibility

## Current Chapter 1 analysis

Use [CURRENT_ANALYSIS.md](CURRENT_ANALYSIS.md). The entry point is `python -m reproducibility.run_current_analysis`. It verifies exact frozen input identities and runs existing models without changing the source cohort, construct definitions, predictors or the frozen test families. The current runner has eight numerical stages; the eighth is the explicitly post-hoc RV estimator-validity sensitivity added during submission audit.

The [2026-09-15 execution receipt](current_replay_execution.json) and [16-file numerical comparison](current_replay_validation.json) record a completed eight-stage replay with `PASS`. The current reference manifest contains 16 checksum-locked files. This closes the current replay-execution gate; it does **not** by itself publish the current v3 archive.

Numerical reproduction begins with frozen measurements, not mutable source photographs. Repeating image acquisition and measurement from scratch is a distinct upstream task: the model, acquisition and measurement provenance is retained under `legacy/ch1_global/v2/` and [actions_artifact_catalog.json](actions_artifact_catalog.json). The numerical package must not be described as a complete, permanent archive of every original photograph or as an independent physical-trait validation. The exact layer boundary is frozen in [UPSTREAM_DATA_BOUNDARY_20260926.md](UPSTREAM_DATA_BOUNDARY_20260926.md).

## Frozen v2 reproduction

The original [v2 runbook](../legacy/v2/ORIGINAL_REPRODUCTION.md) applies to immutable code revision `584af97b050d15701f26ce1facea212d5b648d4d`, paired with Zenodo DOI [10.5281/zenodo.22295791](https://doi.org/10.5281/zenodo.22295791). Do not use its old paths as current entry points. The [path map](legacy_path_map.json) locates historical implementations in the cleaned tree.

The v2 `public_release_manifest.json`, `material_availability.json`, `recovery_inventory.json` and related receipts describe that historical release. Their readiness statements do not certify current v3 public availability. Archived source paths are resolved using the path map rather than changing recorded scientific inputs.

## Current release staging

The [2026-09-12 staging receipt](CURRENT_RELEASE_STAGING_20260912.json), subsequently extended with the taxonomy, 2026-09-15 estimator-validity work and the completed eight-stage replay, records checksum-verified durable copies of the current nine-predictor environment input and current Actions evidence required by the 16-file reference surface. The private staging area now also holds the estimator-validity artifact, the CI-rendered integration figure that exposes the equal-n sensitivity, and the eight-stage/16-file replay artifact. This removes dependence on those expiring Actions artifacts for material recovery.

This staging area is private and is **not** a public release. The eight-stage replay/16-file validation receipt is complete and durably preserved. On 2026-09-26, the repository also built one end-to-end current-v3 Zenodo staging archive from main commit `4ba2304eba83e9a54f616d33cc85fd7ad818caf0`; its inner archive SHA-256 is `489d3ab8b68b646c8eba510e853aa524d36cff138263da30bc3ecd910140b91e`, and the Actions wrapper is durably preserved in the owner Drive archive. The builder reports only `release_metadata` and `figure_document_qa` as unresolved release gaps. See [ZENODO_UPDATE_AUDIT.md](ZENODO_UPDATE_AUDIT.md) for the exact receipt and remaining publication steps.

## Current release bundle builder

`python -m reproducibility.build_current_release_bundle` is the offline, fail-closed packager for the current numerical release. It does not download or publish anything. Supply a directory containing exactly one ZIP for each frozen numerical input artifact ID (`9612943217`, `9633419268`, `8983877726`, `8227254443`). For native status, either pass the exact frozen CSV explicitly or place the checksum-verified WCVP sensitivity artifact `10292140117` in the same directory; the builder then extracts `input/observation_native_status.csv` and verifies the frozen SHA automatically.

A staging bundle can be built before the final manuscript figure surface is frozen:

```bash
python -m reproducibility.build_current_release_bundle \
  --input-dir /path/to/verified-archives \
  --out /path/to/azami_ch1_current_release_staging.zip
```

The builder verifies the four numerical archive SHA-256 values and their required members, verifies the taxonomy archive before extracting native status when that route is used, normalizes the native-status transport only through the same frozen LF/CRLF rule used by the numerical runner, verifies all 16 `current_reference` files, requires a clean Git worktree, snapshots that exact `HEAD` with `git archive`, copies the current replay receipts and metadata, and writes a deterministic outer ZIP plus a `.sha256` sidecar. The numerical dependency file is pinned to the exact completed replay environment, including Biopython 1.88. Figure files may be supplied through the portable manifest generated by `reproducibility.render_current_figure_surface`; the release builder verifies manifest-relative hashes before copying them into the bundle.

For the public release, use `--final`. Final mode refuses to build unless both a frozen figure/provenance manifest and completed release-metadata JSON are supplied:

```bash
python -m reproducibility.build_current_release_bundle \
  --input-dir /path/to/verified-archives \
  --figure-manifest /path/to/final_figure_manifest.json \
  --release-metadata /path/to/zenodo_release_metadata.json \
  --expected-head <FINAL_COMMIT_SHA> \
  --final \
  --out /path/to/azami_ch1_current_release.zip
```

This hard stop is intentional: a durable staging bundle must not silently become a publication-ready claim while the final figure surface or release metadata is still unresolved.

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

The existing Zenodo record is unchanged. A new version is needed for current inputs, code and outputs: see [ZENODO_UPDATE_AUDIT.md](ZENODO_UPDATE_AUDIT.md). An audit, local reconstruction or private durable staging does not mean a new DOI version has been published. Claim credential-free current reproduction only after the new version is downloaded anonymously, its hashes are verified, and the numerical runner passes its reference checks.
