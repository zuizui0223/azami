# Chapter 1 current figure surface — 2026-09-12

## Purpose

Freeze the **roles** of the current GEB figure surface after the Chapter 1 v3 scope freeze. This document does not change any scientific result and is not manuscript text. It prevents a stale v2 figure order from pulling the paper back toward an endpoint-atlas story after the construct-level scale result became the central comparative inference.

The final figure numbers remain subject to document QA, but the scientific roles below are now the preferred current order.

## Proposed Main figure surface

### Main 1 — measurement → biological construct workflow

**Role:** establish what an image-defined phenotype is and how the 22 measured endpoints enter the current construct-level synthesis.

Current recipe:

- `reproducibility/render_layout_revisions.py`
- output stem inherited from the historical Figure 1 renderer.

Must communicate:

- detector localizes the visible capitulum but does not validate biological accuracy;
- deterministic continuous measurements retain endpoint-specific definitions and missingness;
- current v3 uses outcome-blind biological constructs above the frozen endpoint layer;
- image vertical is not gravity and visible colour is not calibrated physical reflectance.

Do not turn this figure into a methods catalogue. Its job is to let the reader trust the transition from public image to defined measurement to biological construct.

### Main 2 — realized global sampling domain

**Role:** show where the broad image-phenomics atlas exists and how uneven the realized sample is.

Current source:

- frozen v2 geographic-sampling-domain figure and its provenance.

Scientific cohort:

- **46,276 spatially thinned observations / 259 source-assigned taxa**.

This is the broad atlas/sampling-domain cohort. It is **not** the direct nine-construct within-vs-among integration cohort.

### Main 3 — scale-dependent construct integration

**Role:** carry the principal current v3 comparative conclusion.

Canonical recipe:

- `reproducibility/render_scale_integration.py`
- generic output stem `Figure_v3_scale_integration` until final numbering is frozen.

Scientific cohort:

- **1,734 observations / 42 taxa**;
- same nine biological constructs and all 36 construct relations at both scales.

Required visual claims:

- within-vs-among matrix alignment `rho = 0.439125`, QAP `P = 0.0041`;
- module organization remains detectable at both scales;
- observed median RV = `0.002238` within vs `0.043212` among;
- **33/36** construct relations stronger among taxa;
- among-minus-within median-RV bootstrap positive in **100%** of 1,000 replicates.

Interpretation boundary:

> The capitulum is partially organized at both scales, but among-taxon visible-phenotype integration is stronger overall while detailed relation geometry is only partly conserved.

Do not convert this pattern into genetic/developmental modularity, causal evolutionary change, or a universal integration law.

### Main 4 — complete construct × environment atlas

**Role:** show the complete environmental context without returning to a cherry-picked endpoint list.

Canonical recipe:

- `reproducibility/render_construct_environment.py`.

The figure displays the complete current nine-construct × nine-predictor family at each biological scale, with scalar signed slopes and joint-construct non-negative magnitudes kept distinct.

Its purpose is contextual: the paper's central result is not that every construct has a climate effect. The environment atlas shows where the two retained ecological anchors sit inside the broader construct-level phenotype.

### Main 5 — robustness of the two retained ecological anchors

**Role:** finish the Results by showing why only two among-taxon associations receive headline ecological interpretation.

Retained anchors:

- lower floral chroma under higher shortwave radiation;
- larger image-referenced presentation angle under higher annual precipitation.

The robustness surface should retain the frozen multiplicity result and the construct-level sampling, broad/residual spatial and 52-placement historical sensitivity sequence. Any caption must state that 52 placement scenarios are sensitivity analyses on the same observational data, not 52 independent confirmations.

## Demotion from the old Main surface

### Old Main Figure 3 — taxon-mean information loss

Move the standalone information-loss figure to Supporting Information.

The result remains useful and should remain in the paper:

- repeated public images retain visible variation that a single taxon mean discards;
- the partition combines biological, photographic and measurement structure and is not genetic variance.

However, after v3 this is supporting methodological evidence, not the paper's highest-level biological conclusion. Keeping it as a full Main figure would over-weight the 'species means lose information' story relative to the stronger scale-dependent integration result.

If space is needed in the Main text, summarize the information-loss result numerically and point to the Supporting figure.

## Evidence-layer separation

The five Main roles deliberately separate four different denominators/questions:

1. **measurement definition** — what the image variables mean;
2. **46,276 / 259 global atlas** — realized geographic and endpoint/environment domain;
3. **1,734 / 42 complete common cohort** — direct within-vs-among integration geometry and strength;
4. **two frozen ecological anchors** — robustness-qualified among-taxon environment associations.

Do not merge these denominators into one apparent sample size or imply that the 46,276-observation cohort is complete for all nine constructs.

## Current release status

This role freeze does **not** yet satisfy the release builder's final figure-manifest gate.

Before `build_current_release_bundle --final` may be used:

1. render every final Main and Supporting export from its declared source;
2. complete visual and document-pagination QA;
3. freeze final numbering and captions;
4. create one checksum manifest containing the exact renderer/source/provenance/export files used by the submitted manuscript;
5. supply that manifest to the fail-closed release builder.

Until then, the existing `reproducibility/figures/` directory remains a labelled frozen v2 reference archive, not the current final manuscript figure set.
