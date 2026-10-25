# Chapter 1 v3 direct integration-scale contrast — 2026-09-10

## Purpose

This result makes the reviewer-highlighted scale dependence explicit rather than inferring it from separate within- and among-taxon significance tables. It reuses the exact complete-18 common cohort and the existing 1,000-replicate taxon bootstrap. No frozen v2 endpoint result, multiplicity family, or manuscript headline candidate is changed.

## Common-cohort observed contrast

The same nine biological constructs and all 36 construct relations are evaluated on 1,734 observations from 42 taxa.

- median within-taxon RV = **0.002238**;
- median among-taxon RV = **0.043212**;
- difference of medians (among - within) = **+0.040974**;
- median pairwise relation delta (among - within) = **+0.031431**;
- **33/36** construct relations are stronger among taxa than within taxa.

This is a contrast in visible-image phenotypic integration. It is not evidence that evolution necessarily increases integration.

## Taxon-bootstrap support

Across 1,000 taxon-bootstrap replicates:

- median difference of median RV (among - within) = **+0.065928**;
- bootstrap 95% interval = **+0.038197 to +0.119048**;
- fraction of bootstrap replicates with median among-taxon RV > median within-taxon RV = **1.000**;
- median number of relations stronger among taxa = **32/36**;
- bootstrap 95% range for that count = **29 to 34 / 36**;
- fraction of replicates in which a majority of relations are stronger among taxa = **1.000**.

Thus the tendency toward stronger among-taxon integration is not driven by one or two taxa in the complete-18 cohort. The magnitude remains a property of this image-defined observational system and should not be interpreted as a causal evolutionary transition.

## Relationship to the existing scale result

This direct strength contrast complements, rather than replaces, the existing matrix-alignment result:

- within-vs-among matrix alignment on the common cohort: rho = **0.439125**, QAP P = **0.0041**;
- 1,000-taxon-bootstrap alignment 95% interval = **0.01708 to 0.53517**, with **97.9%** positive;
- module cohesion is supported within taxa (P = **0.0013**) and among taxa (P = **0.0365**).

Together these results support a more precise statement:

> Capitulum organization is partially conserved across scales, but the among-taxon phenotype is more strongly integrated overall while retaining a reorganized relation geometry.

That is stronger and more informative than saying only that different endpoint-environment associations are significant at different scales.

## Relationship to frozen v2 positive results

The result does not modify the frozen v2 ecological conclusions. The two headline among-taxon patterns remain:

- lower floral chroma under higher shortwave radiation;
- larger image-referenced presentation angle under higher annual precipitation.

Both have already passed biological reaggregation plus the construct-level sampling, spatial/residual and 52-tree historical sensitivity sequence. This scale-contrast layer explains the broader phenotypic context in which those two ecological anchors sit; it does not add a new headline environment association.

## Reproducibility

- workflow run: `34435743710`
- job: `102740315624`
- artifact: `10136229131`
- artifact digest: `sha256:7b4c5703b048ff108cd794921394ae22d1423461c9aaf330b47029c5d1f4908f`
- common cohort: `1,734 observations / 42 taxa`
- bootstrap replicates: `1,000`
- machine-readable summary: `construct_scale_contrast_summary.json` in the workflow artifact

Claim boundary: secondary direct scale contrast of image-defined biological constructs; not functional, developmental or genetic modularity, not proof of plasticity, and not a replacement for frozen v2 endpoint inference.
