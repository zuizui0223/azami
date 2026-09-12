# Chapter 1 GEB Supporting Information numbering plan — 2026-09-12

This is a **presentation/organization plan**, not a scientific reanalysis. It preserves the existing Supporting Information content and separates the frozen full-27 baseline from later v3/revision diagnostics.

## Current state

The current Library supplement is `Azami_Chapter1_Supplement_FINAL.pdf`. It already contains **Appendix S1. Full-27 continuous-trait and full-environment analysis** with detailed methods, audit boundaries and frozen supporting outputs.

Its existing display inventory includes:

- Figures S1–S7;
- Tables S1–S12.

The 2026-09-07 Main manuscript also cites later revision-diagnostic Tables V1–V9 that are not part of the old S1–S12 first-citation sequence. The current GEB submission should not ship two competing supplementary numbering systems.

## Frozen organization

Use two Supporting Information appendices.

### Appendix S1 — Frozen full-27 baseline and original robustness chain

Preserve the current Appendix S1 text and scientific order. Renumber displays mechanically:

| Existing label | GEB submission label |
|---|---|
| Figure S1 | Figure S1.1 |
| Figure S2 | Figure S1.2 |
| Figure S3 | Figure S1.3 |
| Figure S4 | Figure S1.4 |
| Figure S5 | Figure S1.5 |
| Figure S6 | Figure S1.6 |
| Figure S7 | Figure S1.7 |
| Table S1 | Table S1.1 |
| Table S2 | Table S1.2 |
| Table S3 | Table S1.3 |
| Table S4 | Table S1.4 |
| Table S5 | Table S1.5 |
| Table S6 | Table S1.6 |
| Table S7 | Table S1.7 |
| Table S8 | Table S1.8 |
| Table S9 | Table S1.9 |
| Table S10 | Table S1.10 |
| Table S11 | Table S1.11 |
| Table S12 | Table S1.12 |

Do not change the underlying data or inferential status during renumbering.

### Appendix S2 — v3 integration, revision diagnostics and submission-readiness sensitivities

Move the later `V`-series diagnostic surfaces out of the Main numbering namespace and into Appendix S2.

| Working label | GEB submission label | Role |
|---|---|---|
| Table V1 | Table S2.1 | predictor-correlation / interpretation diagnostic |
| Table V2 | Table S2.2 | colour interpretation diagnostic |
| Table V3 | Table S2.3 | recovered orientation / mirror-threshold diagnostic |
| Table V4 | Table S2.4 | conditional candidate models |
| Table V5 | Table S2.5 | precipitation-season comparison |
| Table V6 | Table S2.6 | minimum-observation threshold sensitivity |
| Table V7 | Table S2.7 | individual-taxon slope diagnostic |
| Table V8 | Table S2.8 | linked bounding-box remeasurement sensitivity |
| Table V9 | Table S2.9 | pairwise measured-component coefficients/coverage |
| new WCVP authority sensitivity summary | Table S2.10 | accepted-name / synonym-collapse sensitivity |

The exact captions must be reconciled against the final source tables before export; this map fixes **identity and order**, not wording that has not yet been verified.

## Figure demotion/addition in Appendix S2

The old standalone Main taxon-mean information-loss figure is scientifically valid but no longer a Main display under the v3 story. Move its final export to:

- **Figure S2.1 — Taxon-mean information loss across measured endpoints.**

If any new diagnostic figures are retained for the moved Main post-hoc material, assign them sequentially as Figure S2.2, S2.3, ... only after their final source/provenance is frozen. Do not renumber Appendix S1 figures to accommodate later v3 diagnostics.

The new scale-dependent integration figure stays in the **Main paper** and must not be duplicated as a Supporting figure.

## Main-text citation replacements

During final DOCX synchronization, perform explicit replacement rather than broad search/replace:

- old `Table S1`–`Table S12` citations -> `Table S1.1`–`Table S1.12` respectively;
- old `Figure S1`–`Figure S7` citations -> `Figure S1.1`–`Figure S1.7` respectively;
- `Table V1`–`Table V9` -> `Table S2.1`–`Table S2.9`;
- new authority-taxonomy table -> `Table S2.10`;
- old Main information-loss figure references must be rewritten so the numerical result remains in Main text but the display points to `Figure S2.1`.

Ranges must also be transformed semantically: for example `Tables S6–S8` becomes `Tables S1.6–S1.8`, not `Tables S1.6–S8`.

## Content moves from Main to Appendix S2

To meet the Main-text word budget without discarding evidence, Appendix S2 should absorb the detailed versions of:

1. post-hoc predictor covariance and colour interpretation;
2. recovered-orientation and mirror-threshold checks;
3. conditional/seasonal candidate models;
4. threshold and individual-taxon diagnostics;
5. linked bounding-box remeasurement;
6. pairwise-complete measured-component diagnostic;
7. WCVP accepted-name/synonym-collapse sensitivity;
8. any detailed tables supporting the compact Main statement that physical colour and gravity-referenced orientation remain unresolved.

Main retains only the conclusions from these diagnostics that materially change interpretation.

## Final SI QA gate

Before submission, verify all of the following against the final Supporting Information file:

1. every `S1.x` and `S2.x` display is present exactly once;
2. every Supporting display is cited from Main or the relevant Appendix text;
3. there are no residual `Table V*`, `Figure V*`, bare `Table S1`–`S12`, or bare `Figure S1`–`S7` citations where Appendix numbering is required;
4. Appendix S1 frozen baseline statements are not silently updated with v3 outcomes;
5. Appendix S2 labels all post-hoc diagnostics as post-hoc and the WCVP sensitivity as predeclared before its outcome;
6. display captions retain cohort denominators and claim boundaries;
7. file names, captions, in-text citations and final figure/table manifest agree exactly.
