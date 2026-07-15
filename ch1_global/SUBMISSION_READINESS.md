# Chapter 1 — submission readiness

Review date: 2026-07-15  
Scope: global continuous capitulum traits in *Cirsium*.

## Current verdict

**The core analytical results and manuscript claims are frozen. The paper is not yet
ready to submit because the remaining credibility displays, taxonomic decisions,
spatial diagnostics, manuscript package and durable archive are incomplete.**

The headline analysis no longer relies on unvalidated categorical CLIP states. It
uses continuous deterministic image measurements with explicit QC, missingness,
horizontal-flip repeatability and separate within- and among-species inference.

## Completed

- 6,626-head continuous measurement run from 3,725 primary observations;
- 216 accepted image-analysis taxa;
- strict within-species climate input of 46,276 observations from 259 taxa;
- final two-axis lability cohort of 102 complete taxa;
- QC-retention bias audit by taxon, space and CHELSA environment;
- primary traits joined to environmental metadata;
- high-resolution involucre/spine exploratory supplement;
- hair, mucilage, ploidy and hybrid evidence queues with `unknown` as default;
- separate within-species and between-species models;
- species within-variation and environmental-responsiveness indices;
- minimum-n sensitivity at 5, 10 and 20 observations;
- four-quadrant and module-level summaries;
- bounded Pagel-lambda and Brownian PGLS sensitivity across deterministic and 50
  random grafting trees;
- phylogenetic-signal audit: 636/636 fits successful, 54/216 direct backbone tips;
- live NCBI molecular-database audit;
- final claim registry, validator and CI workflow;
- submission-focused documentation and output-contract tooling.

## Claims currently supported

The manuscript may report:

- global climatic structure in continuous capitulum traits, with clear separation
  of among- and within-species scales;
- a two-axis decomposition of lability into visible within-species variation and
  species-specific climate tracking;
- a moderate negative association between those axes across 102 complete taxa
  (Spearman rho = -0.333);
- module differences: colour shows the strongest median climate tracking, whereas
  shape shows the greatest visible variation;
- zero BH-FDR-significant primary effects in the strict <=10 km cohort;
- four small BH-FDR-significant effects in the expanded pooled coefficient table;
- differential image assessability as a quantified measurement limitation;
- historical/phylogenetic results as grafting-sensitive supplementary evidence;
- exploratory high-resolution involucre/spine patterns.

Every FDR statement must name its cohort. The zero-effect and four-effect summaries
must never be merged into a single unnamed within-species result.

The manuscript must not claim:

- local adaptation or phenotypic plasticity;
- pollinator causation;
- evolutionary rate;
- a resolved *Cirsium* species tree;
- botanical absence from an unmeasurable image;
- categorical bract, spine, hair or mucilage states from continuous proxies;
- independent biological meaning for hue sine and cosine components.

## Remaining gates before submission

### 1. Measurement demonstration

Create a compact, independently inspectable supplementary panel showing:

- representative low, middle and high values for each continuous variable;
- QC failures;
- horizontal-flip technical replicates;
- camera tilt, occlusion and background leakage;
- a small comparison with published descriptions or manual measurements.

This is a face-validity and interpretation audit, not a training-label exercise.

### 2. Detector performance

Report precision, recall and error categories for the frozen capitulum detector on
an independent source-image subset. Trait QC cannot evaluate heads that were never
detected.

### 3. Taxonomic freeze

Review approximate, synonym and multiple-match records. Commit a table containing:

- source name;
- accepted study name;
- external matched name;
- match type;
- decision;
- decision reason.

Do not silently accept an external synonym when it conflicts with the study
taxonomy.

### 4. Spatial robustness

Add the final residual spatial diagnostic and broad-region/block sensitivity
summary. This is especially important for among-species climate associations that
may reflect geographic turnover.

### 5. Manuscript package

Complete:

- title and abstract;
- Introduction, Results, Discussion and Methods;
- main and supplementary figure captions;
- supplementary tables;
- data-availability and code-availability statements;
- author contributions and acknowledgements.

All headline numbers must be copied from `manuscript/final_claims.json`.

### 6. Durable archive and licences

GitHub Actions artifacts expire and are not citable. Deposit the exact final tables,
figures, input manifests, configurations, software environments and checksums in a
permanent repository. Include public-observation identifiers and licence fields
without redistributing images beyond their licence terms.

## Submission freeze rule

The final release must pass:

- `analysis/validate_final_claims.py`;
- `69_validate_submission_outputs.py`;
- `70_build_submission_manifest.py`;
- `Ch1 submission code CI`;
- `Ch1 final manuscript claims`.

Any rerun that changes a primary table, figure or headline number requires a new
analysis version and manifest. Do not overwrite the submitted result.

## Bottom line

The remaining work is validation display, detector reporting, taxonomic and spatial
review, manuscript assembly and durable archiving. It is not another redesign of
the core trait pipeline and not an invitation to add more headline analyses.
