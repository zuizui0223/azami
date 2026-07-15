# Chapter 1 final freeze checklist

## Analysis is frozen

Do not add another headline trait, climate variable or evolutionary model unless a current result is shown to be invalid. The central two-axis lability result, module contrast, minimum-sample sensitivity and phylogenetic robustness decision are frozen.

## Gates that can still affect credibility

### Measurement demonstration

- Export representative low, median and high examples for every continuous endpoint.
- Include QC failures, camera-tilt cases, occlusion and background leakage.
- Include horizontal-flip technical replicates.
- Add a compact manual/published-description face-validity comparison.

### Detector evaluation

- Freeze an independent source-image subset.
- Report precision, recall, missed-head categories and duplicate detections.
- State that trait QC cannot diagnose undetected capitula.

### Taxonomic freeze

- Resolve approximate, synonym and multiple-match names.
- Commit accepted names, original names, decision reason and evidence source.
- Re-run only the name-join and downstream manifest if names change.

### Spatial credibility

- Add residual spatial autocorrelation diagnostics for headline between-species models.
- Summarize broad-region or spatial-block sensitivity.
- Keep mesh and predictor-correlation sensitivity attached to the between-species model family; do not relabel it as sensitivity of species-specific lability slopes.

### Circular hue

- Retain sine and cosine internally.
- Do not interpret them as separate biological traits.
- A joint circular uncertainty test is optional for Supplementary material and is not a submission blocker.

## Manuscript package

- Main text and Supplement use `manuscript/final_claims.json` as the numerical source of truth.
- Figure captions distinguish 216-taxon species summaries, the 259-taxon within-species input and the 102-taxon complete lability cohort.
- Results distinguish zero FDR signals in the strict <=10 km primary cohort from four small pooled signals in the expanded coefficient table.
- Phylogenetic results are described as sensitivity to historical placement, not a resolved species-tree result.
- Methods state the image-derived and observational interpretation boundary.

## Durable release

- Build the submission bundle on the exact final commit.
- Run `69_validate_submission_outputs.py`, `70_build_submission_manifest.py` and `analysis/validate_final_claims.py`.
- Record Python environment, R `sessionInfo()`, workflow run IDs, source artifact IDs, checksums and commit hash.
- Deposit tables, figures, manifests and permitted metadata in a permanent archive.
- Do not rely on expiring GitHub Actions artifacts as the cited data source.

## Stop rule

The analysis is submission-ready when all credibility gates above are documented and the archived bundle passes every validator. New exploratory ideas move to a separate branch or future paper and must not alter the frozen Chapter 1 release.
