# Results source ledger

This ledger records the frozen repository source for every numerical or inferential statement in `03_results.md`. It is deliberately separate from the prose so that manuscript editing cannot silently mix cohorts or artifact versions.

## Dataset scopes

| Result | Frozen source |
|---|---|
| 216 taxa and 6,626 detected capitula | `manuscript/final_claims.json` → `continuous_species_comparison` |
| 46,276 spatially thinned observations and 259 input taxa | `manuscript/final_claims.json` → `strict_within_species` |
| 102 complete taxa; n ≥ 10, ≥ 6 traits, ≥ 3 predictors per trait | `manuscript/final_claims.json` → `complete_lability_cohort` |

## Visible variance and PCA

| Result | Frozen source |
|---|---|
| Within-species visible-variance range 0.817–0.988 | integrated continuous-analysis variance-partition table; manuscript mapping in `FIGURE_TABLE_MAP.md` |
| PC1 32.9%; PC1–2 56.1%; PC1–3 69.3% | species-level PCA summary; manuscript mapping in `FIGURE_TABLE_MAP.md` |

Before submission, copy the exact source CSV names, workflow run ID, artifact digest and commit into this section from the durable bundle.

## Two-axis lability

| Result | Frozen source |
|---|---|
| Spearman rho = -0.333 | `manuscript/final_claims.json` → `primary_results` |
| Quadrants 22 / 29 / 29 / 22 | `manuscript/final_claims.json` → `primary_results.quadrants` |
| Colour medians 0.553 / 0.175 | `manuscript/final_claims.json` → `module_medians.colour` |
| Orientation medians 0.630 / 0.129 | `manuscript/final_claims.json` → `module_medians.orientation` |
| Shape medians 0.665 / 0.153 | `manuscript/final_claims.json` → `module_medians.shape` |
| Species-level output | `species_lability_axes.csv` |
| Module-level output | `module_lability_summary.csv`, `species_module_lability.csv` |

## Scale-specific environmental evidence

| Result | Frozen source |
|---|---|
| Strict ≤10 km cohort: zero BH-supported endpoint–climate effects | `manuscript/final_claims.json` → `strict_le_10km_primary_cohort` |
| Expanded pooled table: four BH-supported small effects | `manuscript/final_claims.json` → `expanded_pooled_coefficient_table` |
| Grouped spatial support concentrated in colour and orientation; no globally BH-supported shape effect | grouped SPDE-INLA fixed-effect and model-status tables; summary in `FINAL_MANUSCRIPT_STRATEGY.md` |
| Orientation has largest environmental niche centroid contrast | trait-extreme niche summary; summary in `FINAL_MANUSCRIPT_STRATEGY.md` |
| High overlap for aspect ratio, circularity and solidity | trait-extreme niche summary; summary in `FINAL_MANUSCRIPT_STRATEGY.md` |
| Climate-only preferred for orientation; soil improves preferred lightness/chroma models | grouped SPDE-INLA WAIC/model-status tables; summary in `FINAL_MANUSCRIPT_STRATEGY.md` |

Exact endpoint–predictor coefficients should not be inserted into the prose until Table 2 has been assembled from the frozen source tables.

## Minimum-sample sensitivity

| Threshold | Variation rank correlation | Responsiveness rank correlation | Source |
|---|---:|---:|---|
| n ≥ 5 | 0.939 | 0.999 | `sensitivity_minimum_sample_size.csv` |
| n ≥ 20 | 0.914 | 0.987 | `sensitivity_minimum_sample_size.csv` |

## Historical sensitivity and molecular coverage

| Result | Frozen source |
|---|---|
| 216 input taxa, 54 direct tips, 50 random trees, 636 fits, zero failures | `manuscript/final_claims.json` → `historical_sensitivity` |
| Zero direct-backbone-supported non-circular endpoints | `manuscript/final_claims.json` |
| Circularity and solidity grafting-sensitive | `manuscript/final_claims.json` |
| Hue sine/cosine requires joint circular test | `manuscript/final_claims.json` |
| Nucleotide 138; ITS 133; plastid 129; plastome 11; SRA 120 | `manuscript/final_claims.json` → `molecular_database_coverage` |

## Prohibited result-language substitutions

Do not replace the manuscript terms with the following stronger causal claims:

- `visible within-species variation` → not “phenotypic plasticity”;
- `environmental responsiveness` or `climate association` → not “local adaptation”;
- `environmental niche contrast` → not “selection” or “adaptive divergence”;
- `historical sensitivity` → not a resolved *Cirsium* species tree;
- orientation–environment structure → not demonstrated rain-protection adaptation;
- visible colour structure → not demonstrated pollinator causation or anthocyanin regulation.

## Remaining verification gate

The prose is now constrained to frozen claims, but the durable submission bundle still needs to record:

1. workflow run IDs;
2. artifact IDs and SHA-256 digests;
3. exact internal CSV paths;
4. final Git commit;
5. software/session information.
