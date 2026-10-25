# Chapter 1 WCVP authority-taxonomy sensitivity — result

## Outcome

**All three predeclared taxonomy-sensitivity gates passed.** The current Chapter 1 headline is robust to replacing uniquely resolved source taxon labels with accepted WCVP keys and collapsing source labels that map to the same accepted key.

This does **not** replace the source-assigned primary analysis. It addresses the narrower synonym / accepted-name grouping threat that had been left open for GEB submission readiness.

## Outcome-blind remapping

The frozen 46,276-observation / 259-source-taxon cohort was joined to the already frozen WCVP name-resolution table. Only rows with `resolution_status=resolved_unique_accepted_key` were retained; `taxon_name` was replaced by `wcvp:{accepted_key}` before the existing v3 methods were rerun. No trait value, environmental value, effect size or significance entered the retention/grouping rule.

- observations retained: **40,785 / 46,276**;
- source taxa uniquely resolved: **245 / 259**;
- accepted WCVP taxa after synonym collapse: **238**;
- accepted taxa receiving >1 source label: **7**;
- complete common cohort for the direct integration test: **1,415 observations / 38 taxa**.

## Gate 1 — cross-scale geometry

Predeclared rule: matrix-alignment rho > 0, one-sided QAP P < 0.05, and taxon-bootstrap rho 95% lower bound > 0.

- rho = **0.293179**;
- QAP P = **0.0385**;
- bootstrap median rho = **0.265251**;
- bootstrap 95% interval = **0.008745 to 0.490631**.

**PASS.** Cross-scale relation geometry remains positively aligned after authority remapping, although the alignment is weaker than in the source-assigned common cohort.

## Gate 2 — stronger integration among taxa

Predeclared rule: bootstrap lower 95% bound of the among-minus-within median-RV difference > 0 and a majority of the 36 observed relations stronger among taxa.

- median within-taxon RV = **0.002594**;
- median among-taxon RV = **0.044057**;
- **32/36** relations stronger among taxa;
- bootstrap difference median = **+0.067758**;
- bootstrap 95% interval = **+0.036905 to +0.121980**;
- fraction of bootstrap replicates with among > within median RV = **1.000**.

**PASS.** The direct scale-strength result is retained after accepted-name resolution and synonym collapse.

## Gate 3 — two ecological anchors

Predeclared rule: each frozen headline construct-environment association must retain its frozen sign and BH q < 0.05 under authority-resolved grouping.

| Anchor | Authority n taxa | beta | BH q | Result |
|---|---:|---:|---:|---|
| floral chroma x shortwave radiation | 131 | -0.327859 | 0.00675 | **PASS** |
| presentation angle x annual precipitation | 129 | +0.307370 | 0.00675 | **PASS** |

**PASS.** Both ecological anchors retain direction and multiplicity-adjusted support.

## Interpretation boundary

The defensible statement is now stronger but still bounded:

> Scale-dependent visible-phenotype integration and the two retained ecological anchors are robust to the predeclared sensitivity that replaces uniquely resolved source labels with accepted WCVP keys and collapses synonyms.

Do **not** turn this into a claim that all source records are correctly identified. The sensitivity does not resolve misidentified observations, disputed species boundaries, hybrid/cytotype structure, ancestry, genetic lineages, adaptation or causal evolutionary mechanism.

## Reproducibility

- workflow run: `34672948337`;
- artifact: `10292140117`;
- artifact SHA-256: `dfb6eec3001e3a984662d5aba06cda5fa80e144b36ccb4af9fdf45973854edc5`;
- machine-readable frozen result: `analysis/v3/authority_taxonomy_sensitivity_summary_20260912.json`.
