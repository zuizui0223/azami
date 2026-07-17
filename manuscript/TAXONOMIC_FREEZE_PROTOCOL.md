# Chapter 1 taxonomic freeze protocol

## Purpose

The manuscript uses several related but non-identical taxon scopes: 216 accepted image-analysis taxa, 259 taxa in the strict within-species layer and 102 taxa in the complete lability cohort. Before Methods, tables and the submission archive are frozen, every name appearing in each source table must have one explicit, dated taxonomic decision.

This protocol does not automatically choose a taxonomy. It separates external name checking from the reproducible validation of the decisions used in the analysis.

## Decision table

Start from `manuscript/templates/taxonomic_decisions_template.csv` and create one row for every input name.

Required fields:

- `input_name`: exact spelling in the frozen analysis table;
- `accepted_name`: manuscript-facing accepted name;
- `decision`: `accepted`, `synonym`, `provisional` or `excluded`;
- `authority_source`: taxonomic authority consulted, such as POWO, World Flora Online, Catalogue of Life or a regional revision;
- `authority_record`: stable record identifier or URL where available;
- `checked_date`: date of the manual decision;
- `notes`: rationale for conflicts, infraspecific treatment, hybrids or regional usage.

A source name can map to one accepted name. Multiple source names may collapse to one accepted name only when the synonym decision is documented. Provisional names prevent the audit from passing.

## Reproducible validation

Run:

```bash
python analysis/audit_taxonomic_freeze.py \
  --taxon-table <continuous-species-table.csv> \
  --taxon-table <strict-within-species-table.csv> \
  --decisions <completed-taxonomic-decisions.csv> \
  --output-dir <taxonomic-audit-output>
```

Use `--name-column` only when the table does not use one of the recognized names: `taxon_name`, `scientific_name`, `species`, `species_name` or `accepted_name`.

The audit fails when:

1. an observed name has no decision;
2. a name has more than one decision row;
3. a required authority or date is missing;
4. a provisional name remains;
5. an excluded name is still present in an analysis input.

Outputs:

- `taxonomic_freeze_decisions.csv`: manuscript- and supplement-facing decision table;
- `taxonomic_freeze_summary.json`: counts, synonym collapses, unused decisions, errors and provenance.

## Scope rule

A taxonomic decision that changes accepted taxa, merges analysis units or removes observations is not documentation-only. It triggers a new analysis version under `analysis/ch1/pipeline.json`. Spelling normalization that leaves analysis units unchanged may be applied at the manuscript layer, but the original input name must remain in the decision table.

## Submission gate

Methods can describe the accepted-name procedure only after:

- every name in all frozen source cohorts is represented;
- no provisional or excluded name remains in an active cohort;
- the decision table is reviewed by a botanist familiar with Cardueae;
- the output table, source files, Git commit and SHA-256 checksums are copied into the durable submission bundle.
