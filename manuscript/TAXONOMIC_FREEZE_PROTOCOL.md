# Chapter 1 operational taxonomic freeze protocol

## Purpose

Chapter 1 is an image-phenomics and macroecology analysis, not a taxonomic revision of *Cirsium*. Public observations arrive with source-platform taxon assignments, while global and regional authorities can legitimately differ in how several *Cirsium* boundaries are lumped or split. Retrospectively forcing one checklist treatment onto photographs cannot adjudicate those boundaries without specimen-level evidence.

The manuscript therefore distinguishes two questions:

1. **What taxonomic unit was assigned to each frozen public observation?** This is the operational analysis unit and must remain traceable exactly.
2. **Would plausible authority-backed lumping change the ecological conclusions?** This is tested explicitly as taxonomic sensitivity rather than hidden inside name cleaning.

The manuscript should use *source-assigned taxonomic unit*, *assigned taxon* or *assigned species* where the distinction matters. It must not imply that the analysis resolved genus-wide species limits.

## WCVP authority audit

`analysis/build_wcvp_taxonomic_review_candidates.py` checks the union of frozen source names against the World Checklist of Vascular Plants (WCVP) through the GBIF checklist API and records raw responses and stable authority identifiers.

Successful run `31152400481` covered 259 source names:

- 243 unique exact canonical-name matches;
- 16 multiple exact canonical-name matches, generally because an accepted record and a homonymous/synonym record share an author-free canonical name;
- 235 routine accepted-name candidates;
- 8 WCVP synonym candidates;
- 24 high-priority nomenclatural-review rows;
- 0 unmatched names.

Automated matching remains candidate generation, not a claim that WCVP must define the operational analysis units.

## Synonym-collapse sensitivity

All eight WCVP synonym candidates map to names that already occur as active source-platform analysis units. They were therefore collapsed simultaneously in an explicit sensitivity analysis rather than silently imposed on the primary data.

Workflow run `31153385658` shows:

- balanced atlas: 216 source-assigned taxa -> 211 collapsed units;
- exhaustive primary: 259 source-assigned taxa -> 251 collapsed units;
- nested minimum below-unit visible-variance fraction: 0.5886 -> 0.5970;
- maximum absolute endpoint change in below-unit fraction: 0.0154;
- primary 36 climate models: 8 BH-supported component rows before and after collapse;
- FDR decisions changed: 0;
- coefficient signs changed: 0;
- maximum absolute standardized-beta change: 0.000385.

Thus the two headline ecological result families are insensitive to simultaneously applying all eight WCVP synonym collapses.

Full provenance is recorded in `manuscript/results/WCVP_TAXONOMIC_SENSITIVITY_2026-08-07.md`.

## Operational freeze rule

For the submitted primary analysis:

- retain the exact source-platform taxon assignment as the observation-level operational unit;
- preserve the WCVP candidate accepted name, status and authority record in the supplementary taxonomic audit table;
- report the simultaneous WCVP synonym-collapse sensitivity;
- do not silently relabel or merge source observations in manuscript-facing data products;
- describe counts as assigned taxonomic units/taxa where a resolved accepted-species interpretation is not required.

This policy prevents a global checklist from retroactively creating unverified specimen-level identifications while still testing whether checklist lumping explains the main ecological results.

## When a new analysis version is required

A new analysis version is required if the authors decide to make a different taxonomy the **primary analysis unit**, because merging, removing or reassigning source units changes cohort counts and species-level products. A spelling/authorship display correction that leaves the operational unit unchanged is documentation-only.

## Supplementary nomenclatural review

The 24 high-priority WCVP rows should remain visible in the supplementary taxonomic table. Regional authorities or a Cardueae specialist may add notes explaining preferred nomenclature, especially for the eight WCVP synonym conflicts. Such notes improve nomenclatural reporting but, under the operational-unit design and successful collapse sensitivity, are no longer a hidden prerequisite for the validity of the two headline ecological analyses.

## Submission gate

The scientific taxonomic-robustness gate is satisfied when:

- all source names are authority-audited and traceable;
- synonym conflicts are explicitly identified;
- a predeclared lumping sensitivity shows whether those conflicts alter headline results;
- the manuscript consistently distinguishes source-assigned analysis units from a resolved species taxonomy;
- authority tables, sensitivity outputs, code, commit identity and checksums enter the durable archive.

Those conditions are met for the frozen Chapter 1 analysis. Final nomenclatural annotations for the 24 high-priority rows remain desirable supplementary metadata, not an untested source of ecological inference.
