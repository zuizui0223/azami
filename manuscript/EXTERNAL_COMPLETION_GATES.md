# External completion gates for Chapter 1

All manuscript and figure work that can be completed from the frozen repository and archived workflow outputs is now represented in `COMPLETE_DRAFT_FOR_SUPERVISOR.md`. The following gates require source material or scientific decisions that are not present in the repository.

## 1. Independent detector audit

The evaluation code and workflow exist, but completed manual source-image annotations and a finalized model card are absent. Precision, recall and F1 must not be inferred from the analysis detections themselves.

Required input: independent source-image annotations, detector predictions, model weights and identifier, training provenance and audit sampling design.

## 2. Authority-backed taxonomic freeze

The audit infrastructure exists, but the accepted-name decision table has not been populated and externally reviewed. Automated matching alone is insufficient for synonyms, approximate matches and taxonomically difficult Cardueae names.

Required input: one citable authority-backed decision per frozen source name.

## 3. Residual spatial and broad-region diagnostics

The audit and preparation workflows exist, but the frozen grouped-SPDE artifact lacks observation-level residuals or fitted values and the reviewed observation-to-region lookup. A new replacement model should not be introduced merely to fill this gate.

Required input: observation-level grouped-SPDE fitted or residual export and reviewed broad-region assignments.

## 4. Authorship and administrative metadata

The final author order, affiliations, CRediT roles, acknowledgements, funding numbers and corresponding-author details require confirmation from the research team.

## 5. Durable release

The final checksummed bundle must be deposited after review revisions are frozen. Add the immutable release tag, final commit and archive DOI to the manuscript only after deposition.
