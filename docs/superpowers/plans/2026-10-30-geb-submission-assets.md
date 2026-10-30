# GEB submission assets implementation plan

> For agentic workers: execute the independently scoped figure, document and aggregate-package tasks with review before integration. The user has authorised continuation of submission preparation.

**Goal:** turn the source-verified V4 manuscript into an editable review package with actual photographic and sampling figures.

**Architecture:** reuse the frozen cohort selection and coefficient-verification functions. Restore protected inputs only in the existing authenticated Actions environment; export pooled 5-degree count cells with at least five records. Assemble manuscript documents and an explicitly limited aggregate-verification package from a file allowlist.

**Tech Stack:** Python, pandas, matplotlib, python-docx, existing GitHub Actions source recovery.

**Spec:** current manuscript, `SCOPE_LOCK_20261025_JA.md`, and GEB author guidelines checked 30 October 2026.

## Global constraints

- Keep original S1/S2/H3c inference, numerical results and test families unchanged.
- Reproduce S1 angle n=168907/taxa=178 and all six S2 coefficient/denominator pairs before labelling coverage source-verified.
- Reuse H3c quality and geographic support gates: 63 taxa/7117 taxon-cell pairs at 1 degree; 31/814 at 5 degrees.
- Never export original observation, uploader or taxon identifiers, exact coordinates or private environmental pixel records.
- Build Figure 2 with a named equal-area projection, graticules, equatorial scale and explicit suppression accounting.
- Retain exact source-photo credits and applicable licences; distinguish measurement examples from analytical cohorts.
- No journal submission, DOI publication, private release visibility change or messages to third parties.

## Review focus

- Cell-edge and polar coordinates: valid indices, no dateline wrap or displaced points.
- Sparse cells: counts below five omitted, omissions counted, no individual identifiers admitted by schema.
- Different analytical cohorts: each panel labelled with its own denominator, overlapping data never summed as independent records.
- Historical assets: verify hashes and remove obsolete V2/V3 construct counts and incorrect orientation overlay.
- Anonymous review: scan document metadata and links; never call an aggregate-only package a complete model refit.

## Tasks

- [x] Write and run failing aggregation/privacy tests, implement `analysis/v4/submission_sampling_coverage.py`, rerun tests.
- [x] Add source-only recovery workflow with explicit aggregate upload paths; run against authenticated original inputs and recover the successful coverage receipt (run 41020812565).
- [x] Render Figure 2 from the verified receipt, inspect the full figure, and preserve a vector master.
- [x] Build Figure 1 from the existing licensed source composite with a verified rights manifest.
- [x] Update figure captions/callouts and availability status only after evidence exists; preserve scientific sections and coefficients.
- [x] Build and visually verify editable main/SI documents and a bounded aggregate/code package.
- [x] Review the complete diff and run the manuscript, bundle and privacy/cohort contracts (22 + 8 + 5 tests).

## Integration gate

Push the reviewed asset/source tree with an expected-head lease. Verify the resulting branch head and its workflows before reporting completion; the live PR checks are the authoritative evidence for that post-commit step. Keep the source-coverage run 41020812565 as the data-recovery receipt, rather than treating an asset-only CI run as a new source reconstruction.
