# Chapter 1 document / figure QA — 2026-09-27

## Purpose

Track the final document-level gate separately from the already frozen numerical and figure-surface checks. This file does not change any scientific result, figure hash or manuscript claim.

## Canonical scientific/editorial source

For release QA, the current repository is authoritative for the active Chapter 1 scope:

- `analysis/v3/GEB_SCOPE_FREEZE_20260912.md`: scale-dependent construct integration is the principal current v3 comparative result; the 46,276-observation / 259-taxon atlas supplies the broad macroecological context; the two retained trait-environment associations remain ecological anchors.
- `reproducibility/CURRENT_FIGURE_SURFACE_20260912.md`: five Main figures are frozen in the order measurement -> global sampling -> scale integration -> construct-environment atlas -> anchor robustness.
- the live manuscript synchronization patch inspected on 2026-09-27 uses the same integration-centred role map.

Earlier manuscript/SI editorial states must not override that current freeze.

## Supporting Information visual QA

Inspected source:

- `Azami_Ch1_GEB_Supplement_v3_FINAL_20260917.pdf`
- 13 rendered pages
- accompanying DOCX exists as `Azami_Ch1_GEB_Supplement_v3_FINAL_20260917.docx`

Every rendered PDF page was visually inspected on 2026-09-27.

### Layout result: PASS

Observed across pages 1–13:

- no obvious page clipping;
- no duplicated rendered page;
- no missing full-page figure/table region;
- headings, tables and figure captions remain inside the page;
- Figure S2.1 is present and legible;
- dense tables remain small but readable at page scale;
- the final Data and code index / References block is present on the last page.

This closes **SI pagination/layout QA only**.

## Content synchronization status

The 2026-09-17 SI was produced during an earlier editorial reframe. The current repository freeze and the later live synchronization patch foreground scale-dependent integration more strongly than that earlier reframe.

Therefore:

- SI binary/layout integrity: **PASS**;
- SI content compatibility with the current frozen analyses: **usable but requires final Main synchronization check**;
- final cross-document labels/captions/primary-vs-secondary wording: **OPEN until the canonical final Main manuscript is available**.

No SI wording is silently rewritten by this receipt.

## Main manuscript status

A canonical final Main DOCX/PDF matching the current repository freeze was not found in the accessible Library or connected Drive during the 2026-09-27 audit.

Available materials include older Main manuscript files plus a later synchronization patch, but a patch is not a rendered submission document. It cannot establish:

- final figure placement;
- caption numbering;
- page breaks;
- clipping;
- duplicate/stale figures;
- Word/PDF pagination;
- identifying metadata removal.

Therefore Main document QA remains **OPEN**.

## Release-gate interpretation

The current figure artifact itself is checksum-frozen and visually QA'd. The document-level gate is narrowed to one missing object:

> the canonical final Main manuscript synchronized to the current five-Main-figure role map.

Do not set `document_pagination_validated=true` in the final figure manifest until that Main manuscript is rendered and inspected together with the already-passed SI.

## What closes the gate

Once the canonical final Main DOCX/PDF is supplied:

1. verify title/abstract/claim hierarchy against the current GEB scope freeze;
2. verify Main Figures 1–5 occur once, in the frozen role order, with matching captions;
3. verify Supporting references/labels against the SI;
4. render every Main page and inspect clipping, pagination, duplication and stale labels;
5. confirm the blinded file does not contain author-identifying metadata or author-identifying repository links;
6. promote the existing figure hashes without changing the figures and set `document_pagination_validated=true`.

Until those steps pass, `figure_document_qa` remains the sole document-side Zenodo release gap.
