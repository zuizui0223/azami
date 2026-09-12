# GEB submission audit — 2026-09-12

This file separates **journal-format / peer-review-access requirements** from scientific analysis. It does not alter any Chapter 1 result. Source checked on 2026-09-12: Wiley, *Global Ecology and Biogeography* Author Guidelines (`https://onlinelibrary.wiley.com/page/journal/14668238/homepage/forauthors.html`).

## Submission category

Target category: **Research Article**.

Current scientific scope remains the frozen Chapter 1 GEB scope in `analysis/v3/GEB_SCOPE_FREEZE_20260912.md`.

## Hard submission gates

| Requirement | GEB requirement | Current Chapter 1 state | Gate |
|---|---|---|---|
| Structured abstract | Aim; Location; Time period; Major taxa studied; Methods; Results; Main conclusions; <=300 words | Synchronization patch revised to 262 words | **PASS in patch; must be verified in final DOCX** |
| Running title | <40 characters | `Scale-dependent capitulum integration` = 37 characters | **PASS in patch** |
| Keywords | 6–10, alphabetical | Eight current keywords; final manuscript ordering must be checked after synchronization | **VERIFY FINAL DOCX** |
| Main length | Research Articles typically ~5000 main-text words | Do not infer from the 2026-09-07 working draft after v3 replacement | **COUNT FINAL SYNCHRONIZED DOCX** |
| Display pieces | Typically 6–8 tables + figures | Frozen role plan is five Main figures; retain only the Main tables needed to keep total in range | **VERIFY FINAL DOCX** |
| Double-anonymous main file | No author-identifying material; separate identifying title page; remove author name from Word properties | Existing public GitHub/Zenodo identities must not be used as the blinded review link | **OPEN** |
| Data/code during peer review | Data and code supporting results accessible during review; GitHub alone is not a stable archive | Current public v2 Zenodo is incomplete for v3 and may reveal authors | **OPEN: anonymous review archive required** |
| Permanent archive | Stable public repository for publication | Final current-release builder and Zenodo metadata gate exist; final v3 public release not yet published | **OPEN** |
| Figure readability | Embedded for review; legends stand alone; panels `(a)`, `(b)`, ... | Current v3 scale-integration renderer passes fixed-width visual QA | **PARTIAL: full Main/SI surface still required** |
| Figure legend content | Stand-alone; state organism and geographic region where applicable | New scale-integration caption must explicitly say global Cirsium dataset | **PATCH BEFORE FINAL DOCX** |
| Supporting Information | Separate file(s); all SI objects cited; Appendix-based numbering (e.g. Fig. S1.1) | Existing historical `Figure S1`, `Table S1`, etc. notation must not be assumed compliant | **OPEN: renumber final SI** |
| Cover letter | Separate PDF; <250-word paragraph explaining reader interest | Not part of scientific repository | **AUTHOR/SUBMISSION TASK** |
| Title page | Author names, affiliations, emails, acknowledgements; only one corresponding author | Deliberately absent from blinded repo surface | **AUTHOR/SUBMISSION TASK** |
| ORCID | Required at submission | Account/author-owned | **AUTHOR/SUBMISSION TASK** |

## Anonymous peer-review archive

GEB requires review access to supporting data/code while operating double-anonymous peer review. The blinded manuscript therefore must **not** point reviewers to an author-identifying repository URL merely because that repository is public.

The journal explicitly documents Dryad's **Private for Peer Review** route, which creates a randomized private review URL. The preferred Chapter 1 submission path is therefore:

1. build the checksum-gated current package from the final frozen code/figure surface;
2. deposit the review package in an anonymous review archive (preferred: Dryad Private for Peer Review, unless an equivalent non-identifying stable review route is approved);
3. verify the review link in a logged-out / non-owner context;
4. inspect the landing page and downloaded package for author names, usernames, e-mail addresses, local paths, ORCID, repository-owner identifiers and identifying metadata;
5. put only that anonymous review link in the blinded Data and Code Availability Statement;
6. retain the permanent public Zenodo release as the publication archive after the double-anonymous review constraint no longer applies.

A GitHub URL is useful for development and audit but **does not satisfy GEB's stable-archive requirement by itself**.

## Blinded Data and Code Availability wording

Do not use this sentence until the anonymous review package exists and has passed the logged-out anonymity check:

> Data and code supporting this study are available to reviewers through an anonymized private archival link supplied with the submission. The review package contains the frozen numerical inputs, analysis code, aggregate reference outputs, figure provenance and numerical replay receipts. A permanent public archival record will be released for publication.

For the final published manuscript, replace the review sentence with the permanent public archive DOI(s), after anonymous redownload/checksum/clean-replay validation.

## Display-piece plan

The current five-Main-figure role map is frozen in `analysis/v3/CURRENT_FIGURE_SURFACE_20260912.md`:

1. measurement -> biological construct workflow;
2. realized global sampling domain;
3. scale-dependent construct integration;
4. complete construct x environment atlas;
5. robustness of the two retained ecological anchors.

The old standalone taxon-mean information-loss figure is Supporting Information. This preserves display space for the paper's current conceptual contribution rather than the superseded v2 hierarchy.

## Final pre-submission sequence

1. Finish the predeclared WCVP taxonomy sensitivity and retain the outcome whether positive or adverse.
2. Synchronize the actual blinded DOCX with the current v3 manuscript patch.
3. Count the final structured abstract and main text; verify running title and alphabetical 6–10 keywords.
4. Freeze all Main/SI figure and table captions, numbering and checksums; convert SI numbering to the journal's Appendix scheme.
5. Create and independently inspect the anonymous peer-review archive/link.
6. Remove identifying metadata from the blinded DOCX and inspect links/file properties.
7. Supply a separate identifying title page and cover letter outside the blinded main file.
8. Only after the final public archive is published, perform credential-free redownload, checksum verification and clean replay before making the permanent availability claim.
