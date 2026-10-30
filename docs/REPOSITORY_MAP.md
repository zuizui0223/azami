# Repository map — active science, manuscript, archives

**Snapshot:** 2026-10-30, working branch of [PR #159](https://github.com/zuizui0223/azami/pull/159). This is a navigation and ownership map, **not a new scientific analysis**.

## Source of truth, in reading order

| Layer | Canonical path | Status and meaning |
| --- | --- | --- |
| **Paper** | [GEB Research Article](../manuscripts/v4_s1_s2_frozen_20261025/GEB_RESEARCH_ARTICLE.md) | Active, revised three-results narrative: 3.1 BIO18–angle; 3.2 colour–environment; 3.3 local climate analogue–joint phenotype resemblance. |
| Anonymous manuscript text | [Blinded version](../manuscripts/v4_s1_s2_frozen_20261025/GEB_RESEARCH_ARTICLE_BLINDED.md) | Working file; confirm anonymisation and final disclosures before submission. |
| Figures/legends | [Figures 1–5](../manuscripts/v4_s1_s2_frozen_20261025/figures/) and [submission asset audit](../manuscripts/v4_s1_s2_frozen_20261025/GEB_SUBMISSION_ASSET_AUDIT_20261030_JA.md) | Actual source-backed PNG/PDF/SVG figures, rights and render manifests; outdated exploratory `figure2_` names outside this directory are not main Figure 2. |
| Source statistics | [Frozen V4 analysis](../analysis/v4/README.md) and [source results](../analysis/v4/results/) | S1/S2 original corrected families and separate exploratory H3c/H3d source-verified aggregates. |
| Interpretation limits | [Claim ledger](../manuscripts/v4_s1_s2_frozen_20261025/CLAIM_LEDGER_AND_SUBMISSION_GATES.md) | Supported claims and excluded causal/phylogenetic/introduced-evolution interpretations. |
| Next submission steps | [GEB checklist](../manuscripts/v4_s1_s2_frozen_20261025/GEB_SUBMISSION_CHECKLIST.md) | Already made: figure and editable reviewer candidate. Outstanding: authorship declarations, permissions as appropriate, permanent rights-cleared data/code access, final anonymous review. |
| Reproducibility | [Runbook](REPRODUCIBILITY_AND_SUBMISSION_RUNBOOK.md) | From zero-dependency navigation tests to checksum-locked archived aggregate reproduction; never publish private intermediates. |

## Directory ownership

| Location | Classification | Move/delete rule |
| --- | --- | --- |
| [`manuscripts/v4_s1_s2_frozen_20261025/`](../manuscripts/v4_s1_s2_frozen_20261025/) | **Current** article text, figures, SI, builders and source-ledgers. Directory name refers to initial S1/S2 freeze, *not* a prohibition on the later H3c exploratory analysis. | Keep paths stable while PR #159 and figure scripts reference them. |
| [`analysis/v4/`](../analysis/v4/) | **Current** original source scripts, contracts and source-verified aggregate results; some exploratory branches are not headline claims. | Do not silently refit or rewrite frozen S1/S2 results. |
| [`tests/`](../tests/) and [`.github/workflows/`](../.github/workflows/) | **Verification** of source locks, manuscript integrity, scientific methods and controlled source recovery. | Preserve until dependency audit. Six active V4 workflow files as of this snapshot. |
| [`legacy/v3/`](../legacy/v3/) | **Historical** V3, including frozen code, figures and provenance. | Read-only for current V4 science; may be referenced for methods/history. |
| [`legacy/v2/`](../legacy/v2/) and [`legacy/ch1_global/`](../legacy/ch1_global/) | **Historical** Chapter 1/V2 work. | Do not infer V4 effects from these artifacts. |
| [`legacy/v4_s2_precipitation_niche_20261025/`](../legacy/v4_s2_precipitation_niche_20261025/) | **Deferred** S2 climate-niche exploration, not current GEB. | Leave original evidence intact. |
| [`analysis/ch1/`](../analysis/ch1/), [`src/azami_ch1/`](../src/azami_ch1/), [`analysis_outputs/`](../analysis_outputs/) | **Older/shared** chapter code and outputs. | Retain unless all imports, tests and audit histories are remapped. |
| [`docs/`](../docs/) | **Navigation and maintenance**. | Updating an index does not change scientific evidence. |

## Keep analysis units and test families distinct

| Result | Source and eligibility | Correct inference |
| --- | --- | --- |
| 3.1 | Original S1 native 80-test family; **168,907 angle-eligible records / 178 taxa** | BIO18–image-frame angle among taxa, not a gravity-based head-angle experiment. |
| 3.2 | Original S2 128-test family; **two species × two range strata** with distinct colour-eligible records | Three introduced-stratum colour associations, not tests of slopes *between* native/introduced ranges. |
| 3.3 | Source S1 photographic records filtered to **63 geographically eligible taxa**, **7,117 taxon-by-1° cells**; 153,843 contributing images prior to aggregation | Greater joint phenotype resemblance at matched climates versus within-taxon vector reassignment; **separate exploratory Holm-2 family**. |
| Auxiliary 8-construct cohort | **156,683 photos / 174 taxa**, 28 descriptive construct pairs | Covariance and observed configurations; supporting context, not a distinct ecological result 3.1 or a new data source. |

H3c 1° effect direction is sensitive to trait selection and geographic coverage; 5° support is more stable. Original S1/S2 model freeze: 2026-10-24, while H3c entered the GEB narrative 2026-10-29. Figures 1–5 and SI show the distinction.

## Workflows: what they do

- [GEB manuscript/source contract](../.github/workflows/v4-geb-main-contract.yml): verifies manuscript and immutable aggregate outputs.
- [Frozen V4 claim check](../.github/workflows/v4-final-bio18-holm-freeze.yml): protects the original S1/S2 estimates.
- [S1 geographical replication](../.github/workflows/v4-s1-bio18-geographic-replication.yml) and [S1 H3 exploratory runner](../.github/workflows/v4-s1-climate-trait-hypervolume.yml): separately tagged validation/exploration.
- [S2 colour spatial-block sensitivity](../.github/workflows/v4-s2-colour-spatial-block-sensitivity.yml) and [source-verified sampling coverage](../.github/workflows/v4-submission-sampling-coverage.yml): may require protected source credentials and archived assets, with public aggregate-only outputs.

**Do not treat a successful document check as detector accuracy or a proof of adaptive convergence.** No old workflow or draft PR is made authoritative merely by keeping it archived.

## Change procedure

Update the current paper and its claim ledger together **only when evidence supports a scientific change**; otherwise put navigation changes in README/docs without changing numerical code or original receipts. Review [the runbook](REPRODUCIBILITY_AND_SUBMISSION_RUNBOOK.md) and [branch inventory](BRANCH_ARCHIVE_TRIAGE_20261030.md). Avoid renaming/moving paths referenced by checksums, Scripts, old PRs, or archived CI results.
