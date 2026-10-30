# GEB V4 runbook — safe checks and journal hand-off

**Working manuscript:** [current GEB Research Article](../manuscripts/v4_s1_s2_frozen_20261025/GEB_RESEARCH_ARTICLE.md). This runbook distinguishes **checking committed aggregate results** from **rerunning private image-level measurements**.

## 1. Verify the source and the paper locally

Run from the repository root in a Python 3.12 environment:

```bash
python -m unittest discover -s tests -p 'test_repository_navigation.py' -v
python -m unittest discover -s tests -p 'test_v4_geb_main_contract.py' -v
python -m unittest discover -s tests -p 'test_v4_review_aggregate_bundle.py' -v
python -m unittest discover -s tests -p 'test_v4_submission_sampling_coverage.py' -v
```

The first three commands use repository files/aggregates; the sampling-coverage tests use synthetic fixture data. Some other scientific tests or workflow jobs require separately installed dependencies and verified source artifacts. Passing these checks establishes **document, code, and frozen-aggregate consistency**, not independent measurement accuracy or a public raw-data rerun.

## 2. Recheck the original S1/S2 statistical freeze

The original numerical source is the GitHub Actions artifact **12547876704** from [run 40062126006](https://github.com/zuizui0223/azami/actions/runs/40062126006). **Expected ZIP SHA-256:** `6e081dc081bb8ee2311a05ac5b7fcddd3f6309be2fa909e87d432789f4e0e205`. If this artifact has expired, *stop*: retrieve the authenticated pinned source from the authorized repository maintainer. Do not substitute arbitrary newer files.

With `gh` authenticated to authorized repository resources, and dependencies from the [GEB verification workflow](../.github/workflows/v4-geb-main-contract.yml) installed:

```bash
mkdir -p work/geb_v4
gh api repos/zuizui0223/azami/actions/artifacts/12547876704/zip > work/geb_v4/source.zip
echo '6e081dc081bb8ee2311a05ac5b7fcddd3f6309be2fa909e87d432789f4e0e205  work/geb_v4/source.zip' | sha256sum -c -
python -m analysis.v4.final_v4_bio18_holm_freeze \
  --artifact-zip work/geb_v4/source.zip --output work/geb_v4/freeze_check.json
python manuscripts/v4_s1_s2_frozen_20261025/build_full_supplement.py \
  --artifact work/geb_v4/source.zip --out work/geb_v4/si
```

The frozen families are **80 S1 + 128 S2** original tests; H3c has a **separate exploratory Holm-2** test family. Never add sample counts across overlapping cohort definitions.

## 3. Reproduce or inspect displayed figures

| Figure | Script or provenance | Important boundary |
| --- | --- | --- |
| 1 | [Photo and trait measurement pipeline](../manuscripts/v4_s1_s2_frozen_20261025/plot_measurement_pipeline_v4.py) · [credit/rights receipt](../manuscripts/v4_s1_s2_frozen_20261025/FIGURE1_PHOTO_CREDITS_20261030.json) | Example photography demonstrates constructs, **not** independently validated detector performance. |
| 2 | [Sampling coverage renderer](../manuscripts/v4_s1_s2_frozen_20261025/plot_sampling_coverage_v4.py) · [coarse-source JSON](../analysis/v4/results/20261030_SAMPLING_COVERAGE_SOURCE_VERIFIED.json) · [render manifest](../manuscripts/v4_s1_s2_frozen_20261025/figure2_render_manifest.json) | Only pooled 5° cell summaries are public; do not publish or attempt to reverse protected record coordinates. |
| 3 | [Original S1 BIO18-angle renderer](../manuscripts/v4_s1_s2_frozen_20261025/plot_bio18_angle_part2.py) | Main Figure 3 is the original Holm-supported among-taxon coefficient. |
| 4 | [Native/introduced S2 colour renderer](../manuscripts/v4_s1_s2_frozen_20261025/plot_s2_native_introduced_colour.py) | Uses six original coefficients (three native counterparts); **not** an interaction test. |
| 5 | [H3 climate-matched multivariate figure](../manuscripts/v4_s1_s2_frozen_20261025/plot_h3_climate_morphospace.py) | Primary finding is the joint distance effect; breadth, KDE and sensitivities are SI. |

Use the scripts' documented CLI options and the source hash checks, rather than rerunning them on unspecified inputs. The existing [figures directory](../manuscripts/v4_s1_s2_frozen_20261025/figures/) contains source-backed SVG/PDF/PNG assets. Do not treat differently named historical figures as replacements.

## 4. Build reviewer copies and aggregate verification materials

The repository provides [`build_submission_documents.py`](../manuscripts/v4_s1_s2_frozen_20261025/build_submission_documents.py) and [`build_review_aggregate_bundle.py`](../manuscripts/v4_s1_s2_frozen_20261025/build_review_aggregate_bundle.py). Consult `--help`, run them with explicit input and output directories, and compare the resulting figure captions, sample sizes, and supplements to the [submission checklist](../manuscripts/v4_s1_s2_frozen_20261025/GEB_SUBMISSION_CHECKLIST.md).

Generated Word/PDF/ZIP review packages live **outside Git** by design. Existing aggregate review candidates do not include private original photographs or a complete repeatable raw-data pipeline.

## 5. Release gates not yet fulfilled

1. Confirm author order, corresponding-author metadata, funding, CRediT, acknowledgements, journal declarations, cover letter and a separate title page.
2. Confirm that every externally reproduced photograph has a compatible licence and proper credit; Figure 1 credits are recorded, but third-party data are **not** automatically covered by the repository's MIT code licence.
3. Publish the **rights-cleared** V4 aggregate dataset and code in a durable archive, record DOI and version/commit; create a **reviewer-appropriate anonymous** access route in accordance with journal rules.
4. Recheck the Word/PDF line numbering, complete figure legends, editable tables, anonymisation, link accessibility and final counts before manuscript submission.
5. Retain the scientific limits (camera-frame angle, uncalibrated colour, possible geographic composition effects, H3c 1° sensitivity). None of the procedural gates are a license to claim demonstrated causal rain protection, convergent evolution, or detector accuracy.

## 6. What never to do automatically

- Never expose tokens, encryption secrets, unlicensed full images, precise source coordinates or uploader IDs through CI artifacts or a public release.
- Never replace the frozen 24 October S1/S2 estimates without a separately documented new analysis and claim revision.
- Never delete archive paths, remote branches or GitHub Actions artifacts solely because a new README no longer links to them.
- Never publish to Zenodo, merge later scientific PRs, or submit to GEB as a side effect of repository maintenance. PR #159 was explicitly merged on 2026-10-30; PR #160 remains separate.

## 7. Zenodo V4 dataset handoff (not yet published)

Use the [V4 Zenodo data-only handoff](ZENODO_V4_RELEASE_HANDOFF.md), [metadata draft](../manuscripts/v4_s1_s2_frozen_20261025/reproducibility_v4/ZENODO_V4_METADATA_WORKING_DRAFT.json), and [read-only validation script](../scripts/validate_zenodo_v4_release.py). Data rights, creators, reviewer anonymity, exact release bytes and explicit approval are release gates. The earlier V2 DOI and the GitHub V3 draft release do **not** identify this dataset.
