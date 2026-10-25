# Azami — V4 S1–S2 | Global Ecology and Biogeography

**Current main scientific baseline (frozen 24 October 2026): V4 S1–S2.** This repository preserves the existing V3 science under [`legacy/v3/`](legacy/v3/README.md) without deleting its original sources. V4 T2 (phylogenetic history) and T3 (other-tribe extension) are **not** part of the present GEB manuscript.

## Start here

1. [GEB Research Article draft](manuscripts/v4_s1_s2_frozen_20261025/GEB_RESEARCH_ARTICLE.md) — macroecological framing, structured abstract, primary Methods/Results/Discussion, tables and planned display pieces.
2. [Blinded reviewer-text working draft](manuscripts/v4_s1_s2_frozen_20261025/GEB_RESEARCH_ARTICLE_BLINDED.md) — identity-sensitive internal provenance hidden; **not yet submission ready** because anonymous reviewer access is not established.
3. [Source and claim ledger](manuscripts/v4_s1_s2_frozen_20261025/CLAIM_LEDGER_AND_SUBMISSION_GATES.md) and [GEB preparation checklist](manuscripts/v4_s1_s2_frozen_20261025/GEB_SUBMISSION_CHECKLIST.md).
4. [Frozen V4 scientific record](analysis/v4/FINAL_V4_BIO18_HOLM_LEGO_FREEZE_20261024_JA.md) and [read-only science checksum workflow](.github/workflows/v4-final-bio18-holm-freeze.yml).
5. [Historical V3 archive](legacy/v3/README.md) and [V2 archive](legacy/v2/ORIGINAL_REPRODUCTION.md).

## Frozen V4 research question

Are climatic associations and covariation among visible *Cirsium* capitulum **presentation, colour and outline** the same at between-taxon and within-species native/introduced scales, or is the visible floral phenotype partly modular and environmentally scale dependent?

**Primary four-variable model:** CHELSA BIO18 (warmest-quarter **mean monthly** precipitation), BIO05, RSDS, plus WorldClim elevation. This is a **post hoc frozen** model, not independently preregistered. S1 is the native-range source census (313,350 records); S2 is the native/introduced wild census of *C. arvense* and *C. vulgare* (297,622). Eight photograph-derived constructs are analysed. Independent Holm families contain **80 S1 and 128 S2** biological comparisons (hue sine/cosine jointly counted as one construct).

**Four Holm-supported observational associations:**
- S1 between taxa: image-frame head angle × BIO18, standardized β +0.045447, Holm P = 0.00920635 (n 168,907).
- S2 introduced *C. arvense*: image colour lightness × elevation, β −0.173732, Holm P = 0.00131233 (n 26,171).
- S2 introduced *C. vulgare*: chroma × RSDS, β +0.134504, Holm P = 0.01456813 (n 27,200).
- S2 introduced *C. vulgare*: lightness × RSDS, β −0.142801, Holm P = 0.01925414 (n 27,200).

S1's 156,683 common photographs / 174 taxa yielded larger between- than within-taxon RV for **24/28** environmentally adjusted pairs; higher same-domain than cross-domain RV is a measured **photographic phenotypic association**, not proof of genetic modules. In S2, the original four species×range strata retain three Holm-supported **colour** associations in introduced ranges, with **no head-angle term** passing the S2 Holm family. S1 native photographs overlap S2 native photographs and cannot be treated as independent replications.

**Current GEB scope exclusion:** no native-versus-introduced **precipitation-niche** or climate-occupancy comparison, no post-freeze 4D centroid/novelty or matched-climate range-effect regression. Those 25 October exploratory products were moved to [deferred historical research](legacy/v4_s2_precipitation_niche_20261025/README.md), outside the frozen primary narrative. See [scope policy](manuscripts/v4_s1_s2_frozen_20261025/SCOPE_LOCK_20261025_JA.md).

**Not established:** causal rainfall selection, true gravity-based nodding, calibrated anthocyanin, heritable LEGO units, evolutionary module reassembly or introduction-driven evolution. Previous BIO12 geographic analyses show sensitivity to taxon composition; photograph-based confounding remains a limitation. Full supplement and source checks are anchored to [frozen artifact 12547876704](https://github.com/zuizui0223/azami/actions/runs/40062126006), original-data success SHA `d2c45b00018a0c04e7c69a29857207dc1ac06d0c`, ZIP SHA-256 `6e081dc081bb8ee2311a05ac5b7fcddd3f6309be2fa909e87d432789f4e0e205`. The scientific snapshot is `7e25e409d39a7e8bb8a1ce1043d6bbfd41b7d159`.

## Layout and reproducibility boundaries

- `analysis/v4/`: **only frozen S1–S2 code, contracts and source-verified summaries**, not the T2/T3 research branches. For full, earlier private-source reconstruction use the source execution commit pinned above; today's `main` is a submission-oriented presentation of that immutable analysis, not a remeasurement.
- `manuscripts/v4_s1_s2_frozen_20261025/`: GEB article and figure set limited to the 24 October frozen S1–S2 **phenotype–environment** science, with a scope-locked claim ledger and checked numerical supplement.
- `legacy/v3/`: archived V3 analytical source, reproducibility package, result receipts, old CI definitions and dedicated tests.
- `legacy/v2/`, `legacy/ch1_global/v2/`, `analysis/ch1/`, `src/azami_ch1/`, `analysis_outputs/`: retained earlier/shared resources. They are **not** the V4 results.

Only a **read-only V4 audit and GEB manuscript integrity** workflow remain active. Historical V3 Zenodo publishing steps are archived and must not auto-run after this change.

## Public-data limitation and journal status

The frozen V4 numerical artifact is verified but GitHub Actions artifacts can expire and private upstream measurements remain access-restricted. A durable deposit of allowable aggregates/code, anonymous review URL, field/photo validation, licences, full references, author metadata and final figure layout are still necessary to submit to *Global Ecology and Biogeography*. No DOI for V4 or published/submitted manuscript is claimed. The previous Zenodo DOI 10.5281/zenodo.22295791 is a **historical V2 numerical archive**, not the V4 dataset.

Original code is MIT licensed; data, images and third-party assets remain subject to their separate licences; see [NOTICE](NOTICE.md).
