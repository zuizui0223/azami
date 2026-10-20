# v4: original observer effort and dominant-taxon sensitivity (2026-10-20)

**Scope:** exploratory, native-first WCVP/TDWG3, original observations, historical photo-derived traits only. **No v3/GEB scientific outputs modified.** The active branch is PR #142 and must stay Draft.

## Why this test differs from a photograph-quality audit

The concern is **human observation and posting effort**: recording more of a few widespread species, repeated observations from prolific accounts, and many uploads by the same account in the same region. Photograph technical eligibility is an independent issue and is **not** the focus of this test. Counts are of **observations**, never photographs-as-independent-plants.

The original pre-detector archived snapshot contains **665,115 observations / 1,122,854 photo rows**; exact WCVP/TDWG3 native-first matching admits **313,350 observations / 306 accepted taxa** and **92,839 original iNaturalist user IDs**. The original metadata ZIP is checksum pinned. All user IDs, observation IDs, exact coordinates and per-user lists are kept offline, **not committed or attached to public CI**.

### Observation intensity / effective representation

| Source metric | Estimate |
|---|---:|
| C. arvense original native observations / unique observers | **78,646 / 26,644** |
| C. vulgare original native observations / unique observers | **62,368 / 22,860** |
| Top two native species' share of **all** observations | **45.002%** |
| Top five species' share | **63.145%** |
| Top ten species' share | **78.110%** |
| Simpson effective species number from *posting share* (not actual biodiversity) | **8.31** of 306 |
| Top 1% of original user IDs | 928 users, **29.661%** of observations |
| Top 10% of original user IDs | **59.950%** of observations |
| Largest single observer | **0.416%** of observations |
| Records inside repeated taxon × 2-degree equal-area cell × observer groups | **60.818%** |

A highly concentrated species histogram **does not itself prove** model distortion. Models below check how much estimates actually change.

### Four weighting designs

1. **Observation pooled:** each observation gets its historical recording share; species and users can dominate.
2. **Equal taxon:** each accepted taxon contributes equal total weight.
3. **Equal taxon and geography:** within each taxon, its occupied approximate 2-degree equal-area cells get equal total weight.
4. **Equal taxon × geography × observer:** within a taxon-cell, each distinct source observer gets equal total weight, with that observer's repeated records split evenly.

For among-taxon regressions, each species is one final unit **even under scheme 1**. Hence observation vs equal-taxon weighting makes no difference there unless the **within-taxon mean or regional balance** changes. These are different estimands from naive pooled observation regression; never compare them as one common causal effect.

## Results — solar radiation and visible chroma

Primary species-level marginal coefficients: complete nine-variable native CHELSA, at least five measured observations per accepted taxon, 177 taxa / 186,313 observation records.

| Calculation | Standardized beta |
|---|---:|
| Naive observation-pooled marginal, all taxa | **−0.0512** |
| Same naive pooled, omit 2 highest-volume species (predefined in original source) | **−0.2275** |
| Same naive pooled, omit 5 highest-volume species | **−0.2973** |
| Observation-pooled with **equal source-taxon** weights | **−0.1869** |
| Taxon-level, source-observation means (taxa equal final weight) | **−0.3466** |
| Taxon-level, **equal cell** within taxon | **−0.3361** |
| Taxon-level, **equal cell × observer** within taxon | **−0.3472** |
| Taxon-level observer-balanced, omit 2 highest-volume species | **−0.3501** |
| Taxon-level observer-balanced, omit 5 highest-volume species | **−0.3552** |

All 177 taxon-level leave-one-taxon-out observer-balanced marginal results range **−0.3743 to −0.3316** (maximum absolute slope change **0.0271**). The named top-two species are not the source of the negative species-level association. Instead, naive **observation-count pooling dilutes that pattern**.

The pooled coefficients above are from response/predictor-only complete cases, so their original n may be two records larger than all-nine-predictor analyses. The comparison intentionally uses different weight estimands; it is an observation-bias diagnostic, not a biological gradient's magnitude comparison.

When nine environments and the eight-column spherical space basis are fitted jointly, the observer-balanced *among-taxon* solar coefficient is **−0.2422**, 95% HC3 interval **[−0.8476,+0.3632]**. It is **not separately identifiable** from other geographically collinear predictors. Within-taxon observer-balanced solar beta is **−0.0214** (marginal; taxon-cluster 95% interval **[−0.0558,+0.0130]**). No within-taxon colour plasticity conclusion.

Observer omission on fixed remaining species support: omit the 928 most prolific source observers, with exactly the same eligible taxa retained on both sides; chroma species-level marginal beta moves **−0.3692 to −0.3970** (161 taxa). This changes magnitude, not sign; if observers are removed without fixing taxa, some low-support species disappear and results cannot be attributed to observer effect alone.

## Results — annual rain and image-referenced presentation

| Calculation | Standardized beta |
|---|---:|
| Naive observation-pooled marginal | **+0.0294** |
| Naive pooled, omit top 2 species | **+0.0757** |
| Naive pooled, omit top 5 species | **+0.0938** |
| Pooled equal source-taxon weights | **+0.0945** |
| Species-level marginal, equal taxon + cells | **+0.2708** |
| Species-level marginal, **equal taxon × cell × observer** | **+0.2256** |
| Observer-balanced, top two taxa omitted | **+0.2248** |
| Observer-balanced, top five taxa omitted | **+0.2255** |

All 178 taxon-level leave-one-out observer-balanced estimates range **+0.2070 to +0.2534**. On identical species support, removal of the most prolific 1% of observers gives **+0.2697 to +0.2775** across 162 retained taxa.

Critically, the **within-taxon** marginal result is **−0.0081** if observations are all counted equally, but **+0.0435** if taxon, region and original observer are balanced. Its taxon-cluster 95% interval is **[−0.0233,+0.1102]**, including zero. **Do not call within-taxon presentation-rain effects supported.**

Observer-balanced among-taxon joint nine-environment plus spherical-space coefficient is **−0.2472 [−0.5964,+0.1020]**. The marginal positive association is geographically entangled and not an independently identified rainfall mechanism.

## Within/among phenotype covariance under observer balancing

Exactly the same complete phenotype cohort, original archived measurement per observation (not 5 photos considered independent).

| Observer/cell weighting | Taxa | Within RV | Among RV | Matched-n taxon bootstrap median delta | Bootstrap 95% interval |
|---|---:|---:|---:|---:|---:|
| Equal taxon × cell | 177 | 0.0062 | 0.0368 | +0.0252 | [−0.0204,+0.0933] |
| Equal taxon × cell × observer | 177 | 0.0068 | 0.0371 | +0.0246 | [−0.0208,+0.0834] |
| Observer balanced; drop top two species | 175 | 0.0066 | 0.0371 | +0.0237 | [−0.0241,+0.0858] |
| Observer balanced; drop top five species | 172 | 0.0066 | 0.0372 | +0.0258 | [−0.0181,+0.0851] |

All observer-balanced 400-replicate intervals include zero, so **a general interspecific strengthening of integration is not demonstrated** in this native-wide test. RV itself is **not a percentage of variance explained**, and the re-sampling proportions are not p-values.

## Results retained vs not identified

**Supported descriptive observation-bias finding:** Recording effort is strongly concentrated in very few source taxa and, to a lesser extent, prolific observers; unbalanced pooled regressions yield radically different slope magnitudes, while equal-taxon species-level marginal radiation–chroma and precipitation–presentation signs persist after top-species and observer controls.

**Not established:** unique/causal solar or rain effects after joint climate/space control; within-taxon plasticity; selective pressure on pigments, rain protection or floral anatomy; genetic covariance; and observational correction for flower images never posted. This analysis includes *existing* photograph phenotypes only; **96,246 eligible native observations still lack legacy trait measurements**.

The observer IDs and exact locations remain PRIVATE; public code and aggregate receipt only. This post hoc sensitivity should not be presented as confirmatory selection-proof evidence.

### Source identity and reproduction boundary

- Original photo metadata SHA256: 2f216652602a0b3c42b64780d9b31a3ecb2e7fc1300b28b303e8c369b31884ae
- Native350 roster SHA256: c4cc1f59afe93dde7c1ad19d7f432fdf79c110a0334b1b137bb198cd4328371c
- Observer private join SHA256: 2cba1a90cd33e4e79285406e4c02a4c85e14739bb04bef55cbdaebed24c7988b
- Historical measured per-observation phenotype SHA256: b49a093c849928c456d6747ec51772b5a6b78debb2ae5eb5dbe2a91b08f66dd7
- Full native climate cache-derived table SHA256: 410e1ae7f71a6283d3bfb57c322aeda698195317fda6dde5391c8ff8d7b400a1

The full-original-source CI that depends on Actions artifact 8066010557 returns **HTTP 410 expired** and must not be called successful. It is now explicitly manual; ordinary PR tests can validate code without obtaining private source materials.
