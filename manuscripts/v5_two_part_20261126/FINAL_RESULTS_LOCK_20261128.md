# 最終採用結果・凍結記録 / Final scientific result lock — 2026-11-28

User approved the final result set after completion of the standard R hypervolume analysis. This freezes the interpretation and display for the conference poster and manuscript, not submission acceptance or data release. No additional outcome-driven model or threshold search; future genuinely independent validation is a separate study.

## Main question
How does climate structure native floral trait diversity, and how is native-observed variation represented in novel climates of introduced ranges?

## Adopted results

| Result | Adopted finding | Evidence and limit |
|---|---|---|
| 1 | Similar climates, similar capitula | Within-taxon support at 1° and 5° (Holm .004); between-taxon 5° support (.012), 1° not significant (.076). Not independent evolutionary convergence. |
| 2 | Rainfall–downward orientation; precipitation seasonality–chroma | Positive two-grain conditional associations; exploratory pointwise taxon-bootstrap intervals. Not measured function or fitness. |
| 3-1 | Novel climate occupancy and broad shared head-trait space coexist | Standard R raw 5D climate ratio (introduced/native) 0.912 and 2.158 for Cirsium arvense and Cirsium vulgare. Background-corrected coarse-grid ratios 0.614 and 1.636 are separate approximations. Standard 3D native/novel-climate trait intersection as fraction of novel volume 98.3% and 81.5%; equal-n sensitivity 72.3–88.2% and 63.5–75.7% (not CIs). Calibrated novelty support separately 96.2% and 92.4% of climate-novel cell-draws; CV regional rank rho .149545, Holm .018; CA rho .062913, Holm .658. Two-trait CV sensitivity Holm .012. No established species difference or robust trait volume expansion/contraction. |
| 3-2 | Context-dependent nonlinear chroma redistribution in Cirsium vulgare | Climate-matched native/introduced comparison: quadratic +.141659, Holm12 .00284. Introduced novel/analogous comparison: quadratic −.165495, Holm24 .02636719, positive fitted tilt not individually significant after correction. Alternative novelty definition does not reproduce negative curvature. Ratios normalized to native median, not fitness. |

Poster uses existing native Part 1 panels; Result 3-1 retains background-corrected climate PCA contours plus R hypervolume 3.1.7 three-dimensional set points from full native vs novel-climate cells; Result 3-2 shows only supported CV continuous curves with pointwise intervals and alternative-definition caveat. Full species names, not CA/CV abbreviations, in displayed panels. Original analyses and full families preserved in SI.

## Discussion only
Trait diversity may reflect long-term environmental responses and provide standing variation for occupation of new climates. Broad trait-space overlap and context-dependent enrichment are compatible with existing variation as a resource; natural selection on standing genetic variation, plasticity and introduction history remain hypotheses. No inference of causal adaptation, fitness optima, evolutionary rate, adaptive radiation, or distinct selection regimes. Do not claim directional selection from the non-significant corrected linear term.

## Supporting and negative results retained
Regional cancellation, variance components, failed reconstruction/prediction, non-significant direct species contrasts, annual climate-slope coefficients (including native Cirsium vulgare negative precipitation-seasonality/chroma association), rarefaction novelty tests, all 108 R comparisons and all 16 grid geometric tests. No selective deletion. See V5_SUPPORTING_INFORMATION_MAP.md, V6_MANUSCRIPT_EN.md, and the audited records.

## Frozen outputs
- V7_MANUSCRIPT_EN.md: revised working article, not submission-ready.
- Local final conference poster: azami_sssb_cho_FINAL.pptx (SHA256 b9103f9899ea092ca38514dac0fb19f5ad99d0b56acbab81c2056b32d4a65622); previous R poster's 181 unchanged shapes, 43 unchanged media; only top JP/EN conclusion and third discussion text changed. Poster binary remains in conversation, not in this repository.
- [Standard R results](../../analysis/v5/results/PART2_R_HYPERVOLUME_RESULTS_20261128.md), runs 46106255971 and 46109779446.
- [Novelty and grid tests](../../analysis/v5/results/PART2_HYPERVOLUME_RESULTS_20261128.md), runs 46086510362 and 46097936236.
- [Continuous enrichment](../../analysis/v5/results/phenotype_continuous_enrichment/RESULTS_JA.md), run 45953204013; climate-matched range run 45972921869.

Main branch and source datasets unchanged. Submission formatting and full SI remain later work, not new analysis required for this lock.
