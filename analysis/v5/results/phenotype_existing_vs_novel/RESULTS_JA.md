# 既存変異 vs 在来未観測変異：最終横断比較 (2026-11-26)

[Run 45773625808](https://github.com/zuizui0223/azami/actions/runs/45773625808) successful. Commit e235e0877955b0f678fa583e356364cd0894b632; artifact 15095455814; downloaded ZIP SHA-256 0ce6381b896e8580944b6c0be678d736b358a3a182728af9198ec1722d0b433a matches published digest. All 80 planned tests estimated and all original matched pairs retained (C. arvense 30/9 blocks and 183/14; C. vulgare 146/14 and 316/15). Independent recalculation of Holm all80 and family20 exactly matches. Partition identity existing_gain + new_gain = TV verified to 1.12e-16; native-to-novel new_gain equals novel_unobserved exactly.

## Direct decomposition, no temporal or genetic modeling
For fixed native-reference phenotype bins j, calculate positive frequency gain g_j=max(0,p_novel,j-p_analogue,j). Split sum(g_j)=TV into E (bins observed in native) and U (bins unobserved in native). Report E/(E+U). Also report the corresponding native-to-novel decomposition, and unobserved-bin mass separately (not interchangeable with U). Primary 95% native-core comparison; density surplus sensitivity. Eight endpoints + two joint spaces (C*, L*; C*, L*, angle), equal-cell-weighted observations. Full native reference (18,505/13,115 joint measurements; hue 18,503/13,106). Rarefaction: 32 disjoint geographic splits with 256 reference cells, one archived observation per cell, and 128 native calibration cells (one observation per cell); same draw used for both introduced groups, introduced profiles retain within-cell observations. No independent replicates claimed.

| Native reference | Climate definition | Existing share E/(E+U), range of 20 specifications |
|---|---|---:|
| full | native core 95% | 99.9104–100% |
| full | density surplus | 97.0383–100% |
| rarefied | native core 95% | 85.2170–100% |
| rarefied | density surplus | 82.7109–99.6847% |

Existing-bin gain is >50% in all 80 specifications (40 full, 40 rarefied); no inferential test of this ratio was performed. Direct native-to-novel positive gain share instead ranges 98.1816–100% (full) and 58.9221–100% (rarefied). In all full-reference cases, five interior marginal bins, all 12 hue sectors, all 25 interior colour combinations and all 125 interior colour-angle combinations were occupied in both species. Thus full-reference near-100% shares partly follow from bin saturation; explicit outside-extrema bins remain evaluable. Finite native reference size sensitivity is essential. In full reference max novel unobserved mass .110666%; in rarefied reference max scalar unobserved mass 1.43554%, joint max 12.9788%. More axes produce more sparse bins; this is not evidence for more genetic novelty.

## Primary rarefied reference, main comparison

All entries percentage points except last column.

| Species | Space | Native held-out unobserved | Analogue unobserved | Novel unobserved | Existing gain | Unobserved gain | Existing share (%) |
|---|---|---:|---:|---:|---:|---:|---:|
| C. arvense | chroma | .6836 | .8504 | .0231 | 14.2096 | .0153 | 99.8926 |
| C. arvense | colour | 1.1475 | .9550 | .0231 | 26.7506 | .0220 | 99.9178 |
| C. arvense | colour-angle | 13.4521 | 15.0626 | 12.8912 | 49.2927 | 8.2602 | 85.6476 |
| C. vulgare | chroma | .8545 | 1.3227 | .9458 | 11.3400 | .1097 | 99.0419 |
| C. vulgare | colour | 1.3916 | 2.0162 | 1.7097 | 17.0448 | .6333 | 96.4176 |
| C. vulgare | colour-angle | 13.3057 | 14.4486 | 12.3501 | 31.0209 | 5.3814 | 85.2170 |

In 10/40 rarefied cases, novel unobserved mass exceeds native held-out benchmark descriptively. No formal test of that difference, no equivalence claim. In 12/40 rarefied and 6/40 full cases the novel unobserved fraction exceeds analogue. Raw p<.05 in only two cases, both *decreases* in C. vulgare (rarefied 95%, solidity p=.01953125, colour-angle p=.02111816). All80 Holm p=1 throughout; family20 min=.390625. No corrected evidence for increased novel-bin mass. This test addresses signed unobserved fraction, not significance of total TV or E/U, not an equivalence test. The earlier 8-endpoint omnibus tests also did not pass correction.

## Final interpretation and scope
The observed shifts are dominated by frequency redistribution within bins represented in the contemporary native phenotype sample, even under equal-size native reference sensitivity. No robust evidence of an increase in native-unobserved phenotype-bin occupancy in novel climates. Cannot infer ancestral alleles, new mutations, fitness advantages, selection, or post-introduction evolution; a new allele can produce a phenotype inside the native envelope and a standing allele can produce a phenotype outside it. Coarse bins, complete-case photo sampling, unequal numbers of observations within cells, region dependence and image measurement error remain limitations. Earlier pooled rarity-chroma association attenuates with regional centering; candidate increases of high-chroma/low-lightness are exploratory. No new temporal, genetic, source-mixture or climate prediction model; user requested the simple cross-sectional comparison. This closes the direct contrast analysis; use these results with the preexisting evidence ledger and do not search new specifications for significance. Main and poster unchanged.
