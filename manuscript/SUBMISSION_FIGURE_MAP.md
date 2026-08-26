# Submission figure map

This concise map mirrors the authoritative `FIGURE_TABLE_MAP.md`. The current
main manuscript contains five figures. Historical-placement sensitivity is in
Appendix S1, and the failed involucral rows are not plotted as supported effects.

| Item | Content | Current source/status |
|---|---|---|
| Figure 1 | Open-licensed photographs, detector boxes and deterministic continuous measurements | frozen provenance; not an independent accuracy test |
| Figure 2 | Geographic sampling, cohort flow and BIO1 × BIO12 domain | frozen cohort tables; Mollweide equal-area maps |
| Figure 3 | Nested taxon → photograph → head visible-variance decomposition | frozen atlas plus one-head and balanced-replication sensitivities |
| Figure 4 | Taxon-level PC1–PC3 trait architecture | 148 complete taxa; frozen PCA |
| Figure 5 | Frozen primary coefficients, grouped-SPDE context and native-only decisions after seasonal/dominant-taxon audits | local regeneration and visual recheck required; immutable release pending |
| Figure S1.8 | Among-taxon environmental separation geometry | frozen permutation design; presentation only |
| Figure S1.9 | Randomized-tree PGLS placement sensitivity | demoted from main; not resolved phylogenetic correction |

## Excluded or withdrawn displays

- Raw absolute-slope RMS lability plots and median-split quadrants remain withdrawn.
- The three unadjusted involucral BIO4 rows failed the locked resolution/sharpness
  audit and are preserved in Table S1.6 rather than shown as supported effects.
- The former main Figure 6 is now Figure S1.9 because direct-backbone coverage is
  54/216, retained Pagel-λ estimates are zero and direct-backbone non-circular
  phylogenetic signal is unsupported.

Figure 5 must be regenerated through
`.github/workflows/ch1-interpretive-figures-ci.yml` after the current changes.
