# V4 S1–S2 paper package (24 Oct 2026 source freeze)

**Status:** manuscript draft + real frozen statistical supplement; not submitted or independently photo-validated.  
**Scientific parent:** \`snapshot/v4-final-bio18-holm-20261024\` at \`7e25e409d39a7e8bb8a1ce1043d6bbfd41b7d159\`.  
**Analysis parent:** original execution commit \`d2c45b00018a0c04e7c69a29857207dc1ac06d0c\`.

## Where to start

1. [MAIN_MANUSCRIPT_DRAFT.md](MAIN_MANUSCRIPT_DRAFT.md): article narrative, Methods and complete primary Results, two principal tables, initial bibliography and figure plan.
2. [CLAIM_LEDGER_AND_SUBMISSION_GATES.md](CLAIM_LEDGER_AND_SUBMISSION_GATES.md): evidence hierarchy, SHA-checked source receipts, non-claims and remaining verification.
3. [table1_holm_associations.csv](table1_holm_associations.csv): **display-only** four associations after Holm, not the full test family.
4. [module_summary.csv](module_summary.csv) and [s2_module_summary.csv](s2_module_summary.csv): descriptive within/across-module RV medians; **not** a replacement for all 28 pairs.
5. [plot_frozen_summary.py](plot_frozen_summary.py) generates Figures 2–3 summary plots from the provided small CSVs, with locked-value assertions.

## Full numerical supplement regenerated from the original archive

The complete original aggregate archive is available in [Actions run 40062126006](https://github.com/zuizui0223/azami/actions/runs/40062126006), **artifact 12547876704**, name \`v4-s1-s2-extreme-lego-40062126006\`. Its ZIP must have SHA-256:

\`\`\`
6e081dc081bb8ee2311a05ac5b7fcddd3f6309be2fa909e87d432789f4e0e205
\`\`\`

Fetch the ZIP through GitHub Actions and run:

\`\`\`bash
python -m pip install numpy scipy statsmodels matplotlib
python manuscripts/v4_s1_s2_frozen_20261025/build_full_supplement.py \
  --artifact /path/to/v4-s1-s2-extreme-lego-40062126006.zip \
  --out /path/to/v4-s1-s2-si
\`\`\`

The source-checking code recreates:

| File | Contains |
|---|---|
| \`Table_S1_full_208_BIO18_tests.csv\` | **All** primary S1 (80) and S2 (128) biological tests, including joint circular-hue Wald tests, uncorrected and BH/Holm P, coefficients, available pointwise CIs, test group, exposure, scale and each test's n |
| \`Table_S2_all_28_RV_pairs.csv\` | All 28 exact common-photo trait pairs, raw and BIO18/BIO05/RSDS/elevation-adjusted within/among-taxon RV and available signed correlations |
| \`Figure_S1_S1_all_28_RV.png\` / \`.svg\` | Four comparable within/among × raw/adjusted 8×8 trait association matrices |
| \`SOURCE_SHA256.txt\` | Pin to checked frozen ZIP source and explicit provenance restrictions |

**Verified against the actual downloaded artifact, 25 October 2026:** SHA-256 matched byte-for-byte; \`statsmodels\` recomputation of all 80 S1 and 128 S2 BH/Holm adjusted P values differed from the original stored values by **0**; Holm retained exactly 1 S1 + 3 S2 associations. Re-extracted 28-pair RV recovered original count **21/28** between>within raw and **24/28** conditional. Pair-domain medians matched all four exact frozen values with 9 same-domain and 19 cross-domain pairs.

To avoid leaking private observation identifiers, precise coordinates or image rights, the generated SI files only include aggregate result rows and taxon-level labels. The script will stop on an unexpected ZIP SHA, cohort change, Holm discrepancy, or RV summary change.

## Scientific integrity checkpoints

- The S1 family includes **8 main constructs + 2 low-coverage bract image proxies** across two scales and four climate exposures = **80** biological tests. Bract proxies enter the multiplicity denominator but remain SI-only ecological interpretation.
- The S2 family includes eight constructs, four exposures and four independent species × range **fit strata**, not independent evolutionary replicate origins = **128** biological tests.
- Hue sine/cosine is a **joint two-degree-of-freedom single trait**, not two discoveries.
- The 4 exposure model was **selected after earlier explorations**. Within-model Holm correction is not correction for the prior search over climate specifications.
- Median same-module RV can reflect common-pixel and shared-mathematical geometry rather than genetically separable functional systems. The fixed 280-label result is a conditional label rank, not 280 evolutionary replicates.
- The photographic angle is to EXIF/image vertical; gravity-validated nodding, rainfall protection, anthocyanin mechanism and causal adaptation remain unsupported.
- S1 source-native observations and S2 native observations overlap; their numerical results should not be presented as independent validation.
- T2 nuclear phylogeny and T3 other-tribe extension remain separately investigated and are not required for this S1–S2 paper.

## Before journal submission

Validate source imaging/measurement, verify spatial/residual bias and taxonomic audit boundaries, finish full figure captions and image-license attribution, specify all model software versions and preprocessing decisions, verify bibliography against primary sources, decide authorship/CRediT and repository data access, then produce journal-specific manuscript and supporting-information PDFs. These are not claimed complete merely because the numerical freeze and first article narrative exist.
