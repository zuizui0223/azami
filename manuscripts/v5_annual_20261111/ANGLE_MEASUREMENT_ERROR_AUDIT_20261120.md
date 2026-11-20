# V5 angle measurement: suspect prevalence vs true error (20 Nov 2026)

## Scope and source verification

The original confidential records and archived measurement tables were restored in GitHub Actions, with only aggregate outputs disclosed. These are **completed successful analyses**, not imputed estimates:

- **S1 exact V5 annual five-climate-axis native joint-phenotype cohort:** [Actions run #44666968101](https://github.com/zuizui0223/azami/actions/runs/44666968101), 157,822 photo-derived **observation/measurement** records. Counts for 1° (153,843 photo-derived records / 63 taxa / 7,117 taxon–cells) and 5° (132,998 / 31 / 814) were independently asserted.
- **S2 exact focal V5 three-trait-complete native/introduced cohort:** [Actions run #44666062718](https://github.com/zuizui0223/azami/actions/runs/44666062718), 79,421 photo-derived records in two species. Source, WCVP range status and CHELSA five annual variables are verified.

Both use the archived angle and rotated projected outline aspect variables. Angle refers to **EXIF-oriented image up**, not gravity, and is computed from the detected head PCA principal axis signed toward the floral-pixel centroid. The existing archived model retains only heads with `orientation_status=usable`, with a minimum automatic confidence gate and left-right flip agreement (see `legacy/ch1_global/v2/52_measure_primary_traits_continuous.py`, `56_run_primary_traits_continuous_v2.py`).

## Source-reported prevalence

| Flag at photo-derived observation grain | S1 all 157,822 | S1 1° qualified 153,843 | S1 5° qualified 132,998 | S2 79,421 |
|---|---:|---:|---:|---:|
| Rotated outline aspect within factor 1.10 of 1, `abs(log(A)) <= log(1.1)` | 96,677 (61.26%) | 94,159 (61.20%) | 81,742 (61.46%) | 50,254 (63.28%) |
| Signed image angle greater than 90° | 48,975 (31.03%) | 47,609 (30.95%) | 42,039 (31.61%) | 25,167 (31.69%) |
| Both preceding flags | 29,851 (18.91%) | 29,010 (18.86%) | 25,565 (19.22%) | 15,702 (19.77%) |
| Image angle >=110° | 29,146 (18.47%) | 28,335 (18.42%) | 25,300 (19.02%) | 15,143 (19.07%) |
| Very near-unity outline aspect within factor 1.05 | 57,016 (36.13%) | 55,510 (36.08%) | 48,240 (36.27%) | 29,908 (37.66%) |

**Do not call these error rates.** The outline aspect is **height/width of the rotated bounding-box silhouette**, not the eigenvalue ratio which controls PCA-axis stability. Values close to 1 by themselves cannot diagnose angle mismeasurement. An angle >90° represents a signed image convention, not a verified defect.

### Near-unity outline is not strongly enriched for angle >90°

- S1 full V5: fraction of signed image angle >90° is **30.88%** for near-unity (factor 1.10) outline vs **31.28%** outside that class.
- S1 5° cohort: **31.28%** vs **32.14%**, respectively.
- S2 two-species source: **31.25%** vs **32.45%**, respectively.

These are descriptive percentages on non-independent biological records, with no inferential statement. They do **not** validate angle measurement accuracy; they show that this proxy alone does not concentrate images with angles >90°.

### Rain gradient check on source-record composition (S1 full)

All five annual BIO12 source-record quintiles contained approximately 31,565 records. For the joint aspect-near-one/angle>90 flag, the fractions (driest to wettest) were **18.36%, 19.10%, 18.82%, 18.54%, 19.75%** (rainfall median 514, 667, 804, 1024, 1471 mm). No large monotonic change in this *proxy* is evident. However, these raw-record summaries are not the original 1°/5° taxon-paired regression units; they do not establish that photography error is independent of climate, space or taxonomic structure.

## What remains unmeasured

These audits do **not** establish: (i) true physical angle or gravity alignment; (ii) human-annotated agreement to PCA-derived head angle; (iii) error distributions by head-axis **eigenvalue anisotropy**; (iv) image camera roll or EXIF reference validity beyond available metadata; (v) whether stricter actual head orientation confidence thresholds change V5 Fig. 2/3.

To estimate an actual error fraction, one needs an **independently annotated, stratified sample of source photographs and head-specific landmarks**. Report mean absolute angular deviation, quantiles, a 90°-orientation confusion matrix and uncertainty by photographed shape, observer, region and climate. In addition, run the exact frozen multivariate climate-analogue and BIO12/angle models on quality-restricted cohorts, holding the statistical comparison families, taxonomic samples and null methodology constant when possible. The original V5 numbers are neither retracted nor validated by the prevalence audit alone.

## Source scripts and results

- `analysis/v5/audit_angle_suspects.py`
- `analysis/v5/run_angle_suspects_s1_recovered.py`
- `analysis/v5/run_angle_suspects_s2.py`
- `.github/workflows/v5-angle-suspects-s1-recovered.yml`
- `.github/workflows/v5-angle-suspects-s2.yml`

No individual image identifiers, observer identities or source coordinates appear in any public audit.
