# GEB v2 final artifact results — 2026-08-26

## Immutable execution provenance

- GitHub Actions run: `32975451732`
- Head SHA: `f4a6fd5e01a2befd4f49174984a99e53856c2330`
- Artifact ID: `9612943217`
- Artifact name: `ch1-continuous-trait-universe-geb_v2_continuous_trait_universe-32975451732`
- Artifact digest: `sha256:101e996b638996a0c5ae79d358bf51293c3585f0e84c4a961b91dcbedf96211e`

This file records the submission-facing interpretation of the completed artifact-backed GEB v2 run. It does not promote candidate image proxies beyond their preregistered evidence tier.

## Endpoint coverage

The frozen continuous-trait contract contains 27 endpoints across eight modules. The final run measured 22 endpoints. Of 19 inferential endpoints, 18 had observation-level measurements: all 9 established primary endpoints and 9 candidate endpoints. The only unexecuted inferential endpoint was `visible_floret_fraction`.

The five unexecuted endpoints were:

- `visible_floret_fraction` — inferential candidate;
- `corolla_purple_pixel_fraction` — descriptive colour composition;
- `corolla_redmagenta_pixel_fraction` — descriptive colour composition;
- `corolla_white_pixel_fraction` — descriptive colour composition;
- `corolla_yellow_pixel_fraction` — descriptive colour composition.

The continuous trait universe contains 46,276 observations from 259 source-assigned taxa and 374,255 analysis-eligible observation × endpoint measurements. Missing or unexecuted measurements remain missing and are never converted to categories or biological zeros.

## Multivariate phenotype geography

The all-inferential taxon morphospace used 127 taxa complete for 18 measured inferential endpoints. Explained variance was 18.49% for PC1, 12.01% for PC2 and 11.78% for PC3 (42.28% cumulative). Thus the expanded phenotype does not collapse onto a single dominant capitulum syndrome; several quantitative axes are needed to describe among-taxon phenotype structure.

Module-level PCA likewise recovered multiple axes. Involucre architecture alone required seven endpoints, with PC1 explaining 37.13%, PC2 15.59% and PC3 14.24% among the 127 complete taxa.

## Environmental evidence hierarchy

The final evidence atlas contains 31 environmental-evidence rows. Evidence is graded rather than reduced to a universal pass/fail rule. Cross-sectional associations describe spatial phenotype–environment structure; they do not demonstrate plasticity, adaptation, selection or mechanism.

### A — robust established primary associations

| Endpoint | Environment | Standardized effect | q | Interpretation |
| --- | --- | ---: | ---: | --- |
| `orientation_image_vertical_angle` | BIO1 annual mean temperature | +0.0171 | 0.0363 | Global association retained through the relevant seasonal, dominant-taxon, hemisphere and native-range sensitivities. |
| `capitulum_outline_aspect_ratio` | BIO4 temperature seasonality | +0.0137 | 0.0207 | Global association retained through the relevant seasonal, dominant-taxon, hemisphere and native-range sensitivities. |

### B — supported but sensitivity-dependent primary associations

`corolla_lab_chroma` increased with BIO12 annual precipitation (β = +0.0393, q = 0.0363) and `capitulum_outline_aspect_ratio` decreased with BIO12 (β = −0.0108, q = 0.0363). Both retained direction under native-range restriction but lost native-only BH support, so they remain global patterns graded as range-sensitive. Circular hue associations are supported but remain colour-calibration dependent and must be interpreted jointly rather than as named colour shifts.

### C — exploratory expanded-trait associations

The registry-wide candidate screen produced three FDR-supported expanded-trait rows. Endpoint-specific image-quality sensitivity then separated them into two quality-robust exploratory signals and one image-sensitive signal:

| Endpoint | Environment | Standardized effect | q | Final grade | Quality interpretation |
| --- | --- | ---: | ---: | --- | --- |
| `involucre_length_width_ratio` | BIO15 precipitation seasonality | +0.054923 | 0.000282 | `C_exploratory_quality_robust` | Same direction after resolution/sharpness adjustment and in all successful fixed resolution strata. Botanical calibration remains pending. |
| `involucre_apical_taper_ratio` | BIO12 annual precipitation | −0.049721 | 0.000821 | `C_exploratory_quality_robust` | Same direction after resolution/sharpness adjustment and in all successful fixed resolution strata. Botanical calibration remains pending. |
| `bract_projection_p95` | BIO15 precipitation seasonality | −0.033673 | 0.014377 | `C_exploratory_image_sensitive` | FDR support remains in the global candidate screen, but resolution-stratum sign instability prevents promotion. |

The quality-adjusted candidate layer fitted 36 climate rows and returned six quality-adjusted FDR signals. That count must not be reported as six independent headline discoveries: only the three rows above also passed the registry-wide candidate screen and enter the final evidence atlas.

## Submission-facing synthesis

The final GEB v2 result is not a universal climate-associated thistle morphology. Public biodiversity photographs recover a multidimensional, continuous capitulum phenotype landscape with substantial repeated variation below taxon means. Environmental structure differs among phenotype modules and analytical scales. The strongest established primary patterns involve orientation, visible colour and selected gross-outline dimensions; the expanded continuous layer adds two quality-robust but still exploratory involucral-geometry signals and one explicitly image-sensitive signal.

The v2 contribution is therefore twofold: it extends global public-image phenomics beyond species means and coarse categories to an explicit 27-endpoint continuous contract, and it makes robustness part of the result by assigning A/B/C evidence grades without confusing pattern detection with causal or functional inference.

## Claim boundary

Do not use this run to claim phenotypic plasticity, local adaptation, pollinator causation, antagonist defence, botanical spine length, phyllary function, secretion function or evolutionary rate. `visible_floret_fraction` was not executed and cannot support a floral-display inference in this version. Candidate involucre/armature endpoints remain image-derived continuous phenotypes until independent botanical reference validation supports more specific naming.