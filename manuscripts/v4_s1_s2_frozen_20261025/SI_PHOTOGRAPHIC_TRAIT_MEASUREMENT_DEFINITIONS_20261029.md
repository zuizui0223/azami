# Table S4. Continuous photographic trait definitions and measurement provenance

## Purpose and observation unit

The workflow retains continuous presentation, colour and outline measurements together with source observation, photograph, taxon and location. This makes the same photographic record available for individual-trait gradients, descriptive trait associations and local climate-analogue comparisons. The table documents the archived measurement procedures used to generate the numerical exports; the current V4 analyses reuse those values and their existing eligibility flags.

Candidate capitula were localised with YOLO. Colour and geometry routines then measured the detected image regions. Thus object detection and continuous trait calculation are distinct stages. The record-level summary is a median of usable head-level measurements within an observation–photo combination. Head eligibility is domain-specific, so different sets of heads can contribute to presentation, colour and outline in one photograph. “Same photographic record” is the supported common unit; an identical individual head across all domains is not guaranteed.

## Table S4. Eight constructs and their archived definitions

| Domain | Construct | Head-level image definition | Information retained |
| --- | --- | --- | --- |
| Presentation | Image-referenced capitulum angle | Direction of the foreground principal axis, oriented towards floral pixels, relative to image vertical after available EXIF correction: image top 0°, horizontal 90°, image bottom 180°. | Continuous projected head presentation. |
| Colour | Lightness | Median CIELAB lightness over the segmented floral-pixel mask. | Variation in photographed floral lightness. |
| Colour | Chroma | Median CIELAB chroma, \(\sqrt{a^{*2}+b^{*2}}\), over the selected floral-pixel mask. When red–magenta is the dominant colour class, its mask is used; otherwise the union of floral masks is used. The OpenCV encoded a and b channels are centred by subtracting 128 before computing chroma. | Variation in photographed colour intensity. |
| Colour | Circular hue | Circular mean of HSV hue over floral pixels whose OpenCV saturation is at least 25. OpenCV H is multiplied by 2 to obtain degrees. Head-level sine and cosine coordinates represent one circular construct. | Continuous hue around the colour circle. |
| Outline | Aspect ratio | Height divided by width of the non-zero mask bounding box after PCA-based image rotation. | Projected outline proportions. |
| Outline | Circularity | \(4\pi A/P^2\), using segmented contour area \(A\) and perimeter \(P\). | Outline compactness. |
| Outline | Solidity | Segmented contour area divided by its convex-hull area. | The fraction of the convex hull occupied by the outline. |
| Outline | Width-profile coefficient of variation | Population standard deviation divided by mean of row-wise foreground-pixel counts after PCA rotation and cropping. Row widths are normalised by the maximum width; only rows with normalised width >0.05 are retained. The implementation uses a denominator floor of \(10^{-6}\). | Variation in width along the rotated outline. |

There are **eight constructs and nine numerical coordinates**, because hue uses sine and cosine. Aspect is the rotated bounding-box height/width, not a ratio of PCA eigenvalues or an assumed minor/major axis ratio. The width-profile calculation uses pixel rows, not a fixed number of longitudinal bins. Hue is HSV-derived, not the CIELAB hue angle \(\operatorname{atan2}(b^*,a^*)\).

Archived head-level routines combine original and horizontally mirrored image measurements, using circular treatment for hue. At photographic-record aggregation, eligible head-level hue sine and cosine values are summarised separately by their medians. This should not be described as a simple circular mean across all heads in the photograph. The full source code defines mask construction, original/mirrored combination and assessability flags; this table does not assert that default script thresholds were rerun on every V4 image.

## How measurements become comparable ecological records

For native-taxon regressions, each endpoint is selected using its original eligibility flag, a finite measurement and the source environmental-completeness gate. The source selector requires all 12 candidate environmental values to be finite, even though the fitted model contains four environmental predictors. Angle is additionally constrained to 0–180°, and the native regression requires at least five eligible observations per taxon.

The eight-construct descriptive cohort contains **156,683 photographs from 174 taxa**. All 28 construct-pair summaries use that common cohort, so differing pairwise missingness does not determine the set of records compared. The native angle model uses its distinct **168,907-record, 178-taxon** cohort.

H3c applies a separate gate for angle, chroma and aspect on the same photographic record, with finite BIO12, BIO01 and RSDS. Its initial eligible source comprises **157,822 records from 242 taxa**. Subsequent requirements of at least 40 records and eight occupied cells retain **63 taxa at 1°** and **31 at 5°**. Cell-level coordinate medians, taxon-weighted standardisation and whole-vector permutations are defined in main Section 2.7 and the [H3c supplement](SI_H3_CLIMATE_ANALOGUE_CONDITIONAL_OVERLAP_20261029.md).

These counts establish the demonstrated coverage of the observation-and-comparison framework. They are not estimates of object-detection precision or measurement error against field ground truth. Pixel colour describes photographed visible appearance; image angle and outline describe a two-dimensional projection. Calibration to gravity, spectra, anatomical dimensions or reproductive performance requires paired reference measurements.

## Source implementation and analysis provenance

Definitions were checked against the repository state at commit `5bc37b22adf4d8411fd7814df1aaa58de3d94fcb` before this editorial revision. The measurement code and numerical results were not changed.

| Stage | Archived source | Relevant content |
| --- | --- | --- |
| Candidate localisation | [07_screen_downloaded_images_with_yolo.py](../../legacy/ch1_global/v2/07_screen_downloaded_images_with_yolo.py) | Detector predictions, candidate regions and source identities. |
| Continuous geometry and head outputs | [52_measure_primary_traits_continuous.py](../../legacy/ch1_global/v2/52_measure_primary_traits_continuous.py) | PCA rotation, contour geometry, row-wise widths, original/mirrored measurements and domain eligibility. |
| Colour refinement and execution | [56_run_primary_traits_continuous_v2.py](../../legacy/ch1_global/v2/56_run_primary_traits_continuous_v2.py) | Floral masks, HSV circular hue, CIELAB lightness and selected-mask chroma. |
| Photographic-record aggregation | [53_aggregate_primary_trait_continuous.py](../../legacy/ch1_global/v2/53_aggregate_primary_trait_continuous.py) | Grouping by observation, photo, taxon and location; trait-domain usable heads and record-level medians. |
| Input and trait-overlap audit | [native_integration_audit.py](../../analysis/v4/native_integration_audit.py) | Source identity, endpoint eligibility and domain overlap. |
| Final common cohort and regression selectors | [s1_s2_vif5_si10_alltraits.py](../../analysis/v4/s1_s2_vif5_si10_alltraits.py) | Final common-cohort selection and 28 construct-pair associations in `s1_joint_28()`; original eligibility flags, source environmental completeness and regression selection. |
| Original V4 evidence | [V4 freeze](../../analysis/v4/FINAL_V4_BIO18_HOLM_LEGO_FREEZE_20261024_JA.md) | Cohort counts and the original regression/covariation results. |
| Climate-analogue evidence | [H3c source-verified results](../../analysis/v4/results/20261029_H3C_CLIMATE_ANALOGUE_TRAIT_OVERLAP_SOURCE_VERIFIED.json) | Matched source, spatial eligibility, frozen distance estimates and permutation results. |

Earlier technical repeatability and botanical-category exercises use different archived cohorts. They are not promoted here as validation of current V4 biological measurement accuracy. The methodological contribution reported in the main article is the demonstrated retention and ecological use of continuous, joint and geographically linked photographic measurements.
