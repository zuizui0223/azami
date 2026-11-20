# Photo-angle exemplar audit for the Species Biology Symposium poster (20 Nov 2026)

## Why the poster needed correction

The three photographic examples in the poster are real and the displayed numeric angle values were transcribed from `legacy/v3/reproducibility/figures/source/Figure_1_real_photo_measurement_provenance.csv`. **Correct transcription is not sufficient for valid biological interpretation.**

| Example | Archived image-based theta | Archived rotated outline aspect | Visual-inference status |
| --- | ---: | ---: | --- |
| *C. pannonicum* (iNaturalist photo 532420148; Bobo-X CC BY 4.0) | 0.732009° | 1.59259 | Elongated example with image-angle near 0°; still not calibrated to gravity |
| *C. mohavense* (photo 4051162; Jesse Rorabaugh CC0) | 111.147764° | 1.01149 | Nearly circular projected aspect: measured image-angle must **not** be presented as a verified downward-facing head |
| *C. subcoriaceum* (photo 475732253; Diego Blanco CC0) | 113.284419° | 1.14035 | Near-circular projected aspect: measured image-angle must **not** be presented as a verified downward-facing head |

The source-matched angular definition is recorded in `manuscripts/v4_s1_s2_frozen_20261025/SI_PHOTOGRAPHIC_TRAIT_MEASUREMENT_DEFINITIONS_20261029.md` and implemented in `legacy/ch1_global/v2/56_run_primary_traits_continuous_v2.py`:

- The foreground head-matrix PCA first axis is **signed using the visible floral-pixel centroid**.
- Its angle is computed against **EXIF-oriented image-up (0°)**, with image horizontal 90° and image-down 180°.
- It is **not** a gravity-referenced measurement, nor the v1 angle between the head and peduncle.
- The rotated outline's bounding-box aspect ratio is not the same as the PCA eigenvalue ratio. Ratios near 1 here **raise a stability concern** but do **not by themselves quantify angle measurement error**.

The existing source licence/figure audit `manuscripts/v4_s1_s2_frozen_20261025/FIGURE1_PHOTO_CREDITS_20261030.json` also records that an older raster label gave **1.1°** for *C. pannonicum*, whereas the archived production table gives **0.732009°**. The corrected figure should use **0.7°** (rounded production value); the legacy 1.1° overlay must not be reused.

## Poster correction

The revised PPTX `SSSB2026_Cirsium_V5_A0縦_角度表示監査修正版.pptx`:
- Retains the source numbers but flags the last two angle values as **reference image measurements whose biological orientation is unverified**.
- Rephrases the measurement as a PCA/image-vertical angle, not a physical gravity or peduncle inclination.
- Separates the statistical observation (BIO12–signed image angle) from the **untested** rain/pollen protection hypothesis.
- Does **not** alter any frozen V5 slope, Holm P, analogue results or original source rows.

## Submission-critical next checks (not yet run)

1. Independently annotate head-axis orientation and peduncle/gravity references on a blinded geographic/taxonomic photo subset; report angular errors and image roll effects.
2. Test V5 angle results after excluding/reweighting weak PCA anisotropy or low floral-end-projection confidence. The source code's anisotropy quality is computed from PCA eigenvalues; report its distribution rather than using aspect ratio as a substitute.
3. Stratify by photographer, region, taxon and flowering stage to check whether image geometry/assessability covary with BIO12/BIO15.
4. Recheck the **joint three-trait climate-analogue result** after angle quality filtering, not only the signed BIO12 regression. Effect sizes and P-values under these sensitivities are unknown and must not be assumed unchanged.

This is an image-phenotype validity issue; it is not evidence that the previously frozen numerical outputs were fabricated or arithmetically mistyped.
