# Chapter 1 Supporting Information

This directory is the submission-facing Supplement source for the current continuous within-taxon image-phenomics manuscript.

- `SUPPLEMENTARY_INFORMATION.md` contains supplementary methods, ecological interpretation boundaries and figure captions through Figure S8.
- `REFERENCES.md` is the self-contained reference list for citations appearing in the Supporting Information.
- `tables/` contains the complete frozen S01–S12 machine-readable summary set tied to the current claim registry.
- `figure_data/` contains the frozen numeric inputs for Figures S1–S5.
- `FIGURE_MANIFEST.md` records source artifacts, interpretation boundaries and SHA-256 hashes for Figures S1–S5.
- `../FIGURE2_GEOGRAPHY_MANIFEST.md` records the frozen source artifacts and interpretation boundary for main Figure 2 and Figures S6–S7.
- `../INTERPRETIVE_FIGURES_MANIFEST.md` records the frozen source artifacts and interpretation boundary for main Figures 1/4/5/6 and Figure S8.
- `analysis/build_supplement_figures.py` regenerates Figures S1–S5 as SVG, PNG and PDF from committed numeric inputs.
- `analysis/build_main_figure2_geography.py` regenerates main Figure 2 and Figures S6–S7 directly from the frozen exhaustive-cohort and enriched-environment artifacts.
- `analysis/build_ch1_interpretive_figures_v3.py` is the final quantitative entry point for main Figures 4–6 and Figure S8; `analysis/augment_figure1_with_analysis_scales.py` appends the multiscale-analysis strip to the validated actual-photo Figure 1.
- `.github/workflows/ch1-supplement-figures-ci.yml` verifies the S1–S5 release hashes.
- `.github/workflows/ch1-main-figure2-geography-ci.yml` validates frozen cohort counts, downloads the frozen observation artifacts, builds the geographic figures and records source/release hashes.
- `.github/workflows/ch1-interpretive-figures-ci.yml` validates the PCA cohort/variance, supported-row counts and stable SPDE directions before regenerating the interpretive figure bundle.

Git stores scientific inputs, provenance and deterministic figure recipes rather than duplicating large observation-level tables or binary release products. Figures S1–S5 use committed numeric inputs with fixed release hashes; Figures S6–S8 remain tied to immutable workflow-artifact IDs and frozen result tables documented in their manifests.

S03 and S06 contain submission-relevant BH-supported coefficient rows; unsupported rows are not reconstructed from rounded manuscript prose. S05 records only SPDE directions classified as stable in the frozen claim registry, while the complete four-group SPDE model comparison used for Figure S3 is separately frozen in `figure_data/` with artifact provenance.

Figures S6–S7 are sampling/measurement diagnostics rather than additional biological tests. Figure S6 shows how geographic sampling density changes through the coordinate and positional-accuracy filters and taxon × 0.25° thinning. Figure S7 shows geographic variation in whether image-derived orientation, colour and outline are assessable; unassessable images remain measurement missingness and are never interpreted as biological trait absence. Figure S8 is a presentation of the observed geometry underlying the already frozen among-taxon niche-permutation test: it makes centroid separation and niche overlap readable across traits but does not create a new P-value, multiple-testing family or causal claim.

Ecological literature is used in the Supplement to motivate candidate functional interpretations and, equally importantly, to define their limits. Orientation references motivate floral microclimate/exposure hypotheses, flower-colour references motivate mixed biotic/abiotic mechanisms, and *Cirsium* pollinator–antagonist literature motivates alternative functions for fine involucral architecture. None of these external studies upgrades the present observational associations to causal adaptation.

The withdrawn raw lability/quadrant result appears only as statistical provenance/QA. It is not a biological supplementary result. Independent detector and continuous-measurement validation remain external scientific gates and are not replaced by these figures or references.
