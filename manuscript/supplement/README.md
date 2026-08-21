# Chapter 1 Supporting Information

This directory is the submission-facing Supplement source for the current continuous within-taxon image-phenomics manuscript.

- `SUPPLEMENTARY_INFORMATION.md` contains supplementary methods, ecological interpretation boundaries and figure captions.
- `REFERENCES.md` is the self-contained reference list for citations appearing in the Supporting Information.
- `tables/` contains the complete frozen S01–S12 machine-readable summary set tied to the current claim registry.
- `figure_data/` contains the frozen numeric inputs for Figures S1–S5.
- `FIGURE_MANIFEST.md` records source artifacts, interpretation boundaries and SHA-256 hashes for all canonical figure outputs.
- `analysis/build_supplement_figures.py` regenerates Figures S1–S5 as SVG, PNG and PDF from the committed numeric inputs.
- `.github/workflows/ch1-supplement-figures-ci.yml` runs the pinned generator, verifies every frozen figure hash and uploads the 15-file submission figure bundle.

Git therefore stores the scientific inputs and deterministic figure recipe rather than duplicating large binary release products. A figure bundle is valid only when its hashes match `FIGURE_MANIFEST.md`.

S03 and S06 contain submission-relevant BH-supported coefficient rows; unsupported rows are not reconstructed from rounded manuscript prose. S05 records only SPDE directions classified as stable in the frozen claim registry, while the complete four-group SPDE model comparison used for Figure S3 is separately frozen in `figure_data/` with artifact provenance.

Ecological literature is used in the Supplement to motivate candidate functional interpretations and, equally importantly, to define their limits. Orientation references motivate floral microclimate/exposure hypotheses, flower-colour references motivate mixed biotic/abiotic mechanisms, and *Cirsium* pollinator–antagonist literature motivates alternative functions for fine involucral architecture. None of these external studies upgrades the present observational associations to causal adaptation.

The withdrawn raw lability/quadrant result appears only as statistical provenance/QA. It is not a biological supplementary result. Independent detector and continuous-measurement validation remain external scientific gates and are not replaced by these figures or references.
