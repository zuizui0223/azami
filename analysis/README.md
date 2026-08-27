# Submission-facing v2 analysis

`analysis/ch1/pipeline.json` is the executable map. The active scientific lane begins with the 27-endpoint contract and ends with the full-environment sampling, broad-space and historical-placement gates.

Primary entry points:

- `run_geb_v2_full27_environment_atlas.py`
- `validate_geb_v2_full27_environment_atlas.py`
- `run_geb_v2_full27_sampling_composition_sensitivity.py`
- `validate_geb_v2_full27_sampling_composition.py`
- `run_geb_v2_full27_spatial_sensitivity.py`
- `run_geb_v2_full27_historical_sensitivity.py`
- `validate_geb_v2_full27_sensitivities.py`
- `validate_v2_submission_artifacts.py` verifies the tracked frozen validation reports, submission manifest, hashes and reviewer ZIP without requiring protected observation-level inputs.
- `validate_manuscript_citations.py` requires a one-to-one match between active author–year citations and the reference list.
- `build_v2_submission_figures.py`
- `build_v2_submission_package.py`

Numbered measurement scripts remain under `ch1_global/v2/` because their paths are part of frozen execution provenance. The old result-selection and submission-bundle scripts were removed from the active tree.
