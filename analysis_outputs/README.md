# Reviewer bias-control outputs

These outputs were generated on 2026-08-26 under the locked rules in
`analysis/ch1/bias_control_contract.json`. They supplement rather than overwrite
the frozen primary analysis.

- `reviewer_bias_control_v1/`: cyclic day-of-year adjustment and outcome-blind
  dominant-taxon omission models for all 36 primary endpoint components.
- `hemisphere_season_sensitivity_v1/`: Southern-Hemisphere phase alignment and
  taxon-by-hemisphere cyclic-season sensitivities under a separately locked
  contract.
- `involucre_resolution_audit_v1/`: multivariable image-resolution and
  sharpness adjustment plus fixed resolution strata.
- `repeat_photo_bias_control_v1/`: outcome-blind repeat-photo cohort and public
  image manifest; image measurement is pending.
- `native_range_sensitivity_v1/`: versioned WCVP name/distribution resolution,
  pinned TDWG level-3 geometry provenance and native-only model outputs.
- `interpretive_figure_preview_v2/` and `supplement_figure_preview_v2/`:
  reproducible review copies of the figures used for manuscript/DOCX QA. The
  corresponding workflows remain the release mechanism.

The strict observation input had SHA-256
`fb6a2bd73937d1b703f355193361abeef4cbe5192c5168b94ba53857607c9469`.
The photo metadata input had SHA-256
`506917aef9a9452dddcf373acfb941e9f3620ba9740ef588cbcf8725557828e4`.

The native-range gate is complete and retained two of four non-circular
primary rows. No output in this directory closes developmental-stage,
colour-calibration, repeat-photo measurement or independent
measurement-validation gates. These files are an auditable PR snapshot, not a
substitute for the planned DOI-backed durable data/code release.
