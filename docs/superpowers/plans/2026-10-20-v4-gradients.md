# V4 native-observation environmental gradients implementation plan

Goal: estimate conditional environmental gradients of colour, orientation and bract traits within and between species using individual observations.
Architecture: retain verified inputs and endpoint-specific eligibility; separate environmental selection, mixed-model estimation, gradient reporting and anatomical measurement coverage. Run all numerical analyses and tests in GitHub Actions.
Spec: user instruction 2026-10-20; this supersedes the common-slope main analysis, whose results remain archived.

- [ ] Restore within/between execute path in analysis/v4/environment_response.py. Keep raw trait rows, all selected environments and species random intercepts. Verify response-specific counts and 48 tests for eight selected predictors in Actions; do not change correction after seeing results.
- [ ] Environmental justification: audit all nine available candidates and six group definitions, temporal/spatial support and derived-variable dependencies. Existing vifstep threshold 10 is retained for the reproducibility run, not asserted to be an optimal causal selection. Inspect within/between/full-design diagnostics before deciding whether final specification needs revision. No significance-based selection.
- [ ] Report gradients in original trait units, with intervals and supported exposure ranges at each level. Preserve circular hue coordinates; PCA optional descriptive figure only, no NMDS needed. Do not extrapolate every species across the pooled exposure range.
- [ ] Measurement coverage: audit the 313350 native source records, distinguish measured, QC failure, absent original image and unprocessed. Preserve all eligible existing trait observations; no cross-module complete-case filter or old spatial-cohort intersection.
- [ ] Anatomical bract measurement: verify original-image workflow 39419419586; review overlays and base/tip definitions, then implement and validate scalable annotation/detection before broad measurement. Four candidate axes in one photo are not validated full-source coverage. Preserve visibility, whorl where identifiable, source identity and failures.
- [ ] Quantify species influence; evaluate spatial and phylogenetic limitations without treating ordinary random intercepts as either correction. Publish traceable effects, uncertainty and conclusions.

Completion requires actual artifacts for gradients, coverage, measurement extension and diagnostics; green CI alone is insufficient. Common-slope run 39418226382 and decomposed run 39416125522 remain separate provenance chains.

Robustness implementation: fixed exclusions of the top 2 and 5 species ranked by eligible observation counts separately for each module; fixed selected exposures, full-native means and scaling. Broad-space comparison adds the existing eight spherical-coordinate terms on identical endpoint rows. These are diagnostic fits, not alternative main-model selection. Record raw effect intervals and sample counts; no significant-result rescue or new primary family.
