# Corrected native-observation environmental gradients

Source: https://github.com/zuizui0223/azami/actions/runs/39422378058
Analysis commit: 93e8d87. These aggregate files are copied from Actions, not recomputed locally.

This input applies all nine verified source scale/offset transforms and masks 1,239 invalid GSP encodings. All 313,350 native rows remain in the source table; 312,001 have all nine environments. Per-endpoint models use their own eligible rows: colour 185,515–185,564 (244 taxa), orientation 169,263 (244), outline 175,223 (243), bract projection SI 1,298 (111). Environmental selection again excludes only mean temperature (BIO1).

All selected environments are decomposed into species means and observation deviations and fitted simultaneously with species random intercepts. Trait observations are not averaged. Main includes colour, orientation and four outline coordinates; bract projection is SI. Joint tests use Holm across 64 comparisons.

Correction changes the scientific conclusions: the previously reported between-species rainfall–orientation association does not survive. Colour and outline retain selected environmental associations; consult the complete main report and intervals rather than transferring uncorrected results. These estimates have not yet passed species-influence, spatial or phylogenetic follow-up checks.

The source JSON SHA256 is 563bead7d30ca59cd00b5c435457a654de53656f5f74e84e1d164b7dc330341a. Physical-unit presentation is rendered separately from this saved fit in Actions without refitting.
