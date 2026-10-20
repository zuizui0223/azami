> Interpretation hold (2026-10-20): uncorrected GSP includes uint32-max sentinel values. Environment-adjusted v4 results are provisional pending input correction and replay. See analysis/v4/ENVIRONMENT_INPUT_INCIDENT.md.

# V4 common-slope environmental response

Verified Actions run: https://github.com/zuizui0223/azami/actions/runs/39418226382

Analysis commit: `d0a588b8c6905fc47ecb7b81405d8c9c48aace9a` on PR143. All 18 numerical tests passed. These files are copied without recomputation from the `v4-observation-model-results` artifact.

Model: observation-level phenotype ~ eight selected environmental variables + (1 | species). No response averaging or within/between decomposition. All nine stored environmental variables were screened with usdm::vifstep before fitting responses; BIO1 alone was excluded. Selection does not eliminate unmeasured confounding or establish causal effects.

- [Main results: whole colour and orientation](environment_response_ja.md)
- [SI: bract contour projection](environment_response_si_ja.md)
- [Complete aggregate output](environment_response.json)
- [Effects and pointwise 95% intervals](environment_response_effects.csv)
- [Environmental selection provenance](environment_selection.json)
- [Measurement coverage](measurement_reuse_inventory.csv)

Whole colour is associated with temperature seasonality, radiation, VPD, wind and growing-season precipitation after Holm correction over 24 module-environment tests. Higher VPD accompanies higher lightness and lower chroma. Orientation has no association surviving this correction, including annual precipitation. Bract projection in SI is associated with NPP (adjusted P=0.04602); it is not a direct measurement of defence.

Colour uses 186,388–186,437 observations from 244 taxa; orientation uses 170,064 from 244 taxa; bract projection uses 1,303 from 111 taxa. The native source population contains 313,350 observations derived from the 665,115-observation upstream source. Measurement coverage does not equal source coverage.

The [earlier within/between fit](../20261020_selected_environment/README.md) remains separately archived. Its between-species rainfall association with orientation must not be transferred to this common-slope fit. Common slopes do not establish that each species has the same response. Spatial and phylogenetic adjustment, pigment chemistry and adaptive mechanisms are not demonstrated by this fit.
