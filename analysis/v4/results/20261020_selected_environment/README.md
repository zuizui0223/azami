> Interpretation hold (2026-10-20): uncorrected GSP includes uint32-max sentinel values. Environment-adjusted v4 results are provisional pending input correction and replay. See analysis/v4/ENVIRONMENT_INPUT_INCIDENT.md.

# V4 selected-environment response results

Verified [Actions run 39416125522](https://github.com/zuizui0223/azami/actions/runs/39416125522)
completed successfully at commit `b8e18fdb39c73326267649fd29dc4000912c5c9d`.
Seventeen numerical/integration tests passed. All three response modules and
48 joint tests were estimated; the selected set was verified unchanged.
These files are copied from that run's aggregate artifact, not locally refitted.

The source is 665,115 observations, of which 313,350 are native-range records.
Selection used 313,240 native records complete for all nine environmental
candidates in six groups. `usdm::vifstep`, threshold 10, no protected variables
and no subsampling excluded BIO1 only. Final raw-environment VIF: 1.10–3.43.
The selected set is BIO4, BIO12, BIO15, radiation, VPD, wind, growing-season
precipitation and NPP. The group definition, stepwise audit, package version
and exact source hashes are included in `environment_selection.json`.

| Response | Observations | Taxa |
|---|---:|---:|
| Colour lightness/chroma | 186,437 each | 244 |
| Colour hue sine/cosine | 186,388 each | 244 |
| Orientation | 170,064 | 244 |
| Bract projection/spread | 1,303 each | 111 |

The colour response jointly comprises lightness, chroma and circular hue.
After Holm correction across 48 tests:

- Between species: colour is associated with precipitation, radiation, VPD,
  growing-season precipitation and NPP.
- Within species: colour is associated with temperature seasonality, radiation,
  VPD, wind and growing-season precipitation.
- Orientation: precipitation has a positive between-species association,
  +2.395 image-reference degrees per 537.774 mm (one native-population SD),
  pointwise 95% interval [1.033, 3.757], Holm P=0.02398.
  The within-species precipitation coefficient is -0.330 degrees
  [-1.609, 0.948]; no within-species orientation association survives correction.
- Bract contour projection: no environmental joint test survives correction.
  This does not establish absence of a biological defence function.

These are conditional observational associations, not independent causal effects.
The model does not include mean temperature after selection; do not describe it
as temperature-adjusted. Within/between decomposition and endpoint restriction
change collinearity: colour between-species radiation VIF is about 7.75, despite
the lower raw-population VIF. Full endpoint design diagnostics are retained.
No spatial or phylogenetic correction is claimed in this selected-model result.
Do not transfer robustness findings from different earlier predictor sets.
The latest user-directed scope excludes further model additions; previous
results have been seen, so this is not an independent preregistered confirmation.

## Files

- `environment_response.json`: complete aggregate models, intervals, diagnostics,
  source receipt, selection receipt and code hashes.
- `environment_response_effects.csv`: all component coefficients and intervals.
- `environment_response_ja.md`: Japanese joint-test and effect tables.
- `environment_selection.json`, `selected_environment_vif.csv`: selection audit.

SHA-256 of original Actions files:

| File | SHA-256 |
|---|---|
| environment_response.json | `707e3cf0041d127e085bfa47ef83a6ff89a1ad84740cd4ec6c09a893bbe1c517` |
| environment_response_effects.csv | `d10cd4036c60a116535c550360e88de569db70bbbdc0d7824dd2a4a2083559e2` |
| environment_response_ja.md | `7f691c01217aea37682c1239412381002c40a3a6acf9f056244da8835ddea51e` |
