# Chapter 1 v3 biological-axis reanalysis — result summary

## Execution

GitHub Actions run: `34417961671`  
Artifact: `10129842814` (`ch1-v3-biological-axes-34417961671`)  
Artifact digest: `sha256:ba016886d8eeaee0fa7b474cb8f3277d1c42b43d45583d12b2762844ef566dbc`  
Status: **success**

Frozen inputs were recovered by immutable artifact identity:

- 22-endpoint v2 trait universe source artifact `9612943217`;
- nine-predictor process environment source artifact `9633419268`;
- five later-restored display / colour-composition endpoints were explicitly excluded.

## Dimensional reduction

Starting measured endpoints: **22**.

They were reduced to **11 biologically named constructs**, of which ten were inferentially testable at both scales. `surface_specularity` was retained descriptively but had no among-taxon min5 variation and was not promoted into an inferential axis.

Each scale therefore contained **90 construct x environment tests**, rather than treating the measurement endpoints as a larger set of separate biological discoveries.

## Among-taxon min5 result

Twelve construct-gradient rows passed the new exploratory BH family.

The central result is preservation of the frozen v2 scientific spine:

- `floral_chroma × chelsa_rsds_mean`: beta = **-0.345372**, p = **0.0001**, q = **0.00225**;
- `presentation_angle × chelsa_bio12`: beta = **+0.304359**, p = **0.0002**, q = **0.00360**;
- `presentation_angle × chelsa_gsp`: beta = **+0.258078**, p = **0.0018**, q = **0.02025**;
- the six frozen v2 hue associations with BIO1, BIO4, BIO15, radiation, wind and NPP all remain FDR-supported.

Thus **9 of the 10 original v2 among-taxon FDR associations are recovered directly in the reduced biological-construct representation**.

The exception is the original `bract_projection_peak_density × VPD` endpoint association. Once projection endpoints are reorganized into biologically broader projection constructs, VPD is not FDR-supported. Instead the two-dimensional `projection_pattern` construct has an exploratory BIO4 association (joint magnitude 0.4745, p = 0.0054, q = 0.0405). This should not replace the frozen v2 result; it indicates that the armature/projection signal is less stable to biological reaggregation than the two main v2 candidates.

Three additional rows enter the smaller v3 BH family: `floral_lightness × GSP`, `floral_chroma × NPP`, and `projection_pattern × BIO4`. They are **v3-only exploratory results**, not retroactive v2 discoveries.

## Within-taxon result

Nine rows pass the reduced-family BH threshold. The main scalar result is:

- `floral_chroma × VPD`: beta = **-0.028962**, q = **0.000463**.

The joint hue construct retains multiple within-taxon environmental associations, and `head_elongation × BIO4` is also supported (beta = +0.01377, q = 0.02499).

These are cross-sectional within-taxon image associations; they are not demonstrated plasticity.

## Full v2-style sensitivity sequence

A second successful run (`34418904597`, artifact `10130210432`) passed the FDR-supported construct rows through the v2 sampling-composition -> broad/residual spatial -> 52-tree historical-placement sequence.

- sampling: 12/12 among-taxon construct-gradient rows were directionally stable in every declared scenario; 18/21 selected rows across both scales were stable in every scenario;
- spatial: 3/12 among-taxon rows passed; 0/9 within-taxon rows passed the full broad-space plus residual-Moran gate;
- historical: all 3 among-taxon spatial passes remained supported on all 52 placement trees.

Most importantly, **both frozen-v2 headline candidates pass the complete chain after biological reaggregation**:

- `floral_chroma × chelsa_rsds_mean`: sampling minimum effect-magnitude ratio = **0.614876**; spatial beta = **-0.712411**, permutation p = **0.005**, residual Moran p = **0.188**; historical placement **52/52** trees with p < 0.05 and lambda = 0 throughout.
- `presentation_angle × chelsa_bio12`: sampling minimum effect-magnitude ratio = **0.780783**; spatial beta = **+0.286086**, permutation p = **0.023**, residual Moran p = **0.072**; historical placement **52/52** trees with p < 0.05, p range **0.000201–0.000231**, lambda range **0–0.053690**.

`floral_chroma × NPP` is the third construct-level row to pass the full chain, but it is retained as **v3-only exploratory evidence** because it becomes FDR-supported in the smaller construct-level multiplicity family. It is not promoted into the frozen-v2 headline conclusion.

`presentation_angle × GSP` remains positive but fails the construct-level broad-space gate (permutation p = 0.171). The hue rows fail the residual-spatial screen despite several spatial permutation signals. The original bract projection-peak-density × VPD result never enters the construct-level chain because VPD is not FDR-supported after projection reaggregation.

Detailed results and the machine-readable receipt are in [`BIOLOGICAL_AXIS_SENSITIVITY_CHAIN_RESULT_20260910.md`](BIOLOGICAL_AXIS_SENSITIVITY_CHAIN_RESULT_20260910.md) and [`biological_axis_sensitivity_chain_receipt_20260910.json`](biological_axis_sensitivity_chain_receipt_20260910.json).

## Interpretation

The reduced analysis strengthens the v2 headline without redefining it:

1. **Chroma-radiation is stable to biological reaggregation and survives the complete v2-style robustness sequence.**
2. **Orientation/presentation-angle–annual-precipitation is stable to biological reaggregation and survives the complete v2-style robustness sequence.**
3. The hue result family remains coherent at the initial construct-level atlas but does not survive the residual-spatial gate.
4. The projection/armature VPD result is less robust to construct-level aggregation and should remain secondary.
5. New v3-only associations are exploratory consequences of the reduced construct family and must not displace the frozen v2 candidate hierarchy.

This makes v3 a biological synthesis layer over v2: fewer biologically meaningful constructs, the same two principal ecological conclusions, and a stronger demonstration that those two patterns are not artifacts of the original endpoint-level parameterization.
