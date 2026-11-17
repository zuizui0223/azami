# 2026-11-07 | Verified original V4 S1 RF exploratory results: 7 vs 23 variables

**These are genuine V4 source-verified empirical RF fits, NOT legacy V2 and NOT simulations.** Original S1 comprised **313,350 native source observations**, authenticated using archived encrypted-source checks in the protected GitHub runners; no raw original data/observer/taxon identifiers were published. Fixed S1/S2/H3c paper is unchanged.

## Reproducibility
- V4 7 native predictor RF: [Actions #42324966272](https://github.com/zuizui0223/azami/actions/runs/42324966272), succeeded, aggregate-only artifact **13554628048**, JSON SHA-256 `7ef9eafbe03a9d6181ecb2ebb53a487f06c894590b151cb650a4a0f291f3a3a7`.
- V4 all-23 predictor RF: [Actions #42325280652](https://github.com/zuizui0223/azami/actions/runs/42325280652), **rerun attempt 2 SUCCESS**, anonymized artifact **13554764797**, JSON SHA-256 `349d72c6f649e5010aa3288745d86dfef0c666cc7d8ebb05d484a624e71060fb`.
- V4 authentic native S1 SHA-256 `342d78e9a8a960e8d9bdcf0bc2cdacd13bee48214a8b962419d5a71d5aef46ed`; protected seven new-layer SHA-256 `3660407a21a99875f720297b24fec0d0bf46d42a91cf576eb2e185cb63899f77`; nine BIO layers SHA-256 `63a6846ffc730496920dd2e471a5ec148dd630f92cdf97038aeb5c10df6e512a`; source receipt SHA-256 `26b80b279fa2b372099ef945972d3e4fdf58127fe9f8fd92530ac08a931d9545`.
- The **23 predictor columns were available** in each of eight traits; V4 head-angle eligibility gave **178 taxa**, other photographic endpoints **177 taxa**. A taxon was eligible with ≥5 photos; no minimum of 8 cells / 40 photos.
- Both models: RF 80 trees, min_samples_leaf=4, max_features=0.75, same fixed seed 20261107, repeated taxon-heldout 5fold×2=10 heldout folds, 2 permutation repeats; training-only median imputation. Geographic sensitivity: leave-one-dominant-longitude-sector-out. Results are descriptive, with no new RF P-values.

## Predictive Skill by trait
Positive Skill means better mean squared prediction error versus training-fold equal-taxon average. See the accompanying machine-readable CSV.
| Trait | N taxa | Native7 Skill | All23 Skill | Native7 3-sector Skill | All23 3-sector Skill | Top individual All23 feature |
|---|---:|---:|---:|---:|---:|---|
|angle|178|+0.00668|−0.004999|−0.12215|−0.05076|chelsa_bio07|
|lightness|177|+0.22880|+0.19463|+0.02440|−0.08771|chelsa_rsds_mean|
|chroma|177|+0.22001|+0.24011|+0.12309|−0.03293|chelsa_bio08|
|hue|177|+0.19191|+0.26886|+0.06882|−0.02082|chelsa_bio08|
|aspect|177|−0.08907|−0.10113|−0.17937|−0.17122|chelsa_bio10|
|circularity|177|−0.00316|−0.05139|−0.02750|−0.12344|chelsa_bio05|
|solidity|177|−0.03025|−0.04441|−0.00929|+0.04083|chelsa_bio05|
|width_profile|177|+0.05384|+0.07520|−0.16236|−0.26868|topo_elevation|

## Key environment signals and their boundaries
- **Flower color:** three color constructs have positive taxon-heldout Skill with all23: L* +0.195, chroma +0.240, joint circular hue +0.269. Chroma improves +0.020 over native7, hue +0.077; L* falls −0.034. In all23, **BIO08** (mean temperature of wettest quarter) is the highest individual permutation importance for chroma (+0.0535, positive 10/10 folds) and hue (+0.1224, positive 10/10 folds); RSDS remains the leading lightness variable (+0.1135). BIO08's importance is **model-reliance, not independent climate causality** and no calibrated floral pigmentation or spectral response was measured.
- **Orientation:** all23 Skill ~−0.005, so broad model predictive improvement is absent. However, jointly permuting the full precipitation family raises normalized heldout loss by **+0.134** and joint BIO12/BIO13/BIO16 correlated cluster by **+0.095** (10 folds). Individual BIO07 (+0.022) ranks first only in the all23 RF; precipitation remained the largest coherent family. This is a **weak exploratory pattern**, not reliable global angle prediction. True gravitational nodding has not been established.
- **Outline:** aspect, circularity and solidity are negative in out-of-taxon skill. Width-profile CV is modest positive +0.075, but regional transfer −0.269. Do not interpret individual BIO05/elevation ranks as independently robust selection on morphology.
- **Geographic transfer:** all23 is negative for 3/3 color traits (lightness −0.088; chroma −0.033; hue −0.021) despite positive taxon-heldout CV; chroma/hue are positive in two of three separate sectors but large negative middle-sector results dominate. Thus increased feature richness does **not** produce globally transferable environmental rules.
- Rank frequencies and environment family/cluster importances are model-condition- and fold-dependent. Correlated-variable permutation is a diagnostic; group ranks with one fold are unstable. The same dataset had been inspected previously, so this RF exploration is not an independent confirmation.
- Legacy V2 is not used in these V4 results. Existing GEB original S1 80 Holm tests, S2 128 Holm tests and H3c are unchanged.

## Planned next model assessment
Document role of BIO08 vs RSDS / elevation under grouped conditional importance and regional sample composition. Do not substitute old 26-taxon within-grid RF results or 44-axis uncollected candidates as if fitted here. Preserve the full anonymous JSON in GitHub Actions artifact or a reviewed aggregate archive; never publish row-level positions, observer IDs or taxon identities.
