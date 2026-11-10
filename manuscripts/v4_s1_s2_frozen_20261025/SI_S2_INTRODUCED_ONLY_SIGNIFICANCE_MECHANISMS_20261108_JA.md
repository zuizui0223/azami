# V4 Results 3.2：「移入域だけ有意」の生態学的・統計学的説明（2026-11-08）

**Position: hypothesis audit, not a new result or a test of rapid adaptation.** Scientific freeze (original 128-Holm tests) unchanged. Literature is newly verified outside the original V4 source; original observations below derive from V4 frozen main manuscript and spatial-block SI.

## Most important correction: three Holm discoveries are not three introduced-only mechanisms

| Original S2 association | Native β / pointwise CI | Introduced β / pointwise CI | Original Holm128 significance |
|---|---|---|---|
| *C. vulgare* RSDS × photographic chroma | +0.136883 / [+0.056600,+0.217166] (n 17,335) | +0.134504 / [+0.066158,+0.202850] (n 27,200) | Native P=.102398; intro P=.014568 |
| *C. vulgare* RSDS × photographic lightness | −0.048808 / [−0.147073,+0.049457] (n 17,335) | −0.142801 / [−0.216718,−0.068884] (n 27,200) | Native P=1; intro P=.019254 |
| *C. arvense* elevation × photographic lightness | +0.056267 / [−0.018662,+0.131195] (n 22,936) | −0.173732 / [−0.250914,−0.096550] (n 26,171) | Native P=1; intro P=.001312 |

All are **four-variable within-range conditional associations** (BIO18,BIO05,RSDS,elevation), standardised **separately within each species×range stratum**. The original models **did not estimate a native–introduced interaction**. Different Holm decisions `significant here / not there` do not themselves imply a difference in slopes. In *C. vulgare* RSDS×chroma, native and introduced β are **nearly identical**, and the native pointwise CI excludes zero. Do not describe this as an introduced-only phenomenon.

* C. arvense* elevation×lightness is a more interesting **descriptive opposing-sign contrast**, but native and introduced SD scales, temperature/elevation distributions, histories, spatial support and image opportunities differ. No direct native–introduced slope difference established. Previous exploratory X-only S2 interaction family 80 did **not** include elevation (or RSDS for *C. vulgare*) and therefore cannot be used to claim these original candidate interaction effects failed.

Spatial block sensitivity of the frozen S2 *C. arvense* introduced elevation coefficient: 399 resamples across only 10 independent-ish 20° source blocks, 95% percentile [−0.22949,−0.05578]; all 5 largest-block exclusions retain negative sign. This increases sampling robustness but not experimental causality or proof of adaptation.

## How could introduced populations show a stronger environmental trait gradient?

**Prioritised hypotheses, separating their evidence type.**

1. **Detection, differing environment and population samples (highest immediate priority).** Holm correction over 128 outcomes, photo-level n vs geographic cells, native and introduced environmental supports, separate z-standardization, covariance/collinearity and photographer/location composition. A stronger apparent effect in a regional pooled introduced dataset can arise without adaptation, and even native positive relationships may fail Holm by a threshold (the *C. vulgare* chroma example proves this risk).
2. **Source-population filtering / preadaptation / repeated introductions and admixture (best species-specific historical mechanism).** *C. arvense* genetic study of 1,522 individuals/58 European/North American populations documented multiple European introductions, bottlenecks, admixture and substantial bioclimatic niche conservatism. A restricted set of introduced successful lineages may form a clearer climate–trait geographical pattern although the native range includes a heterogeneous mosaic. This is a hypothesis for the V4 lightness pattern, **not a demonstrated colour-provenance link**.
3. **Environmentally induced plasticity and floral stage/sex confounding.** A native–introduced *C. arvense* common-environment stress-response transcriptome study (5+5 populations, 3 treatments) found 4,458 ESTs (~7.1%) with treatment×range responses, supporting potential range-dependent response architecture but **not specifically pigmentation**. *C. arvense* is incompletely dioecious, and male vs female heads, flower-stage and mature male-browning versus female white pappus differ appreciably. If corolla-pixel masks or anthesis QC and sex proportions vary across altitude/range, source photographic L* can acquire a false elevation signal. Explicit flowering-stage/sex and mask checks are a crucial inexpensive sensitivity test; their role in V4 **has not yet been measured**.
4. **Rapid evolution/local adaptation after introduction (biologically plausible but unproven for the observed photograph colour).** Formal experiments in *Hypericum perforatum*, *Lythrum salicaria* and recent *Trifolium repens* have found adaptation to climate within introduced ranges; a 2015 review/meta-analysis found evolutionary differentiation in invasive populations but emphasized that adaptation is frequently inferred without confirmation, and source/geography mismatch creates false differences. Rapid evolution cannot be determined from an image-derived gradient without genetic source-population comparisons, common gardens, fitness or reciprocal transplant.
5. **Novel ecological context, enemy release/biotic selection.** Invasion changes herbivores, pollinators, competition and disturbance. These may alter selection or head appearance and can covary with elevation, but V4 lacks species-specific evidence tying these processes to L*; maintain low-priority exploratory status.

**Native populations can be locally adapted yet show no simple global linear cline.** Different optimum mechanisms (temperature, UV, herbivory, pollinator environment), regional populations with different histories, non-monotonicity, or long-standing stabilizing selection can obscure a single across-continent coefficient. A significant photo-climate slope is not a direct local-adaptation test, in either range.

## Decision-grade next tests on V4 original source

**A. Before evolutionary interpretation (first priority, common support and measurement):**
- Use the *same* original four covariates BIO18,BIO05,RSDS,elevation and **the same physical elevation units** (per 100 m or 1 km), and restrict each species to native/introduced shared elevation–climate support; display native and introduced n photographs, equal-cell supports, occupied elevation distributions and blocks before interactions. No extrapolation if support insufficient.
- For *C. arvense*, fit **photo L* ~ elevation × introduced + BIO18 + BIO05 + RSDS + date nuisance**, using common units and pooled environmental centering, observer×cell balancing and geographical-cluster robust or block bootstrap. Explicitly test the original candidate interaction, not *p_native vs p_intro*. Sensitivity: within-georegion-centred elevation vs between-region elevation; quantify available elevation range within each source area.
- For *C. vulgare*, first test RSDS×introduced on chroma on common radiative support; its original two slopes ~+0.1369 vs +0.1345 suggest **no simple introduced-only gradient**, unlike the Holm decision.
- Photo QC for *C. arvense*: blind reviewer subset showing male/female flowering capitula, anthesis vs spent heads, emergence of white female pappus / brown male florets, flower-pixel mask quality and colour-balance/exposure if available. Compare models restricted to clear peak anthesis and sex-classified subset with unidentified category retained. No claim that sex actually caused the published slope absent this test.

**B. If an introduced-only slope persists after A (mechanism discrimination):**
- Molecular provenance/source-lineage matching and climate overlap; test origin and invasion history before interpreting post-introduction evolution. Species-level published history establishes *possibility of admixture*, not proof for photographic colour.
- Common garden or reciprocal transplant of native/introduced genotypes at controlled temperature × radiation to separate plasticity from genetic differences; direct floral spectral reflectance, pigment and temperature are required to argue a thermal-adaptive mechanism.

## Primary publications checked

- **Guggisberg et al. (2012)**, *Invasion history of North American Canada thistle, Cirsium arvense*, Journal of Biogeography, DOI https://doi.org/10.1111/j.1365-2699.2012.02746.x — microsatellites 1,522 individuals, 58 populations, multiple introductions, admixture, conservatism.
- **Guggisberg et al. (2013)**, *Transcriptome divergence between introduced and native populations of Canada thistle*, New Phytologist, DOI https://doi.org/10.1111/nph.12258 — native/introduced 10 populations; 4,458 ESTs treatment×range expression; not colour.
- **Tiley (2010)**, *Biological Flora of the British Isles: Cirsium arvense (L.) Scop.*, Journal of Ecology, DOI https://doi.org/10.1111/j.1365-2745.2010.01678.x — floral dioecy, head-stage morphology, colour changes with maturity.
- **Colautti & Lau (2015)**, *Contemporary evolution during invasion: evidence for differentiation, natural selection, and local adaptation*, Molecular Ecology, DOI https://doi.org/10.1111/mec.13162 — often inferred vs proven selection, geographic sampling and founder effect.
- **Colautti, Maron & Barrett (2009)**, *Common garden comparisons of native and introduced plant populations: latitudinal clines can obscure evolutionary inferences*, Evolutionary Applications, PMID https://pubmed.ncbi.nlm.nih.gov/25567860 — 32 studies / 28 species; latitude changes or reverses range contrast.
- **Innes et al. (2022)**, *Evolution in response to climate in native and introduced ranges of a globally distributed plant*, Evolution, DOI https://doi.org/10.1111/evo.14514 — white clover thermal clines evolving in introduced ranges.
- **Maron et al. (2004)**, *Rapid evolution of an invasive plant*, Ecological Monographs, DOI https://doi.org/10.1890/03-4027 — clines supported by common gardens, not merely observational photography.
- **Fan et al. (2025)**, *Ecological similarities and dissimilarities between donor and recipient regions shape global plant naturalizations*, Nature Communications 16:10485, https://www.nature.com/articles/s41467-025-65455-y — 11,604 naturalized plant species; climate similarity and recipient biota/human modification affect success; evidence for admission/filtering, not automatically local adaptation.

**Decision:** Do not write “introduced-only significant = stronger adaptation after introduction.” The next targeted empirical contrast is the *C. arvense* elevation×range interaction in shared physical units and environmental supports, with parallel region and flower-sex/stage QC. Retain all three original Holm-positive associations and their existing limitations.
