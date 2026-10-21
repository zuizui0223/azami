# Azami capitulum phenotypes

Current analysis and reproducible numerical evidence for the Chapter 1 study prepared for **Global Ecology and Biogeography (GEB)**. Manuscript text, submission files and private reviewer correspondence are kept outside this repository.

Original software is released under the [MIT license](LICENSE); see [NOTICE.md](NOTICE.md) for the separate treatment of images, data and third-party material.

## Reproduction boundary

The current statistical analysis **has executable code**. Its replay starts from frozen image-derived continuous measurements plus fixed environmental, spatial and historical auxiliary inputs; it does not start from the original photographs. Original photographs, detector training/provenance and other external source materials belong to a separate upstream layer and are not all redistributed by the minimum numerical archive. See [upstream-data and reproducibility boundary](reproducibility/UPSTREAM_DATA_BOUNDARY_20260926.md).

## Start here

1. [Current numerical reproduction](reproducibility/CURRENT_ANALYSIS.md): verified inputs, one-command execution, and result validation.
2. [Current results and interpretation](analysis/v3/README.md): biological constructs, environmental associations and whole-capitulum integration.
3. [GEB scope freeze](analysis/v3/GEB_SCOPE_FREEZE_20260912.md): frozen claim boundary, submission-readiness work, analysis stop list, and separate evidence routes for any higher-tier successor.
4. [Current figure surface](reproducibility/CURRENT_FIGURE_SURFACE_20260912.md): Main-figure role freeze aligned to the v3 story and the boundary between the 46,276/259 atlas and 1,734/42 integration cohort.
5. [GEB submission audit](reproducibility/GEB_SUBMISSION_AUDIT_20260912.md): journal-format, double-anonymous review-access, abstract/display-piece and Supporting Information gates.
6. [Reproducibility and archive status](reproducibility/README.md): published Zenodo inputs and their verification.
7. [Historical implementations](legacy/README.md): v2 and superseded entry points. Do not start the current paper analysis from these scripts.

## V4 current research: rain protection and the evolution of capitulum orientation

**Main question:** *Has rain shaped the orientation of thistle capitula?*

The [Japanese S1–S2–T2–T3 research programme](analysis/v4/RAIN_ORIENTATION_SPATIAL_INTRO_EVOLUTION_20261021_JA.md) and its [machine-readable evidence/stop-rule register](analysis/v4/RAIN_ORIENTATION_RESEARCH_DECISIONS_20261021.json) are the **current scientific scope**. This replaces the former broad A–E five-gate presentation as the *primary paper narrative*:

- **S1 Spatial — MEASURED AND ORIGINAL-OBSERVER ROBUSTNESS EXECUTED:** native five-climate BIO12 between-taxon **+0.4382 image-degrees/100 mm** (95% +0.2058 to +0.6706, Holm P=0.00531); [same-source original observer × taxon × cell single-sample sensitivity](analysis/v4/results/20261021_unified_abcde/S1_S2_RAIN_AND_INTRODUCED_SOURCE_RESULTS_JA.md) remains **+0.3885** (interval +0.1486 to +0.6284), while excluding top 1% reporting observers is **+0.3929**. Historical 52 grafting/placement estimates stay positive (+0.4722 to +0.5078) but are **not evolutionary transition counts**. **Rainfall seasonality BIO15 is NOT a second supported main effect**: S1 between-species +0.788 image-degrees per +10 BIO15 units, pointwise 95% interval **−0.266 to +1.842** (within-species −0.153, interval −0.527 to +0.220); both include zero. BIO12 annual rain is the registered primary axis; BIO15 remains exploratory. Camera pose and true head-level photometric/angular calibration remain open.
- **S2 Native → introduced — all 350 original species names screened, two satisfy source adequacy:** [Verified 350→7→2 original source evidence](analysis/v4/results/20261021_unified_abcde/S2_ALL_SPECIES_RANGE_SUPPORT_RESULT_JA.md) (Actions [39544849439](https://github.com/zuizui0223/azami/actions/runs/39544849439)) reports only **7/350 names** with both wild WCVP native and explicit introduced records; exactly **2** pass predeclared photo-orientation × observer × region-cell support. The other five were explicitly screened but lack enough independent introduced-range cells and/or photographs. **This does not mean only seven Cirsium species have introduced populations globally**, and cannot justify generalization from two species. The [executed climate transfer pilot](analysis/v4/results/20261021_unified_abcde/S2_FIRST_CLIMATE_TRANSFER_RESULT_JA.md) and [five-exposure overlap test](analysis/v4/results/20261021_unified_abcde/S2_FIVE_EXPOSURE_COMMON_SUPPORT_RESULT_JA.md) show the two photo-supported cases are not yet evolutionarily identified. The [source-backed climate-occupancy comparison](analysis/v4/results/20261021_unified_abcde/S2_CLIMATE_NICHE_NOVELTY_RESULT_JA.md) adds a further check: introduced photographs outside at least one native climate marginal 5–95% range were **73.4% vs 59.8%** by observer-unique photo records for arvense/vulgare, but only **39.2% vs 48.3% when 5° geographic cells are equally weighted**, close to native self-baseline **35.0% vs 37.2%**. Observational climate novelty is sampling-dependent and does **not** prove niche evolution, heritable angle changes or adaptation.

- **Prior S2 native→introduced pilot limitations:** [the executed climate analysis](analysis/v4/results/20261021_unified_abcde/S2_FIRST_CLIMATE_TRANSFER_RESULT_JA.md) found weak native-trained prediction skill on introduced image-angle observations. Critically, the [verified five-exposure support audit](analysis/v4/results/20261021_unified_abcde/S2_FIVE_EXPOSURE_COMMON_SUPPORT_RESULT_JA.md) (Actions 39543173995) retains only **600 C. arvense** and **1,198 C. vulgare** native records (2.74%/7.37% of observer-cell-balanced native measurements), while range-slope differences still include zero. The introduced-area original-pixel cache is useful mostly for **North America**. No tested evidence for genetic fixation or universal nontransfer follows.
- **S2 direct introduced-continent PHOTO orientation (new real result):** [source-locked continental angle and warm-region report](analysis/v4/results/20261021_unified_abcde/S2_VULGARE_CONTINENT_IMAGE_ANGLE_RESULT_JA.md), [real execution 39550930196](https://github.com/zuizui0223/azami/actions/runs/39550930196). For *C. vulgare*, image-vertical mean angle **66.61° native** vs **66.27° North America**, **66.26° Africa**, **66.03° South America**; cell-bootstrap differences from native are all zero-crossing. North American observations warmer than the native BIO1 q95 had **+3.94°** higher image-angle than other North American observations, but the cell-bootstrap exploratory range **−12.84 to +24.78° includes zero**. The real archived seasonal-rain BIO15 pixels cover only **72/2,376 Africa** and **0/2,075 South America** eligible photo observations, so increased BIO15 seasonal regions could **not** be linked to head angles there. **0° image-up / 90° horizontal / 180° image-down is not gravity-referenced nodding**; no mechanism, local adaptation or rapid evolution is inferred. The [machine receipt](analysis/v4/results/20261021_unified_abcde/S2_VULGARE_ACTUAL_CONTINENT_ORIENTATION_20261021.json) preserves original counts and intervals.
- **S2 direct wetness × actual photo angle — NEW SOURCE-EXECUTED MAIN TEST:** [same-source report and S1/BIO15 evidence boundary](analysis/v4/results/20261021_unified_abcde/S1_S2_BIO12_MAIN_BIO15_SECONDARY_RESULTS_JA.md) · [Actions 39559452415](https://github.com/zuizui0223/azami/actions/runs/39559452415). Frozen 173,737 source-measured photo angles (94,060 original observer×cell and valid date) plus SHA-verified original full CHELSA BIO12/BIO15 for native and introduced ranges. In North American introductions, above the native annual-rain q95, **arvense +0.21°** (576 obs/21 cells; 95% cell-bootstrap interval **−11.41 to +12.77°**) and **vulgare −9.60°** (809/19; **−22.05 to +1.41°**) relative to the other introduced climates. Local upper-vs-lower rain terciles give **+0.28° / +3.90°**, likewise zero-crossing. No clear increased nodding at wet introduced sites. Seasonal BIO15 is secondary and inconsistent across taxa/continents, not confirmed in S1. [SHA-locked receipt](analysis/v4/results/20261021_unified_abcde/S2_BIO12_PRIMARY_BIO15_SECONDARY_VERIFIED_20261021.json).
- **T2 Historical Cirsium transitions — INTAKE + CONDITIONAL Mk CODE IMPLEMENTED; EMPIRICAL NOT RUN:** [T2/T3 historical evidence contract](analysis/v4/T2_T3_HISTORY_CONTRACT_20261021_JA.md) and [audited T2/T3 runner](analysis/v4/rain_history.py) now enforce climate/tree-blind double botanical scoring, exact tip/concept maps, SHA-pinned independently reviewed nuclear trees and conditional symmetric-Mk stochastic mapping. An independent [test suite](tests/test_v4_rain_history.py) and empty [input templates](analysis/v4/rain_history_templates/README.md) are provided. **No real double-scored flora or eligible rooted Newick has been ingested; no actual transition number, event timing, rain/temperature evolutionary coefficient or causation has been inferred.** The existing 52 grafted placements cannot substitute.
- **T3 Asteraceae convergence — COMPARABILITY GATES + CONDITIONAL TRIBE-SPECIFIC Mk CODE IMPLEMENTED; EMPIRICAL NOT RUN:** independently reviewed anthesis head posture and cross-tribe homology must pass before a broader rooted nuclear species tree is analysed. The pilot can separately tally changes within single-tribe subtrees on an eligible tree, **but this is not proof of convergent rain adaptation**. No cross-tribe score/tree has yet met the input conditions.

The native preprocessing already established **313,350 observations / 306 accepted WCVP keys**, including 170,066 orientation-eligible observations (244 keys). **Seven accepted keys combine distinct source taxon names**; never interpret that pooled variation as one species without checking concepts. The response is based on **EXIF-oriented camera vertical**, not actual gravitational head tilt, rain exclusion, pollen protection, reproduction or genotype. Source-observer identities stay on authorised private runners.

**Controls and downstream mechanisms:** repeated images require confirmation of the *same head*, not just the same observation ID. Colour and outline can be image-specificity controls but **cannot by themselves rule out camera-angle bias**. A historical 5% detector bbox shift produced a p95 **54.1° orientation discrepancy**, and outline shows resolution sensitivity; these are from other data, not per-image V4 error models. A frozen 27-endpoint inventory and the selected 40-taxon native [KDE/W2 robustness analysis](analysis/v4/results/20261021_unified_abcde/NATIVE_DENSITY_ROBUSTNESS_RESULT_JA.md) remain documented as **exploratory Supporting Information candidates**, not the paper's main explanation. EAzami's East-Asian history and aza3's experimental mechanisms stay independently sourced.

The current [V4 native integration Actions](.github/workflows/v4-native-integration.yml) and frozen [52-placement](analysis/v4/results/20261020_native_phylogeny/README.md) plus [104 conditional covariance](analysis/v4/results/20261020_native_covariance/README.md) outputs are reusable computational components—not independent scientific papers. The [older five-gate exploratory plan](analysis/v4/UNIFIED_NATIVE_GATE_A_TO_E_JA.md) remains accessible **for audit/provenance only**. The V3/GEB eight-stage, sixteen-reference-output numerical replay and Zenodo source inputs are unchanged. Whether V4 replaces V3 as the publication is **conditional on validating the new claims**.

## Current evidence chain

```text
Frozen image-derived continuous measurements
    22 measured endpoints out of 27 registered
                    ↓
Biological constructs and environmental associations
    same 46,276-observation / 259-taxon source cohort
                    ↓
Sampling, spatial and phylogenetic-placement sensitivities
                    ↓
Whole-capitulum integration on the common complete-18 cohort
    1,734 observations / 42 taxa; 9 constructs / 36 relations
                    ↓
Taxon bootstrap, direct scale contrast and technical audits
```

The source includes native- and introduced-range observations. The current numerical workflow does not download images, retrain YOLO or revive categorical trait classifiers. YOLO localizes visible capitula; deterministic measurements define the continuous image phenotypes. Current construct analyses use the existing measurements, with the endpoint baseline retained separately.

Environmental associations are marginal analyses of nine predictors, not jointly adjusted effects of nine independent environmental causes. Spatial and phylogenetic-placement models are separate sensitivity analyses. Their passage is not a simultaneous spatial–phylogenetic correction or proof of adaptation. Module integration describes visible phenotypes, not genetic or evolutionary modularity.

## Repository boundaries

| Location | Role |
|---|---|
| `analysis/v3/` | Current construct-level numerical entry points and result receipts |
| `analysis/ch1/` | Shared frozen measurement and environmental contracts |
| `src/azami_ch1/` | Shared provenance and tabular utilities |
| `reproducibility/` | Input identities, current runner, validators and archive status |
| `analysis_outputs/` | Labelled frozen endpoint-level reference results |
| `reproducibility/figures/` | Labelled frozen v2 reference figures, not the current Main figure set |
| `legacy/` | Historical acquisition, measurement, endpoint analysis and figure implementations |

The current numerical runner reuses the retained spatial and historical-placement implementations under `legacy/v2/analysis`; those two explicit dependencies preserve the original calculations. Other legacy entry points are not current analysis instructions.

## Public availability

Zenodo DOI [10.5281/zenodo.22703215](https://doi.org/10.5281/zenodo.22703215) archives the **published v3 numerical inputs**. The existing `Reproduce current Chapter 1 analysis` Action and `python -m reproducibility.run_current_analysis --download` reuse this checksum-pinned package without Zenodo credentials or expiring input artifacts. It contains measurements, the nine-predictor environment, spatial lookup, historical trees and native status. Analysis code, runbooks, reference outputs and figures remain in GitHub. The [archive audit](reproducibility/ZENODO_UPDATE_AUDIT.md) records publication evidence.

This archive preserves the existing analysis cohort, not the complete 665,115-observation upstream source or original photographs. These archived numerical inputs neither include nor reprocess the original upstream source. Historical v2 inputs remain at DOI [10.5281/zenodo.22295791](https://doi.org/10.5281/zenodo.22295791).

Exact historical v2 reproduction remains available at commit `584af97b050d15701f26ce1facea212d5b648d4d`. The current scientific reference is main commit `fe25abd46e7235c85e6da191f976e1c8a02d0406`; subsequent layout-only changes must preserve its numerical results.
