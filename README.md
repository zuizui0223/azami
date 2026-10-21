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

- **S1 Spatial — EXISTING PATTERN:** globally native-range, image-vertical capitulum orientation is associated with among-taxon annual rain (five-climate model BIO12 **+0.4382 image-degrees per 100 mm**, 95% **+0.2058 to +0.6706**, Holm P=0.00531). All 52 historical **grafting/placement sensitivity** estimates for rain retain a positive sign (**+0.4722 to +0.5078** on 194 tree-covered taxa), but do not constitute 52 independent origins or a dated molecular reconstruction. Original human-observer effects, camera pose and independent physical-angle calibration remain essential checks.
- **S2 Native → introduced transfer — NOT RUN:** classify all 665,115 original observations using the original WCVP/TDWG3 **explicit introduced** state, separating unresolved or captive records rather than treating the 351,765 non-native-primary records as introductions. Focal candidates *C. vulgare* and *C. arvense* require native/introduced and photographer/region/phenotype support first. A flat slope in introduced range does **not** prove ancestral adaptive evolution or absence of plasticity.
- **T2 Historical Cirsium transitions — NOT RUN:** botanically scored head orientation independent of climate, exact taxon-concept links, and a vetted phylogeny with validated branch lengths. The existing 52 grafted placements are *sensitivity trees*, not sufficient for dating or counting independent state origins. Rain and temperature must be jointly controlled.
- **T3 Asteraceae convergence — NOT RUN:** if T2 and trait comparability hold, expand Cirsium → Cardueae → other Asteraceae lineages including *Cremanthodium*. Independent derived states need molecular/phylogenetic support, not similarity of photos or species labels.

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
