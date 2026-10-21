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

## One unified native A→E research workflow (experimental V4; V3 frozen)

**[Scientific question and A→E execution contract](analysis/v4/UNIFIED_NATIVE_GATE_A_TO_E_JA.md)**

The native branch is **one study**, not independent Gate A, B, C, D and E papers or separately owned PRs. The [single V4 native integration Action](.github/workflows/v4-native-integration.yml) authenticates the same original 313,350-observation source and original human observer metadata, audits endpoint-specific QC and taxon concepts (**A**), compares matched within-/between-taxon covariance (**B**), estimates actual shared-space KDE density overlap and centre/shape-decomposed Wasserstein distances (**C**), compares the resulting distributions along native BIO12/RSDS regimes while reusing the existing fixed five-environment fitted model and spatial/phylogenetic placement sensitivities (**D**), and identifies taxon-concept matches and independently required functional/genomic experiments in EAzami/aza3 (**E**). All public workflow outputs contain **aggregate statistics only**; individual records, photographers, private source metadata and keys stay private.

The original pre-detector native-range roster is **313,350 observations / 306 accepted WCVP keys**, with **157,779 observations / 242 keys** eligible for the same nine primary image outputs. The endpoint inventory retains **27 registered outputs**, of which **22** had historical numeric measurements; those are not 27 independent anatomical or biological traits. The concept audit detected **seven accepted-key source-name collisions**, including the Ryukyu source labels *Cirsium brevicaule* and *C. irumtiense* on the same key. The unified registered sensitivity keeps mixed concepts **out of its balanced A→D cohort** and records them explicitly for future authority-backed resolution.

The new exploratory matched analysis predefines **one source observation per accepted taxon×5° spherical equal-area cell×original observer stratum**, then up to 40 eligible taxa, with shared 50-record/taxon rarefactions. The first unified numerical run is [Actions #39520165155](https://github.com/zuizui0223/azami/actions/runs/39520165155), which completed the connected sampling, covariance, KDE/W2 and native climate-distribution pilot. Rarefaction repetitions and technical camera/pigment validation are not yet publication-grade; photographic overlap does not demonstrate adaptation or genetic modularity. The latest source and final results are identified in its run-specific aggregate `unified_gate_abcde.json`.

The original whole-native joint five-environment results, conditional cell-spatial diagnostics, [52 historical placement sensitivities](analysis/v4/results/20261020_native_phylogeny/README.md) and [104 three-component conditional covariance fits](analysis/v4/results/20261020_native_covariance/README.md) remain part of this **same V4 research question**. The separate long numerical Actions that generated those frozen evidence objects are computational subjobs, **not independent scientific gates**. Their source and estimator limitations are retained. The frozen V3 GEB eight-stage/16-reference-output analysis and main manuscript claims remain untouched.

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
