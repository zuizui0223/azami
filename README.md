# Azami capitulum phenotypes

Current analysis and reproducible numerical evidence for the Chapter 1 study prepared for **Global Ecology and Biogeography (GEB)**. Manuscript text, submission files and private reviewer correspondence are kept outside this repository.

Original software is released under the [MIT license](LICENSE); see [NOTICE.md](NOTICE.md) for the separate treatment of images, data and third-party material.

## Start here

1. [Current numerical reproduction](reproducibility/CURRENT_ANALYSIS.md): verified inputs, one-command execution, and result validation.
2. [Current results and interpretation](analysis/v3/README.md): biological constructs, environmental associations and whole-capitulum integration.
3. [GEB scope freeze](analysis/v3/GEB_SCOPE_FREEZE_20260912.md): frozen claim boundary, submission-readiness work, analysis stop list, and separate evidence routes for any higher-tier successor.
4. [Current figure surface](reproducibility/CURRENT_FIGURE_SURFACE_20260912.md): Main-figure role freeze aligned to the v3 story and the boundary between the 46,276/259 atlas and 1,734/42 integration cohort.
5. [Reproducibility and archive status](reproducibility/README.md): what is publicly archived and what still needs a new Zenodo version.
6. [Historical implementations](legacy/README.md): v2 and superseded entry points. Do not start the current paper analysis from these scripts.

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

Zenodo DOI [10.5281/zenodo.22295791](https://doi.org/10.5281/zenodo.22295791) archives the **v2** minimum numerical input package. It is not a complete archive of the current construct-level analysis. The [archive audit](reproducibility/ZENODO_UPDATE_AUDIT.md) identifies the required new version. GitHub Actions artifacts expire; they are not a permanent substitute.

Exact historical v2 reproduction remains available at commit `584af97b050d15701f26ce1facea212d5b648d4d`. The current scientific reference is main commit `fe25abd46e7235c85e6da191f976e1c8a02d0406`; subsequent layout-only changes must preserve its numerical results.
