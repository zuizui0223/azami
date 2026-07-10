# azami — global image-derived capitulum traits in *Cirsium*

This repository contains code for a multi-chapter research project on the
ecology and evolution of *Cirsium* floral architecture. Chapter 1 is currently
the submission-focused component: a global macroecological analysis of
continuous capitulum traits measured from public photographs.

## Chapter 1 in one sentence

We use a capitulum detector and deterministic image measurements to quantify
corolla colour, capitulum outline and image-referenced head orientation, test
measurement assessability, join the traits to CHELSA climate, and evaluate
species-level associations under alternative historical-tree hypotheses.

The analysis is descriptive and non-causal. Species-level climate associations
are not presented as evidence of local adaptation, plasticity or selection.

## Current executed dataset

The frozen analysis used:

- 6,626 detected capitula;
- 3,725 iNaturalist observations;
- 216 accepted analysis taxa;
- nine primary continuous endpoints;
- three exploratory involucre/spine image proxies.

QC-retained head-level measurements were:

| Trait group | Usable heads | Retention |
|---|---:|---:|
| Corolla colour | 5,777 | 87.2% |
| Capitulum outline | 5,324 | 80.4% |
| Head orientation | 4,585 | 69.2% |

These percentages are assessability/QC-retention rates, not classification
accuracy.

The high-resolution involucre supplement retained 1,443 of 1,819 selected
capitula, representing 1,292 observations and 210 taxa.

## Main interpretation

Several between-species trait–climate associations were stable across 50 random
within-genus grafting trees. In contrast, the positional-accuracy <=10 km
within-species models yielded no FDR-significant effect for the nine primary
endpoints. The paper therefore distinguishes global species turnover and
macroevolutionary/geographic association from within-species environmental
response.

The dated vascular-plant backbone directly contained only 54 of 216 taxa. The
remaining taxa were grafted within *Cirsium*, so PGLS is used as a historical
constraint sensitivity analysis rather than as a uniquely correct species tree.

## Canonical Chapter 1 pipeline

```text
public photographs + metadata
        |
capitulum detection and crop reconstruction
        |
continuous primary traits (52–56)
        |
QC-retention audit + CHELSA join (57)
        |
high-resolution involucre supplement (58–60)
        |
trait integration tables and pre-tree models (61–64)
        |
historical-tree audit + bounded-lambda PGLS (65–68)
        |
submission validation and provenance manifest
```

The authoritative run instructions are in:

- `ch1_global/README.md`
- `ch1_global/v2/SUBMISSION_PIPELINE.md`
- `ch1_global/v2/CONTINUOUS_PRIMARY_TRAITS.md`
- `ch1_global/v2/HISTORICAL_CONSTRAINTS.md`

## Repository structure

```text
ch1_global/v1/       frozen orientation-only baseline
ch1_global/v2/       current data, measurement and analysis scripts
ch1_shared/          shared collection/training utilities
ch1_japan/           regional exploratory machinery, not part of the current
                     Chapter 1 headline analysis
ch3_trait_architecture/
protocols/
models/               detector/model files retained for provenance
.github/workflows/    lightweight CI plus manual research-scale workflows
tests/                deterministic unit and invariant tests
```

## Primary and auxiliary traits

### Primary

- orientation angle relative to EXIF-oriented image vertical;
- corolla Lab lightness and chroma;
- circular hue sine/cosine;
- capitulum aspect ratio, circularity, solidity and width-profile variation.

### Auxiliary

High-resolution image proxies for involucral projection roughness, spread and
spine-like protrusion. These are not categorical botanical claims such as
`recurved`, `unarmed` or `long-spined`.

Hair, mucilage, ploidy and hybrid status are kept in evidence queues. Missing
textual evidence is coded `unknown`, never biological absence.

## Reproducibility policy

- Heavy global workflows are manually dispatched because they consume external
  images and large intermediate artifacts.
- Pull requests run lightweight syntax, unit and submission-contract checks.
- Every submission bundle must include file SHA-256 checksums, table dimensions,
  software versions, the Git commit and analysis configuration.
- GitHub Actions artifacts are temporary working products and are not the final
  citable archive. The accepted-paper release must be deposited in a durable
  repository such as Zenodo.

## Project chapters

| Chapter | Role |
|---|---|
| Ch.1 | Global image-derived continuous traits and macroecological pattern |
| Ch.2 | Molecular phylogeny and ancestral-state questions |
| Ch.3 | Trait architecture and comparative evolutionary hypotheses |
| Ch.4 | Field and experimental tests of pollination, antagonism and abiotic protection |

## Important limitations

Citizen-science photographs are not colour-calibrated, image vertical is only a
proxy for gravity, and automated outline measures remain view dependent.
Hybridization, chloroplast capture and allopolyploid history cannot be resolved
by a bifurcating mega-tree. The code therefore preserves QC failures, separates
within- and between-species inference, and labels auxiliary results as
exploratory.
