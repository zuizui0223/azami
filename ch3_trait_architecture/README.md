# Cirsium floral trait architecture

This chapter extends the `azami` programme from head orientation alone to a
predeclared, flower-centred trait architecture. Head orientation remains the
focal trait. The extension asks whether it is repeatedly associated with other
traits that alter attraction, pollinator access, floral antagonism, or exposure
to rain and wind.

## Scope

The unit of comparison is a *Cirsium* species, population, or individual,
depending on the analysis. Trait roles must not be hard-coded from correlation
alone. A trait can affect more than one route; for example, head orientation
may alter pollinator posture, floral-antagonist access, and rain exposure.

The initial flower-centred trait blocks are:

```text
A  attraction and pollinator access
D  floral barrier or resistance to floral antagonists
X  abiotic exposure and protection
O  measured interaction and fitness outcomes
```

Leaf resource traits, leaf defence, and leaf herbivory are outside this initial
module. They can be added only through an explicit cross-organ hypothesis.

## Data products

- `data/trait_architecture/cirsium_floral_trait_dictionary.csv` defines each
  candidate trait, measurement scale, candidate role, and evidence requirement.
- `data/trait_architecture/cirsium_field_panel_template.csv` defines the
  individual-level field-panel columns required to link traits to interactions
  and reproductive outcomes.
- `protocols/floral_trait_architecture_protocol.md` states the field and
  interpretation rules.

## Relationship to the project

- **Chapter 1** describes the global macroecological pattern of
  image-derived capitulum traits (orientation, colour, involucral cover,
  head shape, corolla:involucre display) from a standardized ROI on
  citizen-science photographs, and their correspondence to climate,
  terrain, and geography. It makes no adaptive or causal claim; see
  `ch1_global/README.md`.
- **Chapter 2** reconstructs gains and losses of orientation and, where data
  permit, correlated floral-trait evolution.
- **Chapter 3** tests whether orientation belongs to recurrent floral trait
  syndromes after accounting for phylogeny and environment.
- **Chapter 4** tests the candidate mechanisms directly: pollination, floral
  antagonism, rain/wind protection, and reproductive consequences.

No trait covariance is interpreted as adaptation without a stated mechanism and
an outcome measurement or an explicitly limited inference claim.

### Traits shared by name with Chapter 1

`involucral_cover` and `flower_colour` are measured in both Chapter 1 and
Chapter 3, but not at the same precision or for the same purpose, and results
from one should not be read as duplicating or validating the other:

- Chapter 1 measures them broadly, from uncontrolled citizen-science photos,
  across many species and no phylogenetic correction, to describe a global
  pattern (image-derived, error-tolerant).
- Chapter 3 remeasures them narrowly, from standardized photography /
  spectrometry under this protocol, at fewer sites, to test trait-syndrome
  structure after accounting for phylogeny.

Similarly, Chapter 1's `head_shape_aspect_ratio` (relative, no scale
reference) is distinct from this chapter's `capitulum_diameter` (calibrated,
absolute mm) and the two must not be merged into a single variable.
