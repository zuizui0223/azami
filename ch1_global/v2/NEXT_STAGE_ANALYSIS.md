# Chapter 1 next-stage analysis

## Order of operations

1. Audit whether primary-trait QC retention varies among species, spatial blocks,
   latitude bands or CHELSA gradients.
2. Join continuous colour, outline and orientation measurements to the existing
   observation-level environment table without converting QC failure into
   biological absence.
3. Measure involucre projection and spine-like protrusion only on high-resolution
   detected heads as exploratory continuous proxies.
4. Keep indumentum and mucilage in a separate source-backed literature/macro-image
   evidence table.
5. Build QC and coverage figures, then fit environmental models.
6. Add historical/phylogenetic constraints only after the measurement system and
   complete-case sensitivity analyses are frozen.

## Primary traits

Main variables remain continuous:

- orientation angle relative to EXIF-oriented image vertical;
- corolla Lab lightness, chroma and circular hue components;
- capitulum aspect ratio, circularity, solidity and width-profile variation.

QC retention is analysed as its own binary outcome. A significant retention–
environment association is reported as differential assessability and triggers
complete-case/retention sensitivity analysis. It is not interpreted as a
biological trait–environment association.

## Auxiliary involucre analysis

The auxiliary image route begins with a minimum 150-pixel detector-box dimension.
The measured variables are contour proxies:

- involucre projection roughness;
- 95th percentile and maximum positive radial projection;
- fraction of involucre contour showing positive projection;
- number and relative maximum length of spine-like contour peaks.

They are deliberately not labelled `appressed`, `recurved`, `unarmed` or
`long_spined`. Horizontal mirroring is used as a technical replicate. Only
high-resolution, sharp, stable, well-segmented crops are retained.

## Hair and mucilage

`involucre_indumentum` and `external_mucilage_visible` are not inferred as absent
from ordinary photographs. The literature queue requires a source ID, page or
figure, and supporting quotation. An LLM may extract evidence and normalize
synonyms, but unsupported completion is prohibited. Mucilage absence requires an
explicit negative statement or reviewed macro observation of the relevant
surface.

## Historical constraint comes later

The world-scale tree is not treated as exact evolutionary history. After the
measurement and QC layers are frozen, the environmental results will be tested
against increasingly strong historical constraints:

1. species identity / held-out-species validation already present in the atlas;
2. broad clade or geographic lineage terms when defensible;
3. PGLS with Pagel's lambda as a sensitivity analysis;
4. alternative trees and exclusion of known hybrids/polyploids where data permit.

This ordering prevents a large but uncertain synthetic or chloroplast-dominated
tree from hiding measurement problems or creating false precision.
