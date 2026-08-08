# WCVP taxonomic review and synonym-collapse robustness

## Authority-backed candidate review

The frozen source-platform taxon names were cross-checked against the World Checklist of Vascular Plants (WCVP) through the GBIF checklist API. The candidate workflow does not silently rewrite source taxonomy.

Successful WCVP review run `31152400481` covered the union of 259 frozen names:

- 243 unique exact canonical-name matches;
- 16 multiple exact canonical-name matches, generally caused by accepted and homonymous/synonym records sharing an author-free canonical name;
- 235 routine accepted-name candidates;
- 8 WCVP synonym candidates;
- 24 high-priority nomenclatural review rows;
- 0 unmatched names.

The eight synonym candidates were:

- `Cirsium chanroenicum` -> `Cirsium setidens`;
- `Cirsium coryletorum` -> `Cirsium vlassovianum`;
- `Cirsium gilense` -> `Cirsium parryi`;
- `Cirsium hupehense` -> `Cirsium lineare`;
- `Cirsium irumtiense` -> `Cirsium brevicaule`;
- `Cirsium kawakamii` -> `Cirsium japonicum`;
- `Cirsium murdockii` -> `Cirsium tweedyi`;
- `Cirsium zhejiangense` -> `Cirsium yezoense`.

Each target already occurred as a separate source-platform analysis unit, so accepting every WCVP synonymy would merge active units.

## Explicit synonym-collapse sensitivity

Rather than silently imposing those merges, workflow run `31153385658` (artifact `8984168407`, digest `sha256:132613f326f929c1c5864aa9e4f9b7ec07f9ce928145108f908b5ced1bb6f712`) collapsed all eight conflicts simultaneously and recomputed the two headline result families most directly affected by lumping.

### Nested visible variance

The balanced atlas changed from 216 source-assigned taxa to 211 collapsed units. Across all nine primary endpoints, the below-unit visible-variance conclusion strengthened or changed only slightly:

- original minimum below-species fraction: 0.5886;
- collapsed minimum: 0.5970;
- largest absolute endpoint change: 0.0154.

No endpoint dropped below 0.5.

### Primary within-species climate coefficients

The exhaustive spatially thinned table changed from 259 source-assigned taxa to 251 collapsed units. Six rows were removed because the synonym collapse created duplicate accepted-taxon x 0.25-degree-cell records. The full 36 endpoint-component climate models were then recomputed from the frozen environmental values.

- original BH-supported component rows: 8;
- WCVP-collapse BH-supported component rows: 8;
- FDR decisions changed: 0;
- coefficient signs changed: 0;
- maximum absolute standardized-beta change: 0.000385.

Thus the headline nested-variance and within-species climate conclusions are insensitive to simultaneously applying all eight WCVP synonym collapses.

## Submission decision

The manuscript should distinguish **source-assigned taxonomic units** from a claim of a resolved genus-wide species taxonomy. This is appropriate because several *Cirsium* boundaries are actively treated differently among global checklists and regional taxonomic systems, and image records cannot retrospectively adjudicate those boundaries.

WCVP provides the authority-backed nomenclatural audit, while the collapse sensitivity demonstrates that the two main ecological conclusions are not artifacts of retaining the eight source-platform splits. If the final author team chooses to adopt WCVP names literally, the displayed taxon counts and downstream species-level products must be regenerated; the present sensitivity shows that the headline biological conclusions are expected to remain unchanged.
