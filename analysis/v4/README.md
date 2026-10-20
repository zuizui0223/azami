# Chapter 1 v4: restart from upstream photographs

**Current scope follows the user's 2026-10-20 correction: begin with the original metadata and photographs, not the precomputed v3 trait table.** Main and frozen v3/GEB results remain untouched. PR142 stays exploratory/Draft.

Start with [UPSTREAM.md](UPSTREAM.md) and [UPSTREAM_RESULT_20261020.md](UPSTREAM_RESULT_20261020.md). The earlier `functional_pilot.py` and its numerical CI remain a **processed-CSV exploratory checkpoint**, not the source of new v4 biological measurements or a final environmental ranking.

## Active sequence

1. Rebuild the pre-detector acquisition frame from SHA-verified historical metadata: `upstream.py`. It does not inherit the old 46,276-observation cohort or filter by predicted colour/orientation.
2. Refresh photo attachment/licence/taxonomy and acquire actual original-size bytes: `run_upstream_acquisition.py`. Preserve failures, source attribution, hashes and ICC/EXIF handling.
3. Re-localize visible heads with the SHA-locked frozen detector: `run_upstream_detection.py`. A box is not an anatomical mask or a precision/recall evaluation.
4. Define floral tissue patches, head/peduncle landmarks and involucre masks; measure original pixels using `measure_upstream_anatomy.py` and `anatomical_geometry.py`. The committed five-head packet is **assistant_visual_proposal**, not independent human ground truth.
5. Build new environmental cohorts only after measurement assessment. Evaluate joint environmental information, widespread-species influence and within/among covariance; phylogeny remains sensitivity.

The parallel `native_scope.py`, `collect_native_sources.py` and native-first workflow are preserved. Their source-backed range rules must be joined before the native-primary ecological analysis. The public-license diagnostic photo batch is not itself a native-primary sample or a full macroecological cohort.

## Original-byte replay

`prepare_upstream_replay.py` verifies acquisition artifact 12241906850, not any v3 phenotype CSV. It preserves a monochrome source for geometry but prevents it from being counted as a white flower. The anatomical replay workflow remeasures five annotated examples from those exact original bytes.

Code success does not establish botanical accuracy, reflectance, anthocyanin quantity, gravity direction, defence function or adaptation. No new macroecological conclusion is asserted from this diagnostic batch.
