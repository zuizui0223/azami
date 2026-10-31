# SoilGrids 28→34 RF: native-window reader (STAGING ONLY)

This optional experiment is **not active in PR #167**. Preserve the running, private source-real scalar-sample workflow [GitHub Actions #41198399573](https://github.com/zuizui0223/azami/actions/runs/41198399573), which must finish or fail before attempting a new 28→34 fit. Never merge or retarget this fallback into its active PR while it could cancel that run.

## What changes / stays identical

- Changes only the GDAL reading implementation: original-pixel 256×256 disjoint windows, sorted by native raster indices. No resampling, smoothing, inferred nearest source, scale change, imputation, changed split, changed taxon gate or photo coordinate selection.
- Keeps original **313,350** S1 source SHA/AESGCM contract, the 16 CHELSA/WorldClim source raster contract, the derived 5 DEM variables, the original SoilGrids six Q0.5 0–5 cm VRT SHA pins, the fixed 1° / 3-photo-per-cell sample, <=14 original target cells per taxon and <=3,000 private positions total, the >=2/3 within-cell six-layer soil complete gate, >=25 taxa target, two unchanged out-of-group evaluation schemes and 80 trees.
- Keeps all original GEB S1/S2 Holm families, H3c, published manuscript and previous 12/21/23/28 results untouched. No source-photo ID, taxon name, individual location, soil pixels, RF predictions or model weights are exported.
- This **does not** solve source tile byte integrity: master VRT hashes are pinned; underlying remote TIFF tile content is not individually byte-archived. Use the previous source authenticity description, not a stronger claim.

## Verification already completed on synthetic data

- 2,301 synthetic original-pixel probes with a masked-source point produce the same values as `rasterio.sample`, NaNs included. Invalid/mismatched geographic arrays are rejected before any network call.
- Local `pytest` suite (existing native-window tests + source-34 tests + optional reader tests): **13/13 passed**.
- A separate 2,300-point tiled GeoTIFF microbenchmark showed identical pixels, scalar sampling 0.109 s vs window sampling 0.008 s (about 13×). **These are in-memory synthetic timings and not remotely measured SoilGrids/WebDAV throughput.** They must not be cited as empirical wall-clock improvement.
- PR #167 public sentinel test in [Actions #41198399615](https://github.com/zuizui0223/azami/actions/runs/41198399615) established exact equality at four fixed public sentinel positions (three finite pH pixels and one masked). This does not establish worldwide S1 coverage.

## Activation gate (not yet passed)

Only when #41198399573 is terminal and its aggregate, failure log or resource-limit issue is inspected, one may add the optional reader and its explicit provenance flag to the protected source workflow in the **original** PR. The new run must keep the original selection seed and validate source pixel equality and sample/cell/taxon support. Any causal inference or strong claim about 34 environmental explanatory variables remains prohibited without held-out predictive skill and uncertainty checks.