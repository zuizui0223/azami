# V4 restarts at the original observations and images

The user's 2026-10-20 instruction supersedes the processed-CSV-only scope of the first PR142 pilot. `functional_pilot.py` remains a downstream diagnostic, NOT the upstream basis of final v4 results. Main/v3 inputs and claims are unchanged.

## Recovered, actually inspected upstream source

- SHA-verified historical pre-detector metadata archive 8066010557: **1,122,854 photo rows / 665,115 observations**, snapshot 2026-07-03; species-rank labels: 350. This is a historical snapshot, not the current full iNaturalist population.
- Its Actions copy is expired; the byte-identical owner Drive backup was retrieved and verified against the repository durable archive contract. Do not publish a signed temporary transfer URL.
- The original detector package was also recovered from its owner backup: ZIP SHA256 `2bc7cc49c2f213d4a5c5dab96e42fcdda7ceda23ffc3130520a5dd44b8c76c05`, frozen weight SHA256 `4078e0510532852681b65ee529cd82237b649ec99b17c4ca5f1da460a62d2bed`. This is a localization model, not independent anatomical truth.
- Some fields named `raw_image_url` actually point to **square thumbnails**. V4 constructs and explicitly records the provider's `original` URL, received dimensions and SHA256; it does not silently substitute large/medium thumbnails.

## Public diagnostic lane, not a new ecological sampling universe

The first replay restricts redistributable photographs to CC0/BY/BY-SA, species-rank noncaptive observations: 131,637 photo rows, 72,602 observations, 214 taxa. Twenty malformed historical image URLs are retained as unsupported rather than guessed. Obscured/private coordinates are removed; images without usable locations remain available for morphology, not environment extraction.

The planner scans the **whole** original metadata archive without referencing detector outcomes, v3 trait values, colour, or the old 46,276 subset. Its geographic-cell interleaving and 10-observation cap per taxon produce **1,409 queued photographs across 214 taxa**. This is a pilot queue, not the final sample cap. Licensing restriction itself changes the sample; it must not silently become the main ecological inclusion rule. The source archive remains the master frame for a broader rights-respecting research lane.

The committed first batch is **64 source photographs / 32 source taxa**, selected before looking at images. It includes C. vulgare/C. arvense as declared leverage sentinels; remaining taxa are selected by stable hash. No replacement is silently made when a photo disappears or its licence changes. Current public metadata, licences and attribution are checked before retrieving original bytes. Downloads are serial, delayed, bounded to 256 MiB per diagnostic batch, and follow iNaturalist recommended API/media limits.

## What is redefined at the measurement stage

1. **Colour:** anatomical floral ROI first; then direct sRGB-to-Lab L*, a*, b*, C* distributions, low-chroma fraction, hue reliability and actual observed RGB swatches. Do NOT reconstruct a*/b* from separately aggregated C*/hue. Do NOT choose floral pixels only because they are purple or saturated. A missing floral mask is missing evidence, never a whole-image colour measurement. ICC conversion and absent-profile assumptions are logged; these are not reflectance, UV or anthocyanin measurements.
2. **Presentation:** identifiable head-base/floral-end/peduncle landmarks yield image-axis and head-relative-to-peduncle angles separately. Gravity remains unknown without an external reference. Camera roll is not fixed by relabelling image vertical.
3. **Involucre:** explicit body mask, base/rim and identifiable phyllary base/tip landmarks are retained separately. No roughness-to-spine, glare-to-mucilage or texture-to-hair conversion is permitted. Anatomical measurements require annotations; masks and uncertainty must precede functional/ecological interpretation.

## Actual entry points and status

- `python -m analysis.v4.upstream --metadata-zip <verified archive> --out-dir <plan>` rebuilds the complete metadata frame and queue. This scan has been run locally.
- `python -m analysis.v4.run_upstream_acquisition --manifest analysis/v4/upstream_first_batch.json --out-dir work/v4-upstream` refreshes public source metadata and fetches original source photos. GitHub CI executes this acquisition, not the frozen v3 measurement replay.
- Seven synthetic tests cover deterministic sampling, licence/privacy separation, original URL routing, white-pixel retention, malformed source URLs and angle semantics.
- Acquisition creates blank anatomical annotation records and real-photo contact sheets. No validated floral masks or botanical accuracy are invented by a green acquisition job.

New environment estimates and within/among covariances must later use the newly measured observation/photo/head identifiers, not inherit v3 complete-case selection. Environmental importance will compare joint, shared and unique out-of-region predictive information; phylogeny remains a sensitivity lane. The current acquisition batch does not estimate macroecological effects.

References: iNaturalist API Recommended Practices (https://www.inaturalist.org/pages/api+recommended+practices); Troscianko & Stevens 2015, doi:10.1111/2041-210X.12439. The latter motivates distinguishing camera colour from calibrated reflectance.
