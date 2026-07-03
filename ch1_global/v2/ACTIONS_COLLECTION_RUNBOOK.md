# Chapter 1: GitHub Actions collection workflow

## Why the original 500-record run stops here

The first Actions run proved that the iNaturalist collector works. Its 500 observations were ordered from the earliest IDs, so they are an API/pipeline test rather than an analysis dataset. Do not build the final ecological sample from that one artifact.

## 1. Accumulate metadata in large non-overlapping chunks

Open **Actions → Ch1 iNaturalist metadata chunk → Run workflow**.

For the next chunk after the completed 500-observation pilot, enter:

| input | value |
|---|---:|
| `start_obs_id` | `412326` |
| `max_observations` | `25000` |
| `chunk_label` | `id_000412326` |

After completion, read `last_obs_id` in the artifact's `collection_summary.csv`. That exact number becomes `start_obs_id` for the following chunk. Give that next run a matching new label such as `id_012345678`.

Each run stores a separate 90-day artifact containing:

- `photo_metadata.csv`
- `collection_summary.csv` and `checkpoint.json` — includes the next cursor
- `collection_provenance.json` — taxon/API parameters and start/end state
- QC tables

Never overwrite a chunk or reuse its label. This makes the global metadata inventory auditable and prevents accidental gaps or duplicate ranges.

## 2. Quantity target before image download

Start with **10–15 chunks** (250,000–375,000 observations), then merge their metadata using `10_merge_metadata_chunks.py` and inspect taxonomic and spatial coverage. Continue toward the complete inventory only if rare species, southern-hemisphere areas, East Asia, Europe, or tropical mountains are still poorly represented.

Do not create final image queues separately for individual ID chunks. Pool and deduplicate the metadata first, then construct one global species × spatially balanced candidate queue.

## 3. Use YOLO for flower/head presence and visible-head count

After downloading a balanced image queue, run detector-only script `07_screen_downloaded_images_with_yolo.py`. Then run:

```powershell
python ch1_global\v2\09_build_yolo_positive_pool.py `
  --queue "...\image_screening_queue.csv" `
  --screening "...\screening_results.csv" `
  --crops "...\yolo_crop_metadata.csv" `
  --out-dir "...\yolo_positive_pool"
```

This makes:

- `yolo_image_summary.csv`: whether each photo has a detected head and how many heads are visible
- `yolo_positive_heads.csv`: one row per detected head crop for blind trait annotation

`visible_head_count_yolo` is photo-level composition data. It is useful for screening and for distinguishing single-head versus multi-head photos, but it is not automatically a whole-plant floral-display estimate because camera framing and occlusion differ across observations.

## 4. Running YOLO in GitHub Actions

GitHub-hosted runners can download the queue images and run YOLO, but the repository does not contain the detector weight file (`best.pt`). Before enabling an Actions YOLO workflow, place an immutable copy of the verified detector weights in a controlled download location such as a GitHub Release asset, and record its SHA-256 and training-label definition. The workflow must use that fixed version, never an untracked local `best.pt`.
