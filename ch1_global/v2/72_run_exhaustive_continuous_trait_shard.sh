#!/usr/bin/env bash
# Run one exhaustive photo shard: download -> YOLO -> corrected v2 continuous traits.
set -euo pipefail

: "${GH_TOKEN:?GH_TOKEN is required}"
: "${GITHUB_REPOSITORY:?GITHUB_REPOSITORY is required}"
: "${QUEUE_SHARDS_ARTIFACT_ID:?QUEUE_SHARDS_ARTIFACT_ID is required}"
: "${DETECTOR_MODEL_ARTIFACT_ID:?DETECTOR_MODEL_ARTIFACT_ID is required}"
: "${SHARD_ID:?SHARD_ID is required}"
: "${OUTPUT_LABEL:?OUTPUT_LABEL is required}"

DOWNLOAD_SLEEP_SEC="${DOWNLOAD_SLEEP_SEC:-0.25}"
mkdir -p source/shards/extracted source/model/extracted
for specification in \
  "${QUEUE_SHARDS_ARTIFACT_ID}:source/shards/shards.zip:source/shards/extracted" \
  "${DETECTOR_MODEL_ARTIFACT_ID}:source/model/model.zip:source/model/extracted"; do
  IFS=':' read -r artifact_id zip_path extract_dir <<< "$specification"
  curl --fail --location --retry 4 --retry-all-errors --retry-delay 2 \
    --header "Authorization: Bearer ${GH_TOKEN}" \
    --header "Accept: application/vnd.github+json" \
    "https://api.github.com/repos/${GITHUB_REPOSITORY}/actions/artifacts/${artifact_id}/zip" \
    --output "$zip_path"
  unzip -tq "$zip_path"
  unzip -q "$zip_path" -d "$extract_dir"
done

WEIGHTS="$(find source/model/extracted -type f -name best.pt -print -quit)"
test -n "$WEIGHTS"

python - <<'PY'
import os
from pathlib import Path
import pandas as pd
root = Path("source/shards/extracted")
manifest = pd.read_csv(root / "global_ai_inference_shard_manifest.csv", dtype=str, keep_default_na=False)
wanted = str(int(os.environ["SHARD_ID"]))
match = manifest.loc[manifest["shard_id"].astype(str).eq(wanted)]
if len(match) != 1:
    raise SystemExit(f"shard_id={wanted} is absent or ambiguous")
queue = pd.read_csv(root / match.iloc[0]["filename"], dtype=str, keep_default_na=False)
if queue.empty or queue["queue_id"].duplicated().any() or queue["photo_id"].duplicated().any():
    raise SystemExit("Selected exhaustive photo shard is invalid")
queue.to_csv("queue_shard.csv", index=False, encoding="utf-8-sig")
print({"shard_id": int(wanted), "n_photos": len(queue), "n_observations": queue["obs_id"].nunique(), "n_species": queue["taxon_name"].nunique()})
PY

python ch1_global/v2/06_download_inat_screening_queue.py \
  --queue queue_shard.csv \
  --out-dir downloads \
  --image-size medium \
  --sleep-sec "$DOWNLOAD_SLEEP_SEC" \
  --timeout-sec 60

python - <<'PY'
import pandas as pd
queue = pd.read_csv("queue_shard.csv", dtype=str, keep_default_na=False)
results = pd.read_csv("downloads/download_results.csv", dtype=str, keep_default_na=False)
results["_order"] = range(len(results))
latest = results.sort_values("_order").groupby("queue_id", as_index=False).tail(1)
if set(latest["queue_id"]) != set(queue["queue_id"]):
    raise SystemExit("Download results do not cover every queued photo")
failed = latest.loc[~latest["download_status"].eq("success")]
if not failed.empty:
    print(failed[["queue_id", "download_status", "http_status", "message"]].head(20).to_dict("records"))
    raise SystemExit("Incomplete image downloads; rerun the shard rather than silently losing photos")
print({"n_download_success": len(latest)})
PY

python ch1_global/v2/07_screen_downloaded_images_with_yolo.py \
  --queue queue_shard.csv \
  --downloads downloads/download_results.csv \
  --images-dir downloads/images \
  --weights "$WEIGHTS" \
  --out-dir detector \
  --conf 0.25 \
  --pad-ratio 0.12 \
  --batch-size 16 \
  --device cpu

test -s detector/screening_results.csv

if [ ! -s detector/yolo_crop_metadata.csv ]; then
  python - <<'PY'
import json
import os
from pathlib import Path
import pandas as pd
queue = pd.read_csv("queue_shard.csv", dtype=str, keep_default_na=False)
screen = pd.read_csv("detector/screening_results.csv", dtype=str, keep_default_na=False)
if len(screen) != len(queue) or not screen["screen_status"].eq("no_detection").all():
    raise SystemExit("No crop table exists, but screening results are not complete no-detection rows")
report = {
    "shard_id": int(os.environ["SHARD_ID"]),
    "n_queue_photos": int(len(queue)),
    "n_queue_observations": int(queue["obs_id"].nunique()),
    "n_detected_photos": 0,
    "n_no_detection_photos": int(len(screen)),
    "n_detected_heads": 0,
    "n_colour_usable": 0,
    "n_shape_usable": 0,
    "n_orientation_usable": 0,
}
Path("shard_run_summary.json").write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
print(json.dumps(report, ensure_ascii=False, indent=2))
PY
  rm -rf downloads/images
  exit 0
fi

python ch1_global/v2/18_build_head_context_crops.py \
  --crops detector/yolo_crop_metadata.csv \
  --images-dir downloads/images \
  --out-dir context \
  --context-pad-ratio 1.5

python ch1_global/v2/28_build_detected_head_trait_packet.py \
  --yolo-crops detector/yolo_crop_metadata.csv \
  --context-metadata context/head_context_crop_metadata.csv \
  --source-images-dir downloads/images \
  --head-crops-dir detector/head_crops \
  --out-dir trait_packet_build \
  --batch-name "global_continuous_${OUTPUT_LABEL}_shard_${SHARD_ID}"

python ch1_global/v2/56_run_primary_traits_continuous_v2.py \
  --packet-manifest trait_packet_build/trait_annotation_packet/blinded_trait_annotation_manifest.csv \
  --packet-root trait_packet_build/trait_annotation_packet \
  --out-dir continuous \
  --colour-confidence-floor 0.55 \
  --shape-confidence-floor 0.50 \
  --orientation-confidence-floor 0.65 \
  --flip-angle-tolerance 20

python - <<'PY'
import json
import os
from pathlib import Path
import pandas as pd
queue = pd.read_csv("queue_shard.csv", dtype=str, keep_default_na=False)
screen = pd.read_csv("detector/screening_results.csv", dtype=str, keep_default_na=False)
traits = pd.read_csv("continuous/primary_trait_continuous_head_measurements.csv", low_memory=False)
if len(screen) != len(queue) or screen["queue_id"].duplicated().any():
    raise SystemExit("Screening does not cover every queued photo exactly once")
if traits["annotation_unit_id"].astype(str).duplicated().any():
    raise SystemExit("Continuous head measurements contain duplicate annotation units")
report = {
    "shard_id": int(os.environ["SHARD_ID"]),
    "n_queue_photos": int(len(queue)),
    "n_queue_observations": int(queue["obs_id"].nunique()),
    "n_detected_photos": int(screen["screen_status"].eq("detected").sum()),
    "n_no_detection_photos": int(screen["screen_status"].eq("no_detection").sum()),
    "n_detected_heads": int(len(traits)),
    "n_colour_usable": int(traits["colour_status"].eq("usable").sum()),
    "n_shape_usable": int(traits["shape_status"].eq("usable").sum()),
    "n_orientation_usable": int(traits["orientation_status"].eq("usable").sum()),
}
Path("shard_run_summary.json").write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
print(json.dumps(report, ensure_ascii=False, indent=2))
PY

rm -rf downloads/images detector/head_crops context/head_context_crops trait_packet_build/trait_annotation_packet/source_images trait_packet_build/trait_annotation_packet/head_crops trait_packet_build/trait_annotation_packet/context_crops
