#!/usr/bin/env bash
# Execute one global automated Cirsium AI-inference shard in the current runner.
# Required environment: GH_TOKEN, GITHUB_REPOSITORY, QUEUE_SHARDS_ARTIFACT_ID,
# DETECTOR_MODEL_ARTIFACT_ID, SHARD_ID, MODEL_IDS_CSV, OUTPUT_LABEL.
set -euo pipefail

: "${GH_TOKEN:?GH_TOKEN is required}"
: "${GITHUB_REPOSITORY:?GITHUB_REPOSITORY is required}"
: "${QUEUE_SHARDS_ARTIFACT_ID:?QUEUE_SHARDS_ARTIFACT_ID is required}"
: "${DETECTOR_MODEL_ARTIFACT_ID:?DETECTOR_MODEL_ARTIFACT_ID is required}"
: "${SHARD_ID:?SHARD_ID is required}"
: "${MODEL_IDS_CSV:?MODEL_IDS_CSV is required}"
: "${OUTPUT_LABEL:?OUTPUT_LABEL is required}"

DOWNLOAD_SLEEP_SEC="${DOWNLOAD_SLEEP_SEC:-0.50}"

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

WEIGHTS="$(find source/model/extracted -type f -name best.pt | head -n 1)"
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
if queue.empty or queue["queue_id"].duplicated().any() or queue["obs_id"].duplicated().any():
    raise SystemExit("Selected queue shard is invalid")
queue["audit_id"] = queue["queue_id"]
queue.to_csv("queue_shard.csv", index=False, encoding="utf-8-sig")
print({"shard_id": int(wanted), "n_observations": len(queue), "n_species": queue["taxon_name"].nunique()})
PY

python ch1_global/v2/06_download_inat_screening_queue.py \
  --queue queue_shard.csv \
  --out-dir downloads \
  --image-size medium \
  --sleep-sec "$DOWNLOAD_SLEEP_SEC" \
  --timeout-sec 60

python - <<'PY'
from pathlib import Path
import pandas as pd
from PIL import Image

results = pd.read_csv("downloads/download_results.csv", dtype=str, keep_default_na=False)
latest = results.assign(_order=range(len(results))).sort_values("_order").groupby("queue_id", as_index=False).tail(1)
success = latest.loc[latest["download_status"].eq("success")].copy()
if len(success) != len(latest):
    failed = latest.loc[~latest["download_status"].eq("success"), ["queue_id", "download_status", "http_status", "message"]].head(10).to_dict("records")
    raise SystemExit(f"Incomplete image download; no partial shard inference is allowed: {failed}")
failures = []
for _, row in success.iterrows():
    path = Path(row["image_local_path"])
    try:
        with Image.open(path) as image:
            image.verify()
    except Exception as error:
        failures.append({"queue_id": row["queue_id"], "error": str(error)})
if failures:
    raise SystemExit(f"Undecodable downloaded images: {failures[:10]}")
print({"n_download_success": len(success)})
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

test -s detector/yolo_crop_metadata.csv
python ch1_global/v2/18_build_head_context_crops.py \
  --crops detector/yolo_crop_metadata.csv \
  --images-dir downloads/images \
  --out-dir context \
  --context-pad-ratio 1.5

python - <<'PY'
import pandas as pd
context = pd.read_csv("context/head_context_crop_metadata.csv", dtype=str, keep_default_na=False)
if context.empty or not context["context_crop_status"].eq("success").all():
    raise SystemExit("Every detector crop must have a successful context crop")
PY

python ch1_global/v2/28_build_detected_head_trait_packet.py \
  --yolo-crops detector/yolo_crop_metadata.csv \
  --context-metadata context/head_context_crop_metadata.csv \
  --source-images-dir downloads/images \
  --head-crops-dir detector/head_crops \
  --out-dir trait_packet_build \
  --batch-name "global_ai_${OUTPUT_LABEL}_shard_${SHARD_ID}"

model_args=()
IFS=',' read -r -a models <<< "$MODEL_IDS_CSV"
for model_id in "${models[@]}"; do
  model_args+=(--model-id "$model_id")
done
python ch1_global/v2/35_score_traits_ai_ensemble.py \
  --packet-manifest trait_packet_build/trait_annotation_packet/blinded_trait_annotation_manifest.csv \
  --packet-root trait_packet_build/trait_annotation_packet \
  --ontology ch1_global/v2/ontology/ch1_trait_ontology.csv \
  --prompt-spec ch1_global/v2/ontology/clip_trait_prompt_spec.json \
  --out-dir trait_ai_ensemble \
  "${model_args[@]}" \
  --device cpu \
  --batch-size 24 \
  --threads 2 \
  --seed 20260705 \
  --strict-margin-percentile 0.50 \
  --majority-vote-fraction 0.75

python - <<'PY'
import json
import os
from pathlib import Path
import pandas as pd

queue = pd.read_csv("queue_shard.csv", dtype=str, keep_default_na=False)
crops = pd.read_csv("detector/yolo_crop_metadata.csv", dtype=str, keep_default_na=False)
long = pd.read_csv("trait_ai_ensemble/ai_ensemble_trait_measurements_long.csv", dtype=str, keep_default_na=False)
variants = pd.read_csv("trait_ai_ensemble/ai_ensemble_variant_predictions.csv", dtype=str, keep_default_na=False)
provenance = json.loads(Path("trait_ai_ensemble/ai_ensemble_measurement_provenance.json").read_text(encoding="utf-8"))
n_heads = crops[["audit_id", "det_index"]].drop_duplicates().shape[0]
if n_heads < 1 or len(long) != n_heads * 8 or long.duplicated(["annotation_unit_id", "trait_id"]).any():
    raise SystemExit("Automated trait table does not cover each detected head and trait exactly once")
if len(variants) != n_heads * 7 * len(provenance["model_locks"]) * 2:
    raise SystemExit("Variant prediction row count is invalid")
report = {
    "shard_id": int(os.environ["SHARD_ID"]),
    "n_queue_observations": int(len(queue)),
    "n_detected_heads": int(n_heads),
    "n_ai_trait_rows": int(len(long)),
    "n_variant_predictions": int(len(variants)),
    "model_locks": provenance["model_locks"],
}
Path("shard_run_summary.json").write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
print(json.dumps(report, ensure_ascii=False, indent=2))
PY
