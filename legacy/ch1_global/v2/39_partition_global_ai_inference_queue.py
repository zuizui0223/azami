#!/usr/bin/env python3
"""Partition a reproducible global photo-inference queue into balanced shards.

Every photo row is preserved exactly once. Multiple photos from the same
observation are allowed and deliberately distributed by stable photo-level keys.
Species are spread cyclically across shards without thinning the source queue.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import pandas as pd

REQUIRED = {"queue_id", "obs_id", "photo_id", "taxon_name"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Split an exhaustive photo queue into reproducible shards.")
    parser.add_argument("--queue", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--n-shards", type=int, default=20)
    parser.add_argument("--seed", default="20260711")
    return parser.parse_args()


def text(value: Any) -> str:
    return "" if pd.isna(value) else str(value).strip()


def stable_hash_int(seed: str, *parts: Any) -> int:
    value = "|".join([seed, *[text(part) for part in parts]])
    return int(hashlib.sha256(value.encode("utf-8")).hexdigest()[:16], 16)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> None:
    args = parse_args()
    if args.n_shards < 1 or args.n_shards > 100:
        raise ValueError("n-shards must be between 1 and 100")
    queue_path = Path(args.queue)
    if not queue_path.is_file():
        raise FileNotFoundError(queue_path)
    queue = pd.read_csv(queue_path, dtype=str, keep_default_na=False, low_memory=False)
    missing = sorted(REQUIRED.difference(queue.columns))
    if missing:
        raise ValueError(f"Queue missing columns: {missing}")
    if queue.empty or queue["queue_id"].duplicated().any() or queue["photo_id"].duplicated().any():
        raise ValueError("Queue must have nonempty unique queue_id and photo_id values")

    block_column = "spatial_block_2deg" if "spatial_block_2deg" in queue.columns else (
        "spatial_block" if "spatial_block" in queue.columns else "obs_id"
    )
    shard_sizes = [0] * args.n_shards
    shard_rows: list[pd.DataFrame] = []
    for species, group in queue.groupby("taxon_name", sort=True):
        work = group.copy()
        work["_hash"] = [
            stable_hash_int(args.seed, species, block, obs_id, photo_id, queue_id)
            for block, obs_id, photo_id, queue_id in zip(
                work[block_column], work["obs_id"], work["photo_id"], work["queue_id"]
            )
        ]
        work = work.sort_values([block_column, "_hash", "queue_id"]).reset_index(drop=True)
        offset = stable_hash_int(args.seed, species) % args.n_shards
        assigned: list[int] = []
        for index in range(len(work)):
            preferred = (offset + index) % args.n_shards
            minimum = min(shard_sizes)
            if shard_sizes[preferred] <= minimum + 1:
                shard = preferred
            else:
                candidates = [shard_id for shard_id, size in enumerate(shard_sizes) if size == minimum]
                shard = candidates[stable_hash_int(args.seed, species, index) % len(candidates)]
            assigned.append(shard + 1)
            shard_sizes[shard] += 1
        work["shard_id"] = assigned
        shard_rows.append(work.drop(columns=["_hash"]))

    partitioned = pd.concat(shard_rows, ignore_index=True)
    if len(partitioned) != len(queue) or partitioned["queue_id"].duplicated().any() or partitioned["photo_id"].duplicated().any():
        raise ValueError("Partitioning failed to preserve each photo row exactly once")

    out_dir = Path(args.out_dir)
    shard_dir = out_dir / "shards"
    shard_dir.mkdir(parents=True, exist_ok=True)
    manifest_rows = []
    for shard_id, group in partitioned.groupby("shard_id", sort=True):
        filename = f"global_ai_inference_queue_shard_{int(shard_id):03d}.csv"
        target = shard_dir / filename
        group.sort_values(["taxon_name", block_column, "obs_id", "photo_id"]).to_csv(
            target, index=False, encoding="utf-8-sig"
        )
        manifest_rows.append({
            "shard_id": int(shard_id),
            "filename": f"shards/{filename}",
            "n_photos": int(len(group)),
            "n_observations": int(group["obs_id"].nunique()),
            "n_species": int(group["taxon_name"].nunique()),
            "n_spatial_blocks": int(group[block_column].nunique()),
        })
    manifest = pd.DataFrame(manifest_rows).sort_values("shard_id")
    manifest.to_csv(out_dir / "global_ai_inference_shard_manifest.csv", index=False, encoding="utf-8-sig")
    provenance = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "queue_sha256": sha256_file(queue_path),
        "n_queue_photos": int(len(queue)),
        "n_queue_observations": int(queue["obs_id"].nunique()),
        "n_shards": int(args.n_shards),
        "seed": args.seed,
        "shard_sizes": manifest.to_dict("records"),
        "allocation_rule": "All photo rows retained; stable cyclic allocation within species with least-filled correction.",
    }
    (out_dir / "global_ai_inference_shard_provenance.json").write_text(
        json.dumps(provenance, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(json.dumps(provenance, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
