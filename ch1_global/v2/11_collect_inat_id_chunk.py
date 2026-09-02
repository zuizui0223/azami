#!/usr/bin/env python3
"""Collect one bounded iNaturalist metadata chunk beginning after a known ID.

GitHub-hosted runners are ephemeral. Each collection remains an immutable,
non-overlapping ID range with its own checkpoint and provenance. Large final
sweeps can retain the raw audit trail as streamed gzip to limit runner disk use.
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Collect one non-overlapping iNaturalist observation-ID chunk.")
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--start-obs-id", type=int, required=True, help="Collect only observations whose ID is greater than this value.")
    parser.add_argument("--max-observations", type=int, default=25000)
    parser.add_argument("--taxon-name", default="Cirsium")
    parser.add_argument("--taxon-id", type=int, default=None)
    parser.add_argument("--quality-grade", choices=["any", "research", "needs_id", "casual"], default="any")
    parser.add_argument("--photo-license", default="any")
    parser.add_argument("--sleep-sec", type=float, default=0.8)
    parser.add_argument("--timeout-sec", type=int, default=90)
    parser.add_argument("--raw-compression", choices=["none", "gzip"], default="none")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.start_obs_id < 0:
        raise ValueError("--start-obs-id must be >= 0")
    if args.max_observations < 1:
        raise ValueError("--max-observations must be at least 1 for a bounded chunk")
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    protected = [
        "checkpoint.json", "photo_metadata.csv", "observation_raw.ndjson",
        "observation_raw.ndjson.gz", "collection_summary.csv", "collection_provenance.json",
    ]
    existing = [name for name in protected if (out_dir / name).exists()]
    if existing:
        raise FileExistsError(f"Chunk output directory must be fresh; found {existing} in {out_dir}")

    initial_state = {
        "last_obs_id": int(args.start_obs_id),
        "observations_seen": 0,
        "photos_written": 0,
        "pages_completed": 0,
    }
    (out_dir / "checkpoint.json").write_text(json.dumps(initial_state, ensure_ascii=False, indent=2), encoding="utf-8")

    collector = Path(__file__).with_name("04_collect_inat_cirsium_metadata.py")
    command = [
        sys.executable, str(collector),
        "--out-dir", str(out_dir),
        "--taxon-name", args.taxon_name,
        "--max-observations", str(args.max_observations),
        "--quality-grade", args.quality_grade,
        "--photo-license", args.photo_license,
        "--sleep-sec", str(args.sleep_sec),
        "--timeout-sec", str(args.timeout_sec),
        "--raw-compression", args.raw_compression,
    ]
    if args.taxon_id is not None:
        command.extend(["--taxon-id", str(args.taxon_id)])
    print("[INFO] starting bounded ID chunk", json.dumps({
        "start_obs_id": args.start_obs_id,
        "max_observations": args.max_observations,
        "raw_compression": args.raw_compression,
    }, ensure_ascii=False))
    subprocess.run(command, check=True)

    provenance_path = out_dir / "collection_provenance.json"
    provenance = json.loads(provenance_path.read_text(encoding="utf-8"))
    provenance["chunk_wrapper_version"] = "1.1.0"
    provenance["chunk_requested_start_obs_id"] = int(args.start_obs_id)
    provenance["resume_state_at_start"] = initial_state
    provenance["chunk_completed_at_utc"] = datetime.now(timezone.utc).isoformat()
    provenance_path.write_text(json.dumps(provenance, ensure_ascii=False, indent=2), encoding="utf-8")
    summary = out_dir / "collection_summary.csv"
    print("[OK] bounded ID chunk finished; read next cursor from", summary)


if __name__ == "__main__":
    main()
