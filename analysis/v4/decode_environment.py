"""Decode the verified raw-value observation archive into physical environments."""
import argparse
import gzip
import hashlib
import io
import json
from pathlib import Path

import numpy as np
import pandas as pd

RAW_INPUT_SHA = "8ceed49e4ac73b31f9a5c726ca84f4198bfc3c51e9ee43303e29524e915a5c29"
GSP_SOURCE_SHA = "342858ce86f8b5c4f088327e4fa68afc7ac9cb64eff46ffaa980577a0f0641b6"
FOUR = ["chelsa_bio01", "chelsa_rsds_mean", "chelsa_vpd_mean", "chelsa_bio12"]


def decode_column(values, metadata, column):
    values = np.asarray(values, dtype=float).copy()
    invalid = ~np.isfinite(values)
    if metadata["nodata"] is not None:
        invalid |= values == metadata["nodata"]
    extra = np.zeros(len(values), dtype=bool)
    if column == "chelsa_gsp":
        if metadata["source_sha256"] != GSP_SOURCE_SHA or metadata["dtype"] != "uint32":
            raise ValueError("GSP source changed: re-audit invalid encoding policy")
        # Source declares int32-max as nodata, but also contains uint32-max.
        # Retain the distinction in the receipt; do not claim both were declared.
        extra = values == 4294967295
        invalid |= extra
    scale, offset = float(metadata["scale"]), float(metadata["offset"])
    if not np.isfinite(scale) or scale <= 0 or not np.isfinite(offset):
        raise ValueError("Invalid environmental decoding metadata")
    values[invalid] = np.nan
    return values * scale + offset, int(extra.sum())


def execute(table, receipt_path, metadata_path, output):
    original_sha = hashlib.sha256(table.read_bytes()).hexdigest()
    receipt = json.loads(receipt_path.read_text())
    if original_sha != RAW_INPUT_SHA or receipt["prepared_table_sha256"] != original_sha:
        raise ValueError("Expected immutable raw archive; refuse repeat or unverified decoding")
    metadata = json.loads(metadata_path.read_text())
    if len(metadata) != 9 or not set(FOUR) <= set(metadata):
        raise ValueError("Expected all nine verified environmental contracts")
    frame = pd.read_csv(table, low_memory=False)
    if len(frame) != receipt["native_observations"] or not frame.obs_id.is_unique:
        raise ValueError("Native population identity changed")
    decoding = {}
    for column, spec in metadata.items():
        raw = frame[column].to_numpy(float)
        decoded, invalid_encoded = decode_column(raw, spec, column)
        frame[column] = decoded
        decoding[column] = {
            "scale": spec["scale"], "offset": spec["offset"],
            "source_sha256": spec["source_sha256"],
            "declared_nodata": spec["nodata"],
            "finite_before": int(np.isfinite(raw).sum()),
            "finite_after": int(np.isfinite(decoded).sum()),
            "additional_invalid_uint32_max": invalid_encoded,
        }
        for suffix in ("__between", "__within"):
            stale = column + suffix
            if stale in frame:
                frame.drop(columns=stale, inplace=True)
    complete_four = np.isfinite(frame[FOUR]).all(axis=1)
    receipt["complete_four_climate"] = int(complete_four.sum())
    for endpoint, counts in receipt["endpoints"].items():
        counts["eligible_with_four_climate"] = int(
            (frame[endpoint+"__eligible"].eq(True) & np.isfinite(frame[endpoint])
             & complete_four).sum())
    output.mkdir(parents=True, exist_ok=True)
    target = output / "observations.csv.gz"
    with target.open("wb") as handle:
        with gzip.GzipFile(fileobj=handle, mode="wb", filename="", mtime=0) as zipped:
            with io.TextIOWrapper(zipped, encoding="utf-8", newline="") as text:
                frame.to_csv(text, index=False)
    receipt["raw_prepared_table_sha256"] = original_sha
    receipt["prepared_table_sha256"] = hashlib.sha256(target.read_bytes()).hexdigest()
    receipt["environment_decoding"] = decoding
    receipt["environment_metadata_sha256"] = hashlib.sha256(metadata_path.read_bytes()).hexdigest()
    receipt["complete_nine_climate"] = int(np.isfinite(frame[list(metadata)]).all(axis=1).sum())
    receipt["environment_decoding_status"] = "physical_units_gsp_invalid_encoding_masked_v1"
    (output / "measurement_receipt.json").write_text(json.dumps(receipt, indent=2)+"\n")


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    for name in ("table", "receipt", "metadata", "output"):
        p.add_argument("--"+name, type=Path, required=True)
    args = p.parse_args()
    execute(args.table, args.receipt, args.metadata, args.output)
