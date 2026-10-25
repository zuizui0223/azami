"""Aggregate-only audit for the native V3-question/V4-workflow consolidation."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

MODULES = {
    "orientation": ["orientation_image_vertical_angle"],
    "colour": ["corolla_lab_lightness", "corolla_lab_chroma",
               "corolla_hue_sin", "corolla_hue_cos"],
    "outline": ["capitulum_outline_aspect_ratio", "capitulum_outline_circularity",
                "capitulum_outline_solidity", "capitulum_width_profile_cv"],
    "bract_si": ["bract_projection_p95", "bract_spread_fraction"],
}


def eligible(frame, endpoints):
    mask = pd.Series(True, index=frame.index)
    for endpoint in endpoints:
        mask &= frame[endpoint+"__eligible"].eq(True) & np.isfinite(frame[endpoint])
    return mask


def environment_levels(frame, variables):
    values = frame[variables].replace([np.inf, -np.inf], np.nan)
    means = values.groupby(frame.accepted_key).transform("mean")
    return {
        "raw_observations": values,
        "between_equal_species": values.groupby(frame.accepted_key).mean(),
        "within_observations": values-means,
    }


def execute(table, receipt_path, output):
    receipt = json.loads(receipt_path.read_text())
    digest = hashlib.sha256(table.read_bytes()).hexdigest()
    if digest != receipt["prepared_table_sha256"]:
        raise ValueError("Decoded input checksum mismatch")
    if receipt.get("environment_decoding_status") != \
            "physical_units_gsp_invalid_encoding_masked_v1":
        raise ValueError("Corrected environmental decoding required")
    frame = pd.read_csv(table, low_memory=False)
    if len(frame) != receipt["native_observations"] or not frame.obs_id.is_unique:
        raise ValueError("Native observation count/identity mismatch")
    env = list(receipt["environment_decoding"])
    counts = {}
    for endpoint in receipt["endpoints"]:
        mask = eligible(frame, [endpoint])
        counts[endpoint] = {
            "finite": int(np.isfinite(frame[endpoint]).sum()),
            "eligible": int(mask.sum()),
            "species": int(frame.loc[mask, "accepted_key"].nunique())}
        if counts[endpoint]["eligible"] != receipt["endpoints"][endpoint]["eligible"]:
            raise ValueError("Endpoint eligibility differs from receipt")
    overlaps = {}
    names = list(MODULES)
    for i, left in enumerate(names):
        for right in names[i:]:
            mask = eligible(frame, list(dict.fromkeys(MODULES[left]+MODULES[right])))
            overlaps[left+"/"+right] = {
                "n": int(mask.sum()),
                "species": int(frame.loc[mask, "accepted_key"].nunique())}
    core = eligible(frame, MODULES["orientation"]+MODULES["colour"]+MODULES["outline"])
    environment = {}
    for level, data in environment_levels(frame, env).items():
        complete = data.dropna()
        corr = complete.corr()
        environment[level] = {
            "n_complete": len(complete), "variables": env,
            "correlation": json.loads(corr.to_json(orient="split")),
            "scope": "Complete nine-candidate rows; diagnostic, not a variable selection."}
    sampling = {}
    for field in ("latitude", "longitude", "observer_id", "user_id",
                  "photographer_id", "region", "tdwg3"):
        sampling[field] = {
            "available": field in frame,
            "nonmissing": int(frame[field].notna().sum()) if field in frame else 0}
    result = {
        "input_sha256": digest, "native_observations": len(frame),
        "native_species": int(frame.accepted_key.nunique()),
        "endpoint_coverage": counts, "module_overlap": overlaps,
        "core_complete": {"n": int(core.sum()),
                          "species": int(frame.loc[core, "accepted_key"].nunique())},
        "environment_diagnostics": environment, "sampling_fields": sampling,
        "scope": "Input gate only. No response-environment fits, variable selection, "
                 "thinning or adaptation tests. Equal-species between correlations "
                 "are not observation-weighted mixed-model design diagnostics."}
    output.mkdir(parents=True, exist_ok=True)
    (output/"native_integration_audit.json").write_text(
        json.dumps(result, indent=2, allow_nan=False)+"\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    for key in ("table", "receipt", "output"):
        parser.add_argument("--"+key, type=Path, required=True)
    args = parser.parse_args()
    execute(args.table, args.receipt, args.output)
