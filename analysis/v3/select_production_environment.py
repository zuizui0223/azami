"""Run the saved outcome-blind candidate selector on a pinned completed matrix.

This is a candidate representation receipt, not a revision of the active model
contract and not permission to inspect or fit trait-environment associations.
"""
import argparse
import json
from pathlib import Path

import pandas as pd

from .production_environment import CONTRACT, choose_process_variables, validate_contract
from .protected_artifacts import new_json, require
from .workflow import canonical_digest, digest, text_digest


def select(matrix: Path, expected: str, out: Path, contract_path: Path = CONTRACT):
    require(not out.exists(), "Preserve previous selection receipt")
    require(len(expected) == 64 and digest(matrix) == expected, "Environment matrix identity differs")
    contract = json.loads(contract_path.read_text(encoding="utf-8"))
    validate_contract(contract)
    columns = ["obs_id", "accepted_key", *contract["acquire"]]
    frame = pd.read_csv(matrix, usecols=columns, dtype={"obs_id": str, "accepted_key": str})
    require(len(frame) == contract["source"]["observations"]
            and frame["accepted_key"].nunique() == contract["source"]["taxa"]
            and frame["obs_id"].notna().all() and frame["obs_id"].is_unique
            and frame["accepted_key"].notna().all(), "Environment source counts or identities differ")
    # The pinned GSP raster contains unmasked uint32 maximum values (scale 0.1),
    # distinct from its advertised nodata=2147483647. Do not let a finite-value
    # check promote these physically invalid amounts into VIF/model selection.
    invalid_gsp = int((frame["GSP"] == 429496729.5).sum())
    unresolved_bio12 = int((frame["BIO12"] == 6553.5).sum())
    if invalid_gsp or unresolved_bio12:
        result = {"status": "ENVIRONMENT_SELECTION_BLOCKED_BY_INVALID_SOURCE_VALUES",
                  "invalid_GSP_uint32_max_rows": invalid_gsp,
                  "unresolved_BIO12_storage_ceiling_rows": unresolved_bio12,
                  "source_rows": len(frame), "source_taxa": frame["accepted_key"].nunique(),
                  "selected": [], "trace": [], "trait_values_read": 0,
                  "ecological_fitting_authorized": False,
                  "reason": "Preserve the original matrix; use a separately documented QC view, not silent repair."}
    else:
        result = choose_process_variables(frame, contract)
    result.update(schema_version=1, matrix_sha256=expected,
                  contract_canonical_sha256=canonical_digest(contract),
                  implementation_sha256_text_lf=text_digest(Path(__file__)),
                  selector_sha256_text_lf=text_digest(Path(__file__).with_name("production_environment.py")),
                  candidate_only=True, active_ecological_contract_changed=False,
                  source_rows_deleted=0, ecological_models_executed=0,
                  limits=["Outcome-blind redundancy screening is not evidence of independent causal effects.",
                          "Actual within/among module designs and nuisance-adjusted identification still need diagnostics before fitting.",
                          "Removed columns remain in the complete saved matrix; their unique variation was not separately tested."])
    new_json(out, result)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--matrix", type=Path, required=True)
    parser.add_argument("--expected-matrix-sha256", required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(select(args.matrix, args.expected_matrix_sha256, args.out)))


if __name__ == "__main__":
    main()
