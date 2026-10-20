"""Outcome-independent multi-level selection under ENVIRONMENT_CONTRACT.md."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

CANDIDATES = ["chelsa_bio01", "chelsa_bio04", "chelsa_bio12", "chelsa_bio15",
              "chelsa_rsds_mean", "chelsa_vpd_mean"]
PROTECTED = ["chelsa_bio12", "chelsa_rsds_mean", "chelsa_bio01"]
REMOVE_FIRST = ["chelsa_vpd_mean", "chelsa_bio04", "chelsa_bio15"]


def make_levels(frame):
    values = frame[CANDIDATES].replace([np.inf, -np.inf], np.nan)
    means = values.groupby(frame.accepted_key).transform("mean")
    complete = values.notna().all(axis=1)
    # All three levels refer to species represented in the fixed complete cohort.
    between = means.loc[complete].groupby(frame.loc[complete, "accepted_key"]).first()
    return {"raw": values.loc[complete], "between": between,
            "within": (values-means).loc[complete]}


def diagnose(levels, selected):
    result = {}
    for name, frame in levels.items():
        values = frame[selected].to_numpy(float)
        sd = values.std(axis=0)
        if len(values) <= len(selected)+1 or not np.isfinite(values).all() or \
                np.any(sd <= 1e-12):
            raise ValueError("Insufficient finite, variable environment support: "+name)
        z = (values-values.mean(axis=0))/sd
        corr = z.T@z/len(z)
        singular = np.linalg.eigvalsh(corr)[0] <= 1e-10
        vif = None if singular else np.diag(np.linalg.inv(corr)).tolist()
        result[name] = {"n": len(values), "correlation": corr.tolist(),
                        "vif": dict(zip(selected, vif)) if vif else None,
                        "singular": bool(singular)}
    return result


def select(levels):
    selected, history = list(CANDIDATES), []
    while True:
        diagnostics = diagnose(levels, selected)
        conflicts = []
        for level, diagnostic in diagnostics.items():
            corr = diagnostic["correlation"]
            for i, left in enumerate(selected):
                for j in range(i+1, len(selected)):
                    if abs(corr[i][j]) > .8:
                        conflicts.append({"level": level, "left": left,
                                          "right": selected[j], "r": corr[i][j]})
        entry = {"variables": selected.copy(), "diagnostics": diagnostics,
                 "correlation_conflicts": conflicts}
        candidates = {v for pair in conflicts for v in (pair["left"], pair["right"])
                      if v not in PROTECTED}
        if candidates:
            removed = next(v for v in REMOVE_FIRST if v in candidates)
            entry.update(removed=removed, reason="absolute_correlation_above_0.8")
        elif conflicts:
            entry["reason"] = "protected_exposures_correlated"
            history.append(entry)
            return {"status": "unresolved", "selected": selected, "history": history}
        elif any(d["singular"] for d in diagnostics.values()):
            entry["reason"] = "singular_design_requires_review"
            history.append(entry)
            return {"status": "unresolved", "selected": selected, "history": history}
        else:
            maximum = {v: max(d["vif"][v] for d in diagnostics.values()) for v in selected}
            if max(maximum.values()) <= 5:
                entry["reason"] = "diagnostic_thresholds_passed"
                history.append(entry)
                return {"status": "selected", "selected": selected, "history": history}
            removable = [v for v in REMOVE_FIRST if v in selected]
            if not removable:
                entry["reason"] = "protected_exposures_vif_above_5"
                history.append(entry)
                return {"status": "unresolved", "selected": selected, "history": history}
            removed = max(removable, key=lambda v: maximum[v])
            entry.update(removed=removed, reason="multilevel_vif_above_5")
        history.append(entry)
        selected.remove(removed)


def execute(table, receipt_path, output):
    receipt = json.loads(receipt_path.read_text())
    digest = hashlib.sha256(table.read_bytes()).hexdigest()
    if digest != receipt["prepared_table_sha256"] or receipt.get(
            "environment_decoding_status") != "physical_units_gsp_invalid_encoding_masked_v1":
        raise ValueError("Verified corrected environmental input required")
    # Do not load response columns or response eligibility into variable selection.
    frame = pd.read_csv(table, usecols=["obs_id", "accepted_key", *CANDIDATES])
    if not frame.obs_id.is_unique or len(frame) != receipt["native_observations"]:
        raise ValueError("Population identity/count mismatch")
    result = select(make_levels(frame))
    result.update(input_sha256=digest, candidates=CANDIDATES, protected=PROTECTED,
                  contract_sha256=hashlib.sha256(
                      Path("analysis/v4/ENVIRONMENT_CONTRACT.md").read_bytes()).hexdigest(),
                  method="predeclared_multilevel_correlation_then_vif",
                  response_access=False)
    result["scales"] = {
        v: {"mean": float(frame[v].mean()), "sd": float(frame[v].std(ddof=0))}
        for v in result["selected"]}
    output.mkdir(parents=True, exist_ok=True)
    (output/"environment_selection.json").write_text(
        json.dumps(result, indent=2, allow_nan=False)+"\n")
    if result["status"] != "selected":
        raise SystemExit("Environment selection unresolved; no trait fits authorized by this gate.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    for key in ("table", "receipt", "output"):
        parser.add_argument("--"+key, type=Path, required=True)
    args = parser.parse_args()
    execute(args.table, args.receipt, args.output)
