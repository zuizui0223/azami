#!/usr/bin/env python3
"""Offline numerical and sample-identity audit of the Chapter 1 revision.

No network, image retrieval, image transformation, detector or workflow execution.
Frozen v2 outputs are read-only. This is numerical verification on the same data,
not independent biological replication. Bbox summaries are never substituted for
linked, per-head measurements. Missing records leave that check NOT_VERIFIED.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import platform
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import norm

PREDICTORS = ["chelsa_bio01", "chelsa_bio04", "chelsa_bio12", "chelsa_bio15",
              "chelsa_rsds_mean", "chelsa_vpd_mean", "chelsa_sfcwind_mean",
              "chelsa_gsp", "chelsa_npp"]
PAIRS = [("corolla_lab_chroma", "chelsa_rsds_mean"),
         ("orientation_image_vertical_angle", "chelsa_bio12")]
HASHES = {
    "traits": "d775794f2bce2dfd0c1f63c5c8e01778c518f6eeb327bf0d9944045143a02344",
    "heads": "779e6f84830da56fb0d561f62042d087b5fff6b3511150c65727f0fff162a733",
    "environment": "e242aa7ce69d12b11937c1335e84b9638799c50b42ef36b95725e77190df98e7",
}
CONDITIONS = {"baseline", "x_minus_5pct", "x_plus_5pct", "y_minus_5pct", "y_plus_5pct"}


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def strict_bool(values: pd.Series) -> pd.Series:
    text = values.astype(str).str.strip().str.lower()
    if not text.isin({"true", "false", "1", "0"}).all():
        raise ValueError("Unknown or missing boolean; not silently treated as False")
    return text.isin({"true", "1"})


def zscore(values: np.ndarray) -> np.ndarray:
    a = np.asarray(values, float)
    sd = a.std(axis=0, ddof=0)
    if not np.isfinite(a).all() or np.any(sd <= 0):
        raise ValueError("Nonfinite or constant data")
    return (a - a.mean(axis=0)) / sd


def ols_hc3(y: np.ndarray, predictors: np.ndarray) -> dict:
    """Explicit HC3 sandwich and normal-reference 95% interval for column 1.

    Each response/predictor is standardized in its specified cohort. These are
    descriptive intervals, neither spatial/phylogenetic nor multiplicity adjusted.
    """
    raw = np.asarray(predictors, float)
    if raw.ndim == 1:
        raw = raw[:, None]
    design = np.column_stack([np.ones(len(raw)), zscore(raw)])
    response = zscore(np.asarray(y, float))
    if len(response) <= design.shape[1] or np.linalg.matrix_rank(design) != design.shape[1]:
        raise ValueError("Insufficient observations or rank-deficient design")
    inverse = np.linalg.pinv(design)
    beta = inverse @ response
    residual = response - design @ beta
    leverage = np.einsum("ij,ji->i", design, inverse)
    if np.any(1 - leverage <= 1e-12):
        raise ValueError("Unit leverage prevents HC3 estimation")
    scaled = inverse * (residual / (1 - leverage))[None, :]
    covariance = scaled @ scaled.T
    se = np.sqrt(np.maximum(np.diag(covariance), 0))
    width = float(norm.ppf(.975)) * float(se[1])
    return {"n_taxa": len(response), "beta": float(beta[1]), "hc3_se": float(se[1]),
            "ci_low": float(beta[1] - width), "ci_high": float(beta[1] + width),
            "design_rank": int(np.linalg.matrix_rank(design)),
            "condition_number": float(np.linalg.cond(design))}


def taxon_table(traits: pd.DataFrame, environment: pd.DataFrame,
                endpoint: str, minimum: int = 5) -> pd.DataFrame:
    g = traits.loc[traits.endpoint_id.eq(endpoint)].groupby("taxon_name").value.agg(["median", "count"])
    return g.loc[g["count"] >= minimum].join(environment, how="inner").sort_index()


def fixed_selection(heads: pd.DataFrame, eligible_taxa: list[str]) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Reconstruct the already stated selection; no image operation is performed."""
    usable = heads.loc[heads.orientation_status.eq("usable") &
                       np.isfinite(heads.orientation_angle_degrees) &
                       heads.taxon_name.isin(eligible_taxa)].copy()
    observations = usable[["taxon_name", "obs_id"]].drop_duplicates().copy()
    observations["selection_hash"] = observations.obs_id.map(
        lambda x: hashlib.sha256(("azami-v3-bbox-20260907:" + str(int(x))).encode()).hexdigest())
    selected = observations.sort_values(["taxon_name", "selection_hash", "obs_id"])
    selected = selected.groupby("taxon_name", sort=True).head(50).reset_index(drop=True)
    selected_heads = usable.loc[usable.obs_id.isin(selected.obs_id)].sort_values("annotation_unit_id")
    return selected, selected_heads.reset_index(drop=True)


def validate_linked_records(records: pd.DataFrame, scheduled: pd.DataFrame) -> pd.DataFrame:
    """Audit an explicit normalized schema; missing conditions or IDs fail closed.

    Failed executions need rows with usable=False and nonfinite angle_deg, not
    dropped records. This validator checks tabular integrity, not image authenticity.
    """
    required = {"annotation_unit_id", "obs_id", "photo_id", "condition", "angle_deg",
                "usable", "head_clipped", "source_sha256"}
    if not required.issubset(records.columns):
        raise ValueError("Missing columns: " + str(sorted(required - set(records.columns))))
    key = "annotation_unit_id"
    if scheduled[key].duplicated().any() or records.duplicated([key, "condition"]).any():
        raise ValueError("Duplicate scheduled head or head-condition")
    if set(records[key]) != set(scheduled[key]):
        raise ValueError("Missing or unexpected head IDs")
    if not records.groupby(key).condition.apply(lambda x: set(x) == CONDITIONS and len(x) == 5).all():
        raise ValueError("Every head needs exactly all five conditions, including failed ones")
    rows = records.copy()
    rows["usable"] = strict_bool(rows.usable)
    rows["head_clipped"] = strict_bool(rows.head_clipped)
    identity = scheduled.set_index(key)
    for column in ["obs_id", "photo_id"]:
        expected = rows[key].map(identity[column])
        if not rows[column].astype(str).eq(expected.astype(str)).all():
            raise ValueError("Source identity mismatch: " + column)
    valid_hash = rows.source_sha256.fillna("").astype(str).str.fullmatch(r"[0-9a-f]{64}")
    if (rows.usable & ~valid_hash).any():
        raise ValueError("Usable record has no SHA256 image identity")
    for _, group in rows.groupby(key):
        hashes = group.loc[group.source_sha256.fillna("").ne(""), "source_sha256"].unique()
        if len(hashes) > 1:
            raise ValueError("One head has multiple source image identities")
    rows["angle_deg"] = pd.to_numeric(rows.angle_deg, errors="raise")
    finite = np.isfinite(rows.angle_deg)
    if (rows.usable & ~finite).any() or (finite & ~rows.angle_deg.between(0, 180)).any():
        raise ValueError("Invalid image-referenced angle")
    angles = rows.pivot(index=key, columns="condition", values="angle_deg")
    qc = identity[["obs_id", "photo_id", "taxon_name", "orientation_angle_degrees"]].copy()
    qc["baseline_error"] = (angles.baseline - qc.orientation_angle_degrees).abs()
    qc["all_five_usable"] = rows.groupby(key).usable.all()
    qc["any_head_clipped"] = rows.groupby(key).head_clipped.any()
    qc["max_change"] = angles.sub(angles.baseline, axis=0).abs().max(axis=1)
    complete = angles.notna().all(axis=1)
    for threshold in (5, 10, 20):
        qc[f"stable_{threshold}"] = (complete & qc.baseline_error.le(1e-6) & qc.all_five_usable &
                                      ~qc.any_head_clipped & qc.max_change.le(threshold))
    return qc.reset_index()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    for role in HASHES:
        parser.add_argument("--" + role, type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--bbox-records", type=Path)
    args = parser.parse_args()
    paths = {role: getattr(args, role) for role in HASHES}
    receipts = {role: digest(path) for role, path in paths.items()}
    for role, expected in HASHES.items():
        if receipts[role] != expected:
            raise ValueError(role + " SHA256 mismatch")
    args.out_dir.mkdir(parents=True, exist_ok=True)
    env = pd.read_csv(args.environment)
    if env.obs_id.duplicated().any() or len(env) != 46276 or env.taxon_name.nunique() != 259:
        raise ValueError("Wrong primary environment cohort")
    env_medians = env.groupby("taxon_name")[PREDICTORS].median()
    columns = ["obs_id", "taxon_name", "endpoint_id", "value", "measurement_available"]
    traits = pd.read_csv(args.traits, usecols=columns)
    traits = traits.loc[strict_bool(traits.measurement_available) & np.isfinite(traits.value) &
                        traits.obs_id.isin(env.obs_id) & traits.endpoint_id.isin([p[0] for p in PAIRS])]
    if traits.duplicated(["obs_id", "endpoint_id"]).any():
        raise ValueError("Duplicate observation-endpoint")
    mapped_taxon = traits.obs_id.map(env.set_index("obs_id").taxon_name)
    if not traits.taxon_name.eq(mapped_taxon).all():
        raise ValueError("Trait/environment taxon mismatch")
    models, thresholds, cohorts, vifs = [], [], {}, []
    for endpoint, predictor in PAIRS:
        cohort = taxon_table(traits, env_medians, endpoint)
        cohorts[endpoint] = cohort
        core = [p for p in ["chelsa_bio01", "chelsa_bio12"] if p != predictor]
        scenarios = [("univariate", []), ("climate_core", core)]
        scenarios += [("plus_" + p, [p]) for p in PREDICTORS if p != predictor]
        scenarios += [("all_nine", [p for p in PREDICTORS if p != predictor])]
        for label, adjustment in scenarios:
            result = ols_hc3(cohort["median"].to_numpy(), cohort[[predictor] + adjustment].to_numpy())
            models.append(dict(endpoint=endpoint, predictor=predictor, model=label,
                               adjustment=";".join(adjustment), **result))
        for minimum in (2, 5, 10, 20, 30):
            g = taxon_table(traits, env_medians, endpoint, minimum)
            thresholds.append(dict(endpoint=endpoint, predictor=predictor, minimum=minimum,
                                   **ols_hc3(g["median"].to_numpy(), g[[predictor]].to_numpy())))
        zx = zscore(cohort[PREDICTORS].to_numpy())
        for i, p in enumerate(PREDICTORS):
            design = np.column_stack([np.ones(len(zx)), np.delete(zx, i, axis=1)])
            residual = zx[:, i] - design @ np.linalg.lstsq(design, zx[:, i], rcond=None)[0]
            vifs.append(dict(endpoint=endpoint, predictor=p,
                             vif=float(np.sum(zx[:, i]**2) / np.sum(residual**2))))
    orientation = cohorts[PAIRS[1][0]]
    head_columns = ["annotation_unit_id", "obs_id", "photo_id", "taxon_name", "orientation_status",
                    "orientation_angle_degrees"]
    heads = pd.read_csv(args.heads, usecols=head_columns)
    if heads.annotation_unit_id.duplicated().any():
        raise ValueError("Duplicate head identifier")
    if not heads.obs_id.isin(env.obs_id).all():
        raise ValueError("Heads outside primary cohort")
    usable = heads.loc[heads.orientation_status.eq("usable") & np.isfinite(heads.orientation_angle_degrees)]
    original = traits.loc[traits.endpoint_id.eq(PAIRS[1][0])].set_index("obs_id").value
    rebuilt = usable.groupby("obs_id").orientation_angle_degrees.median()
    if set(original.index) != set(rebuilt.index):
        raise ValueError("Observation median ID mismatch")
    error = float((original.sort_index() - rebuilt.sort_index()).abs().max())
    if error > 1e-10:
        raise ValueError("Original orientation medians do not reproduce")
    observations, selected_heads = fixed_selection(heads, list(orientation.index))
    counts = {"taxa": int(selected_heads.taxon_name.nunique()), "observations": len(observations),
              "heads": len(selected_heads), "photographs": int(selected_heads.photo_id.nunique())}
    if counts != {"taxa": 142, "observations": 3829, "heads": 10431, "photographs": 6561}:
        raise ValueError("Declared fixed subset does not reproduce: " + str(counts))
    pd.DataFrame(models).to_csv(args.out_dir / "conditional_models.csv", index=False)
    pd.DataFrame(thresholds).to_csv(args.out_dir / "minimum_observations.csv", index=False)
    pd.DataFrame(vifs).to_csv(args.out_dir / "candidate_vif.csv", index=False)
    observations.to_csv(args.out_dir / "reconstructed_selection_observations.csv", index=False)
    selected_heads.to_csv(args.out_dir / "reconstructed_selection_heads.csv", index=False)
    report = {"audit_id": "ch1_v3_offline_consistency_20260907", "status": "VERIFIED_AVAILABLE_SCOPE_ONLY",
              "input_sha256": receipts, "original_orientation": {"usable_heads": len(usable),
              "observations": len(rebuilt), "maximum_absolute_error_degrees": error},
              "fixed_subset": counts, "conditional_models": len(models), "threshold_models": len(thresholds),
              "interval": "HC3 sandwich; normal reference; 95%; exploratory; no multiplicity or spatial adjustment",
              "note_duplicate_model": "orientation climate_core and plus_chelsa_bio01 intentionally coincide",
              "image_operations_performed": False, "bbox_record_audit": {"status": "NOT_VERIFIED_RECORDS_NOT_SUPPLIED"},
              "software": {"python": platform.python_version(), "numpy": np.__version__, "pandas": pd.__version__}}
    if args.bbox_records:
        qc = validate_linked_records(pd.read_csv(args.bbox_records), selected_heads)
        qc.to_csv(args.out_dir / "bbox_record_qc.csv", index=False)
        report["bbox_record_audit"] = {"status": "TABULAR_IDENTITIES_CHECKED_NOT_IMAGE_REMEASUREMENT",
                                      "input_sha256": digest(args.bbox_records),
                                      "stable_head_counts": {str(k): int(qc[f"stable_{k}"].sum()) for k in (5, 10, 20)}}
    (args.out_dir / "audit_report.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
