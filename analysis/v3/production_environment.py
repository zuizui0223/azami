"""Checkpointed local extraction for the full native source; never reads traits.

Only explicit metadata columns are read. Each variable/month is a committed,
hash-checked numerical checkpoint. Remote rasters are sampled by occupied blocks,
not downloaded in full. A completion receipt is impossible until every task has
finished. No ecological regression or original-photo request is implemented here.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import sys

import numpy as np
import pandas as pd
import requests

from .environment_matrix import (DEFAULT_CONTRACT, candidate_specs, load_contract,
                                 month_url, sha256_file, vif_table)
from .workflow import canonical_digest

ROOT = Path(__file__).resolve().parents[2]
CONTRACT = ROOT / "analysis/v3/environment_production_contract.json"
SOURCE_COLUMNS = ["obs_id", "accepted_key", "analysis_latitude", "analysis_longitude",
                  "observation_month", "native_range_status", "sin_doy", "cos_doy",
                  "observed_year", "south_indicator", "south_sin", "south_cos"]


def stamp() -> str:
    return datetime.now(timezone.utc).isoformat()


def write_json(path: Path, value: dict) -> None:
    """Atomic status update inside a previously verified run directory."""
    temporary = path.with_name(path.name + ".tmp")
    with temporary.open("w", encoding="utf-8", newline="\n") as stream:
        json.dump(value, stream, indent=2, ensure_ascii=False, allow_nan=False)
        stream.write("\n")
        stream.flush()
        os.fsync(stream.fileno())
    os.replace(temporary, path)


def validate_contract(contract: dict) -> None:
    groups = contract["processes"]
    eligible = [v for group in groups.values() for v in group]
    all_ids = eligible + list(contract["context_only"])
    if (len(all_ids) != len(set(all_ids)) or set(all_ids) != set(contract["acquire"])
            or len(contract["acquire"]) != 15 or len(groups) != 4
            or contract["selection"]["threshold"] != 10
            or not set(contract["selection"]["protected"]) <= set(eligible)
            or any(not set(values) & set(contract["selection"]["protected"]) for values in groups.values())
            or contract["model"]["primary_tests"] != 36
            or contract["model"]["formulations"] != 1
            or contract["execution"]["ecological_fitting_authorized"] is not False
            or contract["execution"]["production_image_execution_authorized"] is not False):
        raise ValueError("Invalid process-selection contract")


def source_frame(path: Path, contract: dict) -> pd.DataFrame:
    if sha256_file(path) != contract["source"]["enriched_csv_sha256"]:
        raise ValueError("Full native source hash differs")
    data = pd.read_csv(path, usecols=SOURCE_COLUMNS, dtype={"obs_id": str, "accepted_key": str})
    if (data.empty or data["obs_id"].isna().any() or data["obs_id"].duplicated().any()
            or data["accepted_key"].isna().any()
            or len(data) != contract["source"]["observations"]
            or data["accepted_key"].nunique() != contract["source"]["taxa"]
            or not data["native_range_status"].eq("native").all()):
        raise ValueError("Full native source membership differs")
    numeric = data[[c for c in SOURCE_COLUMNS if c not in ("obs_id", "accepted_key", "native_range_status")]]
    if (not np.isfinite(numeric.to_numpy(float)).all()
            or not data["analysis_latitude"].between(-90, 90).all()
            or not data["analysis_longitude"].between(-180, 180).all()
            or not data["observation_month"].between(1, 12).all()
            or not data["observation_month"].eq(data["observation_month"].astype(int)).all()):
        raise ValueError("Invalid source coordinate/calendar")
    data = data.sort_values("obs_id", kind="stable").reset_index(drop=True)
    data["equal_taxon_weight"] = 1.0 / data.groupby("accepted_key")["obs_id"].transform("size")
    return data


def object_identity(url: str) -> dict:
    response = requests.head(url, timeout=(15, 45), allow_redirects=True)
    response.raise_for_status()
    etag = response.headers.get("ETag", "")
    length = response.headers.get("Content-Length", "")
    if not etag or etag.startswith("W/") or not length.isdigit():
        raise ValueError("Remote raster lacks a strong object identity")
    return {"url": url, "resolved_url": response.url, "etag": etag,
            "bytes": int(length), "last_modified": response.headers.get("Last-Modified")}


def sample_blocks(dataset, coordinates: np.ndarray) -> tuple[np.ndarray, dict]:
    """Nearest containing pixel, identical scale/mask semantics to dataset.sample."""
    import rasterio
    from rasterio.warp import transform
    if dataset.crs is None:
        raise ValueError("Raster has no CRS")
    xs, ys = coordinates[:, 0], coordinates[:, 1]
    if str(dataset.crs).upper() not in {"EPSG:4326", "OGC:CRS84"}:
        xs, ys = transform("EPSG:4326", dataset.crs, xs.tolist(), ys.tolist())
    rows, cols = rasterio.transform.rowcol(dataset.transform, xs, ys)
    rows, cols = np.asarray(rows, dtype=int), np.asarray(cols, dtype=int)
    valid = (rows >= 0) & (rows < dataset.height) & (cols >= 0) & (cols < dataset.width)
    values = np.full(len(coordinates), np.nan)
    block_h, block_w = dataset.block_shapes[0]
    positions = np.flatnonzero(valid)
    block_cols = (dataset.width + block_w - 1) // block_w
    keys = (rows[positions] // block_h) * block_cols + cols[positions] // block_w
    order = np.argsort(keys, kind="stable")
    positions, keys = positions[order], keys[order]
    cuts = np.r_[0, np.flatnonzero(np.diff(keys)) + 1, len(keys)] if len(keys) else [0]
    scale, offset = float(dataset.scales[0]), float(dataset.offsets[0])
    for left, right in zip(cuts[:-1], cuts[1:]):
        selected = positions[left:right]
        row0 = rows[selected[0]] // block_h * block_h
        col0 = cols[selected[0]] // block_w * block_w
        window = rasterio.windows.Window(int(col0), int(row0),
                                         min(block_w, dataset.width - col0),
                                         min(block_h, dataset.height - row0))
        block = dataset.read(1, window=window, masked=True).astype(float)
        raw = block[rows[selected] - row0, cols[selected] - col0]
        values[selected] = np.ma.filled(raw, np.nan) * scale + offset
    values[~np.isfinite(values)] = np.nan
    return values, {"crs": str(dataset.crs), "transform": list(dataset.transform),
                    "width": dataset.width, "height": dataset.height,
                    "scale": scale, "offset": offset, "block_shape": [block_h, block_w],
                    "blocks_read": len(cuts) - 1, "out_of_bounds": int((~valid).sum())}


def task_plan(frame: pd.DataFrame, raster_contract: dict, variables: list[str]) -> list[dict]:
    result = []
    specs = {s.output_id: s for s in candidate_specs(raster_contract)}
    for variable in variables:
        spec = specs[variable]
        months = sorted(frame["observation_month"].unique()) if spec.temporal_mode == "monthly" else [0]
        for month in months:
            indices = (np.flatnonzero(frame["observation_month"].eq(month))
                       if month else np.arange(len(frame)))
            result.append({"id": f"{variable}_{int(month):02d}", "variable": variable,
                           "month": int(month), "indices": indices,
                           "url": month_url(raster_contract, spec.variable, int(month)) if month else spec.url})
    return result


def checkpoint(task: dict, frame: pd.DataFrame, directory: Path, run_id: str,
               identity_reader=object_identity, sampler=sample_blocks) -> dict:
    receipt_path = directory / (task["id"] + ".json")
    data_path = directory / (task["id"] + ".npz")
    indices = task["indices"]
    if receipt_path.exists():
        receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
        if (receipt["run_id"] != run_id or receipt["task"] != task["id"]
                or receipt["source"]["url"] != task["url"]
                or sha256_file(data_path) != receipt["numerical_sha256"]):
            raise ValueError("Checkpoint identity/hash differs")
        with np.load(data_path, allow_pickle=False) as saved:
            if (not np.array_equal(saved["indices"], indices)
                    or saved["values"].shape != indices.shape
                    or np.isinf(saved["values"]).any()):
                raise ValueError("Checkpoint membership/values differ")
        return receipt
    # Orphans from an interrupted pre-commit task are preserved for inspection.
    if data_path.exists():
        raise ValueError("Uncommitted checkpoint exists; inspect before resuming")
    import rasterio
    before = identity_reader(task["url"])
    coordinates = frame.iloc[indices][["analysis_longitude", "analysis_latitude"]].to_numpy(float)
    with rasterio.Env(GDAL_DISABLE_READDIR_ON_OPEN="EMPTY_DIR", CPL_VSIL_CURL_ALLOWED_EXTENSIONS=".tif,.TIF",
                      GDAL_HTTP_TIMEOUT="45", GDAL_HTTP_CONNECTTIMEOUT="15", GDAL_HTTP_MAX_RETRY="2",
                      GDAL_HTTP_RETRY_DELAY="1", GDAL_CACHEMAX=128 * 1024 * 1024):
        with rasterio.open(task["url"]) as dataset:
            values, raster = sampler(dataset, coordinates)
    if before != identity_reader(task["url"]):
        raise ValueError("Remote raster changed during sampling")
    with data_path.open("xb") as stream:
        np.savez_compressed(stream, indices=indices, values=values)
        stream.flush()
        os.fsync(stream.fileno())
    receipt = {"task": task["id"], "run_id": run_id, "source": before, "raster": raster,
               "rows": len(indices), "finite_rows": int(np.isfinite(values).sum()),
               "numerical_sha256": sha256_file(data_path), "completed_at": stamp()}
    write_json(receipt_path, receipt)
    return receipt


def choose_process_variables(frame: pd.DataFrame, contract: dict) -> dict:
    validate_contract(contract)
    candidates = [v for group in contract["processes"].values() for v in group]
    values = frame[candidates].apply(pd.to_numeric, errors="raise")
    finite = np.isfinite(values.to_numpy(float)).all(axis=1)
    data = frame.loc[finite, ["accepted_key", *candidates]].copy()
    complete_taxa = data["accepted_key"].nunique()
    report = {"status": "ENVIRONMENT_SELECTION_NOT_ESTIMABLE", "source_rows": len(frame),
              "complete_rows": len(data), "source_taxa": frame["accepted_key"].nunique(),
              "complete_taxa": complete_taxa, "complete_case_scope_fixed": True,
              "source_rows_excluded_from_selection_only": int((~finite).sum()),
              "selected": [], "trace": [], "context_only": contract["context_only"],
              "trait_values_read": 0, "ecological_fitting_authorized": False}
    if len(data) < 3 or complete_taxa < 2:
        report["reason"] = "Insufficient fixed complete-case support"
        return report
    constants = [v for v in candidates if data[v].nunique() < 2]
    if constants:
        report.update(reason="Constant process candidates", constant_variables=constants)
        return report
    data["equal_taxon_weight"] = 1.0 / data.groupby("accepted_key")["accepted_key"].transform("size")
    selected = candidates.copy()
    protected = set(contract["selection"]["protected"])
    threshold = contract["selection"]["threshold"]
    while True:
        table, matrix = vif_table(data, selected, "equal_taxon_weight")
        vif = {row.variable: float("inf") if row.vif == "infinite" else float(row.vif)
               for row in table.itertuples(index=False)}
        removable = [v for v in selected if v not in protected and vif[v] >= threshold]
        removed = sorted(removable, key=lambda v: (-vif[v], v))[0] if removable else None
        report["trace"].append({"variables": selected.copy(), "vif": table.to_dict(orient="records"),
                                "matrix": matrix, "removed": removed})
        if removed is None:
            break
        selected.remove(removed)
    if matrix["matrix_rank"] != len(selected) or any(v >= threshold for v in vif.values()):
        report["reason"] = "Protected process representation remains rank deficient or VIF >= 10"
        return report
    w = data["equal_taxon_weight"].to_numpy(float)
    centers, scales = {}, {}
    for variable in selected:
        x = data[variable].to_numpy(float)
        centers[variable] = float(np.average(x, weights=w))
        scales[variable] = float(np.sqrt(np.average((x - centers[variable]) ** 2, weights=w)))
    report.update(status="FULL_NATIVE_ENVIRONMENT_REPRESENTATION_SELECTED_NO_ECOLOGY", selected=selected,
                  processes={key: [v for v in group if v in selected] for key, group in contract["processes"].items()},
                  centers=centers, scales=scales,
                  removed=[v for v in candidates if v not in selected],
                  primary_test_slots=contract["model"]["primary_tests"])
    return report


def run(source: Path, out: Path, workers: int = 2, resume: bool = False) -> dict:
    if not 1 <= workers <= 2:
        raise ValueError("Use one or two bounded raster workers")
    contract = json.loads(CONTRACT.read_text(encoding="utf-8"))
    validate_contract(contract)
    raster_contract = load_contract(DEFAULT_CONTRACT)
    data = source_frame(source, contract)
    identity = {"source_sha256": contract["source"]["enriched_csv_sha256"],
                "selection_contract": contract, "raster_contract": raster_contract,
                "implementation_sha256_text_lf": hashlib.sha256(Path(__file__).read_text(encoding="utf-8").encode()).hexdigest(),
                "helper_sha256_text_lf": hashlib.sha256((ROOT / "analysis/v3/environment_matrix.py").read_text(encoding="utf-8").encode()).hexdigest(),
                "runtime": {name: importlib.metadata.version(name) for name in ("numpy", "pandas", "rasterio", "requests")},
                "python": sys.version}
    run_id = canonical_digest(identity)
    if out.exists():
        if not resume or json.loads((out / "run_identity.json").read_text())["identity"] != identity:
            raise ValueError("Existing run requires --resume and exact frozen input/code/runtime identity")
    else:
        out.mkdir(parents=True, exist_ok=False)
        (out / "checkpoints").mkdir()
        write_json(out / "run_identity.json", {"run_id": run_id, "started_at": stamp(), "identity": identity})
    tasks = task_plan(data, raster_contract, contract["acquire"])
    progress = {"status": "FULL_NATIVE_ENVIRONMENT_EXTRACTION_RUNNING", "run_id": run_id,
                "observations": len(data), "taxa": data["accepted_key"].nunique(),
                "tasks_total": len(tasks), "tasks_completed": 0, "tasks_failed": [],
                "trait_values_read": 0, "image_requests_executed": 0, "ecological_models_executed": 0,
                "ecological_fitting_authorized": False, "off_device_restore_verified": False}
    def emit():
        progress["updated_at"] = stamp()
        write_json(out / "progress.json", progress)
        print(json.dumps(progress), flush=True)
    emit()
    with ThreadPoolExecutor(max_workers=workers) as pool:
        futures = {pool.submit(checkpoint, task, data, out / "checkpoints", run_id): task for task in tasks}
        for future in as_completed(futures):
            task = futures[future]
            try:
                future.result()
                progress["tasks_completed"] += 1
            except Exception as exc:
                progress["tasks_failed"].append({"task": task["id"], "error_type": type(exc).__name__,
                                                  "error": str(exc)[:600]})
            emit()
    if progress["tasks_failed"]:
        progress["status"] = "INCOMPLETE_PROVIDER_OR_CHECKPOINT_FAILURE_NO_SELECTION"
        emit()
        return progress
    for variable in contract["acquire"]:
        data[variable] = np.nan
    for task in tasks:
        with np.load(out / "checkpoints" / (task["id"] + ".npz"), allow_pickle=False) as saved:
            data.loc[saved["indices"], task["variable"]] = saved["values"]
    matrix_path = out / "environment_candidate_matrix_private.csv"
    if not matrix_path.exists():
        data.to_csv(matrix_path, index=False, lineterminator="\n", mode="x")
    else:
        old = pd.read_csv(matrix_path, dtype={"obs_id": str, "accepted_key": str})
        pd.testing.assert_frame_equal(old, data, check_dtype=False, rtol=1e-12, atol=1e-12)
    from .environment_matrix import diagnose_environment
    diagnostic, tables = diagnose_environment(data, contract["acquire"], 0.8, "equal_taxon_weight")
    for key, table in tables.items():
        table.to_csv(out / ("environment_" + key + ".csv"), index=False)
    write_json(out / "environment_diagnostics.json", diagnostic)
    progress.update(status="FULL_NATIVE_ENVIRONMENT_EXTRACTION_COMPLETE_NO_ECOLOGY",
                    selection_status="NOT_EXECUTED_ACQUISITION_ONLY",
                    matrix_sha256=sha256_file(matrix_path))
    emit()
    return progress


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=2)
    parser.add_argument("--resume", action="store_true")
    args = parser.parse_args()
    result = run(args.source, args.out_dir, args.workers, args.resume)
    return 0 if result["status"] == "FULL_NATIVE_ENVIRONMENT_EXTRACTION_COMPLETE_NO_ECOLOGY" else 2


if __name__ == "__main__":
    raise SystemExit(main())
