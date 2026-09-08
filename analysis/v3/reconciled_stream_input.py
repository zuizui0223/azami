"""Location-blind, whole-component pilot input from a pinned reconciled schedule.

This reads no ecological columns and fetches no images. It never falls back to
the historical merged photo table. A pilot is the hash-ordered component prefix
that fits its observation budget, not a selection for successful image yield.
"""
from __future__ import annotations

import json
from pathlib import Path
import sqlite3

from .reconciled_photo_schedule import STATUS, component_score
from .workflow import ROOT, canonical_digest, digest, text_digest


def pilot_input(path: Path, expected_sha256: str, maximum_observations: int = 128, *, component_ids: list[str] | None = None) -> dict:
    if not 1 <= maximum_observations <= 128:
        raise ValueError("Reconciled pilot size must be 1..128 observations")
    if not expected_sha256 or digest(path) != expected_sha256:
        raise ValueError("Reconciled schedule identity mismatch")
    contract = json.loads((ROOT / "analysis/v3/reconciled_photo_schedule_contract.json").read_text(encoding="utf-8"))
    with sqlite3.connect(path.resolve().as_uri() + "?mode=ro", uri=True) as db:
        db.row_factory = sqlite3.Row
        execution = {row[0]: json.loads(row[1]) for row in db.execute("SELECT key,value_json FROM execution")}
        if (execution.get("status") != STATUS or execution.get("production_image_execution_authorized") is not False
                or execution.get("contract_canonical_sha256") != canonical_digest(contract)):
            raise ValueError("Reconciled schedule execution contract mismatch")
        groups = [(component_score(row[0], contract["ordering_salt"]), row[0], row[1]) for row in
                  db.execute("SELECT component_id,COUNT(*) FROM native_observations GROUP BY component_id")]
        groups.sort()
        selected_groups, actual = [], 0
        if component_ids is not None:
            if not component_ids or len(component_ids) != len(set(component_ids)):
                raise ValueError("Chunk components must be nonempty and unique")
            requested = set(component_ids)
            groups = [g for g in groups if g[1] in requested]
            if {g[1] for g in groups} != requested or sum(g[2] for g in groups) > maximum_observations:
                raise ValueError("Chunk components are absent or exceed the observation bound")
        for score, component, count in groups:
            if actual + count > maximum_observations:
                break
            selected_groups.append((component, score))
            actual += count
        if not selected_groups:
            raise ValueError("First complete component exceeds the bounded pilot budget")
        db.execute("CREATE TEMP TABLE selected_component(component_id TEXT PRIMARY KEY,score TEXT)")
        db.executemany("INSERT INTO selected_component VALUES (?,?)", selected_groups)
        scores = dict(db.execute("SELECT obs_id,score FROM native_observations JOIN selected_component USING(component_id)"))
        selected = sorted(scores, key=lambda obs: (scores[obs], obs))
        if len(selected) != actual:
            raise ValueError("Selected component membership mismatch")
        jobs = [dict(row) for row in db.execute("SELECT j.* FROM photo_jobs j JOIN selected_component s USING(component_id) ORDER BY s.score,j.photo_id")]
        queue, links, states = [], [], {}
        components = dict(selected_groups)
        allowed_states = {"source_photo_version_missing", "license_conflict", "license_unavailable",
                          "url_unavailable_or_invalid", "url_conflict", "request_candidate_not_authorized"}
        for job in jobs:
            photo, state = job["photo_id"], job["state"]
            if state not in allowed_states or job["component_score"] != components[job["component_id"]]:
                raise ValueError("Unknown schedule state or component score mismatch")
            native = [row[0] for row in db.execute("SELECT obs_id FROM native_links WHERE photo_id=? ORDER BY obs_id", (photo,))]
            if not native or not set(native).issubset(scores):
                raise ValueError("Photo crosses selected component boundary")
            known = [row[0] for row in db.execute("SELECT obs_id FROM known_photo_links WHERE photo_id=? ORDER BY obs_id", (photo,))]
            versions = []
            for row in db.execute("SELECT * FROM photo_versions WHERE photo_id=? ORDER BY kind,origin,source_row,photo_index", (photo,)):
                value = dict(row)
                value["photo_fields"] = json.loads(value.pop("photo_fields_json"))
                versions.append(value)
            if len(versions) != job["source_version_count"] or set(known) != {v["obs_id"] for v in versions}:
                raise ValueError("Source version count or known link mismatch")
            states[state] = states.get(state, 0) + 1
            if state == "request_candidate_not_authorized":
                if not job["original_url"] or job["license_code"] not in contract["license_codes"]:
                    raise ValueError("Request candidate lacks a permitted source")
                queue.append({"photo_id": photo, "component_id": job["component_id"], "obs_ids": native,
                              "known_metadata_obs_ids": known, "original_url": job["original_url"],
                              "license_code": job["license_code"], "source_versions": versions})
            for obs in native:
                links.append({"source_row": "", "photo_id": photo, "obs_id": obs, "status": state,
                              "license_code": job["license_code"], "source_url": "", "original_url": job["original_url"],
                              "source_metadata_json": json.dumps({"component_id": job["component_id"], "verified_versions": versions}, sort_keys=True),
                              "known_photo_versions_json": json.dumps(sorted({(v["license_code"], v["original_url"]) for v in versions})),
                              "known_photo_observation_ids_json": json.dumps(known)})
        expected_links = db.execute("SELECT SUM(expected_photo_count) FROM native_observations JOIN selected_component USING(component_id)").fetchone()[0]
        if len(links) != expected_links or {link["obs_id"] for link in links} != set(selected):
            raise ValueError("Selected observations or photo links were lost")
    return {"selected": selected, "selection_scores": scores, "queue": queue, "links": links,
            "input_sha256": execution["input_sha256"], "schedule_sha256": expected_sha256,
            "report": {"mode": "reconciled_whole_component_chunk" if component_ids is not None else "reconciled_whole_component_pilot", "ordering_salt": contract["ordering_salt"],
                       "maximum_observations": maximum_observations, "selected_observations": actual,
                       "selected_components": len(selected_groups), "selected_photo_links": len(links),
                       "selected_photo_jobs": len(jobs), "request_candidates": len(queue), "photo_states": states,
                       "source_reconciliation_verified": True, "production_image_execution_authorized": False,
                       "images_fetched": 0, "trait_values_read": 0,
                       "adapter_sha256_text_lf": text_digest(Path(__file__))}}
