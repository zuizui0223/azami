"""Retain all source records while grouping known dependence before evaluation.

Connected components join observations sharing a photo identity, exact image
bytes or EXIF-oriented decoded pixels. These are dependence groups, not inferred
biological individuals. No taxon, coordinates, trait values or human labels enter.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path, PurePosixPath
import sqlite3

from .build_image_workspace import readable
from .workflow import ROOT, digest, text_digest

DEVELOPMENT_POOL = "historical_detector_development_pool"
TRAINING_SHA = "9baf5bfe477586e2e4c3c35b15af1ec34424cb2187869903a2a1f5f32bdc3871"
SPECIFICATION = {
    "version": "v3_source_dependence_components_v1", "source_retention": "all photos and observation-photo links",
    "component_edges": ["shared_photo_id", "exact_byte_identity", "exact_oriented_pixel_identity"],
    "component_id": "obs: followed by the lexicographically smallest original observation ID in the component",
    "development_exposure": "propagate historical usage pools and recorded training/validation membership to every connected image/observation",
    "fold_rule": "SHA-256 of v3-dependence-fold:20260907: plus component ID, interpreted as a big-endian integer modulo 5",
    "fold_meaning": "operational grouping only; not a new independent evaluation split or evidence that historical training was group-separated",
    "unknown_dependence": "different encodings/crops of the same scene without exact pixel identity, photographer/site dependence and same-head identity remain unassessed",
}


class Components:
    def __init__(self):
        self.parents = {}

    def find(self, key):
        if not key:
            raise ValueError("Missing observation identity")
        parent = self.parents.setdefault(key, key)
        while parent != self.parents[parent]:
            parent = self.parents[parent]
        while key != parent:
            following = self.parents[key]
            self.parents[key] = parent
            key = following
        return parent

    def union(self, left, right):
        left, right = self.find(left), self.find(right)
        if left == right:
            return False
        lower, upper = sorted((left, right))
        self.parents[upper] = lower
        return True


def fold(component):
    return int.from_bytes(hashlib.sha256(("v3-dependence-fold:20260907:"+component).encode()).digest(), "big") % 5


def training_membership(source, manifest):
    """Exact filename join to the archived development cache, not ID guessing."""
    candidates = {}
    for photo, path in source.execute("SELECT photo_id,source_path FROM cache_records WHERE pool=?", (DEVELOPMENT_POOL,)):
        basename = PurePosixPath(path.replace("\\", "/")).name
        if basename in candidates:
            raise ValueError("Ambiguous development-cache filename")
        candidates[basename] = photo
    output, seen = [], set()
    with readable(manifest).open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            if None in row or any(v is None for v in row.values()):
                raise ValueError("Malformed training provenance row")
            if row["split"] not in ("train", "val") or row["queue_id"] in seen:
                raise ValueError("Invalid or repeated historical training entry")
            seen.add(row["queue_id"])
            basename = PurePosixPath(row["source_image"].replace("\\", "/")).name
            if basename not in candidates:
                raise ValueError("Historical training image is not linked to the development cache")
            output.append((row["queue_id"], candidates[basename], row["split"]))
    return output


def build(workspace, manifest, out, expected_training_sha=TRAINING_SHA):
    workspace, out = workspace.resolve(), out.resolve()
    if out == ROOT or (ROOT in out.parents and not any(out.is_relative_to(ROOT/p) for p in ("local_data", "outputs"))):
        raise ValueError("Use an external or ignored local output directory")
    if out.exists():
        raise ValueError("Preserve earlier dependence assignments")
    if digest(readable(manifest)) != expected_training_sha:
        raise ValueError("Historical training-manifest identity mismatch")
    report = json.loads((workspace/"image_workspace_report.json").read_text(encoding="utf-8"))
    input_db = workspace/"image_workspace.sqlite"
    source_sha = digest(input_db)
    if source_sha != report["workspace_database_sha256"]:
        raise ValueError("Image-workspace identity mismatch")
    out.mkdir(parents=True)
    try:
        with sqlite3.connect(input_db.as_uri()+"?mode=ro", uri=True) as source, sqlite3.connect(out/"dependence_groups.sqlite") as db:
            db.executescript("""
                CREATE TABLE observations(obs_id TEXT PRIMARY KEY,component_id TEXT,fold INTEGER);
                CREATE TABLE photos(photo_id TEXT PRIMARY KEY,component_id TEXT);
                CREATE TABLE source_links(obs_id TEXT,photo_id TEXT,PRIMARY KEY(obs_id,photo_id));
                CREATE TABLE objects(sha256 TEXT PRIMARY KEY,component_id TEXT,decoded_rgb_sha256 TEXT);
                CREATE TABLE photo_versions(photo_id TEXT,sha256 TEXT,PRIMARY KEY(photo_id,sha256));
                CREATE TABLE evidence(kind TEXT,identity_key TEXT,left_obs TEXT,right_obs TEXT,merged_components INTEGER);
                CREATE TABLE historical_usage(photo_id TEXT,pool TEXT);
                CREATE TABLE historical_training(queue_id TEXT PRIMARY KEY,photo_id TEXT,split TEXT);
                CREATE TABLE components(component_id TEXT PRIMARY KEY,n_observations INTEGER,n_photos INTEGER,n_cached_objects INTEGER,fold INTEGER,development_pool_exposed INTEGER,training_exposed INTEGER,validation_exposed INTEGER,pools_json TEXT);
            """)
            groups = Components()
            owners = {}
            batch = []
            for obs, photo in source.execute("SELECT obs_id,photo_id FROM source_links ORDER BY photo_id,obs_id"):
                groups.find(obs)
                if photo in owners:
                    left = owners[photo]
                    db.execute("INSERT INTO evidence VALUES (?,?,?,?,?)", ("shared_photo_id", photo, left, obs, int(groups.union(left, obs))))
                else:
                    owners[photo] = obs
                batch.append((obs, photo))
                if len(batch) == 20000:
                    db.executemany("INSERT INTO source_links VALUES (?,?)", batch)
                    batch.clear()
            db.executemany("INSERT INTO source_links VALUES (?,?)", batch)
            n_photos = 0
            for (photo,) in source.execute("SELECT photo_id FROM photos"):
                if photo not in owners:
                    raise ValueError("Photo universe and source-link coverage differ")
                n_photos += 1
            if n_photos != len(owners):
                raise ValueError("Photo universe and source-link coverage differ")
            versions = list(source.execute("SELECT photo_id,sha256 FROM photo_versions ORDER BY photo_id,sha256"))
            object_pixels = dict(source.execute("SELECT sha256,decoded_rgb_sha256 FROM objects"))
            by_byte, by_pixel = {}, {}
            for photo, sha in versions:
                if photo not in owners or sha not in object_pixels:
                    raise ValueError("Orphan cached photo version")
                obs = owners[photo]
                for kind, key, mapping in (("exact_byte_identity", sha, by_byte), ("exact_oriented_pixel_identity", object_pixels[sha], by_pixel)):
                    if key in mapping:
                        left = mapping[key]
                        db.execute("INSERT INTO evidence VALUES (?,?,?,?,?)", (kind, key, left, obs, int(groups.union(left, obs))))
                    else:
                        mapping[key] = obs
            if set(by_byte) != set(object_pixels):
                raise ValueError("Cached object lacks an original photo link")
            assigned = {obs: "obs:"+groups.find(obs) for obs in groups.parents}
            db.executemany("INSERT INTO observations VALUES (?,?,?)", ((o, assigned[o], fold(assigned[o])) for o in sorted(assigned)))
            db.executemany("INSERT INTO photos VALUES (?,?)", ((p, assigned[owners[p]]) for p in sorted(owners)))
            db.executemany("INSERT INTO photo_versions VALUES (?,?)", versions)
            db.executemany("INSERT INTO objects VALUES (?,?,?)", ((s, assigned[by_byte[s]], object_pixels[s]) for s in sorted(object_pixels)))
            usage = list(source.execute("SELECT photo_id,pool FROM cache_records ORDER BY record_id"))
            if any(p not in owners for p, _ in usage):
                raise ValueError("Usage provenance outside the source universe")
            training = training_membership(source, manifest)
            db.executemany("INSERT INTO historical_usage VALUES (?,?)", usage)
            db.executemany("INSERT INTO historical_training VALUES (?,?,?)", training)
            summary, pools, splits = {}, {}, {}
            for obs, component in assigned.items():
                summary.setdefault(component, [0, 0, 0])[0] += 1
            for p, owner in owners.items():
                summary[assigned[owner]][1] += 1
            for sha, obs in by_byte.items():
                summary[assigned[obs]][2] += 1
            for p, pool in usage:
                pools.setdefault(assigned[owners[p]], set()).add(pool)
            for _, p, split in training:
                splits.setdefault(assigned[owners[p]], set()).add(split)
            db.executemany("INSERT INTO components VALUES (?,?,?,?,?,?,?,?,?)", (
                (c, v[0], v[1], v[2], fold(c), int(DEVELOPMENT_POOL in pools.get(c, ())), int("train" in splits.get(c, ())), int("val" in splits.get(c, ())), json.dumps(sorted(pools.get(c, ())))) for c, v in sorted(summary.items())))
            db.executescript("CREATE INDEX observation_component ON observations(component_id); CREATE INDEX photo_component ON photos(component_id); CREATE INDEX object_component ON objects(component_id);")
            scalar = lambda query: db.execute(query).fetchone()[0]
            checks = {
                "source_links_crossing_components": scalar("SELECT COUNT(*) FROM source_links l JOIN observations o USING(obs_id) JOIN photos p USING(photo_id) WHERE o.component_id!=p.component_id"),
                "photo_versions_crossing_components": scalar("SELECT COUNT(*) FROM photo_versions v JOIN photos p USING(photo_id) JOIN objects o USING(sha256) WHERE p.component_id!=o.component_id"),
                "identical_pixels_crossing_components": scalar("SELECT COUNT(*) FROM (SELECT decoded_rgb_sha256 FROM objects GROUP BY decoded_rgb_sha256 HAVING COUNT(DISTINCT component_id)>1)"),
                "folds_crossing_components": scalar("SELECT COUNT(*) FROM observations o JOIN components c USING(component_id) WHERE o.fold!=c.fold"),
            }
            if any(checks.values()):
                raise ValueError("Dependence closure or fold integrity failed")
            counts = {"observations_retained": len(assigned), "photos_retained": len(owners),
                      "source_links_retained": scalar("SELECT COUNT(*) FROM source_links"), "components": len(summary),
                      "components_with_multiple_observations": scalar("SELECT COUNT(*) FROM components WHERE n_observations>1"),
                      "maximum_observations_in_component": scalar("SELECT MAX(n_observations) FROM components"),
                      "cached_objects_retained": len(object_pixels), "components_with_cached_objects": scalar("SELECT COUNT(*) FROM components WHERE n_cached_objects>0"),
                      "components_with_multiple_cached_objects": scalar("SELECT COUNT(*) FROM components WHERE n_cached_objects>1"),
                      "historical_training_images": len(training), "historical_train_images": sum(s=="train" for _,_,s in training),
                      "historical_validation_images": sum(s=="val" for _,_,s in training),
                      "components_shared_by_historical_train_and_validation": scalar("SELECT COUNT(*) FROM components WHERE training_exposed=1 AND validation_exposed=1"),
                      "cached_objects_in_development_exposed_components": scalar("SELECT SUM(n_cached_objects) FROM components WHERE development_pool_exposed=1"),
                      "components_with_multiple_historical_usage_pools": sum(len(v)>1 for v in pools.values()),
                      "source_rows_deleted": 0}
            if counts["photos_retained"] != report["counts"]["source_photo_ids"] or counts["observations_retained"] != report["counts"]["source_observation_ids"] or counts["source_links_retained"] != report["counts"]["source_observation_photo_links"]:
                raise ValueError("Source denominators changed")
        result = {"status": "FULL_SOURCE_KNOWN_DEPENDENCE_GROUPS_BUILT", "specification": SPECIFICATION,
                  "source_workspace_sha256": source_sha, "historical_training_manifest_sha256": digest(readable(manifest)),
                  "implementation_sha256_text_lf": text_digest(Path(__file__)), "counts": counts, "integrity_checks": checks,
                  "output_database_sha256": digest(out/"dependence_groups.sqlite"),
                  "limits": ["Exact-identity closure does not find every duplicated scene or identify the same biological individual.",
                             "Operational folds do not retroactively repair historical training leakage or create fresh independent validation after outcome inspection.",
                             "Component grouping preserves source rows; it does not make all photos within a component identical biological measurements.",
                             "Historical development exposure is propagated conservatively; absence of that record does not establish absence of all model exposure."]}
        (out/"dependence_groups_report.json").write_text(json.dumps(result, indent=2)+"\n", encoding="utf-8", newline="\n")
        print(json.dumps(result, indent=2))
        return result
    except (Exception, KeyboardInterrupt) as error:
        (out/"incomplete_run.json").write_text(json.dumps({"status": "INCOMPLETE_DO_NOT_USE", "error_type": type(error).__name__, "reason": str(error)}, indent=2), encoding="utf-8")
        raise


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workspace", type=Path, required=True)
    parser.add_argument("--training-manifest", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    args = parser.parse_args()
    build(args.workspace, args.training_manifest, args.out_dir)


if __name__ == "__main__":
    main()
