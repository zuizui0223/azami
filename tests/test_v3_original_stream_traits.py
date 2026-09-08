import csv
import json
from pathlib import Path

from analysis.v3.stream_original_traits import (
    photo_schedule,
    qualified_endpoints,
    select_observation_ids,
)


ROOT = Path(__file__).resolve().parents[1]
DECISION = ROOT / "analysis" / "v3" / "measurement_qualification_decision_20260908.json"


def test_decision_exposes_exactly_14_original_stream_endpoints():
    decision = json.loads(DECISION.read_text(encoding="utf-8"))
    qualified, bbox = qualified_endpoints(decision)
    assert len(qualified) == 14
    assert "orientation_image_vertical_angle" in qualified
    assert "corolla_lab_chroma" in qualified
    assert "visible_floret_fraction" in qualified
    assert set(bbox) == {
        "orientation_image_vertical_angle",
        "capitulum_outline_aspect_ratio",
        "capitulum_outline_circularity",
        "capitulum_outline_solidity",
        "capitulum_width_profile_cv",
    }


def test_observation_selection_is_deterministic_and_sharded(tmp_path):
    cohort = tmp_path / "cohort.csv"
    with cohort.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["obs_id"])
        writer.writeheader()
        for obs in range(1, 101):
            writer.writerow({"obs_id": str(obs)})
    first = select_observation_ids(cohort, n=12, shard_index=0, shard_count=1)
    second = select_observation_ids(cohort, n=12, shard_index=0, shard_count=1)
    assert first == second and len(first) == 12
    shards = [set(select_observation_ids(cohort, n=0, shard_index=i, shard_count=4)) for i in range(4)]
    assert len(set.union(*shards)) == 100
    assert all(not (shards[i] & shards[j]) for i in range(4) for j in range(i+1, 4))


def test_photo_schedule_keeps_all_licensed_photos_for_selected_observations(tmp_path):
    metadata = tmp_path / "metadata.csv"
    fields = ["obs_id", "photo_id", "photo_license_code", "medium_image_url", "large_image_url"]
    rows = [
        {"obs_id":"1","photo_id":"11","photo_license_code":"cc-by","medium_image_url":"https://static.inaturalist.org/photos/11/medium.jpg","large_image_url":"https://static.inaturalist.org/photos/11/large.jpg"},
        {"obs_id":"1","photo_id":"12","photo_license_code":"cc0","medium_image_url":"https://static.inaturalist.org/photos/12/medium.jpg","large_image_url":"https://static.inaturalist.org/photos/12/large.jpg"},
        {"obs_id":"2","photo_id":"21","photo_license_code":"","medium_image_url":"https://static.inaturalist.org/photos/21/medium.jpg","large_image_url":"https://static.inaturalist.org/photos/21/large.jpg"},
        {"obs_id":"3","photo_id":"31","photo_license_code":"cc-by","medium_image_url":"https://static.inaturalist.org/photos/31/medium.jpg","large_image_url":"https://static.inaturalist.org/photos/31/large.jpg"},
    ]
    with metadata.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader(); writer.writerows(rows)
    schedule, report = photo_schedule(metadata, {"1", "2"})
    assert {row["photo_id"] for row in schedule} == {"11", "12"}
    assert report["selected_observations"] == 2
    assert report["scheduled_unique_photos"] == 2
    assert report["unavailable_states"]["license_unavailable"] == 1
    assert all("/original.jpg" in row["original_url"] for row in schedule)
