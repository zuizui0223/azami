from analysis.v3 import resolution_stream_gate as gate


def row(obs, photo, license="cc-by", url=None):
    url = url or f"https://inaturalist-open-data.s3.amazonaws.com/photos/{photo}/large.jpg"
    return {
        "obs_id": str(obs),
        "photo_id": str(photo),
        "photo_license_code": license,
        "medium_image_url": url.replace("/large.", "/medium."),
        "large_image_url": url,
        "raw_image_url": "",
    }


def test_size_url_preserves_extension_and_query():
    assert gate.sized_url("https://x/photos/4/large.JPG?foo=1", "original") == "https://x/photos/4/original.JPG?foo=1"


def test_selection_is_deterministic_and_one_photo_per_observation():
    rows = [row(1, 11), row(1, 12), row(2, 21), row(3, 31), row(4, 41), row(5, 51)]
    a = gate.select_rows(iter(rows), 3)
    b = gate.select_rows(iter(rows), 3)
    assert [(r["obs_id"], r["photo_id"]) for r in a] == [(r["obs_id"], r["photo_id"]) for r in b]
    assert len({r["obs_id"] for r in a}) == 3
    assert not any(r["photo_id"] == "12" for r in a)


def test_unlicensed_rows_are_not_admitted():
    assert not gate.eligible_source_row(row(1, 1, license=""))
    assert gate.eligible_source_row(row(1, 1, license="cc0"))


def test_normalized_box_matching_handles_resolution_change():
    medium = [[0, 0, 50, 50], [50, 50, 100, 100]]
    original = [[0, 0, 100, 100], [100, 100, 200, 200]]
    assert gate.greedy_matches(medium, original, (100, 100), (200, 200), .5) == [(0, 0, 1.0), (1, 1, 1.0)]


def test_summary_flags_original_gain_when_predeclared_threshold_crossed():
    rows = []
    for i in range(20):
        rows.append({
            "medium_status": "usable" if i < 16 else "low_resolution",
            "original_status": "usable",
            "medium_value": float(i) if i < 16 else None,
            "original_value": float(i),
        })
    result = gate.summarize_endpoint(rows, min_common=15, gain_threshold=.10)
    assert result["status"] == "ORIGINAL_MAY_ADD_INFORMATION"
    assert "original_eligibility_gain" in result["flags"]


def test_small_overlap_stays_insufficient():
    rows = [
        {"medium_status": "usable", "original_status": "usable", "medium_value": 1.0, "original_value": 2.0}
        for _ in range(5)
    ]
    assert gate.summarize_endpoint(rows, min_common=15)["status"] == "INSUFFICIENT_COMMON_USABLE"
