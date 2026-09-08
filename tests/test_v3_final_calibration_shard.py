import hashlib

import numpy as np

from analysis.v3.dependence_resampling import cohort_partition, crossed_draw, source_partition
from analysis.v3.multicoordinate_calibration import generate
from analysis.v3.run_final_calibration_shard import definition, outer_seed


def test_final_calibration_seed_and_shard_geometry_are_frozen():
    spec, prelim = definition()
    assert spec["seed"] == 2026090893
    assert spec["bootstrap_replicates"] == 999
    assert spec["module_dimensions"] == {"orientation": 1, "visible_colour": 8, "gross_shape": 4}
    assert spec["predictors"] == 9
    assert spec["nuisance_columns"] == 13
    assert spec["family_slots"] == 36
    assert outer_seed(spec, "iid_null", 0) == 2026090893
    assert outer_seed(spec, "crossed_spatial_null", 0) == 2026190893
    assert outer_seed(spec, "scale_difference", 399) == 2026391292
    assert list(prelim["process_indices"]) == ["wetting_moisture", "radiation", "heat_drying", "mechanical"]
    assert not spec["ecological_fitting_authorized"]


def test_deterministic_crossed_draw_sequence_is_module_independent():
    spec, _ = definition()
    seed = outer_seed(spec, "crossed_spatial_null", 3)
    data = generate("crossed_spatial_null", seed)
    source = source_partition(
        np.arange(len(data["taxa"])), data["taxa"], data["components"],
        data["latitude"], data["longitude"], grid_degrees=2,
    )
    positions, _, _ = cohort_partition(source, np.arange(len(data["taxa"])))

    def fingerprint():
        digest = hashlib.sha256()
        for replicate in range(5):
            indices, copy_taxa = crossed_draw(source, positions, seed=seed, replicate=replicate)
            digest.update(np.asarray(indices, dtype=np.int64).tobytes())
            digest.update(b"\0")
            digest.update("\n".join(copy_taxa.tolist()).encode("utf-8"))
            digest.update(b"\xff")
        return digest.hexdigest()

    assert fingerprint() == fingerprint()
