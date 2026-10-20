import numpy as np
import pandas as pd
from analysis.v4.native_integration_audit import eligible, environment_levels


def test_overlap_requires_only_requested_traits_and_their_qc():
    frame = pd.DataFrame({"a": [1, 2, np.nan], "a__eligible": [True]*3,
                          "b": [4, 5, 6], "b__eligible": [False, True, True]})
    assert eligible(frame, ["a"]).tolist() == [True, True, False]
    assert eligible(frame, ["a", "b"]).tolist() == [False, True, False]


def test_environment_decomposition_is_population_based_and_species_equal():
    frame = pd.DataFrame({"accepted_key": [1, 1, 1, 2],
                          "rain": [100., 200., 600., 400.]})
    levels = environment_levels(frame, ["rain"])
    assert levels["between_equal_species"]["rain"].tolist() == [300., 400.]
    assert levels["within_observations"]["rain"].tolist() == [-200., -100., 300., 0.]
