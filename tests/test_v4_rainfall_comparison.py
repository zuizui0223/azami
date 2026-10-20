import numpy as np
import pandas as pd

from analysis.v4.rainfall_comparison import RAIN, TRAIT, comparison_frame


def test_same_cohort_keeps_population_means_and_raw_responses():
    frame = pd.DataFrame({
        "obs_id": [1, 2, 3, 4], "accepted_key": [1, 1, 1, 2],
        RAIN: [100., 200., 600., 400.], "other": [1., np.nan, 3., 4.],
        TRAIT: [17., 28., 39., 50.], TRAIT+"__eligible": [True, True, False, True]})
    sample = comparison_frame(frame, [RAIN, "other"])
    assert sample.obs_id.tolist() == [1, 4]
    assert sample[TRAIT].tolist() == [17., 50.]
    # Population mean includes observations excluded from this trait's fit.
    assert sample[RAIN+"__between"].tolist() == [300., 400.]
    assert sample[RAIN+"__within"].tolist() == [-200., 0.]
    assert RAIN+"__between" not in frame
