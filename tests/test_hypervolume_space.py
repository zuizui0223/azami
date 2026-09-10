import sys
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'analysis/v3'))
from plot_hypervolume_space import region


def test_region_probability_and_determinism():
    sample=np.random.default_rng(8).normal(size=(50,3))
    cutoff,mass=region(sample,.57,'test')
    assert np.isfinite(cutoff) and .87<mass<.93
    assert (cutoff,mass)==region(sample,.57,'test')
