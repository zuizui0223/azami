#!/usr/bin/env python3
"""Source-verified descriptive Figure 2: photographic colour, angle and outline."""
import argparse
import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

FROZEN = {
    ("angle", "chroma"): (-0.22512839486857492, -0.028067356918921597),
    ("angle", "aspect"): (-0.24694620634514822, -0.0693870326972044),
    ("chroma", "aspect"): (0.01993972231831634, -0.01824691397495542),
}
