"""Check in-plane rotation of the legacy projected head outline measurement."""
from pathlib import Path
import importlib.util
import sys

import cv2
import numpy as np
import pandas as pd

def measure(mask):
    y, x = np.nonzero(mask)
    return (np.ptp(y)+1)/(np.ptp(x)+1)

def audit():
    root = Path(__file__).resolve().parents[2]
    path = root/'legacy/ch1_global/v2/52_measure_primary_traits_continuous.py'
    spec = importlib.util.spec_from_file_location('archived_shape',path)
    source = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(source)
    rows = []
    for angle in range(0,181,15):
        mask = np.zeros((301,301),np.uint8)
        cv2.ellipse(mask,(150,150),(72,36),angle,0,360,1,-1)
        rotated, _ = source.rotate_mask_to_major_axis(mask)
        y,x = np.nonzero(mask)
        points = np.column_stack([x,y]).astype(np.float32)
        _, vec = cv2.PCACompute(points, mean=None, maxComponents=2)
        major = vec[0]
        # Comparison control: direct projection, without image resampling.
        ext = np.ptp(points @ vec.T,axis=0)+1
        rows.append(dict(input_angle=angle,archived_aspect=measure(rotated),geometric_long_short=ext[0]/ext[1]))
    return pd.DataFrame(rows)

if __name__=='__main__':
    data=audit()
    assert data.archived_aspect.max()/data.archived_aspect.min()>3.5
    assert data.geometric_long_short.max()/data.geometric_long_short.min()<1.06
    data.to_csv(sys.argv[1],index=False)
