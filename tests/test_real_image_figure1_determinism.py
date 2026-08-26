import importlib.util
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "analysis" / "build_real_image_figure1.py"
SPEC = importlib.util.spec_from_file_location("build_real_image_figure1", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def test_svg_export_is_byte_deterministic(tmp_path: Path) -> None:
    fig, ax = plt.subplots()
    ax.imshow(np.arange(16, dtype=np.uint8).reshape(4, 4))
    first = tmp_path / "first.svg"
    second = tmp_path / "second.svg"

    MODULE.save_svg(fig, first)
    MODULE.save_svg(fig, second)
    plt.close(fig)

    assert first.read_bytes() == second.read_bytes()
