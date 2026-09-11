"""Render two presentation-only revisions without changing frozen v2 figures.

The historical renderer and scientific tables remain unchanged. Exact source
substitutions fail if the renderer changes, rather than silently losing a fix.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import inspect
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def replace_once(source: str, old: str, new: str) -> str:
    if source.count(old) != 1:
        raise ValueError(f"Expected exactly one layout anchor: {old!r}")
    return source.replace(old, new)


def render(output: Path) -> None:
    output = output.resolve()
    frozen = ROOT / "reproducibility/figures"
    if output == frozen or frozen in output.parents:
        raise ValueError("Output must not overwrite the frozen figure archive")
    output.mkdir(parents=True, exist_ok=True)
    source_path = ROOT / "legacy/v2/figures/_build_v2_figures_impl.py"
    spec = importlib.util.spec_from_file_location("frozen_figure_layout", source_path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    source = inspect.getsource(module.figure_measurement)
    replacements = [
        ("source_image[2190:2860, 480:1000]", "source_image[1300:2010, 100:610]"),
        ('"(i) Orientation relative to image vertical"', '"(i) Image-vertical reference"'),
        ('orientation.axis("off")', 'orientation.axis("off")\n    orientation.annotate("", xy=(-0.10, 0.9), xytext=(-0.10, 0.1), xycoords="axes fraction", arrowprops=dict(arrowstyle="->", color="black", lw=1), annotation_clip=False)'),
        ("finite v2 outputs", "continuous measurements"),
        ('("Detector = crop", "does not score traits"', '("11 constructs", "biological interpretation"'),
        ("Three open-licensed photographs are presentation/source provenance only; the v2 cohort and results come exclusively from the prespecified full-27 lane.", "Source photographs illustrate measurement only. Arrow indicates image vertical, not measured gravity. Values are frozen production outputs."),
        ("photo 532420148; CC BY.", "photo 532420148; Bobo-X, CC BY."),
    ]
    for old, new in replacements:
        source = replace_once(source, old, new)
    exec(compile(source, str(source_path), "exec"), vars(module))
    module.figure_measurement(output)
    source = replace_once(inspect.getsource(module.figure_s3_spatial_diagnostics), ">= 2.5", ">= 2.2")
    exec(compile(source, str(source_path), "exec"), vars(module))
    module.figure_s3_spatial_diagnostics(module.DEFAULT_INPUT, output)
    stems = [module.MAIN_STEMS[0], module.SUPP_STEMS[2]]
    files = [output / (stem + suffix) for stem in stems for suffix in (".png", ".pdf")]
    receipt = {
        "data_refitted": False,
        "historical_source_commit": "1857f0ad398b8a5c41d6080443ef6146a675c185",
        "figure_1_angle": "Production CSV 0.732009 degrees rounded to 0.7; historical 1.1-degree overlay is not reused",
        "figure_s5_change": "Long label aligned inward; points, axes and P values unchanged",
        "source_renderer_sha256": hashlib.sha256(source_path.read_bytes()).hexdigest(),
        "outputs": {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in files},
        "document_pagination_validated": False,
    }
    (output / "layout_receipt.json").write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=ROOT / "work/layout-revisions")
    render(parser.parse_args().output)
