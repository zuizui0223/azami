"""Render the complete core construct atlas without refitting or selecting rows."""
from __future__ import annotations

import argparse
import ast
import hashlib
import json
import math
from pathlib import Path

import pandas as pd
from PIL import Image, ImageDraw, ImageFont

ROOT = Path(__file__).resolve().parents[1]
ENVS = ['chelsa_bio01', 'chelsa_bio04', 'chelsa_bio12', 'chelsa_bio15',
        'chelsa_rsds_mean', 'chelsa_vpd_mean', 'chelsa_sfcwind_mean', 'chelsa_gsp', 'chelsa_npp']
LABELS = ['BIO1', 'BIO4', 'BIO12', 'BIO15', 'Radiation', 'VPD', 'Wind', 'GSP', 'NPP']


def registry():
    tree = ast.parse((ROOT / 'analysis/v3/run_biological_axis_reanalysis.py').read_text(encoding='utf-8-sig'))
    expr = next(n.value for n in tree.body if isinstance(n, ast.Assign)
                and any(isinstance(t, ast.Name) and t.id == 'CONSTRUCTS' for t in n.targets))
    # Read the literal registry; never import or execute the analysis pipeline.
    return eval(compile(ast.Expression(expr), '<registry>', 'eval'), {'__builtins__': {}, 'dict': dict})


def validate(frame, constructs):
    expected = {(c, e) for c in constructs if c != 'surface_specularity' for e in ENVS}
    actual = list(zip(frame.construct_id, frame.predictor))
    if len(actual) != len(expected) or set(actual) != expected:
        raise ValueError('Expected the complete unique 90-row construct family')
    for row in frame.itertuples():
        value = row.effect_magnitude if constructs[row.construct_id]['kind'] == 'joint' else row.beta_std
        if not math.isfinite(float(value)) or not math.isfinite(float(row.q_bh)) or not 0 <= row.q_bh <= 1:
            raise ValueError('Non-finite effect or invalid q value in the input atlas')


def render(axis_dir: Path, output: Path, font_path: Path):
    output = output.resolve()
    frozen = ROOT / 'reproducibility/figures'
    if output == frozen or frozen in output.parents:
        raise ValueError('Output must not overwrite the frozen figure archive')
    constructs = registry()
    core = [c for c in constructs if not c.startswith('surface_')]
    assert len(constructs) == 11 and len(core) == 9
    inputs = [axis_dir / 'biological_axes_within.csv', axis_dir / 'biological_axes_among_min5.csv']
    frames = [pd.read_csv(p) for p in inputs]
    for frame in frames:
        validate(frame, constructs)
    image = Image.new('RGB', (2160, 1850), 'white')
    draw = ImageDraw.Draw(image)
    font = lambda n: ImageFont.truetype(str(font_path), n)
    for pi, (name, frame) in enumerate(zip(['Within taxa', 'Among taxa'], frames)):
        top = 70 + pi * 850
        draw.text((40, top), f'({chr(97 + pi)}) {name}', font=font(48), fill='black')
        for j, label in enumerate(LABELS):
            draw.text((610 + j * 165 + 82, top + 88), label, font=font(30), fill='black', anchor='mm')
        for i, c in enumerate(core):
            y = top + 130 + i * 72
            joint = constructs[c]['kind'] == 'joint'
            label = c.replace('_', ' ').capitalize() + (' [joint]' if joint else '')
            draw.text((35, y + 35), label, font=font(37), fill='black', anchor='lm')
            for j, e in enumerate(ENVS):
                row = frame[(frame.construct_id == c) & (frame.predictor == e)].iloc[0]
                value = float(row.effect_magnitude if joint else row.beta_std)
                background = '#eeeeee' if joint else ('#e8f0f7' if value >= 0 else '#faeee2')
                x = 610 + j * 165
                draw.rectangle((x, y, x + 160, y + 67), fill=background)
                text = f'{value:.3f}' + ('*' if row.q_bh < .05 else '')
                draw.text((x + 80, y + 34), text, font=font(32), fill='black', anchor='mm')
        draw.text((610, top + 795), '* BH q < 0.05 in the 90-test family; before sensitivity checks', font=font(30), fill='black')
    draw.text((40, 1800), 'Scalar cells: signed standardized slopes. Joint cells: non-negative vector magnitude.', font=font(31), fill='black')
    output.mkdir(parents=True, exist_ok=True)
    destination = output / 'Figure_4_construct_environment.png'
    image.save(destination, dpi=(330, 330))
    sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
    receipt = {'refitted': False, 'rows_per_family': 90, 'displayed_core_constructs': 9,
               'inputs': {p.name: sha(p) for p in inputs}, 'font_sha256': sha(font_path),
               'output_sha256': sha(destination), 'document_pagination_validated': False}
    (output / 'construct_environment_receipt.json').write_text(json.dumps(receipt, indent=2) + '\n', encoding='utf-8')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--axis-dir', type=Path, default=ROOT / 'work/current/results/axes')
    parser.add_argument('--output', type=Path, default=ROOT / 'work/layout-revisions')
    parser.add_argument('--font', type=Path, required=True, help='Path to licensed Arial TTF used by the manuscript; not redistributed')
    args = parser.parse_args()
    render(args.axis_dir, args.output, args.font)
