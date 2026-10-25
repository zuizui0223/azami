import pandas as pd
import pytest

from reproducibility.render_construct_environment import ENVS, default_font_path, registry, validate


def test_complete_family_preserves_all_core_and_diagnostic_rows():
    constructs = registry()
    frame = pd.DataFrame([{'construct_id': c, 'predictor': e, 'effect_magnitude': 0.1,
                           'beta_std': -0.1, 'q_bh': 0.8}
                          for c in constructs if c != 'surface_specularity' for e in ENVS])
    validate(frame, constructs)
    assert len(frame) == 90
    for invalid in [frame.iloc[:-1], pd.concat([frame.iloc[:-1], frame.iloc[[0]]])]:
        with pytest.raises(ValueError, match='complete unique'):
            validate(invalid, constructs)
    frame.loc[0, 'q_bh'] = float('nan')
    with pytest.raises(ValueError, match='Non-finite'):
        validate(frame, constructs)

def test_default_construct_figure_font_is_resolvable():
    path = default_font_path()
    assert path.is_file()
    assert path.suffix.lower() in {".ttf", ".otf", ".ttc"}

