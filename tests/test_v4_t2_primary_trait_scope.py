import json
from pathlib import Path
from analysis.v4 import s1s2_multitrait_branch_core as core
from analysis.v4.native_integration import FEATURES

def test_primary_scope_is_frozen_eight_constructs_not_inventory22():
    inventory=json.loads(Path('analysis/v4/trait27_functional_crosswalk_v1.json').read_text())
    reg=core.primary_trait_registry(inventory)
    units=core.primary_trait_tests(reg)
    assert len(units)==8
    assert set(reg['groups'])=={'presentation','flower_colour','whole_head_shape'}
    assert len(reg['measured'])==9
    assert set(reg['measured'])=={c for cols in FEATURES.values() for c in cols}
    assert len(units['construct:hue'])==2
    assert not any('bract' in c or 'involucre' in c for cols in units.values() for c in cols)
    assert 'block:flower_colour_joint' not in units

def test_primary_runner_has_only_three_domain_pairs():
    from analysis.v4.run_s1s2_multitrait_nuclear_branch import PAIRS
    assert len(PAIRS)==3
    assert set(PAIRS)=={('presentation','flower_colour'),('presentation','whole_head_shape'),('flower_colour','whole_head_shape')}
