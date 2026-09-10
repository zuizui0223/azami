"""Join existing evidence for two questions; never select pairs by significance."""
import argparse
import hashlib
import json
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from analysis.v3.build_environment_first_synthesis import MODULES


def build(source, environment_map, out):
    report = json.loads((source / 'construct_scale_upgrade_report.json').read_text())
    pairs = pd.read_csv(source / 'complete18_construct_pairwise.csv')
    expected = {frozenset(p) for p in combinations(MODULES, 2)}
    actual = [frozenset(p) for p in pairs[['construct_left', 'construct_right']].itertuples(index=False, name=None)]
    assert len(pairs) == 36 and len(set(actual)) == 36 and set(actual) == expected
    pairs['module_left'] = pairs.construct_left.map(MODULES)
    pairs['module_right'] = pairs.construct_right.map(MODULES)
    pairs['same_module'] = pairs.module_left.eq(pairs.module_right)
    for scale in ['within_taxon', 'among_taxon']:
        suffix = scale.split('_')[0]
        similarity = pd.read_csv(source / f'environment_similarity_{suffix}.csv', index_col=0)
        assert set(similarity.index) == set(MODULES) == set(similarity.columns)
        assert np.isfinite(similarity.to_numpy()).all()
        assert np.allclose(similarity, similarity.T)
        pairs[f'{scale}_environment_magnitude_similarity'] = [
            similarity.loc[l, r] for l, r in pairs[['construct_left', 'construct_right']].itertuples(index=False, name=None)]
        values = pairs[f'{scale}_rv']
        assert np.isfinite(values).all() and values.between(0, 1).all()
        for same, key in [(True, 'mean_within_module_rv'), (False, 'mean_between_module_rv')]:
            assert np.isclose(values[pairs.same_module == same].mean(), report['module_cohesion'][scale][key])
        rho = spearmanr(values, pairs[f'{scale}_environment_magnitude_similarity']).statistic
        assert np.isclose(rho, report['integration_environment_coupling'][scale]['rho'])
    alignment = spearmanr(pairs.within_taxon_rv, pairs.among_taxon_rv).statistic
    assert np.isclose(alignment, report['common_cohort_matrix_alignment']['rho'])
    assert np.allclose(pairs.among_taxon_rv - pairs.within_taxon_rv, pairs.delta_among_minus_within)
    atlas = pd.read_csv(environment_map)
    assert len(atlas) == 162 and not atlas.duplicated(['scale', 'construct_id', 'predictor']).any()
    summary = {
        'analysis_type': 'replay and synthesis of existing PR93 evidence, no new tests',
        'upstream_run': 34432661967, 'upstream_artifact': 10135139012,
        'environment_question': 'Which biological traits covary with which gradients after declared sensitivity checks?',
        'environment_full_chain_rows': atlas.loc[atlas.robustness_status.eq('full_declared_chain_pass'), ['construct_id', 'predictor', 'scale', 'beta_std', 'q_bh']].to_dict('records'),
        'environment_interpretation': 'Inherited associations, not independent replication; chroma-NPP remains exploratory, not a frozen-v2 headline',
        'organization_question': 'How strongly are traits integrated within and across biological modules?',
        'common_cohort': report['common_cohort'],
        'all_pairs': len(pairs), 'same_module_pairs': int(pairs.same_module.sum()),
        'module_cohesion': report['module_cohesion'],
        'cross_scale_alignment': report['common_cohort_matrix_alignment'],
        'cross_scale_bootstrap': report['taxon_bootstrap'],
        'exploratory_bridge': report['integration_environment_coupling'],
        'supplementary_omnibus': 'separate question; retain unsupported results; not a test of modularity',
        'limits': ['No evolutionary reconstruction or demonstrated adaptation', 'Presentation is a singleton module', 'Magnitude-profile similarity does not test signed response alignment', 'Image and mathematical covariance may contribute', 'Cohorts differ between trait-specific environment tests and complete-case integration; the bridge uses these existing profiles, not a common-cohort environmental refit'],
        'source_sha256': {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(source.iterdir()) if p.is_file()},
        'environment_map_sha256': hashlib.sha256(environment_map.read_bytes()).hexdigest(),
    }
    out.mkdir(parents=True, exist_ok=True)
    pairs.to_csv(out / 'organization_environment_pairs.csv', index=False)
    (out / 'summary.json').write_text(json.dumps(summary, indent=2, allow_nan=False) + '\n')
    return summary


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--source', type=Path, required=True)
    parser.add_argument('--environment-map', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(build(args.source, args.environment_map, args.out), indent=2))
