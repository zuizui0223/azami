from pathlib import Path
import pandas as pd
from analysis.v3.build_capitulum_organization_synthesis import build


def test_two_questions_preserve_all_pairs_and_existing_results(tmp_path):
    root = Path(__file__).resolve().parents[1]
    summary = build(root / 'analysis_outputs/capitulum_organization_20260910/source',
                    root / 'analysis_outputs/environment_first_20260910/trait_environment_map.csv', tmp_path)
    pairs = pd.read_csv(tmp_path / 'organization_environment_pairs.csv')
    assert summary['all_pairs'] == 36
    assert summary['same_module_pairs'] == 7
    assert len(summary['environment_full_chain_rows']) == 3
    assert summary['common_cohort'] == dict(observations=1734, taxa=42, minimum_complete_observations_per_taxon=5)
    assert len(pairs) == 36
    assert summary['exploratory_bridge']['among_taxon']['qap_p_two_sided'] > .05
    assert summary['exploratory_bridge']['within_taxon']['qap_p_two_sided'] > .05
