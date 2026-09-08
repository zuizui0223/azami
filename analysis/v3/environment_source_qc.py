"""Preserve raw environmental acquisition and make an explicit source-QC view."""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from .protected_artifacts import new_json, require
from .workflow import ROOT, canonical_digest, digest, text_digest

CONTRACT = ROOT / 'analysis/v3/environment_source_qc_contract_20260908.json'


def apply_qc(matrix: Path, out: Path, contract_path: Path = CONTRACT):
    require(not out.exists(), 'Preserve previous environmental QC output')
    contract = json.loads(contract_path.read_text(encoding='utf-8'))
    require(contract['status'] == 'source_value_qc_before_trait_join_not_model_selection'
            and contract['ecological_fitting_authorized'] is False, 'Unknown source-QC authority')
    require(digest(matrix) == contract['input_matrix_sha256'], 'Raw environment matrix identity differs')
    raw = pd.read_csv(matrix, dtype={'obs_id':str, 'accepted_key':str}, float_precision='round_trip')
    require(len(raw) == contract['observations'] and raw['obs_id'].is_unique
            and raw['obs_id'].notna().all() and raw['accepted_key'].nunique() == contract['taxa'], 'Source identity/counts differ')
    data, counts = raw.copy(), {}
    for rule in contract['rules']:
        variable = rule['variable']
        require(variable in {'GSP','BIO12'} and variable not in counts, 'Unexpected or duplicate QC variable')
        mask = raw[variable].eq(rule['exact_scaled_value'])
        require(int(mask.sum()) == rule['expected_rows'], 'Source anomaly counts differ')
        require(not any(variable+s in data for s in ('_source_raw','_source_qc')), 'QC columns already present')
        data[variable+'_source_raw'] = raw[variable]
        data[variable+'_source_qc'] = np.where(mask, rule['flag'], 'no_rule_flag')
        data.loc[mask, variable] = np.nan
        counts[variable] = int(mask.sum())
    unchanged = [c for c in raw if c not in counts]
    pd.testing.assert_frame_equal(raw[unchanged], data[unchanged])
    out.mkdir(parents=True)
    target = out / 'environment_candidate_matrix_qc_private.csv'
    data.to_csv(target, index=False, mode='x', lineterminator='\n')
    reread = pd.read_csv(target, dtype={'obs_id':str,'accepted_key':str}, float_precision='round_trip')
    pd.testing.assert_frame_equal(data, reread, check_exact=False, rtol=1e-15, atol=1e-15)
    report = {'schema_version':1, 'status':'ENVIRONMENT_SOURCE_QC_VIEW_VERIFIED_NO_ECOLOGY',
              'raw_matrix_sha256':contract['input_matrix_sha256'], 'qc_matrix_sha256':digest(target),
              'contract_canonical_sha256':canonical_digest(contract),
              'implementation_sha256_text_lf':text_digest(Path(__file__)),
              'source_rows':len(data), 'source_taxa':data['accepted_key'].nunique(),
              'flagged_working_values_set_missing':counts,
              'raw_values_retained':True, 'source_rows_deleted':0, 'ecological_models_executed':0,
              'ecological_fitting_authorized':False,
              'limits':['This is invalid/ceiling-value QC, not provider-confirmed missingness, interpolation or ecological inference.']}
    new_json(out/'source_qc_report.json', report)
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--matrix', type=Path, required=True)
    parser.add_argument('--out-dir', type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(apply_qc(args.matrix, args.out_dir)))


if __name__ == '__main__':
    main()
