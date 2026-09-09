import json
from pathlib import Path
from analysis.v3.workflow import canonical_digest,text_digest

ROOT = Path(__file__).resolve().parents[1]


def read(relative):
    return json.loads((ROOT/relative).read_text(encoding='utf-8'))


def test_current_solver_is_bound_to_returned_final_geometry_smoke_not_old_replay():
    proof = read('reproducibility/v3_final_geometry_smoke_current_code_20260909.json')
    smoke = read(proof['source_receipt'])
    assert canonical_digest(smoke)==proof['source_receipt_canonical_sha256']
    assert proof['source_workflow_run_id']==smoke['source_workflow_run_id']==34248192547
    assert proof['source_head_sha']==smoke['source_head_sha']
    for path,sha in proof['implementation_sha256_text_lf'].items():
        assert sha==text_digest(ROOT/path)
    assert proof['implementation_sha256_text_lf']['analysis/v3/joint_partial_pooling.py']=='c506afa529f2a3ac125bbd72e0334cf7e125b774ef80ae50fb91b0f764c2a15c'
    assert proof['independent_generated_datasets']==4 and proof['grid_cases']==8
    assert proof['verified_shared_draw_records']==1592
    returned = {r['scenario']:r for r in proof['records']}
    assert set(returned)==set(smoke['scenarios'])
    for entry in smoke['input_manifest']:
        record = returned[entry['scenario']]
        assert record['public_report_sha256']==entry['sha256']
        assert record['artifact_path']==entry['artifact_path']
        assert {r['grid_degrees'] for r in record['cases']}=={2,5}
        assert all(r['draws_verified']==199 for r in record['cases'])
    assert proof['ecological_fitting_authorized'] is False
    assert proof['ecological_models_executed']==proof['empirical_trait_environment_values_read']==0
