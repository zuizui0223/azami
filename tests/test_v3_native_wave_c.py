import json
import re

from analysis.v3.cloud_measurement_chunks import RUNTIME
from analysis.v3.measurement_chunks import contract
from analysis.v3.workflow import ROOT, canonical_digest


def read(path):
    return json.loads((ROOT/path).read_text(encoding='utf-8'))


def test_wave_b_completed_receipt_reconciles_cloud_and_local_all27():
    done = read('reproducibility/v3_native_raw_wave_b_20260908.json')
    assert done['run_id'] == 34204622177 and done['run_attempt'] == 2
    assert done['run_conclusion'] == 'success'
    assert len(done['records']) == 16
    assert done['batch_contract_canonical_sha256'] == canonical_digest(read(done['batch_contract']))
    for r in done['records']:
        c, v = r['cloud'], r['local_restoration']
        assert c['output_asset'] == v['asset']
        assert c['detected_heads'] == v['detected_heads']
        assert c['all27_raw_slots'] == v['raw_endpoint_slots'] == 27*v['detected_heads']
        assert c['bbox_slots'] == v['bbox_slots'] == 25*v['detected_heads']
        assert c['request_slots'] == v['verified_photo_units']
        assert c['new_photo_units_executed'] + c['units_restored_without_requests'] == c['request_slots']
        assert c['environment_values_read'] == c['ecological_models_executed'] == v['local_image_requests'] == 0
    assert done['aggregate']['verified_photo_units'] == sum(r['cloud']['request_slots'] for r in done['records']) == 1000
    assert done['aggregate']['selected_observations'] == sum(r['cloud']['selected_observations'] for r in done['records']) == 711
    assert done['aggregate']['detected_heads'] == sum(r['cloud']['detected_heads'] for r in done['records']) == 1684
    assert done['recovery']['completed_units_restored_without_requests'] == 250
    assert done['recovery']['new_requests_for_those_completed_units'] == 0


def test_wave_c_extends_verified_source_order_without_new_measurement_policy():
    batch = read('analysis/v3/native_measurement_wave_20260908_c.json')
    previous = read(batch['previous_completion_receipt'])
    assert batch['previous_completion_receipt_canonical_sha256'] == canonical_digest(previous)
    assert batch['plan_id'] == previous['plan_id']
    assert list(batch['chunks']) == [f'c{i:06d}' for i in range(18,146)]
    assert not set(batch['chunks']) & {r['cloud']['chunk_id'] for r in previous['records']}
    assert batch['measurement_contract_canonical_sha256'] == canonical_digest(contract())
    assert batch['runtime_contract_canonical_sha256'] == canonical_digest(json.loads(RUNTIME.read_text()))
    assert batch['maximum_parallel_chunks'] == 4
    assert batch['local_input_roundtrip_verified'] is True
    assert batch['input_asset']['files'] == 129 and not batch['input_asset']['source_images_included']
    for key in ('observations','requests','photo_jobs'):
        assert batch['aggregate'][key] == sum(row[key] for row in batch['chunks'].values())
    assert batch['aggregate'] == {'observations':5784,'requests':8086,'photo_jobs':9944}
    assert all(row['observations'] <= 128 and row['requests'] <= 64 for row in batch['chunks'].values())
    workflow = (ROOT/'.github/workflows/ch1-v3-native-raw-wave-c.yml').read_text()
    matrix = re.search(r'chunk: \[([^]]+)\]', workflow).group(1)
    assert [x.strip() for x in matrix.split(',')] == list(batch['chunks'])
    assert 'max-parallel: 4' in workflow and 'cancel-in-progress: false' in workflow
    assert '--batch analysis/v3/native_measurement_wave_20260908_c.json' in workflow
    assert 'analysis/v3/native_measurement_wave_20260908_b.json' not in workflow
