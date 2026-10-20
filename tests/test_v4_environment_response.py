import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

from analysis.v4.environment_candidate_audit import CANDIDATES
from analysis.v4.environment_response import execute
from analysis.v4.observation_models import MODULES


def test_fixed_workflow_keeps_response_specific_missingness(tmp_path):
    rng=np.random.default_rng(384)
    g=np.repeat(np.arange(35),10)
    frame=pd.DataFrame({'obs_id':np.arange(len(g)),'accepted_key':g})
    for v in CANDIDATES:
        frame[v]=rng.normal(size=35)[g]+rng.normal(size=len(g))
    crosswalk=json.loads((Path(__file__).parents[1]/'analysis/v4/trait27_functional_crosswalk_v1.json').read_text())
    primary=set(sum(MODULES.values(),[]))
    for entry in crosswalk['endpoints']:
        name=entry['id']
        frame[name]=rng.normal(size=len(g)) if name in primary else np.nan
        frame[name+'__eligible']=name in primary
    frame.loc[0,'orientation_image_vertical_angle']=np.nan
    frame.loc[0,'orientation_image_vertical_angle__eligible']=False
    receipt={'native_observations':len(frame),'endpoints':{}}
    for entry in crosswalk['endpoints']:
        name=entry['id']
        n=int(frame[name+'__eligible'].sum())
        receipt['endpoints'][name]={'finite':n,'eligible':n,
                                   'eligible_with_four_climate':n,'source':'synthetic'}
    table=tmp_path/'input.csv.gz'
    frame.to_csv(table,index=False,compression='gzip')
    receipt['prepared_table_sha256']=hashlib.sha256(table.read_bytes()).hexdigest()
    receipt_path=tmp_path/'receipt.json'
    receipt_path.write_text(json.dumps(receipt))
    selected=[v for v in CANDIDATES if v!='chelsa_vpd_mean']
    selection_path=tmp_path/'selection.json'
    selection_path.write_text(json.dumps({
        'status':'selected_before_trait_models','selected':selected,
        'input':{'input_sha256':receipt['prepared_table_sha256'],'candidates':CANDIDATES}}))
    result=execute(table,receipt_path,tmp_path/'output',selection_path)
    assert result['inference_status']=='estimated'
    assert set(result['models'])==set(MODULES)
    assert result['predictors']==selected
    assert sum(len(m['tests']) for m in result['models'].values())==30
    assert result['models']['orientation']['components']['orientation_image_vertical_angle']['n']==349
    assert result['models']['colour']['components']['corolla_lab_lightness']['n']==350
    assert result['models']['bract_projection']['components']['bract_projection_p95']['n']==350
    assert (tmp_path/'output/environment_response_ja.md').exists()
