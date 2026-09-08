import json

import numpy as np
import pandas as pd
import pytest

from analysis.v3 import environment_model as model


def test_source_selected_model_is_one_joint_four_process_design():
    spec=model.definition()
    assert spec['formulations']==1 and len(spec['variables'])==9
    assert spec['processes']['wetting_moisture']==['pr_month','BIO12','BIO18','GSP']
    assert spec['processes']['heat_drying']==['vpd_month','tas_month','BIO1']
    assert spec['transformation_source_rows']==317986
    assert spec['transformation_source_taxa']==354
    assert not spec['ecological_fitting_authorized']
    assert len(model.test_family())==len(set(model.test_family()))==36
    assert sum(m=='orientation' for m,p,q in model.test_family())==12


def test_saved_transform_is_identical_on_subsets_and_scales():
    spec=model.definition()
    data=pd.DataFrame({v:[spec['centers'][v],spec['centers'][v]+spec['scales'][v]] for v in spec['variables']})
    data['trait_must_not_be_read']=['unreadable','unreadable']
    full=model.transform(data)
    assert np.allclose(full,np.array([[0]*9,[1]*9]))
    assert np.array_equal(full[1:],model.transform(data.iloc[1:]))
    data.loc[0,'pr_month']=np.nan
    with pytest.raises(ValueError,match='no silent row deletion'):
        model.transform(data)


def test_missing_test_probabilities_keep_full_denominator_without_invented_values():
    family=model.test_family()
    values=dict.fromkeys(family)
    values[family[0]]=0.001
    values[family[1]]=0.002
    result=model.holm_complete_family(values)
    assert result[family[0]]==pytest.approx(0.036)
    assert result[family[1]]==pytest.approx(0.070)
    assert result[family[2]] is None
    values.pop(family[-1])
    with pytest.raises(ValueError,match='complete planned family'):
        model.holm_complete_family(values)


def test_changed_candidate_receipt_is_rejected(tmp_path):
    for file in (model.RECEIPT,model.SELECTION):
        path=tmp_path/file
        path.parent.mkdir(parents=True,exist_ok=True)
        path.write_text((model.ROOT/file).read_text())
    path=tmp_path/model.RECEIPT
    receipt=json.loads(path.read_text())
    receipt['candidate_selection']['selected'].append('v2_favourable_variable')
    path.write_text(json.dumps(receipt))
    with pytest.raises(ValueError,match='selection changed'):
        model.definition(tmp_path)
