"""Export environment-only native rows for selection before phenotype modelling."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd


def execute(table, receipt_path, output):
    receipt=json.loads(receipt_path.read_text())
    digest=hashlib.sha256(table.read_bytes()).hexdigest()
    if digest!=receipt['prepared_table_sha256']:
        raise ValueError('Prepared input mismatch')
    contract=Path(__file__).parents[1]/'ch1/capitulum_environment_blocks_contract.json'
    blocks=json.loads(contract.read_text())['blocks']
    columns=[v for block in blocks for v in block['predictors']]
    if len(blocks)!=6 or len(columns)!=9 or len(set(columns))!=9:
        raise ValueError('Expected six groups and all nine stored candidates')
    frame=pd.read_csv(table,usecols=columns)
    if frame.chelsa_gsp.eq(4294967295).any():
        raise ValueError('Undecoded GSP sentinel: correct environmental input before selection')
    if len(frame)!=receipt['native_observations']:
        raise ValueError('Native population changed')
    mask=np.isfinite(frame).all(axis=1)
    output.mkdir(parents=True,exist_ok=True)
    frame.loc[mask,columns].to_csv(output/'environment_only.csv',index=False)
    metadata={'input_sha256':digest,'group_definition_sha256':hashlib.sha256(contract.read_bytes()).hexdigest(),
              'groups':blocks,'candidates':columns,'native_n':len(frame),
              'selection_n':int(mask.sum()),'missing_any_candidate_n':int((~mask).sum()),
              'finite_counts':{v:int(np.isfinite(frame[v]).sum()) for v in columns},
              'value_audit':{
                  v:{
                      'quantiles':{
                          str(q):float(frame.loc[np.isfinite(frame[v]),v].quantile(q))
                          for q in (0,.001,.01,.05,.5,.95,.99,.999,1)},
                      'most_frequent_values':[
                          {'value':float(value),'count':int(count)}
                          for value,count in frame[v].value_counts().head(12).items()],
                      'integer_sentinel_counts':{
                          str(value):int(frame[v].eq(value).sum())
                          for value in (-9999,32767,65535,2147483647,4294967295)},
                  } for v in columns},
              'threshold':10,'keep':[],'sampling':'all complete native observation rows',
              'scope':'Environment-only selection; no trait values or eligibility used',
              'environment_csv_sha256':hashlib.sha256((output/'environment_only.csv').read_bytes()).hexdigest()}
    (output/'selection_input.json').write_text(json.dumps(metadata,indent=2)+'\n')


if __name__=='__main__':
    p=argparse.ArgumentParser()
    for name in ('table','receipt','output'):
        p.add_argument('--'+name,type=Path,required=True)
    a=p.parse_args()
    execute(a.table,a.receipt,a.output)
