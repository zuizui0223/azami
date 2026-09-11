"""Compare every archived aggregate reference against a complete numerical replay."""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

REFERENCE = Path(__file__).resolve().parent/'current_reference'


def compare(expected, actual, path='root'):
    if isinstance(expected, dict):
        if not isinstance(actual, dict) or expected.keys() != actual.keys():
            raise AssertionError(f'{path}: object keys differ')
        for key in expected:
            compare(expected[key], actual[key], f'{path}.{key}')
    elif isinstance(expected, list):
        if not isinstance(actual, list) or len(expected) != len(actual):
            raise AssertionError(f'{path}: array length differs')
        for i,(a,b) in enumerate(zip(expected,actual)):
            compare(a,b,f'{path}[{i}]')
    elif isinstance(expected,bool) or expected is None:
        if expected != actual:
            raise AssertionError(f'{path}: {actual!r} != {expected!r}')
    elif isinstance(expected,(int,float)):
        if not isinstance(actual,(int,float)) or not math.isclose(expected,actual,rel_tol=1e-8,abs_tol=1e-10):
            raise AssertionError(f'{path}: {actual!r} != {expected!r}')
    elif expected != actual:
        # CSV numeric fields and JSON-encoded component vectors retain tolerances.
        try:
            a,b=json.loads(expected),json.loads(actual)
        except (ValueError,TypeError):
            raise AssertionError(f'{path}: {actual!r} != {expected!r}') from None
        compare(a,b,path)


def load(path):
    if path.suffix=='.json':
        return json.loads(path.read_text(encoding='utf-8'))
    with path.open(encoding='utf-8',newline='') as handle:
        return list(csv.reader(handle))


def validate(results):
    manifest=json.loads((REFERENCE/'manifest.json').read_text(encoding='utf-8'))
    checked=[]
    for row in manifest['files']:
        ref=REFERENCE/row['path']
        if hashlib.sha256(ref.read_bytes()).hexdigest()!=row['sha256']:
            raise AssertionError(f'Reference bytes changed: {row["path"]}')
        actual=results/row['path']
        if not actual.is_file():
            raise FileNotFoundError(actual)
        compare(load(ref),load(actual),row['path'])
        checked.append(row['path'])
    return {'status':'PASS','reference_commit':manifest['scientific_reference_commit'],
            'aggregate_files_compared':len(checked),'files':checked,
            'numerical_tolerance':{'relative':1e-8,'absolute':1e-10},
            'scope':'Numerical replay of existing models; not independent biological validation or public archive certification'}


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--results',required=True,type=Path)
    args=p.parse_args()
    report=validate(args.results)
    (args.results/'validation.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(report,indent=2))


if __name__=='__main__': main()
