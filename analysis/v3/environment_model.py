"""One pinned process-block model definition, without reading any trait values.

This adopts the verified source-QC selection, not the rejected raw-value screen.
It supplies shared predictor transformations and the complete probability family;
neither source VIF nor these definitions authorizes ecological fitting.
"""
from __future__ import annotations

import itertools
import json
from pathlib import Path

import numpy as np

from .workflow import ROOT, canonical_digest

RECEIPT='reproducibility/v3_full_native_environment_qc_selection_20260908.json'
RECEIPT_SHA='1dbddaa1cfab0075e89fa3fed1b1934e7aaa5e25041678513611542bc346e2d8'
SELECTION='analysis/v3/environment_production_contract.json'
MODULES=('orientation','visible_colour','gross_shape')
PROCESSES=('wetting_moisture','radiation','heat_drying','mechanical')
QUESTIONS=('within','among','among_minus_within')


def definition(root: Path=ROOT):
    receipt=json.loads((root/RECEIPT).read_text(encoding='utf-8'))
    rule=json.loads((root/SELECTION).read_text(encoding='utf-8'))
    if canonical_digest(receipt)!=RECEIPT_SHA:
        raise ValueError('Pinned source-QC environment selection changed')
    selected=receipt['candidate_selection']
    if (selected['contract_canonical_sha256']!=canonical_digest(rule)
            or selected['status']!='FULL_NATIVE_ENVIRONMENT_REPRESENTATION_SELECTED_NO_ECOLOGY'
            or selected['matrix_sha256']!=receipt['source_qc']['qc_matrix_sha256']
            or selected['trait_values_read']!=0 or selected['ecological_models_executed']!=0):
        raise ValueError('Environment source, selection rule or outcome-blind boundary differs')
    variables=selected['selected']
    if ([v for process in PROCESSES for v in selected['processes'][process]]!=variables
            or len(set(variables))!=len(variables)
            or any(not selected['processes'][p] for p in PROCESSES)
            or not set(rule['selection']['protected'])<=set(variables)):
        raise ValueError('Retained process partition or representative differs')
    if (set(selected['centers'])!=set(variables) or set(selected['scales'])!=set(variables)
            or any(not np.isfinite(selected['centers'][v]) or not np.isfinite(selected['scales'][v])
                   or selected['scales'][v]<=0 for v in variables)):
        raise ValueError('Common predictor transformation is undefined')
    final=selected['trace'][-1]
    if (final['variables']!=variables or final['matrix']['matrix_rank']!=len(variables)
            or any(not isinstance(r['vif'],(int,float)) or r['vif']>=rule['selection']['threshold'] for r in final['vif'])):
        raise ValueError('Selected source design did not pass the declared redundancy screen')
    return {'status':'ONE_PROCESS_BLOCK_ENVIRONMENT_DEFINITION_VERIFIED_NOT_FITTED',
            'receipt':RECEIPT,'receipt_canonical_sha256':RECEIPT_SHA,
            'variables':variables,'processes':selected['processes'],
            'centers':selected['centers'],'scales':selected['scales'],
            'transformation_source_rows':selected['complete_rows'],
            'transformation_source_taxa':selected['complete_taxa'],
            'source_qc_matrix_sha256':selected['matrix_sha256'],
            'formulations':1,'primary_probability_slots':36,
            'ecological_fitting_authorized':False}


def transform(frame, *, root: Path=ROOT):
    """Read predictors only; caller must declare the actual model membership.

    Do not drop missing rows, refit centers by module/taxon, or impute values.
    Output columns follow the process order in the frozen source selection.
    """
    spec=definition(root)
    values=frame[spec['variables']].to_numpy(dtype=float)
    if not np.isfinite(values).all():
        raise ValueError('Declare complete model support before transformation; no silent row deletion')
    centers=np.array([spec['centers'][v] for v in spec['variables']])
    scales=np.array([spec['scales'][v] for v in spec['variables']])
    return (values-centers)/scales


def test_family():
    return tuple((m,p,q) for m,p,q in itertools.product(MODULES,PROCESSES,QUESTIONS))


def holm_complete_family(probabilities):
    """Reserve all 36 slots; unavailable probabilities stay explicitly None."""
    family=test_family()
    if set(probabilities)!=set(family):
        raise ValueError('Report exactly the complete planned family, including non-estimable slots')
    known=[]
    for key,value in probabilities.items():
        if value is None:
            continue
        if isinstance(value,bool) or not np.isfinite(value) or not 0<=value<=1:
            raise ValueError('Invalid primary probability')
        known.append((float(value),key))
    result={key:None for key in family}
    previous=0.0
    for rank,(value,key) in enumerate(sorted(known)):
        previous=max(previous,min(1.0,(len(family)-rank)*value))
        result[key]=previous
    return result
