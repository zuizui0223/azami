"""One native A→B→C→D→E handoff, not five unrelated result families.

The same verified 313350-observation input, same nine measured coordinates and
taxon×cell×observer rarefaction support ALL descriptive sampling analyses.
Existing five-environment model and spatial constraints are reused unchanged.
All outputs are aggregate-only; the source image/observer roster never leaves
the authorised private runner. Experimental function and causal evolution remain
independent validation requirements, NOT inferred from photographs.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

from analysis.v4.native_trait_models import prepare
from analysis.v4.native_gates_ae import authority_mapping, gate_a, gate_e, CORE9
from analysis.v4.native_gates_bcd import gate_b, gate_c, gate_d_same_sample

EXPECTED_SHA="342d78e9a8a960e8d9bdcf0bc2cdacd13bee48214a8b962419d5a71d5aef46ed"
FIXED_PREDICTORS=[
    "chelsa_bio01","chelsa_bio04","chelsa_bio12",
    "chelsa_bio15","chelsa_rsds_mean"
]


def validate_inputs(
    table: Path, receipt_path: Path, models_path: Path,
    integration_path: Path, spatial_path: Path
) -> tuple[dict,dict,dict,dict]:
    input_sha=hashlib.sha256(table.read_bytes()).hexdigest()
    receipt=json.loads(receipt_path.read_text(encoding="utf-8"))
    model=json.loads(models_path.read_text(encoding="utf-8"))
    integration=json.loads(integration_path.read_text(encoding="utf-8"))
    spatial=json.loads(spatial_path.read_text(encoding="utf-8"))
    if input_sha != EXPECTED_SHA:
        raise ValueError("The corrected native source bytes do not match frozen gate input")
    if any(x!=input_sha for x in (
        receipt.get("prepared_table_sha256"),model.get("input_sha256"),
        integration.get("input_sha256"),spatial.get("input_sha256")
    )):
        raise ValueError("A-E stages do not share exactly the same decoded native table")
    model_sha=hashlib.sha256(models_path.read_bytes()).hexdigest()
    if integration.get("model_sha256")!=model_sha or spatial.get("model_sha256")!=model_sha:
        raise ValueError("Models, integration and spatial diagnostics from different runs")
    if receipt.get("environment_decoding_status") != "physical_units_gsp_invalid_encoding_masked_v1":
        raise ValueError("Decoded native source lacks audited physical environmental units")
    if model.get("selection",{}).get("selected")!=FIXED_PREDICTORS:
        raise ValueError("Environmental predictors differ from frozen native five-variable contract")
    if integration.get("common_cohort",{}).get("n")!=157654 or integration["common_cohort"]["species"]!=174:
        raise ValueError("Original B common cohort changed")
    if receipt.get("native_observations")!=313350 or len(receipt.get("endpoints",{}))!=27:
        raise ValueError("Native source or 27-endpoint inventory changed")
    return receipt,model,integration,spatial


def execute(
    table:Path,receipt_path:Path,models_path:Path,integration_path:Path,
    spatial_path:Path,mapping_path:Path,output:Path,
    observer_path:Path | None = None,*,b_bootstrap:int=200,
    c_replicates:int=12,max_taxa:int=40
)->dict:
    receipt,models,integration,spatial=validate_inputs(
        table,receipt_path,models_path,integration_path,spatial_path)
    frame=pd.read_csv(table,dtype={"obs_id":str,"accepted_key":str},low_memory=False)
    if len(frame)!=313350 or not frame.obs_id.is_unique:
        raise ValueError("Expected frozen 313350 unique source observations")
    observer = (pd.read_csv(observer_path,dtype={"obs_id":str,"user_id":str})
                if observer_path else None)
    mapping=authority_mapping(mapping_path)
    full=prepare(frame,FIXED_PREDICTORS)
    # Strict one-way information flow: taxon/observer/QC first; no climate
    # effect, KDE or bootstrap outcome used to decide which rows qualify.
    gateA,candidate=gate_a(full,mapping,observer,max_taxa=max_taxa)
    selected=set(candidate.accepted_key.unique())
    gateE=gate_e(mapping,currently_selected_keys=selected,native_frame=frame)
    output_dict={
        "schema_version":"unified_native_capitulum_abcde_20261021_v1",
        "status":"scientific_gates_partially_open_independent_field_validation_needed",
        "native_source_sha256":EXPECTED_SHA,
        "model_sha256":hashlib.sha256(models_path.read_bytes()).hexdigest(),
        "integration_sha256":hashlib.sha256(integration_path.read_bytes()).hexdigest(),
        "spatial_sha256":hashlib.sha256(spatial_path.read_bytes()).hexdigest(),
        "taxonomy_git_blob_sha1":"b662f9687f11d4069965531d5214962b07f6df8c",
        "source_denominators":{"native_observations":313350,"native_taxa":306,
                              "joint_nine_output":157779,"joint_nine_output_taxa":242},
        "shared_sampling_contract":{
            "same_acceptance":"all nine measured coordinates from SAME source observation",
            "joint_feature_concepts":8,
            "biological_axes":["presentation","visible_floret_colour","whole_head_outline"],
            "rare_taxa":"predeclared retained after excluding WCVP mixed-source concepts",
            "taxon_cell_observer_rarefaction":True,
            "no_new_pixels_or_environmental_reselection":True,
        },
        "A":gateA,
        "B":{"status":"not_run_until_A_observer_selection_verified"},
        "C":{"status":"not_run_until_A_observer_selection_verified"},
        "D":{
            "status":"PR144_full_native_models_already_executed",
            "five_selected_climate_variables":FIXED_PREDICTORS,
            "full_reference_model_cohort":{"observations":157654,"taxa":174},
            "full_reference_spatial_diagnostics":"conditional cell Moran I and fixed balanced sampling, not causal",
            "environmental_overlap_same_sample":{"status":"waiting_for_A_selection"},
            "phylogenetic_evidence":"52 frozen alternative placements, not independent trees",
            "covariance_evidence":"104 PSD-constrained estimates; boundary sensitivity, not genetic covariance",
        },
        "E":gateE,
        "limitations":[
            "Biological adaptive functions, genetic modularity and repeated gains require independent experiments,",
            "Photographic colour, shape, head angle and genus names are subject to measurement and classification error.",
            "Original human observers are not independent plant identities.",
        ],
    }
    if gateA["status"]!="observer_cell_balanced_selection_constructed":
        return _save(output,output_dict)
    if len(selected)<5:
        output_dict["B"]={"status":"blocked_insufficient_independently_supported_taxa"}
        output_dict["C"]={"status":"blocked_insufficient_independently_supported_taxa"}
        output_dict["D"]["environmental_overlap_same_sample"]={
            "status":"blocked_insufficient_independently_supported_taxa"}
        return _save(output,output_dict)

    # Gate B, C and D use the same frozen candidate rows and no competing
    # taxon/observer/region filters. Gate B significance floor >=20 taxa.
    if len(selected)>=20:
        output_dict["B"]=gate_b(candidate,models,integration,replicates=b_bootstrap)
    else:
        output_dict["B"]={
            "status":"insufficient_taxa_for_registered_20_taxon_multivariate_B",
            "available_taxa":len(selected)
        }
    output_dict["C"]=gate_c(candidate,replicates=c_replicates,max_taxa=max_taxa)
    output_dict["D"]["environmental_overlap_same_sample"]=gate_d_same_sample(
        candidate,full,models
    )
    output_dict["status"]="numerical_native_ABCDE_screen_complete_without_field_functional_gate"
    return _save(output,output_dict)


def _save(output:Path,result:dict)->dict:
    output.mkdir(parents=True,exist_ok=True)
    target=output/"unified_gate_abcde.json"
    target.write_text(json.dumps(result,indent=2,ensure_ascii=False,allow_nan=False)+"\n",
                      encoding="utf-8")
    print(json.dumps({
        "status":result["status"],
        "matched_nine_output":result["A"]["matched_nine_output_observations"],
        "selected_taxa":result["A"]["predeclared_selection"]["selected_taxa"],
        "B":result["B"]["status"],"C":result["C"]["status"],
        "E_collision_keys":result["E"]["collision_keys"],
    },ensure_ascii=False),flush=True)
    return result


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ("table","receipt","models","integration","spatial","mapping","output"):
        p.add_argument("--"+name,type=Path,required=True)
    p.add_argument("--observers",type=Path)
    p.add_argument("--b-bootstrap",type=int,default=200)
    p.add_argument("--c-replicates",type=int,default=12)
    p.add_argument("--max-taxa",type=int,default=40)
    a=p.parse_args()
    execute(a.table,a.receipt,a.models,a.integration,a.spatial,a.mapping,a.output,
            a.observers,b_bootstrap=a.b_bootstrap,c_replicates=a.c_replicates,
            max_taxa=a.max_taxa)


if __name__=="__main__":
    main()
