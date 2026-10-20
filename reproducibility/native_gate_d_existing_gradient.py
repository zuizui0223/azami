"""Re-audit corrected PR143 native *gradient* aggregate, not Gate D density.

Input: actual GitHub Actions 39422378058 'environment_response.json'
(SHA-256 locked). Environmental affine decode and invalid GSP masking were
already performed by PR143. Does NOT refit models or test climate-associated
changes of covariance/TPD/W2. This is not a new azami GEB V3 conclusion.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import math
from pathlib import Path

ORIGINAL_RUN=39422378058
ORIGINAL_ARTIFACT=12263957798
SHA256="563bead7d30ca59cd00b5c435457a654de53656f5f74e84e1d164b7dc330341a"
CORRECTED_INPUT_SHA256="342d78e9a8a960e8d9bdcf0bc2cdacd13bee48214a8b962419d5a71d5aef46ed"
MODELS=("colour","orientation","head_outline","bract_projection")
CLIMATE=("chelsa_bio04","chelsa_bio12","chelsa_bio15","chelsa_rsds_mean",
         "chelsa_vpd_mean","chelsa_sfcwind_mean","chelsa_gsp","chelsa_npp")
EXPECTED_JOINT_FAMILY=64


def verify_corrected_report(path:Path)->dict:
    raw=path.read_bytes()
    if hashlib.sha256(raw).hexdigest()!=SHA256:
        raise ValueError("Wrong or modified corrected PR143 numerical artifact")
    value=json.loads(raw)
    if value.get("input_sha256")!=CORRECTED_INPUT_SHA256:
        raise ValueError("Unknown environmental source input")
    if str(value.get("run_id"))!=str(ORIGINAL_RUN):
        raise ValueError("Unexpected workflow run")
    return value


def summarize_existing_gradients(data:dict)->dict:
    if data.get("input_sha256")!=CORRECTED_INPUT_SHA256:
        raise ValueError("Not source-corrected native input")
    if data.get("family_size")!=EXPECTED_JOINT_FAMILY:
        raise ValueError("Original Holm family changed")
    if data.get("inference_status")!="estimated":
        raise ValueError("Not an executed model report")
    if set(data.get("predictors",[]))!=set(CLIMATE) or len(data["predictors"])!=8:
        raise ValueError("Environmental candidates were silently changed")
    if set(data.get("models",{}))!=set(MODELS):
        raise ValueError("Model family changed")

    expected={name+"__"+scale for name in CLIMATE for scale in ("within","between")}
    found={};joint={}
    for model in MODELS:
        tests=data["models"][model]["tests"]
        if set(tests)!=expected:
            raise ValueError("Missing or substituted predictor/scale joint test")
        count={"within":0,"between":0}
        found[model]=count
        selected=[]
        for name,result in tests.items():
            if result.get("status")!="estimated":
                raise ValueError("Non-estimated joint contrast")
            p=result.get("p");holm=result.get("holm_p")
            if not isinstance(p,(int,float)) or not isinstance(holm,(int,float)):
                raise ValueError("Missing numerical p/adjustment")
            if not (math.isfinite(p) and math.isfinite(holm) and 0<=p<=1 and 0<=holm<=1):
                raise ValueError("Invalid adjusted p")
            if holm+1e-12<p:
                raise ValueError("Holm adjustment unexpectedly less than raw p")
            if holm<.05:
                scale=name.rsplit("__",1)[1]
                count[scale]+=1
                selected.append({"predictor_scale":name,"holm_p":float(holm)})
        joint[model]=sorted(selected,key=lambda x:x["predictor_scale"])
    if sum(sum(v.values()) for v in found.values())>EXPECTED_JOINT_FAMILY:
        raise ValueError("Impossible joint test count")
    orientation=data["models"]["orientation"]["tests"]
    if orientation["chelsa_bio12__between"]["holm_p"] < .05 or orientation["chelsa_bio12__within"]["holm_p"] < .05:
        raise ValueError("Disproved historical source claim about corrected rain and angle")
    return {
        "status":"EXECUTED_CORRECTED_PR143_GRADIENT_AUDIT_NOT_NEW_GATE_D",
        "workflow_run":ORIGINAL_RUN,
        "artifact_id":ORIGINAL_ARTIFACT,
        "input_sha256":CORRECTED_INPUT_SHA256,
        "report_sha256":SHA256,
        "joint_holm_tests":EXPECTED_JOINT_FAMILY,
        "scale_specific_supported_counts":found,
        "total_supported":sum(sum(v.values()) for v in found.values()),
        "individual_supported_joint_tests":joint,
        "between_species_rain_orientation_holm_p":float(orientation["chelsa_bio12__between"]["holm_p"]),
        "within_species_rain_orientation_holm_p":float(orientation["chelsa_bio12__within"]["holm_p"]),
        "inference_limitations":[
            "This retrospective gradient report is a different estimand from Gate D covariance/TPD/W2 environmental reorganization.",
            "Species-random-intercept observational fit; geographic and observer holdouts incomplete.",
            "Fixed downloaded photo phenotypes are image-relative, uncalibrated and viewpoint dependent.",
            "Multivariate joint tests do not identify individual trait coefficients as independently supported.",
            "Bract-contour projection is not botanically validated spine function.",
            "No phylogenetic correction, time-direction, causal adaptation or genetic inference.",
            "No new ecological fit performed; original GEB V3 frozen numerical outputs unchanged.",
        ]
    }


def main()->None:
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--report-json",required=True,type=Path)
    p.add_argument("--out",required=True,type=Path)
    a=p.parse_args()
    report=summarize_existing_gradients(verify_corrected_report(a.report_json))
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(report,indent=2,ensure_ascii=False,allow_nan=False)+"\n",encoding="utf-8")
    print(json.dumps({"status":report["status"],
                      "scale_specific_supported_counts":report["scale_specific_supported_counts"]}))


if __name__=="__main__":
    main()
