"""Choose REAL external Asteraceae photo QC images for eventual S1 continuous-angle.

Phenotype-/climate-BLIND sampling: all species must already pass original GBIF
>=10 photo and >=3 geographic-cell SOURCE gate. Take <=5 source-public image
references per species with deterministic seed, no U/D, no angle, no climate
columns, no secret coordinates, no re-hosting copyrighted image bytes.

This packet is for actual original Cirsium YOLO transfer plus independent
capitulum/pedicel QC. A photo record with no original validated angle is never
converted to a phenotype by this script.
"""
from __future__ import annotations
import argparse,csv,hashlib,json,random
from pathlib import Path

SHA_JSON="ffb90a009b031522480f3b5618db2c7ad70050b993e6eee70f0a428e05cd20d2"
N_FULL=131
N_READY=23
N_PER_SPECIES=5
SEED=20261024

def prepare(result_path,photo_csv,output):
    raw=Path(result_path).read_bytes()
    if hashlib.sha256(raw).hexdigest()!=SHA_JSON:
        raise ValueError("Cannot change original T3 DNA-tip and photo four-climate sources")
    r=json.loads(raw)
    if r["original_genbank_ITS_biological_taxa"]!=N_FULL or r["n_species_qualify_for_MANUAL_photo_angle_audit"]!=N_READY:
        raise ValueError("Unexpected T3 all131 original source support")
    eligible={x["original_GenBank_ITS_taxon"] for x in r["by_species"]
              if x["candidate_for_manual_S1_angle_protocol_audit"]}
    if len(eligible)!=N_READY:
        raise ValueError("Original true GBIF + four-climate source quality filter changed")
    with Path(photo_csv).open(encoding="utf-8",newline="") as f:
        data=list(csv.DictReader(f))
    by={}
    seen=set()
    for row in data:
        n=row["source_taxon"]
        if n not in eligible:continue
        k=row["gbif_occurrence_key"]
        if not k or k in seen:
            raise ValueError("GBIF original photo key duplicated or empty")
        if not row["source_image_url"].startswith("https://") or not row["image_license_url"]:
            raise ValueError("Unlicensed or inaccessible original source image metadata")
        seen.add(k);by.setdefault(n,[]).append(row)
    if set(by)!=eligible:
        raise ValueError("Missing source image references from validated T3 quality-eligible taxa")
    rng=random.Random(SEED)
    qc=[]
    for taxon in sorted(eligible):
        rows=sorted(by[taxon],key=lambda x:x["gbif_occurrence_key"])
        rng.shuffle(rows)
        if len(rows)<N_PER_SPECIES:raise ValueError("T3 quality-gated source taxa lack five original images")
        for ordinal,row in enumerate(rows[:N_PER_SPECIES],1):
            qc.append({
                "blind_queue_id":f"T3_PHOTO_{len(qc)+1:04}",
                "DNA_tip_original_source_taxon":taxon,
                "gbif_occurrence_key":row["gbif_occurrence_key"],
                "source_public_image_url":row["source_image_url"],
                "source_license_url":row["image_license_url"],
                "original_rights_holder":row.get("source_rightsHolder",""),
                "head_present": "",
                "pedicel_or_stem_visible": "",
                "whole_head_unobstructed": "",
                "image_axis_correct_EXIF": "",
                "camera_tilt_unresolved": "",
                "anthesis_stage_qc": "",
                "original_S1_YOLO_head_detector_confidence": "",
                "original_S1_angle_measurement_status": "",
                "original_S1_image_vertical_0_to_180_angle_degrees": "",
                "reviewer_1": "",
                "reviewer_2": "",
                "QC_note": "",
            })
    if len(qc)!=N_READY*N_PER_SPECIES:
        raise ValueError("T3 sample size drift")
    root=Path(output);root.mkdir(parents=True,exist_ok=True)
    with (root/"T3_OTHER_TRIBE_23_SPECIES_115_BLIND_S1_ANGLE_QC_QUEUE.csv").open(
        "w",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=list(qc[0]))
        w.writeheader();w.writerows(qc)
    scheme={
        "schema":"azami_V4_T3_23_species_115_photo_blind_angle_QC_packet_v1",
        "status":"READY_FOR_INDEPENDENT_EXTERNAL_IMAGE_ANNOTATION_NOT_A_TRAIT_RESULT",
        "original_GenBank_ITS_source_taxa":N_FULL,
        "n_QC_eligible_independent_source_taxa":N_READY,
        "n_selected_original_GBIF_image_references":len(qc),
        "n_image_angle_measurements_in_packet":0,
        "phenotype_and_BIO18_BIO05_RSDS_elevation_hidden_when_selecting_images":True,
        "study_common_4axes":["chelsa_bio18","chelsa_bio05","chelsa_rsds_mean","elevation_m"],
        "original_S1_phenotype":"orientation_image_vertical_angle 0..180°; NO gravity/plant-nodding inference",
        "source_S1_implementation":"legacy/ch1_global/v2/52_measure_primary_traits_continuous.py",
        "original_S1_capitulum_model_sha256":"4078e0510532852681b65ee529cd82237b649ec99b17c4ca5f1da460a62d2bed",
        "minimum_head_width_px_original_S1":96,
        "original_S1_orientation_confidence_floor":0.65,
        "original_S1_flip_angle_tolerance_degrees":20,
        "next_step":"retrieve original weighted detector from long-term archive in analysis environment, run original YOLO ROI + pedicel angle calculation and manually QC blinded T3 images; do not claim domain transfer until validated against dual independent reviewers.",
        "manually_review_before_accepting_numeric_angle":True,
        "two_independent_reviewers_per_image_preferred":True,
        "image_downloads_or_reuploads_performed_here":0,
        "source_public_image_license_and_attribution_retained":True,
        "source_lat_lon_and_climate_values_NOT_in_annotation_packet":True,
        "T3_image_angle_BIO18_effect_test_status":"BLOCKED_AWAITING_VALIDATED_EXTERNAL_ANGLE_MEASUREMENTS",
        "source_131_taxa_fourclimate_JSON_SHA256":SHA_JSON,
        "image_selection_seed":SEED,
    }
    (root/"T3_BLIND_ORIENTATION_ANGLE_QC_REQUIREMENTS.json").write_text(
        json.dumps(scheme,indent=2,sort_keys=True,ensure_ascii=False)+"\n")
    print("T3_S1_ANGLE_BLIND_TRANSFER_QUEUE="+json.dumps({
        "eligible_source_taxa":N_READY,
        "photo_URLs_for_quality_review":len(qc),
        "already_measured_S1_equivalent_photo_angles":0,
        "no_flora_U_D_substitution":True},sort_keys=True))
    return scheme,qc

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--source-json",type=Path,required=True)
    p.add_argument("--source-photo-csv",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    prepare(a.source_json,a.source_photo_csv,a.out)
if __name__=="__main__":main()
