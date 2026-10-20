"""Synthetic privacy and missingness gates for independently blinded native review."""
import hashlib
import io
import json
import zipfile

import pandas as pd
import pytest
from PIL import Image

from analysis.v4 import blinded_native_anatomy_review as q


def source(tmp_path):
    pic=io.BytesIO()
    Image.new("RGB",(30,40),(200,90,170)).save(pic,format="JPEG")
    rgb=pic.getvalue()
    ledger=[]
    for obs,photo in [(10,1000),(11,1001)]:
        ledger.append({
            "obs_id":obs,"photo_id":photo,"status":"downloaded_original",
            "original_path":f"originals/{photo}.jpg",
            "sha256":hashlib.sha256(rgb).hexdigest(),
            "license":"cc0","attribution":"fixture source",
            "width":30,"height":40,
        })
    original=tmp_path/"source.zip"
    with zipfile.ZipFile(original,"w") as z:
        z.writestr("acquisition_ledger.json",json.dumps(ledger))
        for _,photo in [(10,1000),(11,1001)]:
            z.writestr(f"originals/{photo}.jpg",rgb)
    return original


def test_native_membership_by_obs_id_not_assigned_photo(tmp_path):
    original=source(tmp_path)
    n=tmp_path/"native.csv"
    pd.DataFrame({
        "obs_id":[10,11],"photo_id":[991,998],
        "native_primary":[True,False],"accepted_key":[111,222],
    }).to_csv(n,index=False)
    out=tmp_path/"out"
    result=q.build(original,n,out,hashlib.sha256(original.read_bytes()).hexdigest())
    assert result["source_photos"]==2
    assert result["source_native_observations"]==1
    assert result["native_original_photos_recovered"]==1
    assert result["native_accepted_taxa"]==1
    assert result["independent_labels_completed"]==0
    for reviewer in "AB":
        with zipfile.ZipFile(out/f"blinded_review_{reviewer}.zip") as z:
            rows=json.loads(z.read("annotations_to_complete.json"))
            assert len(rows)==1
            assert rows[0]["anatomical_phyllaries"] is None
            assert rows[0]["review_status"]=="unreviewed"
            assert "accepted_key" not in json.dumps(rows)
            assert "original_observation_id" not in json.dumps(rows)
    a=json.loads(zipfile.ZipFile(out/"blinded_review_A.zip").read("annotations_to_complete.json"))
    b=json.loads(zipfile.ZipFile(out/"blinded_review_B.zip").read("annotations_to_complete.json"))
    assert a[0]["image_token"] != b[0]["image_token"]


def test_no_native_status_not_inferred_from_species_name(tmp_path):
    source=tmp_path/"input.csv"
    pd.DataFrame({"obs_id":[1],"accepted_key":[22]}).to_csv(source,index=False)
    with pytest.raises(ValueError,match="explicitly supplied"):
        q.native_roster(source)


def test_empty_annotation_is_not_phyllary_absence():
    template=q.blank_review("A-1","img.jpg")
    assert template["anatomical_phyllaries"] is None
    assert template["floral_polygons"] is None
    assert q.polygon_mask(None,(100,100)) is None
    assert q.polygon_mask([], (100,100)).sum()==0


def test_no_review_means_no_agreement_estimate(tmp_path):
    book=tmp_path/"codebook.json"
    book.write_text(json.dumps([{"token_A":"A-1","token_B":"B-1"}]))
    a=q.blank_review("A-1","a.jpg")
    b=q.blank_review("B-1","b.jpg")
    left=tmp_path/"A.json"
    right=tmp_path/"B.json"
    left.write_text(json.dumps([a]))
    right.write_text(json.dumps([b]))
    result=q.evaluate(book,left,right)
    assert result["both_raters_attempted"]==0
    assert result["median_polygon_iou"]["floral_polygons"] is None
    assert result["assessment"]=="PENDING_INDEPENDENT_REVIEW"
    assert result["independent_anatomical_validation_complete"] is False


def test_different_raters_geometry_agrees_when_median_iou_is_one(tmp_path):
    book=tmp_path/"book.json"
    book.write_text(json.dumps([{"token_A":"A-1","token_B":"B-1"}]))
    a=q.blank_review("A-1","a.jpg")
    b=q.blank_review("B-1","b.jpg")
    for obj,origin in ((a,[.2,.2]),(b,[.3,.2])):
        obj.update(review_status="completed",flowering_stage="anthesis",
                   floral_polygons=[[[0,0],[.5,0],[.5,.5],[0,.5]]],
                   head_base=origin,anatomical_phyllaries=[])
    left=tmp_path/"a.json"
    right=tmp_path/"b.json"
    left.write_text(json.dumps([a]))
    right.write_text(json.dumps([b]))
    result=q.evaluate(book,left,right)
    assert result["both_raters_attempted"]==1
    assert result["median_polygon_iou"]["floral_polygons"]==pytest.approx(1.)
    assert result["median_normalized_landmark_disagreement"]["head_base"]==pytest.approx(.1)
    assert result["median_disagreement_in_annotated_visible_phyllary_count"]==0


def test_unscorable_without_reason_and_bad_whorl_fail_closed():
    a=q.blank_review("A-1","a.jpg")
    a["review_status"]="unscorable"
    with pytest.raises(ValueError,match="reason"):
        q.validate_review([a],["A-1"])
    a["review_status"]="completed"
    a["anatomical_phyllaries"]=[{"base":[.1,.2],"tip":[.3,.4],"whorl":"spuriously-appressed"}]
    with pytest.raises(ValueError,match="whorl"):
        q.validate_review([a],["A-1"])


def test_source_zip_hash_fails_before_creating_blind_pack(tmp_path):
    original=source(tmp_path)
    roster=tmp_path/"native.csv"
    pd.DataFrame({"obs_id":[10],"native_primary":[True],"accepted_key":[10]}).to_csv(roster,index=False)
    with pytest.raises(ValueError,match="checksum"):
        q.build(original,roster,tmp_path/"out",source_sha256="0"*64)
    assert not (tmp_path/"out"/"blinded_review_A.zip").exists()
