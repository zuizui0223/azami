import importlib
import json
from pathlib import Path
import sys

import cv2
import numpy as np
from PIL import Image
import pytest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
perturb=importlib.import_module("analysis.v3.perturb_cached_heads")


def image():
    array=np.full((160,200,3),(20,150,20),np.uint8)
    cv2.ellipse(array,(100,85),(22,40),0,0,360,(60,90,50),-1)
    cv2.ellipse(array,(100,50),(20,15),0,0,360,(180,50,200),-1)
    return array


def conditions():
    return json.loads(perturb.CONTRACT.read_text())["conditions"]


def test_full_condition_grid():
    grid=conditions()
    assert len(grid)==len({c["id"] for c in grid})==14
    assert grid[0]=={"id":"baseline","kind":"identity"}


def test_identity_exact_pixels_and_crop_rounding():
    original=image()
    h,c,r=perturb.transform(original,[75,30,125,130],[],conditions()[0])
    assert r["head_box"]==(69,18,131,142)
    assert np.array_equal(h,original[18:142,69:131])
    assert r["transformed_bgr_pixel_sha256"]==perturb.pixel_hash(original)


def test_shift_is_relative_to_unpadded_box_and_moves_only_focal():
    _,_,r=perturb.transform(image(),[75,30,125,130],[[2,2,12,12]],{"id":"shift","kind":"bbox_shift","dx":.05,"dy":0.})
    assert r["head_box"]==(71,18,134,142)
    assert r["all_excluded_head_boxes"][1]==(0,0,14,14)


def test_half_resolution_coordinates_and_no_upsampling():
    _,_,r=perturb.transform(image(),[75,30,125,130],[],{"id":"half","kind":"resolution","scale":.5})
    assert r["actual_scale_x"]==r["actual_scale_y"]==.5
    assert r["transformed_image_width"]==100
    with pytest.raises(ValueError,match="upsample"):
        perturb.transform(image(),[75,30,125,130],[],{"id":"bad","kind":"resolution","scale":2})


def test_brightening_reports_clipping_without_modifying_source():
    original=np.full((10,10,3),250,np.uint8)
    h,_,r=perturb.transform(original,[2,2,8,8],[],{"id":"bright","kind":"intensity","multiplier":1.2})
    assert np.all(h==255) and np.all(original==250)
    assert r["out_of_range_channel_fraction_before_clipping"]==1


def test_boundary_clipping_is_explicit():
    _,_,r=perturb.transform(image(),[0,0,20,20],[],conditions()[1])
    assert r["head_clipped"] and r["context_clipped"]


def test_unrecognized_condition_rejected():
    with pytest.raises(ValueError,match="Unknown"):
        perturb.transform(image(),[75,30,125,130],[],{"id":"bad","kind":"automatic_repair"})


def test_dimension_aware_pixel_hash():
    assert perturb.pixel_hash(np.zeros((10,20),np.uint8))!=perturb.pixel_hash(np.zeros((20,10),np.uint8))


def task(tmp_path):
    original=image()
    path=tmp_path/"image.png"
    Image.fromarray(original[:,:,::-1]).save(path)
    # Cache pixel identity is width then height and RGB, unlike named BGR hashes.
    rgb=original[:,:,::-1].copy()
    sha=perturb.hashlib.sha256((200).to_bytes(8,"big")+(160).to_bytes(8,"big")+rgb.tobytes()).hexdigest()
    row=(perturb.digest(path),path.name,200,160,sha)
    box=[75,30,125,130]
    h,c,recipe=perturb.transform(original,box,[],conditions()[0])
    baseline=perturb.features.measure(h,c,recipe["context_box"],recipe["all_excluded_head_boxes"])[0]["endpoints"]
    return (str(tmp_path),row,"h1",box,[],baseline,conditions(),"measured")


def test_all_conditions_endpoints_and_replay(tmp_path):
    result=perturb.run_head(task(tmp_path))
    assert result["status"]=="completed"
    assert len(result["records"])==14
    assert all(len(r["result"]["endpoints"])==27 for r in result["records"])
    assert all(r["recipe"] for r in result["records"])
    json.dumps(result,allow_nan=False)


def test_baseline_mismatch_prevents_further_measurement(tmp_path):
    values=list(task(tmp_path))
    values[5][0]["value"]=999
    result=perturb.run_head(values)
    assert result["status"]=="baseline_replay_mismatch"
    assert all(r["status"]=="not_run_baseline_mismatch" for r in result["records"][1:])


def test_failed_source_retains_every_slot(tmp_path):
    values=list(task(tmp_path))
    values[-1]="invalid_roi"
    result=perturb.run_head(values)
    assert result["status"]=="source_not_evaluable"
    assert len(result["records"])*27==378


def test_changed_source_bytes_are_errors_not_qc(tmp_path):
    values=list(task(tmp_path))
    values[1]=("0"*64,*values[1][1:])
    result=perturb.run_head(values)
    assert result["status"]=="source_error"
    assert all("error" in r for r in result["records"])


def test_worker_lock_is_released_on_exit(tmp_path):
    with perturb.runner_lock(tmp_path):
        with pytest.raises((OSError,BlockingIOError)):
            with perturb.runner_lock(tmp_path):
                pass
    with perturb.runner_lock(tmp_path):
        pass
