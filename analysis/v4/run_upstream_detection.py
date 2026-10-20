"""Replay the frozen visible-capitulum localizer on recovered ORIGINAL bytes."""
import argparse,json,os
from pathlib import Path
from PIL import ImageDraw
from analysis.v4.upstream import digest
from analysis.v4.run_upstream_acquisition import decode
WEIGHT_SHA='4078e0510532852681b65ee529cd82237b649ec99b17c4ca5f1da460a62d2bed'

def run(source,weights,out):
    if digest(weights)!=WEIGHT_SHA:raise ValueError('Frozen detector weight mismatch')
    import torch
    from ultralytics import YOLO, __version__
    torch.set_num_threads(2)
    model=YOLO(str(weights))
    if list(model.names.values())!=['visible_capitulum']:raise ValueError('Not the visible-head localizer')
    source=Path(source);out=Path(out);out.mkdir(parents=True,exist_ok=True)
    ledger=json.loads((source/'acquisition_ledger.json').read_text());rows=[]
    for row in ledger:
        if row['status']!='downloaded_original':
            rows.append({'photo_id':row['photo_id'],'status':'source_unavailable','detections':[]});continue
        path=source/row['original_path']
        if digest(path)!=row['sha256']:raise ValueError('Original photograph checksum changed')
        image,profile=decode(path)
        pred=model.predict(image,imgsz=960,conf=.25,iou=.7,device='cpu',verbose=False)[0]
        boxes=[]
        for box in pred.boxes:
            boxes.append({'xyxy':box.xyxy[0].tolist(),'confidence':float(box.conf[0])})
        boxes.sort(key=lambda b:(b['xyxy'][1],b['xyxy'][0]))
        rows.append({'photo_id':row['photo_id'],'source_sha256':row['sha256'],'status':'detected' if boxes else 'no_detection',
                     'width':image.width,'height':image.height,'detections':boxes,'botanically_validated':False})
        overlay=image.copy();d=ImageDraw.Draw(overlay)
        for i,b in enumerate(boxes):
            d.rectangle(b['xyxy'],outline='red',width=max(2,image.width//300));d.text(tuple(b['xyxy'][:2]),str(i+1),fill='yellow')
        overlay.thumbnail((800,800));overlay.save(out/f"{row['photo_id']}_detections.jpg")
        print(row['photo_id'],len(boxes),flush=True)
    receipt={'weight_sha256':WEIGHT_SHA,'ultralytics':__version__,'torch':torch.__version__,
             'execution_environment':('GitHub_Actions_CPU' if os.environ.get('GITHUB_ACTIONS')=='true' else 'local_CPU'),
             'confidence':.25,'imgsz':960,'iou':.7,
             'original_photos':sum(r['status']!='source_unavailable' for r in rows),'detected_photos':sum(r['status']=='detected' for r in rows),
             'heads_detected':sum(len(r['detections']) for r in rows),'v3_detections_reused':False,'anatomical_segmentation_completed':False,'rows':rows}
    (out/'detector_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    return receipt
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--source-dir',required=True);p.add_argument('--weights',required=True);p.add_argument('--out-dir',required=True)
    a=p.parse_args();r=run(a.source_dir,a.weights,a.out_dir);print(json.dumps({k:v for k,v in r.items() if k!='rows'},indent=2))
