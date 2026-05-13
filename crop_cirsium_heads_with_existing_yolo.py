from pathlib import Path
import pandas as pd
from PIL import Image
from ultralytics import YOLO

# =========================
# paths
# =========================
BASE = Path(r"C:\Users\zuizui\cirsium_inat")

# 今回はこのフォルダだけ使う
IMAGE_DIR = BASE / "inat_cirsium_japan"

# 昔のYOLO学習結果の候補
YOLO_CANDIDATES = [
    Path(r"C:\Users\zuizui\mc\runs\cirsium_detect\weights\best.pt"),
    Path(r"C:\Users\zuizui\mc\runs\cirsium_detect\best.pt"),
    Path(r"C:\Users\zuizui\mc\runs\cirsium_detect\train\weights\best.pt"),
    Path(r"C:\Users\zuizui\mc\runs\cirsium_detect\detect\train\weights\best.pt"),
]

# 出力先：前回の flowering_any_species と混ざらないように別名
OUTDIR = BASE / "inat_cirsium_japan_all_headcrops"
CROP_DIR = OUTDIR / "crops"
META_CSV = OUTDIR / "head_crop_metadata.csv"

CONF_THRES = 0.25
IMG_EXTS = {".jpg", ".jpeg", ".png", ".webp"}

# cropを少し広めに取る倍率
PAD_RATIO = 0.12


def find_existing_path(candidates, what):
    for p in candidates:
        if p.exists():
            print(f"[INFO] found {what}: {p}")
            return p

    raise FileNotFoundError(
        f"[ERROR] {what} not found. Checked:\n"
        + "\n".join(str(p) for p in candidates)
    )


def list_images(image_dir):
    imgs = []
    for ext in IMG_EXTS:
        imgs.extend(image_dir.rglob(f"*{ext}"))
        imgs.extend(image_dir.rglob(f"*{ext.upper()}"))
    return sorted(set(imgs))


def safe_crop_box(x1, y1, x2, y2, w, h, pad_ratio=0.12):
    bw = x2 - x1
    bh = y2 - y1

    pad_x = bw * pad_ratio
    pad_y = bh * pad_ratio

    x1 = max(0, int(x1 - pad_x))
    y1 = max(0, int(y1 - pad_y))
    x2 = min(w, int(x2 + pad_x))
    y2 = min(h, int(y2 + pad_y))

    return x1, y1, x2, y2


def main():
    if not IMAGE_DIR.exists():
        raise FileNotFoundError(f"[ERROR] image folder not found: {IMAGE_DIR}")

    OUTDIR.mkdir(parents=True, exist_ok=True)
    CROP_DIR.mkdir(parents=True, exist_ok=True)

    yolo_path = find_existing_path(YOLO_CANDIDATES, "YOLO model")

    images = list_images(IMAGE_DIR)

    print(f"[INFO] image_dir: {IMAGE_DIR}")
    print(f"[INFO] images found: {len(images)}")
    print(f"[INFO] output dir: {OUTDIR}")

    if len(images) == 0:
        raise RuntimeError(f"No images found in {IMAGE_DIR}")

    model = YOLO(str(yolo_path))

    rows = []
    crop_count = 0
    no_det_count = 0
    error_count = 0

    for i, img_path in enumerate(images, start=1):
        if i % 100 == 0:
            print(
                f"[INFO] processing {i}/{len(images)} "
                f"| crops={crop_count} | no_det={no_det_count}"
            )

        try:
            im = Image.open(img_path).convert("RGB")
            W, H = im.size

            results = model.predict(
                source=str(img_path),
                conf=CONF_THRES,
                verbose=False
            )

            if (
                len(results) == 0
                or results[0].boxes is None
                or len(results[0].boxes) == 0
            ):
                no_det_count += 1
                continue

            boxes = results[0].boxes

            for j, box in enumerate(boxes):
                xyxy = box.xyxy[0].detach().cpu().numpy()
                conf = (
                    float(box.conf[0].detach().cpu().numpy())
                    if box.conf is not None
                    else None
                )
                cls = (
                    int(box.cls[0].detach().cpu().numpy())
                    if box.cls is not None
                    else -1
                )

                x1, y1, x2, y2 = map(float, xyxy)
                cx1, cy1, cx2, cy2 = safe_crop_box(
                    x1, y1, x2, y2, W, H, PAD_RATIO
                )

                crop = im.crop((cx1, cy1, cx2, cy2))

                stem = img_path.stem
                crop_name = f"{stem}_det{j:02d}_conf{conf:.3f}.jpg"
                crop_path = CROP_DIR / crop_name

                crop.save(crop_path, quality=95)

                rows.append({
                    "source_image": str(img_path),
                    "source_name": img_path.name,
                    "crop_path": str(crop_path),
                    "crop_name": crop_name,
                    "det_index": j,
                    "yolo_conf": conf,
                    "yolo_class": cls,
                    "image_width": W,
                    "image_height": H,
                    "x1": x1,
                    "y1": y1,
                    "x2": x2,
                    "y2": y2,
                    "crop_x1": cx1,
                    "crop_y1": cy1,
                    "crop_x2": cx2,
                    "crop_y2": cy2,
                })

                crop_count += 1

        except Exception as e:
            error_count += 1
            print(f"[WARN] failed: {img_path} | {e}")

    meta = pd.DataFrame(rows)
    meta.to_csv(META_CSV, index=False, encoding="utf-8-sig")

    print("\n=== DONE ===")
    print(f"images total : {len(images)}")
    print(f"crops saved  : {crop_count}")
    print(f"no detection : {no_det_count}")
    print(f"errors       : {error_count}")
    print(f"crop dir     : {CROP_DIR}")
    print(f"metadata csv : {META_CSV}")


if __name__ == "__main__":
    main()