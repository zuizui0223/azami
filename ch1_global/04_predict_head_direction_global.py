from pathlib import Path
import numpy as np
import pandas as pd
from PIL import Image, ImageDraw, ImageFont
import tensorflow as tf
from tensorflow import keras


# ============================================================
# paths
# ============================================================
BASE = Path(r"C:\Users\zuizui\cirsium_inat")

CROP_META = BASE / "inat_cirsium_japan_all_headcrops" / "head_crop_metadata.csv"
MODEL_PATH = BASE / "cirsium_head_direction_model_yolo_crops_fixed" / "cirsium_head_direction_yolo_crop_model.keras"

OUTDIR = BASE / "inat_cirsium_japan_all_head_direction_predictions_fixed"
OUTDIR.mkdir(parents=True, exist_ok=True)

CROP_PRED_OUT = OUTDIR / "all_crop_head_direction_predictions.csv"
IMAGE_SUMMARY_OUT = OUTDIR / "source_image_head_direction_summary.csv"
HIGH_CONF_OUT = OUTDIR / "high_confidence_crop_predictions.csv"
UNCERTAIN_OUT = OUTDIR / "uncertain_crop_predictions.csv"

# contact sheet outputs
SHEET_NODDING = OUTDIR / "contact_sheet_high_conf_nodding.jpg"
SHEET_UPWARD = OUTDIR / "contact_sheet_high_conf_upward.jpg"
SHEET_UNCERTAIN = OUTDIR / "contact_sheet_uncertain.jpg"
SHEET_FALSELIKE = OUTDIR / "contact_sheet_low_yolo_or_uncertain_check.jpg"


# ============================================================
# settings
# ============================================================
IMG_SIZE = 224
BATCH_SIZE = 32

# 解析用の安全閾値
NODDING_HIGH = 0.70
UPWARD_HIGH = 0.30

# YOLO cropそのものの信頼度
MIN_YOLO_CONF_FOR_MAIN = 0.25

# contact sheet
THUMB = 160
N_PER_SHEET = 100


# ============================================================
# image loading
# ============================================================
def load_images_numpy(paths):
    """
    MobileNetV3Small(include_preprocessing=True)で学習したので、
    ここでも /255.0 せず、0-255 float32で渡す。
    """
    imgs = []

    for p in paths:
        im = Image.open(p).convert("RGB").resize((IMG_SIZE, IMG_SIZE))
        arr = np.asarray(im).astype("float32")
        imgs.append(arr)

    return np.stack(imgs, axis=0)


def predict_in_batches(model, paths):
    probs = []

    for i in range(0, len(paths), BATCH_SIZE):
        batch_paths = paths[i:i + BATCH_SIZE]
        X = load_images_numpy(batch_paths)
        p = model.predict(X, verbose=0).reshape(-1)
        probs.append(p)

        if i == 0 or (i // BATCH_SIZE) % 10 == 0:
            print(f"[INFO] predicted {min(i + BATCH_SIZE, len(paths))}/{len(paths)}")

    return np.concatenate(probs)


# ============================================================
# classification
# ============================================================
def assign_class(prob):
    if prob >= NODDING_HIGH:
        return "nodding_high"
    elif prob <= UPWARD_HIGH:
        return "upward_high"
    else:
        return "uncertain"


def binary_label(prob):
    return "nodding" if prob >= 0.5 else "upward"


def confidence_from_prob(prob):
    return max(prob, 1.0 - prob)


# ============================================================
# source image aggregation
# ============================================================
def summarize_by_source(df):
    """
    1元画像から複数cropが出た場合の代表判定。
    方針：
    - 高信頼nodding cropが1つでもあれば、その画像は nodding_high
    - そうでなく、高信頼upward cropだけなら upward_high
    - 高信頼が混在したら mixed_high
    - 高信頼がなければ uncertain
    """
    rows = []

    for source_name, g in df.groupby("source_name", dropna=False):
        g_main = g[g["yolo_conf"] >= MIN_YOLO_CONF_FOR_MAIN].copy()

        if len(g_main) == 0:
            rows.append({
                "source_name": source_name,
                "n_crops": len(g),
                "n_main_crops": 0,
                "source_decision": "no_reliable_crop",
                "best_crop_path": np.nan,
                "best_prob_nodding": np.nan,
                "max_prob_nodding": np.nan,
                "min_prob_nodding": np.nan,
                "mean_prob_nodding": np.nan,
                "n_nodding_high": 0,
                "n_upward_high": 0,
                "n_uncertain": 0,
            })
            continue

        n_nod = int((g_main["decision"] == "nodding_high").sum())
        n_up = int((g_main["decision"] == "upward_high").sum())
        n_unc = int((g_main["decision"] == "uncertain").sum())

        # 一番確信度が高いcrop
        best_idx = g_main["confidence"].idxmax()
        best = g_main.loc[best_idx]

        if n_nod > 0 and n_up > 0:
            source_decision = "mixed_high"
        elif n_nod > 0:
            source_decision = "nodding_high"
        elif n_up > 0 and n_unc == 0:
            source_decision = "upward_high"
        elif n_up > 0 and n_unc > 0:
            # 高信頼upwardがあり、不確実もあるがnodding高信頼はない
            source_decision = "upward_high_with_uncertain"
        else:
            source_decision = "uncertain"

        rows.append({
            "source_name": source_name,
            "n_crops": len(g),
            "n_main_crops": len(g_main),
            "source_decision": source_decision,
            "best_crop_path": best["crop_path"],
            "best_prob_nodding": best["pred_prob_nodding"],
            "best_binary_label_050": best["pred_label_050"],
            "best_confidence_050": best["confidence"],
            "max_prob_nodding": g_main["pred_prob_nodding"].max(),
            "min_prob_nodding": g_main["pred_prob_nodding"].min(),
            "mean_prob_nodding": g_main["pred_prob_nodding"].mean(),
            "n_nodding_high": n_nod,
            "n_upward_high": n_up,
            "n_uncertain": n_unc,
        })

    return pd.DataFrame(rows)


# ============================================================
# contact sheets
# ============================================================
def get_font():
    # Windowsで使えればメイリオ、だめならデフォルト
    candidates = [
        r"C:\Windows\Fonts\meiryo.ttc",
        r"C:\Windows\Fonts\arial.ttf",
    ]
    for p in candidates:
        if Path(p).exists():
            try:
                return ImageFont.truetype(p, 14)
            except Exception:
                pass
    return ImageFont.load_default()


def make_contact_sheet(df, out_path, title, max_n=N_PER_SHEET):
    df = df.head(max_n).copy()

    if len(df) == 0:
        print(f"[WARN] no rows for contact sheet: {title}")
        return

    cols = 5
    rows = int(np.ceil(len(df) / cols))

    label_h = 54
    W = cols * THUMB
    H = rows * (THUMB + label_h)

    canvas = Image.new("RGB", (W, H), "white")
    draw = ImageDraw.Draw(canvas)
    font = get_font()

    for i, row in enumerate(df.itertuples(index=False)):
        x = (i % cols) * THUMB
        y = (i // cols) * (THUMB + label_h)

        try:
            img = Image.open(row.crop_path).convert("RGB")
            img.thumbnail((THUMB, THUMB))
            canvas.paste(img, (x, y))
        except Exception:
            continue

        prob = getattr(row, "pred_prob_nodding", np.nan)
        yconf = getattr(row, "yolo_conf", np.nan)
        dec = getattr(row, "decision", "")
        src = getattr(row, "source_name", "")

        text1 = f"{dec}"
        text2 = f"P(nod)={prob:.2f} YOLO={yconf:.2f}"
        text3 = str(src)[:24]

        draw.text((x + 3, y + THUMB + 2), text1, fill=(0, 0, 0), font=font)
        draw.text((x + 3, y + THUMB + 18), text2, fill=(0, 0, 0), font=font)
        draw.text((x + 3, y + THUMB + 34), text3, fill=(0, 0, 0), font=font)

    canvas.save(out_path, quality=95)
    print(f"[INFO] saved contact sheet: {out_path}")


# ============================================================
# main
# ============================================================
def main():
    if not CROP_META.exists():
        raise FileNotFoundError(f"crop metadata not found: {CROP_META}")
    if not MODEL_PATH.exists():
        raise FileNotFoundError(f"model not found: {MODEL_PATH}")

    print(f"[INFO] crop metadata: {CROP_META}")
    print(f"[INFO] model: {MODEL_PATH}")

    meta = pd.read_csv(CROP_META)

    if "crop_path" not in meta.columns:
        raise ValueError("metadata needs crop_path column")
    if "source_name" not in meta.columns:
        raise ValueError("metadata needs source_name column")
    if "yolo_conf" not in meta.columns:
        raise ValueError("metadata needs yolo_conf column")

    meta["crop_path"] = meta["crop_path"].astype(str)
    meta["crop_exists"] = meta["crop_path"].apply(lambda x: Path(x).exists())
    meta = meta[meta["crop_exists"]].copy()

    print(f"[INFO] existing crops: {len(meta)}")

    model = keras.models.load_model(MODEL_PATH)

    crop_paths = [Path(x) for x in meta["crop_path"].tolist()]
    prob = predict_in_batches(model, crop_paths)

    meta["pred_prob_nodding"] = prob
    meta["pred_label_050"] = np.where(meta["pred_prob_nodding"] >= 0.5, "nodding", "upward")
    meta["confidence"] = meta["pred_prob_nodding"].apply(confidence_from_prob)
    meta["decision"] = meta["pred_prob_nodding"].apply(assign_class)
    meta["is_main_yolo_conf"] = meta["yolo_conf"] >= MIN_YOLO_CONF_FOR_MAIN

    meta.to_csv(CROP_PRED_OUT, index=False, encoding="utf-8-sig")

    high = meta[meta["decision"].isin(["nodding_high", "upward_high"])].copy()
    uncertain = meta[meta["decision"] == "uncertain"].copy()

    high.to_csv(HIGH_CONF_OUT, index=False, encoding="utf-8-sig")
    uncertain.to_csv(UNCERTAIN_OUT, index=False, encoding="utf-8-sig")

    summary = summarize_by_source(meta)
    summary.to_csv(IMAGE_SUMMARY_OUT, index=False, encoding="utf-8-sig")

    print("\n=== crop prediction summary ===")
    print(meta["decision"].value_counts(dropna=False))
    print("\n=== binary threshold=0.5 summary ===")
    print(meta["pred_label_050"].value_counts(dropna=False))
    print("\n=== source image summary ===")
    print(summary["source_decision"].value_counts(dropna=False))

    print("\n=== probability summary ===")
    print(meta["pred_prob_nodding"].describe())

    # contact sheets
    nodding_df = (
        meta[meta["decision"] == "nodding_high"]
        .sort_values("pred_prob_nodding", ascending=False)
    )
    upward_df = (
        meta[meta["decision"] == "upward_high"]
        .sort_values("pred_prob_nodding", ascending=True)
    )
    uncertain_df = (
        meta[meta["decision"] == "uncertain"]
        .assign(dist_to_05=lambda d: (d["pred_prob_nodding"] - 0.5).abs())
        .sort_values("dist_to_05", ascending=True)
    )
    check_df = (
        meta[(meta["decision"] == "uncertain") | (meta["yolo_conf"] < 0.40)]
        .sort_values(["decision", "yolo_conf"], ascending=[True, True])
    )

    make_contact_sheet(nodding_df, SHEET_NODDING, "high confidence nodding")
    make_contact_sheet(upward_df, SHEET_UPWARD, "high confidence upward")
    make_contact_sheet(uncertain_df, SHEET_UNCERTAIN, "uncertain")
    make_contact_sheet(check_df, SHEET_FALSELIKE, "low yolo or uncertain check")

    print("\nSaved:")
    print(f"  {CROP_PRED_OUT}")
    print(f"  {IMAGE_SUMMARY_OUT}")
    print(f"  {HIGH_CONF_OUT}")
    print(f"  {UNCERTAIN_OUT}")
    print(f"  {SHEET_NODDING}")
    print(f"  {SHEET_UPWARD}")
    print(f"  {SHEET_UNCERTAIN}")
    print(f"  {SHEET_FALSELIKE}")


if __name__ == "__main__":
    main()
