from pathlib import Path
import random
import re
import numpy as np
import pandas as pd
from PIL import Image

import tensorflow as tf
from tensorflow import keras
from tensorflow.keras import layers

from sklearn.model_selection import train_test_split
from sklearn.metrics import classification_report, confusion_matrix
from sklearn.utils.class_weight import compute_class_weight


# ============================================================
# paths
# ============================================================
BASE = Path(r"C:\Users\zuizui\cirsium_inat")

# YOLOで作ったcrop metadata
CROP_META = BASE / "inat_cirsium_japan_all_headcrops" / "head_crop_metadata.csv"

# 既存ラベルCSV候補
# 優先順に探す
LABEL_CANDIDATES = [
    BASE / "cirsium_head_direction_training_from_kahaku_old1379.csv",
    BASE / "cirsium_head_direction_training_from_kahaku_flowering_any_species.csv",
    BASE / "cirsium_head_direction_training_from_kahaku.csv",
    BASE / "cirsium_head_direction_image_training_from_species.csv",
    BASE / "cirsium_head_direction_binary_rechecked.csv",
]

# 今回は fixed 版として別フォルダに保存
OUTDIR = BASE / "cirsium_head_direction_model_yolo_crops_fixed"
OUTDIR.mkdir(parents=True, exist_ok=True)

TRAIN_TABLE_OUT = OUTDIR / "training_table_yolo_crops.csv"
MODEL_OUT = OUTDIR / "cirsium_head_direction_yolo_crop_model.keras"
LABELS_OUT = OUTDIR / "labels.csv"
REPORT_OUT = OUTDIR / "validation_report.txt"
VAL_PRED_OUT = OUTDIR / "validation_predictions.csv"


# ============================================================
# settings
# ============================================================
IMG_SIZE = 224
BATCH_SIZE = 32
EPOCHS_HEAD = 30
EPOCHS_FINETUNE = 15
SEED = 42

MIN_YOLO_CONF = 0.25

CLASS_NAMES = ["upward", "nodding"]


# ============================================================
# utilities
# ============================================================
def set_seed(seed=42):
    random.seed(seed)
    np.random.seed(seed)
    tf.random.set_seed(seed)


def find_existing(candidates, name):
    for p in candidates:
        if p.exists():
            print(f"[INFO] found {name}: {p}")
            return p

    raise FileNotFoundError(
        f"[ERROR] {name} not found. Checked:\n"
        + "\n".join(str(p) for p in candidates)
    )


def normalize_label(x):
    """
    head_direction / orientation_label などの表記ゆれを
    upward / nodding に正規化する。
    """
    s = str(x).strip().lower()

    if s in ["", "nan", "none", "na"]:
        return np.nan

    nodding_keys = [
        "nodding",
        "down",
        "downward",
        "pendant",
        "pendulous",
        "drooping",
        "点頭",
        "下向",
        "下向き",
        "垂",
        "垂れる",
        "うつむ",
        "俯",
    ]

    upward_keys = [
        "upward",
        "up",
        "upright",
        "erect",
        "erected",
        "ascending",
        "直立",
        "上向",
        "上向き",
        "斜上",
        "立つ",
        "立ち上",
    ]

    for k in nodding_keys:
        if k in s:
            return "nodding"

    for k in upward_keys:
        if k in s:
            return "upward"

    return np.nan


def guess_label_column(df):
    """
    ラベル列を推定。
    今回は head_direction を優先。
    """
    candidates = [
        "head_direction",
        "orientation_label",
        "orientation",
        "direction",
        "label",
        "binary_label",
        "class",
        "y",
    ]

    for c in candidates:
        if c in df.columns:
            vals = df[c].astype(str).apply(normalize_label)
            hit = vals.notna().sum()
            if hit > 0:
                return c

    best_col = None
    best_hit = 0

    for c in df.columns:
        vals = df[c].astype(str).apply(normalize_label)
        hit = vals.notna().sum()
        if hit > best_hit:
            best_hit = hit
            best_col = c

    if best_col is None or best_hit == 0:
        raise ValueError(f"Could not guess label column. columns={list(df.columns)}")

    return best_col


def make_join_keys_from_string(x):
    """
    local_path / source_name / image_url / photo_id などから
    対応に使えそうなキーを複数作る。
    """
    s = str(x)

    if s.lower() in ["nan", "none", ""]:
        return []

    s2 = s.split("?")[0].replace("\\", "/")
    stem = Path(s2).stem
    name = Path(s2).name

    keys = set()

    if stem and stem.lower() not in ["nan", "none"]:
        keys.add(stem)

    if name and name.lower() not in ["nan", "none"]:
        keys.add(name)
        keys.add(Path(name).stem)

    # iNaturalistのphoto_idっぽい数字も拾う
    nums = re.findall(r"\d{5,}", s)
    for n in nums:
        keys.add(n)

    return list(keys)


def make_crop_long_table(crop):
    """
    crop metadata 側からjoin keyを作る。
    """
    crop = crop.copy()

    crop["crop_key_list"] = crop["source_name"].apply(make_join_keys_from_string)

    if "source_image" in crop.columns:
        extra = crop["source_image"].apply(make_join_keys_from_string)
        crop["crop_key_list"] = [
            list(set(a) | set(b))
            for a, b in zip(crop["crop_key_list"], extra)
        ]

    if "crop_name" in crop.columns:
        extra = crop["crop_name"].apply(make_join_keys_from_string)
        crop["crop_key_list"] = [
            list(set(a) | set(b))
            for a, b in zip(crop["crop_key_list"], extra)
        ]

    crop_long = crop.explode("crop_key_list").rename(
        columns={"crop_key_list": "join_key"}
    )
    crop_long["join_key"] = crop_long["join_key"].astype(str)
    crop_long = crop_long[crop_long["join_key"].str.len() > 0].copy()

    return crop_long


def make_label_long_table(lab):
    """
    label CSV 側からjoin keyを作る。
    local_path, image_url, photo_id など複数列を使う。
    """
    key_columns = []

    for c in [
        "local_path",
        "source_image",
        "image_path",
        "path",
        "filename",
        "image_name",
        "source_name",
        "file",
        "photo_file",
        "photo_id",
        "image_url",
        "obs_id",
    ]:
        if c in lab.columns:
            key_columns.append(c)

    if len(key_columns) == 0:
        raise ValueError("No usable image key columns found in label CSV")

    print(f"[INFO] label key columns used: {key_columns}")

    lab_parts = []

    for c in key_columns:
        tmp = lab[[c, "label_norm"]].copy()
        tmp["join_key_list"] = tmp[c].apply(make_join_keys_from_string)
        tmp = tmp.explode("join_key_list").rename(
            columns={"join_key_list": "join_key"}
        )
        tmp["join_key"] = tmp["join_key"].astype(str)
        tmp = tmp[["join_key", "label_norm"]]
        lab_parts.append(tmp)

    lab_long = pd.concat(lab_parts, ignore_index=True)
    lab_long = lab_long.dropna(subset=["join_key"])
    lab_long = lab_long[lab_long["join_key"].astype(str).str.len() > 0].copy()

    # 同じjoin_keyに複数ラベルがある場合は最頻ラベル
    lab2 = (
        lab_long
        .groupby("join_key")["label_norm"]
        .agg(lambda x: x.value_counts().idxmax())
        .reset_index()
    )

    return lab_long, lab2


# ============================================================
# training table
# ============================================================
def load_and_join_training_table():
    if not CROP_META.exists():
        raise FileNotFoundError(f"crop metadata not found: {CROP_META}")

    crop = pd.read_csv(CROP_META)

    required = ["crop_path", "source_name", "yolo_conf"]
    for c in required:
        if c not in crop.columns:
            raise ValueError(f"crop metadata needs column: {c}")

    crop["crop_path"] = crop["crop_path"].astype(str)
    crop["crop_exists"] = crop["crop_path"].apply(lambda x: Path(x).exists())

    crop = crop[crop["crop_exists"]].copy()
    crop = crop[crop["yolo_conf"] >= MIN_YOLO_CONF].copy()

    print("\n=== crop metadata check ===")
    print(f"crop rows after exists/conf filter: {len(crop)}")
    cols_show = [c for c in ["source_name", "crop_name", "yolo_conf"] if c in crop.columns]
    print(crop[cols_show].head(5).to_string(index=False))

    crop_long = make_crop_long_table(crop)

    label_path = find_existing(LABEL_CANDIDATES, "label CSV")
    lab = pd.read_csv(label_path)

    print(f"\n[INFO] label CSV columns: {list(lab.columns)}")

    label_col = guess_label_column(lab)
    print(f"[INFO] guessed label column: {label_col}")

    lab["label_norm"] = lab[label_col].apply(normalize_label)

    print("\n=== raw label check ===")
    print(lab[label_col].value_counts(dropna=False).head(30))

    print("\n=== normalized label check ===")
    print(lab["label_norm"].value_counts(dropna=False))

    lab = lab.dropna(subset=["label_norm"]).copy()

    lab_long, lab2 = make_label_long_table(lab)

    train = crop_long.merge(lab2, on="join_key", how="inner")

    # 同じcropが複数キーで同じラベルにmatchした場合を整理
    train = train.drop_duplicates(subset=["crop_path", "label_norm"]).copy()

    # 同じcrop_pathに複数ラベルが付いた場合は除外
    conflict = train.groupby("crop_path")["label_norm"].nunique()
    conflict_paths = conflict[conflict > 1].index

    if len(conflict_paths) > 0:
        print(f"[WARN] conflicting labels found and removed: {len(conflict_paths)} crops")
        train = train[~train["crop_path"].isin(conflict_paths)].copy()

    train = train[train["label_norm"].isin(CLASS_NAMES)].copy()
    train["label_id"] = train["label_norm"].map(
        {"upward": 0, "nodding": 1}
    ).astype(int)

    print("\n=== join summary ===")
    print(f"crop rows        : {len(crop)}")
    print(f"crop join rows   : {len(crop_long)}")
    print(f"label usable rows: {len(lab)}")
    print(f"label join keys  : {len(lab2)}")
    print(f"matched crops    : {len(train)}")
    print(train["label_norm"].value_counts(dropna=False))

    if len(train) == 0:
        print("\n[DEBUG] crop source_name examples:")
        print(crop["source_name"].head(20).to_string(index=False))

        print("\n[DEBUG] crop join_key examples:")
        print(crop_long[["source_name", "join_key"]].head(30).to_string(index=False))

        print("\n[DEBUG] label join_key examples:")
        print(lab_long.head(30).to_string(index=False))

        raise RuntimeError(
            "No crop matched with label CSV. "
            "source_name / local_path / photo_id の対応がまだ合っていません。"
        )

    if train["label_id"].nunique() < 2:
        raise RuntimeError(
            "Only one class found after join. "
            "upward/nodding の両方が必要です。"
        )

    train.to_csv(TRAIN_TABLE_OUT, index=False, encoding="utf-8-sig")
    print(f"\n[INFO] saved training table: {TRAIN_TABLE_OUT}")

    return train


# ============================================================
# dataset
# ============================================================
def load_image_for_tf(path, label):
    """
    MobileNetV3Small(include_preprocessing=True) に合わせて
    0-255 float32 のまま渡す。
    ここで /255.0 しない。
    """
    img = tf.io.read_file(path)
    img = tf.image.decode_image(img, channels=3, expand_animations=False)
    img = tf.image.resize(img, [IMG_SIZE, IMG_SIZE])
    img = tf.cast(img, tf.float32)
    return img, label


def make_dataset(df, training=True):
    paths = df["crop_path"].astype(str).values
    labels = df["label_id"].astype(np.int32).values

    ds = tf.data.Dataset.from_tensor_slices((paths, labels))
    ds = ds.map(load_image_for_tf, num_parallel_calls=tf.data.AUTOTUNE)

    if training:
        ds = ds.shuffle(buffer_size=len(df), seed=SEED)

        aug = keras.Sequential([
            layers.RandomFlip("horizontal"),
            layers.RandomRotation(0.08),
            layers.RandomZoom(0.10),
            layers.RandomContrast(0.15),
        ])

        def augment(x, y):
            return aug(x, training=True), y

        ds = ds.map(augment, num_parallel_calls=tf.data.AUTOTUNE)

    ds = ds.batch(BATCH_SIZE).prefetch(tf.data.AUTOTUNE)
    return ds


# ============================================================
# model
# ============================================================
def build_model():
    """
    MobileNetV3Smallの内部preprocessingを使う。
    入力画像は0-255のfloat32。
    """
    base = keras.applications.MobileNetV3Small(
        input_shape=(IMG_SIZE, IMG_SIZE, 3),
        include_top=False,
        weights="imagenet",
        pooling="avg",
        include_preprocessing=True
    )

    base.trainable = False

    inputs = keras.Input(shape=(IMG_SIZE, IMG_SIZE, 3))
    x = base(inputs, training=False)
    x = layers.Dropout(0.30)(x)
    outputs = layers.Dense(1, activation="sigmoid")(x)

    model = keras.Model(inputs, outputs)

    model.compile(
        optimizer=keras.optimizers.Adam(1e-3),
        loss="binary_crossentropy",
        metrics=[
            "accuracy",
            keras.metrics.AUC(name="auc"),
            keras.metrics.Precision(name="precision"),
            keras.metrics.Recall(name="recall"),
        ]
    )

    return model, base


def load_images_numpy(paths):
    """
    validation prediction用。
    ここも /255.0 しない。
    """
    imgs = []

    for p in paths:
        im = Image.open(p).convert("RGB").resize((IMG_SIZE, IMG_SIZE))
        arr = np.asarray(im).astype("float32")
        imgs.append(arr)

    return np.stack(imgs, axis=0)


def make_class_weight(train_df):
    classes = np.array([0, 1])
    cw = compute_class_weight(
        class_weight="balanced",
        classes=classes,
        y=train_df["label_id"].values
    )
    class_weight = {0: float(cw[0]), 1: float(cw[1])}

    print("\n=== class weight ===")
    print(f"0 upward : {class_weight[0]:.4f}")
    print(f"1 nodding: {class_weight[1]:.4f}")

    return class_weight


# ============================================================
# main
# ============================================================
def main():
    set_seed(SEED)

    train_table = load_and_join_training_table()

    train_df, val_df = train_test_split(
        train_table,
        test_size=0.2,
        random_state=SEED,
        stratify=train_table["label_id"]
    )

    print("\n=== split ===")
    print("train:")
    print(train_df["label_norm"].value_counts())
    print("val:")
    print(val_df["label_norm"].value_counts())

    class_weight = make_class_weight(train_df)

    train_ds = make_dataset(train_df, training=True)
    val_ds = make_dataset(val_df, training=False)

    model, base = build_model()

    callbacks = [
        keras.callbacks.ModelCheckpoint(
            MODEL_OUT,
            monitor="val_auc",
            save_best_only=True,
            mode="max"
        ),
        keras.callbacks.EarlyStopping(
            monitor="val_auc",
            patience=8,
            restore_best_weights=True,
            mode="max"
        ),
        keras.callbacks.ReduceLROnPlateau(
            monitor="val_loss",
            factor=0.3,
            patience=3,
            min_lr=1e-6
        )
    ]

    print("\n=== phase 1: train classifier head ===")
    model.fit(
        train_ds,
        validation_data=val_ds,
        epochs=EPOCHS_HEAD,
        callbacks=callbacks,
        class_weight=class_weight
    )

    print("\n=== phase 2: fine tune last layers ===")
    base.trainable = True

    for layer in base.layers[:-30]:
        layer.trainable = False

    model.compile(
        optimizer=keras.optimizers.Adam(1e-5),
        loss="binary_crossentropy",
        metrics=[
            "accuracy",
            keras.metrics.AUC(name="auc"),
            keras.metrics.Precision(name="precision"),
            keras.metrics.Recall(name="recall"),
        ]
    )

    model.fit(
        train_ds,
        validation_data=val_ds,
        epochs=EPOCHS_FINETUNE,
        callbacks=callbacks,
        class_weight=class_weight
    )

    # best modelを読み直す
    model = keras.models.load_model(MODEL_OUT)

    val_paths = val_df["crop_path"].astype(str).tolist()
    y_true = val_df["label_id"].astype(int).values

    X_val = load_images_numpy(val_paths)
    prob = model.predict(X_val, batch_size=BATCH_SIZE).reshape(-1)
    y_pred = (prob >= 0.5).astype(int)

    cm = confusion_matrix(y_true, y_pred)
    rep = classification_report(
        y_true,
        y_pred,
        target_names=CLASS_NAMES,
        digits=3
    )

    print("\n=== Validation report threshold=0.5 ===")
    print(cm)
    print(rep)

    val_out = val_df.copy()
    val_out["pred_prob_nodding"] = prob
    val_out["pred_id_050"] = y_pred
    val_out["pred_label_050"] = np.where(y_pred == 1, "nodding", "upward")
    val_out["correct_050"] = val_out["label_id"] == val_out["pred_id_050"]
    val_out.to_csv(VAL_PRED_OUT, index=False, encoding="utf-8-sig")

    with open(REPORT_OUT, "w", encoding="utf-8") as f:
        f.write("=== Validation report threshold=0.5 ===\n")
        f.write(str(cm))
        f.write("\n\n")
        f.write(rep)
        f.write("\n\n")
        f.write("=== class weight ===\n")
        f.write(str(class_weight))
        f.write("\n\n")
        f.write(f"training_table: {TRAIN_TABLE_OUT}\n")
        f.write(f"validation_predictions: {VAL_PRED_OUT}\n")
        f.write(f"model: {MODEL_OUT}\n")

    pd.DataFrame({"label": CLASS_NAMES}).to_csv(
        LABELS_OUT,
        index=False,
        encoding="utf-8-sig"
    )

    print("\nSaved:")
    print(f"  {MODEL_OUT}")
    print(f"  {LABELS_OUT}")
    print(f"  {TRAIN_TABLE_OUT}")
    print(f"  {VAL_PRED_OUT}")
    print(f"  {REPORT_OUT}")


if __name__ == "__main__":
    main()
