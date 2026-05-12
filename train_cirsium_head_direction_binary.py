#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
train_cirsium_head_direction_binary.py

cirsium_head_direction_binary_rechecked.csv を使って、
アザミ頭花の向き binary classifier を学習する。

想定:
- CSV内に画像パス列がある: local_path / image_path / image_file / path / filename など
- CSV内にラベル列がある: label / head_direction / direction / orientation / orientation_label など
- ラベルは upward / up / erect / 0 / 1 / nodding / down / downward などを自動で解釈する

出力:
  cirsium_head_direction_model/
    cirsium_head_direction_binary.keras
    labels.csv
    training_history.csv
    validation_predictions.csv
"""

import argparse
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd
from PIL import Image

from sklearn.model_selection import train_test_split
from sklearn.metrics import classification_report, confusion_matrix

import tensorflow as tf
from tensorflow import keras
from tensorflow.keras import layers


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--csv", default="cirsium_head_direction_binary_rechecked.csv")
    p.add_argument("--out-dir", default="cirsium_head_direction_model")
    p.add_argument("--img-size", type=int, default=224)
    p.add_argument("--batch-size", type=int, default=16)
    p.add_argument("--epochs-head", type=int, default=15)
    p.add_argument("--epochs-finetune", type=int, default=8)
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--path-col", default="", help="画像パス列名を手動指定したい場合")
    p.add_argument("--label-col", default="", help="ラベル列名を手動指定したい場合")
    p.add_argument("--base-dir", default="", help="画像パスが相対パスの場合の基準フォルダ。空ならCSVの親フォルダ")
    return p.parse_args()


def find_column(df, candidates):
    lower_map = {c.lower(): c for c in df.columns}
    for cand in candidates:
        if cand.lower() in lower_map:
            return lower_map[cand.lower()]
    # partial match
    for c in df.columns:
        cl = c.lower()
        for cand in candidates:
            if cand.lower() in cl:
                return c
    return None


def normalize_label(x):
    if pd.isna(x):
        return None
    s = str(x).strip().lower()

    # remove spaces
    s2 = s.replace(" ", "_").replace("-", "_")

    upward_terms = {
        "up", "upward", "upright", "erect", "ascending", "上", "上向き",
        "直立", "斜上", "上向", "0_up", "upward_or_erect"
    }
    nodding_terms = {
        "down", "downward", "nodding", "pendulous", "drooping", "declined",
        "下", "下向き", "点頭", "垂下", "下垂", "下向", "1_down",
        "nodding_or_downward"
    }

    if s2 in upward_terms:
        return "upward"
    if s2 in nodding_terms:
        return "nodding"

    # numeric binary: 0/1 の意味が不明な場合もあるが、ここでは一般的に 0=upward, 1=nodding とする
    if s2 in {"0", "0.0"}:
        return "upward"
    if s2 in {"1", "1.0"}:
        return "nodding"

    # Japanese partial
    if ("上" in s) or ("直立" in s) or ("斜上" in s):
        return "upward"
    if ("下" in s) or ("点頭" in s) or ("垂" in s):
        return "nodding"

    return None


def resolve_path(p, base_dir):
    if pd.isna(p):
        return None
    path = Path(str(p).strip().strip('"'))
    if path.exists():
        return path
    if base_dir:
        cand = Path(base_dir) / path
        if cand.exists():
            return cand
    return None


def check_image(path):
    try:
        with Image.open(path) as im:
            im.verify()
        return True
    except Exception:
        return False


def make_dataset(df, img_col, label_col, img_size, batch_size, shuffle):
    paths = df[img_col].astype(str).values
    labels = df[label_col].astype(int).values

    def load_image(path, label):
        img = tf.io.read_file(path)
        img = tf.image.decode_image(img, channels=3, expand_animations=False)
        img = tf.image.resize(img, (img_size, img_size))
        img = tf.cast(img, tf.float32) / 255.0
        return img, label

    ds = tf.data.Dataset.from_tensor_slices((paths, labels))
    ds = ds.map(load_image, num_parallel_calls=tf.data.AUTOTUNE)
    if shuffle:
        ds = ds.shuffle(buffer_size=len(df), seed=42)
    ds = ds.batch(batch_size).prefetch(tf.data.AUTOTUNE)
    return ds


def main():
    args = parse_args()

    np.random.seed(args.seed)
    tf.random.set_seed(args.seed)

    csv_path = Path(args.csv)
    if not csv_path.exists():
        raise FileNotFoundError(f"CSVが見つかりません: {csv_path}")

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(csv_path)
    print("Loaded CSV:", csv_path)
    print("Rows:", len(df))
    print("Columns:", list(df.columns))

    path_col = args.path_col or find_column(
        df,
        ["local_path", "image_path", "image_file", "label_image_path", "path", "filepath", "file", "filename"]
    )
    label_col = args.label_col or find_column(
        df,
        ["head_direction", "direction", "orientation_label", "orientation", "label", "binary_label", "class"]
    )

    if path_col is None:
        raise ValueError("画像パス列が見つかりません。--path-col で指定してください。")
    if label_col is None:
        raise ValueError("ラベル列が見つかりません。--label-col で指定してください。")

    print("Using path_col:", path_col)
    print("Using label_col:", label_col)

    base_dir = args.base_dir or str(csv_path.parent)
    df["_resolved_path"] = df[path_col].apply(lambda x: resolve_path(x, base_dir))

    before = len(df)
    df = df[df["_resolved_path"].notna()].copy()
    print(f"Existing image paths: {len(df)} / {before}")

    df["_label_name"] = df[label_col].apply(normalize_label)
    print("Raw normalized label counts:")
    print(df["_label_name"].value_counts(dropna=False))

    df = df[df["_label_name"].isin(["upward", "nodding"])].copy()
    if len(df) < 50:
        raise ValueError("有効なラベル付き画像が少なすぎます。CSVの列名・ラベル表記を確認してください。")

    # Optional image validation
    print("Checking image files...")
    ok = []
    for p in df["_resolved_path"]:
        ok.append(check_image(p))
    df = df[np.array(ok)].copy()
    print("Valid images:", len(df))

    label_to_id = {"upward": 0, "nodding": 1}
    id_to_label = {0: "upward", 1: "nodding"}
    df["_label_id"] = df["_label_name"].map(label_to_id).astype(int)
    df["_resolved_path"] = df["_resolved_path"].astype(str)

    print("Final label counts:")
    print(df["_label_name"].value_counts())

    # save cleaned training table
    df.to_csv(out_dir / "training_table_cleaned.csv", index=False, encoding="utf-8-sig")

    train_df, val_df = train_test_split(
        df,
        test_size=0.2,
        random_state=args.seed,
        stratify=df["_label_id"]
    )

    train_ds = make_dataset(train_df, "_resolved_path", "_label_id", args.img_size, args.batch_size, shuffle=True)
    val_ds = make_dataset(val_df, "_resolved_path", "_label_id", args.img_size, args.batch_size, shuffle=False)

    data_aug = keras.Sequential([
        layers.RandomFlip("horizontal"),
        layers.RandomRotation(0.06),
        layers.RandomZoom(0.10),
        layers.RandomContrast(0.12),
    ])

    base = keras.applications.MobileNetV3Small(
        input_shape=(args.img_size, args.img_size, 3),
        include_top=False,
        weights="imagenet"
    )
    base.trainable = False

    inputs = keras.Input(shape=(args.img_size, args.img_size, 3))
    x = data_aug(inputs)
    x = keras.applications.mobilenet_v3.preprocess_input(x * 255.0)
    x = base(x, training=False)
    x = layers.GlobalAveragePooling2D()(x)
    x = layers.Dropout(0.35)(x)
    outputs = layers.Dense(2, activation="softmax")(x)

    model = keras.Model(inputs, outputs)

    callbacks = [
        keras.callbacks.ModelCheckpoint(
            filepath=str(out_dir / "best_head.keras"),
            monitor="val_accuracy",
            save_best_only=True,
            mode="max"
        ),
        keras.callbacks.EarlyStopping(
            monitor="val_loss",
            patience=5,
            restore_best_weights=True
        )
    ]

    model.compile(
        optimizer=keras.optimizers.Adam(1e-3),
        loss="sparse_categorical_crossentropy",
        metrics=["accuracy"]
    )

    print("\n=== Training classifier head ===")
    hist1 = model.fit(
        train_ds,
        validation_data=val_ds,
        epochs=args.epochs_head,
        callbacks=callbacks
    )

    print("\n=== Fine-tuning last layers ===")
    base.trainable = True
    for layer in base.layers[:-25]:
        layer.trainable = False

    model.compile(
        optimizer=keras.optimizers.Adam(1e-5),
        loss="sparse_categorical_crossentropy",
        metrics=["accuracy"]
    )

    hist2 = model.fit(
        train_ds,
        validation_data=val_ds,
        epochs=args.epochs_finetune,
        callbacks=callbacks
    )

    model_path = out_dir / "cirsium_head_direction_binary.keras"
    model.save(model_path)

    labels_df = pd.DataFrame({
        "label_id": [0, 1],
        "label": ["upward", "nodding"]
    })
    labels_df.to_csv(out_dir / "labels.csv", index=False, encoding="utf-8-sig")

    # Save history
    h1 = pd.DataFrame(hist1.history)
    h1["phase"] = "head"
    h1["epoch"] = np.arange(1, len(h1) + 1)
    h2 = pd.DataFrame(hist2.history)
    h2["phase"] = "finetune"
    h2["epoch"] = np.arange(1, len(h2) + 1)
    hist = pd.concat([h1, h2], ignore_index=True)
    hist.to_csv(out_dir / "training_history.csv", index=False, encoding="utf-8-sig")

    # Validation predictions
    val_paths = val_df["_resolved_path"].astype(str).values
    val_pred = model.predict(val_ds)
    val_pred_id = np.argmax(val_pred, axis=1)
    val_prob = np.max(val_pred, axis=1)

    out_val = val_df.copy()
    out_val["pred_id"] = val_pred_id
    out_val["pred_label"] = [id_to_label[int(i)] for i in val_pred_id]
    out_val["pred_prob"] = val_prob
    out_val.to_csv(out_dir / "validation_predictions.csv", index=False, encoding="utf-8-sig")

    print("\n=== Validation report ===")
    print(confusion_matrix(val_df["_label_id"], val_pred_id))
    print(classification_report(
        val_df["_label_id"],
        val_pred_id,
        target_names=["upward", "nodding"]
    ))

    print("\nSaved:")
    print(" ", model_path)
    print(" ", out_dir / "labels.csv")
    print(" ", out_dir / "training_table_cleaned.csv")
    print(" ", out_dir / "validation_predictions.csv")
    print(" ", out_dir / "training_history.csv")


if __name__ == "__main__":
    main()
