#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
predict_cirsium_head_direction_all.py

学習済みモデルを使って、metadata.csv など全画像に頭花向き予測を付与する。

基本:
  python predict_cirsium_head_direction_all.py

例:
  python predict_cirsium_head_direction_all.py --csv inat_cirsium_japan/metadata.csv --path-col local_path
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import tensorflow as tf
from tensorflow import keras


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--csv", default="inat_cirsium_japan/metadata.csv")
    p.add_argument("--model", default="cirsium_head_direction_model/cirsium_head_direction_binary.keras")
    p.add_argument("--labels", default="cirsium_head_direction_model/labels.csv")
    p.add_argument("--out-csv", default="inat_cirsium_japan/metadata_with_head_direction_predictions.csv")
    p.add_argument("--path-col", default="local_path")
    p.add_argument("--img-size", type=int, default=224)
    p.add_argument("--batch-size", type=int, default=32)
    p.add_argument("--review-threshold", type=float, default=0.65)
    return p.parse_args()


def load_image(path, img_size):
    img = tf.io.read_file(path)
    img = tf.image.decode_image(img, channels=3, expand_animations=False)
    img = tf.image.resize(img, (img_size, img_size))
    img = tf.cast(img, tf.float32) / 255.0
    return img


def main():
    args = parse_args()

    df = pd.read_csv(args.csv)
    if args.path_col not in df.columns:
        raise ValueError(f"path column not found: {args.path_col}. Columns: {list(df.columns)}")

    df[args.path_col] = df[args.path_col].astype(str)
    df = df[df[args.path_col].apply(lambda p: Path(p).exists())].copy()

    model = keras.models.load_model(args.model)
    labels = pd.read_csv(args.labels)
    id_to_label = dict(zip(labels["label_id"], labels["label"]))

    paths = df[args.path_col].astype(str).values

    ds = tf.data.Dataset.from_tensor_slices(paths)
    ds = ds.map(lambda p: load_image(p, args.img_size), num_parallel_calls=tf.data.AUTOTUNE)
    ds = ds.batch(args.batch_size).prefetch(tf.data.AUTOTUNE)

    pred = model.predict(ds)
    pred_id = np.argmax(pred, axis=1)
    pred_prob = np.max(pred, axis=1)

    df["head_direction_pred"] = [id_to_label[int(i)] for i in pred_id]
    df["head_direction_prob"] = pred_prob
    df["head_direction_final"] = df["head_direction_pred"]
    df.loc[df["head_direction_prob"] < args.review_threshold, "head_direction_final"] = "review"

    out = Path(args.out_csv)
    out.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out, index=False, encoding="utf-8-sig")

    print("Saved:", out)
    print(df["head_direction_final"].value_counts())


if __name__ == "__main__":
    main()
