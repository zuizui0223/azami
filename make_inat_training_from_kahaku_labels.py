#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
make_inat_training_from_kahaku_labels.py

科博「日本のアザミ」から取得した全記録CSVを、iNaturalist metadata.csv に種名で結合し、
機械学習用の画像単位ラベルCSVを作る。

入力:
  kahaku_azami_all_records_selenium.csv
  inat_cirsium_japan/metadata.csv

出力:
  cirsium_head_direction_training_from_kahaku.csv
  kahaku_inat_label_join_summary.csv

使い方:
  python make_inat_training_from_kahaku_labels.py

その後:
  python train_cirsium_head_direction_binary.py ^
    --csv cirsium_head_direction_training_from_kahaku.csv ^
    --path-col local_path ^
    --label-col head_direction
"""

import argparse
import re
from pathlib import Path

import pandas as pd


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--kahaku-csv", default="kahaku_azami_all_records_selenium.csv")
    p.add_argument("--metadata-csv", default="inat_cirsium_japan/metadata.csv")
    p.add_argument("--out-csv", default="cirsium_head_direction_training_from_kahaku.csv")
    p.add_argument("--summary-csv", default="kahaku_inat_label_join_summary.csv")
    p.add_argument("--metadata-taxon-col", default="taxon_name")
    p.add_argument("--metadata-path-col", default="local_path")
    p.add_argument("--min-confidence", default="medium", choices=["low", "medium", "high"])
    p.add_argument("--include-mixed", action="store_true")
    return p.parse_args()


def canonical_species(x):
    if pd.isna(x):
        return None
    s = str(x).strip().replace("_", " ")
    s = re.sub(r"\s+", " ", s)
    m = re.search(r"\b(Cirsium\s+[A-Za-z][A-Za-z-]+)\b", s)
    if not m:
        return None
    genus, epithet = m.group(1).split()[:2]
    return f"{genus[0].upper()}{genus[1:].lower()} {epithet.lower()}"


def normalize_orientation(x):
    if pd.isna(x):
        return None
    s = str(x).strip().lower()
    if s in {"upward", "up", "upright", "erect", "ascending"}:
        return "upward"
    if s in {"nodding", "downward", "down", "pendulous", "drooping"}:
        return "nodding"
    return None


def confidence_rank(x):
    if pd.isna(x):
        return 0
    s = str(x).strip().lower()
    return {"low": 1, "medium": 2, "high": 3}.get(s, 0)


def main():
    args = parse_args()

    kahaku_path = Path(args.kahaku_csv)
    meta_path = Path(args.metadata_csv)

    if not kahaku_path.exists():
        raise FileNotFoundError(f"科博CSVが見つかりません: {kahaku_path}")
    if not meta_path.exists():
        raise FileNotFoundError(f"metadata.csv が見つかりません: {meta_path}")

    kahaku = pd.read_csv(kahaku_path)
    meta = pd.read_csv(meta_path)

    print("Loaded Kahaku:", kahaku_path, len(kahaku), "rows")
    print("Kahaku columns:", list(kahaku.columns))
    print("Loaded metadata:", meta_path, len(meta), "rows")
    print("metadata columns:", list(meta.columns))

    # 科博側の種名列をなるべく自動検出
    species_candidates = [
        "species_extracted", "species", "種名", "scientific_name", "scientificName"
    ]
    species_col = None
    for c in species_candidates:
        if c in kahaku.columns:
            species_col = c
            break
    if species_col is None:
        raise ValueError("科博CSVに種名列が見つかりません。候補: " + ", ".join(species_candidates))

    # 向き列
    orient_candidates = [
        "head_orientation_binary", "manual_orientation_binary", "orientation", "head_direction"
    ]
    orient_col = None
    # manual列が埋まっているなら優先
    if "manual_orientation_binary" in kahaku.columns and kahaku["manual_orientation_binary"].notna().any():
        tmp = kahaku["manual_orientation_binary"].astype(str).str.strip()
        if (tmp != "").any():
            orient_col = "manual_orientation_binary"
    if orient_col is None:
        for c in orient_candidates:
            if c in kahaku.columns:
                orient_col = c
                break
    if orient_col is None:
        raise ValueError("科博CSVに頭花向き列が見つかりません。")

    # confidence列がなければ medium 扱い
    if "confidence" not in kahaku.columns:
        kahaku["confidence"] = "medium"

    # evidence列
    evidence_col = None
    for c in ["orientation_phrase", "evidence_phrase", "キャッチフレーズ", "ノート", "記載"]:
        if c in kahaku.columns:
            evidence_col = c
            break

    kahaku["_species_key"] = kahaku[species_col].apply(canonical_species)
    kahaku["_label"] = kahaku[orient_col].apply(normalize_orientation)
    kahaku["_conf_rank"] = kahaku["confidence"].apply(confidence_rank)

    min_rank = {"low": 1, "medium": 2, "high": 3}[args.min_confidence]

    use = kahaku.copy()
    if not args.include_mixed:
        use = use[use["_label"].isin(["upward", "nodding"])].copy()
    use = use[use["_conf_rank"] >= min_rank].copy()
    use = use[use["_species_key"].notna()].copy()

    # 同じ種でラベル衝突がある場合は除外
    conflict = (
        use.groupby("_species_key")["_label"]
        .nunique()
        .reset_index(name="n_label")
    )
    conflict_species = conflict.loc[conflict["n_label"] > 1, "_species_key"].tolist()
    if conflict_species:
        print("\nWARNING: conflicting labels found; excluded:")
        print(conflict_species)
        use = use[~use["_species_key"].isin(conflict_species)].copy()

    label_cols = ["_species_key", "_label", "confidence"]
    if evidence_col is not None:
        label_cols.append(evidence_col)
    if "japanese_name_extracted" in use.columns:
        label_cols.append("japanese_name_extracted")
    elif "japanese_name" in use.columns:
        label_cols.append("japanese_name")

    label_map = use[label_cols].drop_duplicates("_species_key").copy()

    meta["_species_key"] = meta[args.metadata_taxon_col].apply(canonical_species)

    joined = meta.merge(label_map, on="_species_key", how="inner")
    joined = joined.rename(columns={"_label": "head_direction", "_species_key": "species_matched"})

    # 画像ファイルが存在するものだけ
    joined["image_exists"] = joined[args.metadata_path_col].astype(str).apply(lambda p: Path(p).exists())
    joined = joined[joined["image_exists"]].copy()

    joined.to_csv(args.out_csv, index=False, encoding="utf-8-sig")

    summary = (
        joined.groupby(["species_matched", "head_direction"])
        .size()
        .reset_index(name="n_images")
        .sort_values(["head_direction", "species_matched"])
    )
    summary.to_csv(args.summary_csv, index=False, encoding="utf-8-sig")

    print("\nSaved:", args.out_csv)
    print("Rows:", len(joined))
    print("Label counts:")
    print(joined["head_direction"].value_counts(dropna=False))
    print("\nSpecies matched:", joined["species_matched"].nunique())
    print("Saved summary:", args.summary_csv)
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
