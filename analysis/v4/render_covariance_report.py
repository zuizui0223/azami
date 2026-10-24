"""Report all covariance components and conditional uncertainty."""
import argparse
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from analysis.v4.native_integration import FEATURES


def render(source, output):
    receipt = json.loads((source/"receipt.json").read_text())
    matrices = json.loads((source/"matrices.json").read_text())
    failures = json.loads((source/"failures.json").read_text())
    data = pd.read_csv(source/"covariance_pairs.csv")
    if failures or len(matrices)!=104 or len(data)!=52*2*3*28:
        raise ValueError("Incomplete covariance decomposition")
    output.mkdir(parents=True, exist_ok=True)
    rows = []
    for (version, level, left, right), group in data.groupby(["version", "level", "left", "right"]):
        if len(group)!=52:
            raise ValueError("Expected 52 placements per pair")
        rows.append(dict(version=version, level=level, left=left, right=right,
            rv_median=group.rv.median(), rv_min=group.rv.min(), rv_max=group.rv.max(),
            r_median=group.correlation.median(), r_min=group.correlation.min(),
            r_max=group.correlation.max(), interval_lower=group.correlation_lower.min(),
            interval_upper=group.correlation_upper.max(),
            finite_intervals=int(group.correlation_lower.notna().sum()),
            positive_intervals=int(group.correlation_lower.gt(0).sum()),
            negative_intervals=int(group.correlation_upper.lt(0).sum())))
    summary = pd.DataFrame(rows)
    summary.to_csv(output/"pair_summary.csv", index=False)
    fractions = []
    for key, item in matrices.items():
        tree, version = key.rsplit("/", 1)
        for k, level in enumerate(("phylogenetic", "independent_species", "within")):
            for j, endpoint in enumerate(receipt["columns"]):
                fractions.append(dict(tree=tree, version=version, level=level,
                    endpoint=endpoint, fraction=item["variance_fractions"][k][j],
                    lower=item["variance_fraction_intervals"][0][k][j],
                    upper=item["variance_fraction_intervals"][1][k][j]))
    f = pd.DataFrame(fractions)
    f.to_csv(output/"variance_fractions.csv", index=False)
    lines = ["# 頭花形質の系統・非系統種間・種内共分散", "",
        f"数値Actions: {receipt['run_id']}。コード: `{receipt['git_sha']}`。",
        f"元の共通集団 {receipt['original_common_n']:,}観測・{receipt['original_common_species']}種から、"
        f"樹に対応した {receipt['n']:,}観測・{receipt['species']}種を使用。",
        "色・向き・輪郭の9測定座標を同時に扱い、色相はsin/cosの2座標を組として比較。",
        "", "## 方法と解釈範囲", "",
        "- 3成分の共分散行列をPSD制約付き二次モーメント法で推定。REMLやMCMCではない。",
        "- 種平均と種内交差積を併用し、種ごとの観測数による平均値の標本誤差を補正。",
        "- 環境調整は既存の5環境要因の種内・種間係数を差し引く。種効果は差し引かない。",
        "- 52樹×調整前後×200反復のGaussianパラメトリックブートストラップ。"
        "環境係数・分布仮定に条件付きの点別区間で、境界では過小被覆の可能性がある。",
        "- 非系統種間成分は遺伝的成分、系統成分は適応・種分化の証明とは同義ではない。",
        "- 種内成分は観測数でプールした推定。以前の種均等重みの記述相関と直接同一視しない。",
        "- 相関が大きくても、その成分の分散が微小なら強い進化的結合と解釈しない。",
        "- 52本は代替配置で、独立した証拠ではない。範囲の外包は統合95%信頼区間ではない。",
        "", "## 環境調整後の分散割合", "",
        "|測定座標|系統成分の範囲|非系統種間の範囲|種内の範囲|",
        "|---|---:|---:|---:|"]
    for endpoint in receipt["columns"]:
        cells = []
        for level in ("phylogenetic", "independent_species", "within"):
            subset = f[f.version.eq("adjusted") & f.endpoint.eq(endpoint) & f.level.eq(level)]
            cells.append(f"{100*subset.fraction.min():.2f}–{100*subset.fraction.max():.2f}%")
        lines.append("|"+endpoint+"|"+"|".join(cells)+"|")
    lines += ["", "## 代表的な形質間相関（環境調整後）", "",
              "区間は樹ごとの区間の外包。正/負の列は区間がゼロをまたがない樹の数。",
              "", "|組み合わせ|成分|rの範囲|区間外包|正/負（52樹）|",
              "|---|---|---:|---:|---:|"]
    chosen = {("angle", "chroma"), ("angle", "aspect"), ("lightness", "chroma"),
              ("circularity", "solidity"), ("circularity", "width_profile")}
    for _, row in summary[summary.version.eq("adjusted")].iterrows():
        if (row.left, row.right) not in chosen:
            continue
        lines.append(f"|{row.left} × {row.right}|{row.level}|"
                     f"{row.r_min:.3f}–{row.r_max:.3f}|"
                     f"[{row.interval_lower:.3f}, {row.interval_upper:.3f}]|"
                     f"{row.positive_intervals}/{row.negative_intervals}|")
    conditions = [v["gram_condition"] for v in matrices.values()]
    bounded = sum(min(v["eigenvalues"][0])<1e-8 for v in matrices.values())
    lines += ["", "## 推定診断", "",
              f"104分解すべて完了。カーネルGram条件数: {min(conditions):.3g}–{max(conditions):.3g}。",
              f"系統共分散が少なくとも一方向で境界にある分解: {bounded}/104。",
              "完全な行列・制約前推定・固有値はmatrices.json、全28組はpair_summary.csv参照。",
              "境界推定の相関や点別区間だけから、特定の進化的モジュールを確定しない。", ""]
    (output/"covariance_report_ja.md").write_text("\n".join(lines), encoding="utf-8")
    names = list(FEATURES)
    fig, axes = plt.subplots(2, 3, figsize=(14, 9), constrained_layout=True)
    for row, version in enumerate(("raw", "adjusted")):
        for col, level in enumerate(("phylogenetic", "independent_species", "within")):
            ax = axes[row, col]
            values = np.full((8, 8), np.nan)
            for _, entry in summary[summary.version.eq(version)&summary.level.eq(level)].iterrows():
                i, j = names.index(entry.left), names.index(entry.right)
                values[i, j] = values[j, i] = entry.rv_median
            im = ax.imshow(values, vmin=0, vmax=1, cmap="viridis")
            ax.set_xticks(range(8), names, rotation=60, ha="right", fontsize=8)
            ax.set_yticks(range(8), names, fontsize=8)
            title = f"{version}: {level.replace('_', ' ')}"
            if level == "phylogenetic":
                title += "\nBoundary estimates: inspect uncertainty"
            ax.set_title(title, fontsize=10)
    fig.colorbar(im, ax=axes.ravel().tolist(), label="Median RV across placement trees")
    fig.suptitle("Conditional covariance components — PSD moment estimates\n"
                 "Large RV alone is not evidence: small component variances can produce extreme correlations.",
                 fontsize=12)
    fig.savefig(output/"covariance_components.png", dpi=180, bbox_inches="tight", pad_inches=.2)
    plt.close(fig)
    manifest = {p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                for p in output.iterdir() if p.is_file() and p.name!="manifest.json"}
    (output/"manifest.json").write_text(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--source", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    render(**vars(p.parse_args()))
