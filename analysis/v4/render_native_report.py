"""Render saved aggregate evidence without refitting or selecting hypotheses."""
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

ENV = {"chelsa_bio01": "Temperature (deg C)", "chelsa_bio04": "Thermal seasonality (BIO4)",
       "chelsa_bio12": "Annual rainfall (mm/year)", "chelsa_bio15": "Rainfall seasonality (CV pp)",
       "chelsa_rsds_mean": "Radiation (MJ/m2/day)"}
MODULE = {"orientation": "向き", "colour": "色全体", "outline": "頭花輪郭", "bract_si": "総苞片 SI"}
LABEL = {
    "orientation_image_vertical_angle": "Image-vertical angle (degrees)",
    "corolla_lab_lightness": "Floral lightness (L*)",
    "corolla_lab_chroma": "Floral chroma (C*)",
    "corolla_hue_sin": "Floral hue: sine", "corolla_hue_cos": "Floral hue: cosine",
    "capitulum_outline_aspect_ratio": "Head aspect ratio",
    "capitulum_outline_circularity": "Head circularity",
    "capitulum_outline_solidity": "Head solidity",
    "capitulum_width_profile_cv": "Head width-profile CV",
    "bract_projection_p95": "Bract contour projection p95",
    "bract_spread_fraction": "Bract contour spread fraction"}


def save(fig, output, name):
    fig.savefig(output/(name+".png"), dpi=180, bbox_inches="tight")
    fig.savefig(output/(name+".pdf"), bbox_inches="tight")
    plt.close(fig)


def read_sources(source):
    names = ["native_trait_models", "native_integration", "native_spatial", "native_integration_audit"]
    data = {name: json.loads((source/(name+".json")).read_text()) for name in names}
    digest = hashlib.sha256((source/"native_trait_models.json").read_bytes()).hexdigest()
    for name in ("native_integration", "native_spatial"):
        if data[name]["model_sha256"] != digest:
            raise ValueError("Integration/space results do not use these fitted models")
    if len({value["input_sha256"] for value in data.values()}) != 1:
        raise ValueError("Aggregate stages have different inputs")
    return data


def execute(source, output):
    data = read_sources(source)
    models, integration, space, audit = [data[k] for k in (
        "native_trait_models", "native_integration", "native_spatial", "native_integration_audit")]
    output.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False})
    predictors = models["selection"]["selected"]
    for module in models["models"].values():
        for endpoint, component in module["components"].items():
            fig, axes = plt.subplots(len(predictors), 2, figsize=(11, 2.25*len(predictors)))
            for row, variable in enumerate(predictors):
                for col, level in enumerate(("between", "within")):
                    ax = axes[row, col]
                    points = component["gradients"][variable+"__"+level]
                    ax.axhline(0, color="#999999", linewidth=.7)
                    ax.fill_between(points["x"], points["lower"], points["upper"],
                                    color="#4379a5", alpha=.2)
                    ax.plot(points["x"], points["change"], color="#235789")
                    ax.set_xlabel(ENV[variable]+(" deviation" if level == "within" else ""))
                    ax.set_ylabel("Trait change")
                    if row == 0:
                        ax.set_title("Between species" if level == "between" else "Within species")
            fig.suptitle(f'{LABEL[endpoint]} | n={component["n"]:,}; species={component["n_species"]}',
                         y=1.01, fontsize=13)
            fig.text(.5, -.015, "Conditional change from exposure median; pointwise 95% intervals.\n"
                     "Marginal 5–95% exposure range; not joint covariate overlap or a prediction interval.",
                     ha="center", fontsize=9)
            fig.tight_layout()
            save(fig, output, endpoint+"_gradients")
    core = integration["common_cohort"]
    names = list(FEATURES)
    fig, axes = plt.subplots(2, 2, figsize=(12, 11), constrained_layout=True)
    for row, version in enumerate(("raw", "adjusted")):
        for col, level in enumerate(("within", "between")):
            matrix = np.full((len(names), len(names)), np.nan)
            for pair in core["pairs"]:
                i, j = names.index(pair["left"]), names.index(pair["right"])
                matrix[i, j] = matrix[j, i] = pair["estimates"][version+"/"+level]["rv"]
            ax = axes[row, col]
            im = ax.imshow(matrix, vmin=0, vmax=1, cmap="viridis")
            ax.set_xticks(range(len(names)), names, rotation=50, ha="right")
            ax.set_yticks(range(len(names)), names)
            ax.set_title(version.capitalize()+" / "+level+" species")
    fig.colorbar(im, ax=axes, label="RV (unsigned association; not a signed effect)", shrink=.7)
    fig.suptitle(f'Common cohort: {core["n"]:,} observations / {core["species"]} species\n'
                 "Raw RV has sample-size bias; interpret with equal-n diagnostic", fontsize=13)
    save(fig, output, "integration_matrices")
    fig, axes = plt.subplots(1, 2, figsize=(10, 4), constrained_layout=True)
    rain = space["endpoints"]["orientation_image_vertical_angle"]
    factor = 100/models["selection"]["scales"]["chelsa_bio12"]["sd"]
    for ax, level in zip(axes, ("between", "within")):
        effects = [rain["main_effects"]["chelsa_bio12__"+level],
                   rain["balanced_fit"]["components"]["orientation_image_vertical_angle"]["effects"]["chelsa_bio12__"+level]]
        for y, effect in enumerate(effects):
            point = effect["raw_units_per_environment_sd"]*factor
            lo, hi = np.array(effect["raw_ci95"])*factor
            ax.errorbar(point, y, xerr=[[point-lo], [hi-point]], fmt="o", color="#235789", capsize=4)
        ax.axvline(0, color="#999999", linewidth=.8)
        ax.set_yticks([0, 1], ["All observations", "One per species x cell"])
        ax.set_title(level.capitalize()+" species")
        ax.set_xlabel("Image-angle change per 100 mm rain (degrees)")
        ax.set_ylim(-.6, 1.6)
    fig.suptitle("Rainfall–orientation: fixed spatial concentration comparison\n"
                 "Intervals cluster by species; thinning does not certify spatial independence", fontsize=11)
    save(fig, output, "rainfall_spatial_comparison")
    lines = [
        "# V4：自生域での頭花形質の組み合わせと環境対応",
        "",
        f'数値実行: https://github.com/zuizui0223/azami/actions/runs/{models["run_id"]}',
        f'コード: `{models["git_sha"]}`。入力SHA-256: `{models["input_sha256"]}`。',
        "",
        "V3の問いを維持し、自生域の観測に限定して、種内・種間の形質間関連と環境対応を評価した。",
        "これは改訂解析であり、既に見た結果から独立した事前登録検証ではない。",
        "",
        f'母集団は{audit["native_observations"]:,}観測／{audit["native_species"]}分類群。',
        f'共変動の共通標本は{core["n"]:,}観測／{core["species"]}種（各種5観測以上）。',
        "形質別の環境モデルは共通標本に制限せず、それぞれの利用可能観測を使用した。",
        "",
        "環境は温度・温度季節性・年間降水量・降水季節性・日射。"
        "VPDは温度との重複で除外。VPDの機能が否定されたわけではない。",
        "主解析の共同検定30件とSI10件をそれぞれHolm補正。成分の区間は点別95%。",
        "",
        "## 環境との対応（共同検定）",
        "",
        "|形質群|環境|種間 Holm P|種内 Holm P|",
        "|---|---|---:|---:|"]
    for module, model in models["models"].items():
        for variable in predictors:
            tests = model["tests"]
            lines.append(f'|{MODULE[module]}|{ENV[variable]}|'
                         f'{tests[variable+"__between"].get("holm_p", float("nan")):.4g}|'
                         f'{tests[variable+"__within"].get("holm_p", float("nan")):.4g}|')
    lines += ["", "## 向きと降水", "",
              "角度の増加は画像上でより下向き。重力基準の姿勢を直接検証したものではない。",
              "", "|尺度|全観測：100 mmあたり度［95%区間］|空間均等化：同単位［95%区間］|",
              "|---|---:|---:|"]
    for level, label in (("between", "種間"), ("within", "種内")):
        rendered = []
        for effect in (rain["main_effects"]["chelsa_bio12__"+level],
                       rain["balanced_fit"]["components"]["orientation_image_vertical_angle"]["effects"]["chelsa_bio12__"+level]):
            point = effect["raw_units_per_environment_sd"]*factor
            lo, hi = np.array(effect["raw_ci95"])*factor
            rendered.append(f"{point:+.3f} [{lo:+.3f}, {hi:+.3f}]")
        lines.append("|"+label+"|"+"|".join(rendered)+"|")
    lines += ["", f'向きの空間均等化で観測数は{rain["main_n"]:,}から{rain["balanced_n"]:,}へ変わる。',
              "区間は引き続き種クラスタに基づく。間引きによる標本変更は空間誤差モデルの代用ではない。",
              "", "## 種内・種間の共変動", "",
              "共通標本の全28関係を表示した。調整前後は同じ観測で比較し、種のランダム切片は差し引かない。",
              "調整後の区間は環境係数を固定した種ブートストラップであり、その係数の推定誤差は含まない。",
              "", "|条件|等標本数診断：種間−種内RVの関係間中央値（2.5%, 50%, 97.5%）|",
              "|---|---|"]
    for key, value in core["equal_n_diagnostic"].items():
        lines.append("|"+key+"|"+", ".join(f"{v:.5f}" for v in value["median_pair_delta_q025_q50_q975"])+"|")
    lines += ["", "この範囲は再標本化診断であり、信頼区間ではない。生のRVだけで種間が強いと結論しない。",
              "", "## 空間診断と間引きの判断", "",
              "|形質|全観測|均等化後|セル数|残差 Moran I|診断 P|",
              "|---|---:|---:|---:|---:|---:|"]
    for endpoint, value in space["endpoints"].items():
        diagnostic = value["residual_spatial_diagnostic"]
        lines.append(f'|{LABEL[endpoint]}|{value["main_n"]}|{value["balanced_n"]}|'
                     f'{diagnostic["cells"]}|{diagnostic.get("I", float("nan")):.4g}|'
                     f'{diagnostic.get("permutation_p", float("nan")):.4g}|')
    lines += ["", "自生域限定は対象集団の定義であり、空間依存を解消しない。"
              "この固定比較だけで間引き不要・空間独立とは断定しない。"
              "残差依存が残る場合、P値の空間的頑健性は未確認として扱う。",
              "", "## 解釈の限界", "",
              "写真の色はアントシアニン濃度ではなく、輪郭は三次元形態・防御機能の直接測定ではない。"
              "環境との条件付き対応や共変動は、適応・遺伝的統合・因果的な環境説明率を証明しない。"
              "種ランダム切片は系統補正ではない。撮影者IDは準備済み表にないため撮影者効果は未調整。",
              "", "全成分の実単位効果は `native_trait_effects.csv`、"
              "共変動の区間とペア別標本は `integration_pairs.csv` を参照。"]
    (output/"report_ja.md").write_text("\n".join(lines)+"\n", encoding="utf-8")
    (output/"native_trait_effects.csv").write_bytes((source/"native_trait_effects.csv").read_bytes())
    rows = []
    for scope, result in [("common", core)]+list(integration["pair_specific"].items()):
        for pair in result.get("pairs", []):
            for estimate, value in pair["estimates"].items():
                rows.append(dict(scope=scope, n=result["n"], species=result["species"],
                                 left=pair["left"], right=pair["right"], estimate=estimate,
                                 details=json.dumps(value)))
    pd.DataFrame(rows).to_csv(output/"integration_pairs.csv", index=False)
    manifest = {"source_run": models["run_id"], "source_commit": models["git_sha"],
                "sources": {p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                            for p in source.glob("*.json")},
                "outputs": {p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                            for p in output.iterdir() if p.is_file() and p.name != "manifest.json"}}
    (output/"manifest.json").write_text(json.dumps(manifest, indent=2)+"\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    execute(args.source, args.output)
