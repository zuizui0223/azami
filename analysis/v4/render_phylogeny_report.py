"""Summarise completed phylogenetic sensitivity without refitting."""
import argparse
import json
from pathlib import Path
import pandas as pd


def render(source, output):
    receipt = json.loads((source/"receipt.json").read_text())
    coverage = json.loads((source/"coverage.json").read_text())
    data = pd.read_csv(source/"coefficients.csv")
    output.mkdir(parents=True, exist_ok=True)
    rows = []
    for (endpoint, term), group in data.groupby(["endpoint", "term"]):
        base = group[group.model.eq("same_cohort_independent")]
        phy = group[group.model.eq("phylogenetic")]
        if len(base) != 1 or len(phy) != 52:
            raise ValueError(f"Incomplete result: {endpoint} {term}")
        increment = {"chelsa_bio01": 1, "chelsa_bio04": 100, "chelsa_bio12": 100,
                     "chelsa_bio15": 10, "chelsa_rsds_mean": 1}[term.rsplit("__", 1)[0]]
        rows.append(dict(endpoint=endpoint, term=term, increment=increment,
            n=int(base.iloc[0]["n"]), species=int(base.iloc[0].species),
            independent=base.iloc[0].coefficient*increment,
            independent_lower=base.iloc[0].lower95*increment,
            independent_upper=base.iloc[0].upper95*increment,
            phy_min=phy.coefficient.min()*increment, phy_max=phy.coefficient.max()*increment,
            lower_envelope=phy.lower95.min()*increment,
            upper_envelope=phy.upper95.max()*increment,
            positive_intervals=int(phy.lower95.gt(0).sum()),
            negative_intervals=int(phy.upper95.lt(0).sum()),
            zero_phy_variance=int(phy.phylogenetic_variance.eq(0).sum()),
            max_phy_variance=float(phy.phylogenetic_variance.max())))
    summary = pd.DataFrame(rows)
    summary.to_csv(output/"phylogenetic_summary.csv", index=False)
    text = [
        "# V4 系統配置感度検証",
        "",
        f"Actions run: {receipt['run_id']}; code: `{receipt['git_sha']}`.",
        "",
        "自生域の観測値を保持し、既存の5環境要因の種内・種間成分を同時に投入した。",
        "系統に沿う種効果と独立な種効果を分けた観測レベルのGaussian REMLを使用。",
        "比較対象は、同じ樹掲載種・同じ観測だけを使う独立種効果モデル。",
        "52本は既存の代替配置樹であり、独立した52回の実証ではない。",
        "区間は分散成分を推定値に固定したモデルベース近似95%区間であり、"
        "以前の種クラスター区間とは異なる。全樹の区間の外包は統合信頼区間ではない。",
        "個別成分の区間は多重比較補正前。色全体・輪郭全体の合同検定、"
        "系統的形質間共分散、適応・種分化の検証とは区別する。",
        "",
        "## 樹の網羅率", "",
        "|応答|元の観測数|樹掲載観測数|元の種数|樹掲載種数|",
        "|---|---:|---:|---:|---:|"]
    for key, item in coverage.items():
        text.append(f"|{key}|{item['full_n']}|{item['covered_n']}|"
                    f"{item['full_species']}|{item['covered_species']}|")
    text += ["", "## 向きと降水", ""]
    for _, r in summary[summary.endpoint.eq("orientation_image_vertical_angle") &
                         summary.term.str.startswith("chelsa_bio12")].iterrows():
        text.append(
            f"- {r.term}: 降水100 mm当たり、独立種モデル {r.independent:.4f} "
            f"[{r.independent_lower:.4f}, {r.independent_upper:.4f}] 度。"
            f"系統モデルの係数範囲 {r.phy_min:.4f}〜{r.phy_max:.4f} 度、"
            f"区間の外包 [{r.lower_envelope:.4f}, {r.upper_envelope:.4f}]。"
            f"正の区間 {r.positive_intervals}/52、負の区間 {r.negative_intervals}/52。")
    text += ["", "## 全成分", "",
             "係数は表記した環境増分当たりの元の形質単位。",
             "", "|応答・項|環境増分|系統係数範囲|区間の外包|正/負区間の樹数|",
             "|---|---:|---:|---:|---:|"]
    for _, r in summary.iterrows():
        text.append(f"|{r.endpoint} / {r.term}|{r.increment}|"
                    f"{r.phy_min:.5g}〜{r.phy_max:.5g}|"
                    f"[{r.lower_envelope:.5g}, {r.upper_envelope:.5g}]|"
                    f"{r.positive_intervals}/{r.negative_intervals}|")
    text += ["", f"失敗記録: {json.dumps(receipt['failures'], ensure_ascii=False)}", ""]
    (output/"phylogeny_report_ja.md").write_text("\n".join(text), encoding="utf-8")


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--source", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    render(**vars(p.parse_args()))
