#!/usr/bin/env python3
from __future__ import annotations

import importlib.util
from pathlib import Path

MODULE_PATH = Path(__file__).with_name("build_ch1_interpretive_figures.py")
spec = importlib.util.spec_from_file_location("ch1_interpretive_core", MODULE_PATH)
if spec is None or spec.loader is None:
    raise RuntimeError(f"Cannot load {MODULE_PATH}")
core = importlib.util.module_from_spec(spec)
spec.loader.exec_module(core)


def fixed_plot_figure6(primary, pgls, out):
    predictors = ["BIO1", "BIO4", "BIO12", "BIO15"]
    endpoint_order = ["Orientation", "Lightness", "Chroma", "Hue sine", "Hue cosine", "Aspect ratio", "Circularity", "Solidity", "Width-profile CV"]
    fig, axes = core.plt.subplots(1, 2, figsize=(13.2, 6.0), gridspec_kw={"width_ratios": [1.0, 1.15]})

    ax = axes[0]
    ax.set_xlim(-0.5, len(predictors)-0.5)
    ax.set_ylim(-0.5, len(endpoint_order)-0.5)
    ax.set_xticks(range(len(predictors)), predictors)
    ax.set_yticks(range(len(endpoint_order)), endpoint_order)
    ax.invert_yaxis()
    ax.grid(True, linewidth=0.4)
    ax.set_title("a  Supported climate axes differ across ecological scales", fontsize=9)
    ax.text(0.01, -0.12, "W = within-taxon primary; H = among-taxon historical sensitivity. Signs are coefficient directions.", transform=ax.transAxes, fontsize=7.5)

    endpoint_map = {"orientation_angle": "Orientation", "corolla_chroma": "Chroma", "hue_sin": "Hue sine", "hue_cos": "Hue cosine", "shape_aspect_ratio": "Aspect ratio"}
    for _, row in primary.iterrows():
        ep = endpoint_map.get(str(row["endpoint"]))
        pred = str(row["predictor"])
        if ep is None or pred not in predictors:
            continue
        x = predictors.index(pred); y = endpoint_order.index(ep)
        sign = "+" if float(row["standardized_beta"]) > 0 else "−"
        ax.text(x-0.15, y, f"W{sign}", ha="center", va="center", fontsize=8, fontweight="bold")

    for _, row in pgls.iterrows():
        ep = str(row["Endpoint"]).replace("Orientation angle", "Orientation")
        pred = str(row["Predictor"]).split()[0]
        if ep not in endpoint_order or pred not in predictors:
            continue
        beta = float(str(row["Median β"]).replace("−", "-").replace("+", ""))
        x = predictors.index(pred); y = endpoint_order.index(ep)
        sign = "+" if beta > 0 else "−"
        ax.text(x+0.15, y, f"H{sign}", ha="center", va="center", fontsize=8)

    ax = axes[1]
    work = pgls.copy()
    betas = core.pd.to_numeric(work["Median β"].astype(str).str.replace("−", "-", regex=False).str.replace("+", "", regex=False), errors="raise").to_numpy(float)
    ranges = [core.parse_range(v) for v in work["2.5–97.5% across trees"]]
    lo = core.np.array([x[0] for x in ranges]); hi = core.np.array([x[1] for x in ranges])
    labels = [f"{row['Endpoint']} · {str(row['Predictor']).replace(' annual ', ' ').replace(' mean ', ' ')}" for _, row in work.iterrows()]
    y = core.np.arange(len(work))[::-1]
    ax.errorbar(betas, y, xerr=core.np.vstack([betas-lo, hi-betas]), fmt="o", capsize=2.5, linewidth=1.0)
    ax.axvline(0, linewidth=0.8, linestyle="--")
    ax.set_yticks(y, labels)
    ax.set_xlabel("Median standardized coefficient across 50 randomized trees", fontsize=8)
    ax.set_title("b  Associations retained in all 50 alternative PGLS trees", fontsize=9)
    ax.tick_params(axis="both", labelsize=7)
    ax.text(0.01, -0.12, "Historical sensitivity only: 54/216 atlas taxa are direct dated-backbone tips.", transform=ax.transAxes, fontsize=7.5)
    fig.suptitle("Within-taxon and among-taxon environmental structure are not interchangeable", fontsize=13)
    fig.tight_layout()
    core.save(fig, out, "Figure_S1_9_historical_placement_sensitivity")


core.plot_figure6 = fixed_plot_figure6

if __name__ == "__main__":
    core.main()
