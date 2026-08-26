#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

plt.rcParams["svg.hashsalt"] = "azami-ch1-interpretive-figures"

TRAITS = [
    "orientation_angle_degrees_median",
    "corolla_lab_lightness_median",
    "corolla_lab_chroma_median",
    "corolla_hue_sin_median",
    "corolla_hue_cos_median",
    "shape_aspect_ratio_median",
    "shape_circularity_median",
    "shape_solidity_median",
    "shape_width_cv_median",
]

FRIENDLY = {
    "orientation_angle_degrees_median": "Orientation",
    "corolla_lab_lightness_median": "Lightness",
    "corolla_lab_chroma_median": "Chroma",
    "corolla_hue_sin_median": "Hue sine",
    "corolla_hue_cos_median": "Hue cosine",
    "shape_aspect_ratio_median": "Aspect ratio",
    "shape_circularity_median": "Circularity",
    "shape_solidity_median": "Solidity",
    "shape_width_cv_median": "Width-profile CV",
    "orientation_angle": "Orientation",
    "corolla_chroma": "Chroma",
    "hue_sin": "Hue sine",
    "hue_cos": "Hue cosine",
    "shape_aspect_ratio": "Aspect ratio",
    "involucre_projection_roughness": "Involucre roughness",
    "involucre_spread_fraction": "Involucre spread",
    "spine_relative_length_max_proxy": "Spine-like projection",
}

PRED_LABEL = {
    "BIO1": "BIO1 mean temperature",
    "BIO4": "BIO4 temperature seasonality",
    "BIO12": "BIO12 annual precipitation",
    "BIO15": "BIO15 precipitation seasonality",
    "soil_pH": "Soil pH",
    "chelsa_bio01": "BIO1",
    "chelsa_bio04": "BIO4",
    "chelsa_bio12": "BIO12",
    "chelsa_bio15": "BIO15",
    "soil_phh2o_0_30cm": "Soil pH",
}

STABLE_SPDE = [
    ("orientation_angle_degrees_median", "chelsa_bio01", "positive"),
    ("corolla_lab_lightness_median", "chelsa_bio01", "positive"),
    ("corolla_lab_lightness_median", "soil_phh2o_0_30cm", "negative"),
    ("corolla_lab_chroma_median", "chelsa_bio01", "negative"),
    ("corolla_lab_chroma_median", "chelsa_bio12", "positive"),
    ("corolla_hue_sin_median", "chelsa_bio04", "positive"),
    ("corolla_hue_cos_median", "chelsa_bio01", "negative"),
    ("corolla_hue_cos_median", "soil_phh2o_0_30cm", "positive"),
]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--environment", required=True)
    p.add_argument("--spde-fixed", required=True)
    p.add_argument("--primary-supported", required=True)
    p.add_argument("--bias-control", required=True)
    p.add_argument("--niche-summary", required=True)
    p.add_argument("--pgls-supported", required=True)
    p.add_argument("--out-dir", required=True)
    return p.parse_args()


def save(fig: plt.Figure, out: Path, stem: str) -> None:
    fig.savefig(out / f"{stem}.svg", bbox_inches="tight", metadata={"Date": None})
    fig.savefig(out / f"{stem}.png", dpi=300, bbox_inches="tight", metadata={"Software": "azami"})
    fig.savefig(out / f"{stem}.pdf", bbox_inches="tight", metadata={"Creator": "azami", "CreationDate": None, "ModDate": None})
    plt.close(fig)


def pca_products(env: pd.DataFrame, out: Path) -> tuple[pd.DataFrame, pd.DataFrame, np.ndarray]:
    counts = env.groupby("taxon_name").size()
    keep = counts[counts >= 5].index
    med = env.loc[env["taxon_name"].isin(keep)].groupby("taxon_name")[TRAITS].median(numeric_only=True)
    complete = med.dropna()
    if len(complete) != 148:
        raise ValueError(f"PCA complete taxa changed: expected 148, got {len(complete)}")
    z = StandardScaler().fit_transform(complete)
    pca = PCA().fit(z)
    expected = np.array([0.32933522, 0.23159018, 0.13242861])
    if not np.allclose(pca.explained_variance_ratio_[:3], expected, atol=2e-7):
        raise ValueError(f"PCA variance changed: {pca.explained_variance_ratio_[:3]}")
    scores = pd.DataFrame(pca.transform(z)[:, :3], index=complete.index, columns=["PC1", "PC2", "PC3"]).reset_index()
    load = pd.DataFrame(pca.components_[:3].T, index=TRAITS, columns=["PC1", "PC2", "PC3"]).reset_index(names="trait")
    scores.to_csv(out / "Figure4_taxon_pca_scores.csv", index=False)
    load.to_csv(out / "Figure4_trait_pca_loadings.csv", index=False)
    return scores, load, pca.explained_variance_ratio_


def plot_figure4(scores: pd.DataFrame, load: pd.DataFrame, variance: np.ndarray, out: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(12.3, 5.9), gridspec_kw={"width_ratios": [1.25, 0.9]})
    ax = axes[0]
    ax.scatter(scores["PC1"], scores["PC2"], s=18, alpha=0.42, linewidths=0)
    ax.axhline(0, linewidth=0.6, linestyle="--")
    ax.axvline(0, linewidth=0.6, linestyle="--")
    sx = max(abs(scores["PC1"].min()), abs(scores["PC1"].max()))
    sy = max(abs(scores["PC2"].min()), abs(scores["PC2"].max()))
    scale = min(sx, sy) * 0.72 / max(np.abs(load[["PC1", "PC2"]].to_numpy()).max(), 1e-9)
    short = {
        "orientation_angle_degrees_median": "Orientation",
        "corolla_lab_lightness_median": "Lightness",
        "corolla_lab_chroma_median": "Chroma",
        "corolla_hue_sin_median": "Hue sin",
        "corolla_hue_cos_median": "Hue cos",
        "shape_aspect_ratio_median": "Aspect ratio",
        "shape_circularity_median": "Circularity",
        "shape_solidity_median": "Solidity",
        "shape_width_cv_median": "Width CV",
    }
    for _, row in load.iterrows():
        x, y = row["PC1"] * scale, row["PC2"] * scale
        ax.annotate("", xy=(x, y), xytext=(0, 0), arrowprops=dict(arrowstyle="->", linewidth=1.0))
        ax.text(x * 1.05, y * 1.05, short[row["trait"]], fontsize=7.5, ha="center", va="center")
    ax.set_xlabel(f"PC1 ({variance[0]*100:.1f}%)")
    ax.set_ylabel(f"PC2 ({variance[1]*100:.1f}%)")
    ax.set_title("a  Taxon scores and trait directions\n148 taxa complete for all nine endpoints", fontsize=10)

    ax = axes[1]
    work = load.copy()
    work["label"] = work["trait"].map(FRIENDLY)
    work = work.iloc[::-1]
    y = np.arange(len(work))
    ax.barh(y, work["PC3"])
    ax.axvline(0, linewidth=0.7)
    ax.set_yticks(y, work["label"])
    ax.set_xlabel(f"PC3 loading ({variance[2]*100:.1f}% variance)")
    ax.set_title("b  Third independent trait dimension", fontsize=10)
    ax.text(0.02, 0.02, f"PC1–PC3 cumulative variance = {variance[:3].sum()*100:.1f}%", transform=ax.transAxes, fontsize=8)
    fig.suptitle("Taxon-level capitulum trait architecture is multidimensional", fontsize=13)
    fig.tight_layout()
    save(fig, out, "Figure_4_taxon_trait_architecture")


def forest(ax, df: pd.DataFrame, value: str, low: str, high: str, labels: list[str], title: str, xlabel: str) -> None:
    y = np.arange(len(df))[::-1]
    vals = pd.to_numeric(df[value], errors="coerce").to_numpy(float)
    lo = pd.to_numeric(df[low], errors="coerce").to_numpy(float)
    hi = pd.to_numeric(df[high], errors="coerce").to_numpy(float)
    ax.errorbar(vals, y, xerr=np.vstack([vals - lo, hi - vals]), fmt="o", capsize=2.5, linewidth=1.0)
    ax.axvline(0, linewidth=0.8, linestyle="--")
    ax.set_yticks(y, labels)
    ax.set_xlabel(xlabel, fontsize=8)
    ax.set_title(title, fontsize=9)
    ax.tick_params(axis="both", labelsize=7)


def plot_figure5(primary: pd.DataFrame, bias: pd.DataFrame, fx: pd.DataFrame, out: Path) -> pd.DataFrame:
    if len(primary) != 8 or len(bias) != 4:
        raise ValueError(f"Expected primary=8 and bias-control=4 rows; got {len(primary)} and {len(bias)}")
    required_native = {
        "frozen_beta", "native_only_beta", "native_only_ci_low", "native_only_ci_high",
        "native_only_q", "native_range_robust",
    }
    missing_native = required_native.difference(bias.columns)
    if missing_native:
        raise ValueError(f"Native-range audit columns missing: {sorted(missing_native)}")
    retained = bias["native_range_robust"].astype(str).str.lower().eq("true")
    if int(retained.sum()) != 2:
        raise ValueError(f"Expected 2/4 native-range-retained rows; got {int(retained.sum())}")
    if not (
        pd.to_numeric(bias["frozen_beta"], errors="raise")
        * pd.to_numeric(bias["native_only_beta"], errors="raise")
        > 0
    ).all():
        raise ValueError("A native-only coefficient changed sign relative to the frozen primary row")
    fig, axes = plt.subplots(1, 3, figsize=(14.7, 6.5), gridspec_kw={"width_ratios": [1.05, 1.35, 1.05]})

    p = primary.copy()
    labels = [f"{FRIENDLY.get(r.endpoint, r.endpoint)} · {r.predictor}" for r in p.itertuples()]
    forest(axes[0], p, "standardized_beta", "ci95_low", "ci95_high", labels,
           "a  Exhaustive within-taxon primary\n8/36 BH-supported components", "Standardized β (95% CI)")
    axes[0].text(0.01, -0.11, "Hue sine/cosine are circular components and require joint interpretation.", transform=axes[0].transAxes, fontsize=7)

    stable_rows = []
    group_order = ["climate", "climate_topography", "climate_soil", "full"]
    marker = {"climate": "o", "climate_topography": "s", "climate_soil": "^", "full": "D"}
    for trait, term, direction in STABLE_SPDE:
        sub = fx.loc[(fx["trait"] == trait) & (fx["term"] == term)].copy()
        if sub.empty:
            raise ValueError(f"Missing SPDE stable pair: {trait} {term}")
        means = pd.to_numeric(sub["mean"], errors="coerce")
        sign = "positive" if (means > 0).all() else "negative" if (means < 0).all() else "mixed"
        if sign != direction:
            raise ValueError(f"SPDE sign changed for {trait} {term}: {sign}")
        stable_rows.append({"trait": trait, "term": term, "direction": direction, "groups": len(sub)})
    rows = pd.DataFrame(stable_rows)
    y_base = np.arange(len(rows))[::-1]
    offsets = dict(zip(group_order, np.linspace(-0.22, 0.22, len(group_order))))
    for g in group_order:
        for i, r in rows.iterrows():
            sub = fx.loc[(fx["trait"] == r.trait) & (fx["term"] == r.term) & (fx["model_group"] == g)]
            if sub.empty:
                continue
            row = sub.iloc[0]
            m = float(row["mean"]); lo = float(row["0.025quant"]); hi = float(row["0.975quant"])
            supported = float(row["q_bh_posterior_tail_global"]) < 0.05
            axes[1].errorbar(m, y_base[i] + offsets[g], xerr=[[m-lo], [hi-m]], fmt=marker[g], markersize=4.5,
                             fillstyle="full" if supported else "none", capsize=1.8, linewidth=0.8,
                             label=g.replace("_", " + ") if i == 0 else None)
    axes[1].axvline(0, linewidth=0.8, linestyle="--")
    spde_labels = [f"{FRIENDLY[t]} · {PRED_LABEL[term]}" for t, term, _ in STABLE_SPDE]
    axes[1].set_yticks(y_base, spde_labels)
    axes[1].set_xlabel("SPDE posterior coefficient (95% credible interval)", fontsize=8)
    axes[1].set_title("b  Grouped SPDE-INLA stable patterns\nfilled marker = global BH q < 0.05 in that model group", fontsize=9)
    axes[1].tick_params(axis="both", labelsize=7)
    axes[1].legend(fontsize=6.5, loc="lower right")

    b = bias.copy()
    blabels = [
        f"{FRIENDLY.get(r.endpoint, r.endpoint)} · {PRED_LABEL.get(r.predictor, r.predictor)}"
        for r in b.itertuples()
    ]
    y_positions = np.arange(len(b))[::-1]
    vals = pd.to_numeric(b["native_only_beta"], errors="raise").to_numpy(float)
    lo = pd.to_numeric(b["native_only_ci_low"], errors="raise").to_numpy(float)
    hi = pd.to_numeric(b["native_only_ci_high"], errors="raise").to_numpy(float)
    kept = b["native_range_robust"].astype(str).str.lower().eq("true").to_numpy()
    for value, lower, upper, y, is_kept in zip(vals, lo, hi, y_positions, kept):
        axes[2].errorbar(
            value, y, xerr=[[value - lower], [upper - value]], fmt="o",
            fillstyle="full" if is_kept else "none", markersize=5.5,
            capsize=2.5, linewidth=1.0,
        )
    axes[2].axvline(0, linewidth=0.8, linestyle="--")
    axes[2].set_yticks(y_positions, blabels)
    axes[2].set_xlabel("Native-only standardized β (95% CI)", fontsize=8)
    axes[2].set_title(
        "c  Native-range restriction after timing/taxon audits\n2/4 non-circular rows meet the locked rule",
        fontsize=9,
    )
    axes[2].tick_params(axis="both", labelsize=7)
    axes[2].text(
        0.01, -0.15,
        "Filled: same sign and native-only BH q < 0.05. Open: same sign but BH failed. Stage and colour gates remain open.",
        transform=axes[2].transAxes, fontsize=6.8,
    )

    fig.suptitle("Effect sizes make the trait-specific environmental structure explicit", fontsize=13)
    fig.tight_layout()
    save(fig, out, "Figure_5_environmental_effect_sizes")
    rows.to_csv(out / "Figure5_SPDE_stable_pairs.csv", index=False)
    return rows


def parse_range(text: str) -> tuple[float, float]:
    clean = str(text).replace("−", "-").replace("+", "")
    left, right = [v.strip() for v in clean.split(" to ")]
    return float(left), float(right)


def plot_figure6(primary: pd.DataFrame, pgls: pd.DataFrame, out: Path) -> None:
    predictors = ["BIO1", "BIO4", "BIO12", "BIO15"]
    endpoint_order = ["Orientation", "Lightness", "Chroma", "Hue sine", "Hue cosine", "Aspect ratio", "Circularity", "Solidity", "Width-profile CV"]
    fig, axes = plt.subplots(1, 2, figsize=(13.2, 6.0), gridspec_kw={"width_ratios": [1.0, 1.15]})

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
    for r in primary.itertuples():
        ep = endpoint_map.get(r.endpoint)
        if ep is None or r.predictor not in predictors:
            continue
        x = predictors.index(r.predictor); y = endpoint_order.index(ep)
        sign = "+" if float(r.standardized_beta) > 0 else "−"
        ax.text(x-0.15, y, f"W{sign}", ha="center", va="center", fontsize=8, fontweight="bold")

    for _, r in pgls.iterrows():
        ep = str(r["Endpoint"]).replace("Orientation angle", "Orientation")
        pred = str(r["Predictor"]).split()[0]
        if ep not in endpoint_order or pred not in predictors:
            continue
        beta = float(str(r["Median β"]).replace("−", "-").replace("+", ""))
        x = predictors.index(pred); y = endpoint_order.index(ep)
        sign = "+" if beta > 0 else "−"
        ax.text(x+0.15, y, f"H{sign}", ha="center", va="center", fontsize=8)

    ax = axes[1]
    work = pgls.copy()
    betas = pd.to_numeric(work["Median β"].astype(str).str.replace("−", "-", regex=False).str.replace("+", "", regex=False), errors="raise").to_numpy(float)
    ranges = [parse_range(v) for v in work["2.5–97.5% across trees"]]
    lo = np.array([x[0] for x in ranges]); hi = np.array([x[1] for x in ranges])
    labels = [f"{r.Endpoint} · {r.Predictor.replace(' annual ', ' ').replace(' mean ', ' ')}" for r in work.itertuples()]
    y = np.arange(len(work))[::-1]
    ax.errorbar(betas, y, xerr=np.vstack([betas-lo, hi-betas]), fmt="o", capsize=2.5, linewidth=1.0)
    ax.axvline(0, linewidth=0.8, linestyle="--")
    ax.set_yticks(y, labels)
    ax.set_xlabel("Median standardized coefficient across 50 randomized trees", fontsize=8)
    ax.set_title("b  Associations retained in all 50 alternative PGLS trees", fontsize=9)
    ax.tick_params(axis="both", labelsize=7)
    ax.text(0.01, -0.12, "Historical sensitivity only: 54/216 atlas taxa are direct dated-backbone tips.", transform=ax.transAxes, fontsize=7.5)
    fig.suptitle("Within-taxon and among-taxon environmental structure are not interchangeable", fontsize=13)
    fig.tight_layout()
    save(fig, out, "Figure_S1_9_historical_placement_sensitivity")


def safe_cov(x: np.ndarray) -> np.ndarray:
    c = np.cov(x, rowvar=False)
    if np.ndim(c) == 0:
        c = np.array([[float(c)]])
    return c + np.eye(c.shape[0]) * 1e-8


def overlap(a: np.ndarray, b: np.ndarray) -> float:
    m1, m2 = a.mean(0), b.mean(0)
    c1, c2 = safe_cov(a), safe_cov(b)
    c = (c1+c2)/2
    diff = (m1-m2)[:, None]
    term1 = float(np.asarray((diff.T @ np.linalg.pinv(c) @ diff)/8).item())
    _, ld = np.linalg.slogdet(c); _, ld1 = np.linalg.slogdet(c1); _, ld2 = np.linalg.slogdet(c2)
    db = term1 + 0.5*(float(ld)-0.5*(float(ld1)+float(ld2)))
    return float(np.exp(-max(db, 0.0)))


def niche_metrics(env: pd.DataFrame) -> pd.DataFrame:
    env_cols = [c for c in [
        "chelsa_bio01", "chelsa_bio04", "chelsa_bio12", "chelsa_bio15", "topo_elevation", "topo_slope", "topo_roughness",
        "soil_bdod_0_30cm", "soil_cec_0_30cm", "soil_cfvo_0_30cm", "soil_clay_0_30cm", "soil_sand_0_30cm", "soil_silt_0_30cm",
        "soil_nitrogen_0_30cm", "soil_phh2o_0_30cm", "soil_soc_0_30cm", "soil_ocd_0_30cm"] if c in env.columns]
    counts = env.groupby("taxon_name").size(); keep = counts[counts >= 5].index
    sp = env.loc[env["taxon_name"].isin(keep)].groupby("taxon_name")[TRAITS+env_cols].median(numeric_only=True).dropna()
    if len(sp) != 148:
        raise ValueError(f"Niche complete taxa changed: {len(sp)}")
    epca = PCA(n_components=3).fit_transform(StandardScaler().fit_transform(sp[env_cols]))
    rows=[]
    for trait in TRAITS:
        vals=sp[trait].to_numpy(float); q1,q3=np.quantile(vals,[0.25,0.75])
        low=epca[vals<=q1]; high=epca[vals>=q3]
        rows.append({"trait":trait,"centroid_distance":float(np.linalg.norm(low.mean(0)-high.mean(0))),"overlap":overlap(low,high)})
    return pd.DataFrame(rows)


def plot_s8(env: pd.DataFrame, niche_summary: pd.DataFrame, out: Path) -> None:
    metrics=niche_metrics(env)
    support={}
    direct={}
    for r in niche_summary.itertuples(index=False):
        ep=str(r[0])
        support[ep]=(float(r[1])<0.05) or (float(r[2])<0.05)
        direct[ep]=(float(r[3])<0.05) or (float(r[4])<0.05)
    reverse={"orientation_angle_degrees_median":"Orientation","corolla_lab_lightness_median":"Lightness","corolla_lab_chroma_median":"Chroma","corolla_hue_sin_median":"Hue sine","corolla_hue_cos_median":"Hue cosine","shape_aspect_ratio_median":"Aspect ratio","shape_circularity_median":"Circularity","shape_solidity_median":"Solidity","shape_width_cv_median":"Width-profile CV"}
    fig,ax=plt.subplots(figsize=(8.4,6.2))
    for _,r in metrics.iterrows():
        label=reverse[r.trait]; sup=support.get(label,False); ds=direct.get(label,False)
        ax.scatter(r.overlap,r.centroid_distance,s=70 if sup else 38,alpha=0.9 if sup else 0.45)
        if ds:
            ax.scatter(r.overlap,r.centroid_distance,s=145,facecolors="none",edgecolors="black",linewidths=1.0)
        ax.annotate(label,(r.overlap,r.centroid_distance),xytext=(4,4),textcoords="offset points",fontsize=8)
    ax.set_xlabel("Gaussian Bhattacharyya overlap (lower = stronger separation)")
    ax.set_ylabel("Environmental centroid distance (higher = stronger separation)")
    ax.set_title("Environmental sorting of low- versus high-trait taxa")
    ax.text(0.01,-0.13,"Larger points: at least one all-taxa metric BH q < 0.05. Black ring: direct-backbone subset has at least one supported metric.",transform=ax.transAxes,fontsize=7.5)
    fig.tight_layout()
    save(fig,out,"Figure_S8_environmental_sorting")
    metrics.to_csv(out/"FigureS8_environmental_sorting_metrics.csv",index=False)


def main() -> None:
    a=parse_args(); out=Path(a.out_dir); out.mkdir(parents=True,exist_ok=True)
    env=pd.read_csv(a.environment,low_memory=False)
    fx=pd.read_csv(a.spde_fixed,low_memory=False)
    primary=pd.read_csv(a.primary_supported)
    bias=pd.read_csv(a.bias_control)
    niche=pd.read_csv(a.niche_summary)
    pgls=pd.read_csv(a.pgls_supported)
    scores,load,var=pca_products(env,out)
    plot_figure4(scores,load,var,out)
    stable=plot_figure5(primary,bias,fx,out)
    plot_figure6(primary,pgls,out)
    plot_s8(env,niche,out)
    products=sorted(p.name for p in out.iterdir() if p.suffix in {".png",".svg",".pdf"})
    if len(products)!=12:
        raise ValueError(f"Expected 12 rendered files, got {len(products)}")
    report={"pca_complete_taxa":len(scores),"pca_variance_pc1_pc3":[float(v) for v in var[:3]],"primary_supported":len(primary),"spde_stable_pairs":len(stable),"bias_control_rows":len(bias),"historical_figure_status":"Supporting Information only","rendered_products":products,"scope":"presentation of frozen or predeclared audit results; no post-outcome model selection"}
    (out/"interpretive_figure_report.json").write_text(json.dumps(report,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(report,indent=2))

if __name__=="__main__":
    main()
