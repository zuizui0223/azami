"""Expanded Supporting Information figures for the frozen GEB release."""
from __future__ import annotations

from geb_submission_release_common import *
from geb_submission_release_main_figures import figure_5, plot_loading_biplot

def figure_s110(
    primary: PCAResult,
    expanded_scores: pd.DataFrame,
    expanded_loadings: pd.DataFrame,
    expanded_explained: Sequence[float],
    module_map: dict[str, str],
    out: Path,
) -> list[Path]:
    fig = plt.figure(figsize=(11.8, 8.6), constrained_layout=True)
    gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 0.72])
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[1, :])
    plot_loading_biplot(ax_a, primary.scores, primary.loadings, 1, 2, primary.explained, module_map)
    ax_a.set_title(f"Established nine endpoints (n = {len(primary.scores)} taxa)")
    plot_loading_biplot(ax_b, expanded_scores, expanded_loadings, 1, 2, expanded_explained, module_map)
    ax_b.set_title(f"Expanded 18 endpoints (n = {len(expanded_scores)} taxa)")
    x = np.arange(3)
    width = 0.34
    ax_c.bar(x - width / 2, np.asarray(primary.explained[:3]) * 100, width, label="Nine endpoints", color="#777777")
    ax_c.bar(x + width / 2, np.asarray(expanded_explained[:3]) * 100, width, label="18 endpoints", color="#0072B2")
    ax_c.set_xticks(x, ["PC1", "PC2", "PC3"])
    ax_c.set_ylabel("Explained variance (%)")
    ax_c.set_title(
        f"Cumulative PC1–PC3: {100*np.sum(primary.explained[:3]):.1f}% vs "
        f"{100*np.sum(expanded_explained[:3]):.1f}%"
    )
    ax_c.legend(frameon=False)
    for label, ax in zip("abc", [ax_a, ax_b, ax_c]):
        add_panel_label(ax, label)
    fig.suptitle("Figure S1.10. Established and expanded taxon morphospaces", fontsize=13, fontweight="bold")
    return save_figure(fig, out / "Figure_S1_10_primary_vs_expanded_PCA")


def figure_s111(
    scores_all: pd.DataFrame,
    report: dict[str, object],
    out: Path,
) -> list[Path]:
    scopes = [
        ("module:colour", "Visible colour"),
        ("module:shape", "Outline shape"),
        ("module:involucre_architecture", "Involucre architecture"),
        ("module:armature", "Armature"),
    ]
    report_map = {str(item["scope"]): item for item in report["pca_scopes"]}
    fig, axes = plt.subplots(2, 2, figsize=(10.8, 8.6), constrained_layout=True)
    for label, ax, (scope, title) in zip("abcd", axes.ravel(), scopes):
        df = scores_all[scores_all["scope"].eq(scope)]
        if df.empty:
            raise ValueError(f"Missing module PCA scope: {scope}")
        rep = report_map[scope]
        ev = rep["explained_variance"]
        ax.scatter(df["PC1"], df["PC2"], s=15, alpha=0.65, color=MODULE_COLOUR[scope.split(":", 1)[1]])
        ax.axhline(0, color="#BBBBBB", lw=0.6)
        ax.axvline(0, color="#BBBBBB", lw=0.6)
        ax.set_xlabel(f"PC1 ({100*ev['PC1']:.1f}%)")
        ax.set_ylabel(f"PC2 ({100*ev['PC2']:.1f}%)")
        ax.set_title(f"{title} (n = {len(df)} taxa)")
        add_panel_label(ax, label)
    fig.suptitle("Figure S1.11. Module-specific taxon morphospaces", fontsize=13, fontweight="bold")
    return save_figure(fig, out / "Figure_S1_11_module_specific_PCA")


def matrix_from_long(
    frame: pd.DataFrame,
    rows: Sequence[str],
    columns: Sequence[str],
    value: str,
    row_key: str = "endpoint_id",
    col_key: str = "predictor",
) -> np.ndarray:
    lookup = {(getattr(r, row_key), getattr(r, col_key)): getattr(r, value) for r in frame.itertuples(index=False)}
    out = np.full((len(rows), len(columns)), np.nan, dtype=float)
    for i, row in enumerate(rows):
        for j, column in enumerate(columns):
            if (row, column) in lookup:
                out[i, j] = float(lookup[(row, column)])
    return out


def figure_s112(quality: pd.DataFrame, contract: pd.DataFrame, out: Path) -> list[Path]:
    module_map = endpoint_module_map(contract)
    endpoints = sorted_units(quality["endpoint_id"].unique(), module_map)
    beta = matrix_from_long(quality, endpoints, PREDICTOR_ORDER, "quality_adjusted_beta")
    category_order = [
        "quality_adjusted_not_fdr_supported",
        "image_sensitive_or_uncertain",
        "quality_robust",
    ]
    category_code = {value: i for i, value in enumerate(category_order)}
    cat = np.full((len(endpoints), 4), np.nan)
    q = matrix_from_long(quality, endpoints, PREDICTOR_ORDER, "quality_adjusted_q")
    lookup = {(r.endpoint_id, r.predictor): category_code[r.v2_interpretation] for r in quality.itertuples(index=False)}
    for i, endpoint in enumerate(endpoints):
        for j, predictor in enumerate(PREDICTOR_ORDER):
            cat[i, j] = lookup[(endpoint, predictor)]
    fig, axes = plt.subplots(1, 2, figsize=(11.8, 7.5), constrained_layout=True)
    limit = max(0.08, float(np.nanmax(np.abs(beta))))
    im = axes[0].imshow(beta, cmap="RdBu_r", norm=TwoSlopeNorm(vmin=-limit, vcenter=0, vmax=limit), aspect="auto")
    axes[0].set_xticks(range(4), [PREDICTOR_LABEL[p] for p in PREDICTOR_ORDER])
    axes[0].set_yticks(range(len(endpoints)), [endpoint_label(e, True) for e in endpoints])
    for i in range(beta.shape[0]):
        for j in range(beta.shape[1]):
            star = "*" if q[i, j] < 0.05 else ""
            axes[0].text(j, i, f"{beta[i,j]:.3f}{star}", ha="center", va="center", fontsize=6,
                         color="white" if abs(beta[i,j]) > 0.55*limit else "black")
    axes[0].set_title("Quality-adjusted standardized coefficients\n* candidate-quality BH q < 0.05")
    plt.colorbar(im, ax=axes[0], fraction=0.046, pad=0.02, label="Adjusted beta")
    cmap = ListedColormap(["#F2F2F2", "#D55E00", "#009E73"])
    axes[1].imshow(cat, cmap=cmap, norm=BoundaryNorm([-0.5,0.5,1.5,2.5], cmap.N), aspect="auto")
    axes[1].set_xticks(range(4), [PREDICTOR_LABEL[p] for p in PREDICTOR_ORDER])
    axes[1].set_yticks(range(len(endpoints)), [endpoint_label(e, True) for e in endpoints])
    symbols = {0: "not FDR", 1: "sensitive", 2: "robust"}
    for i in range(cat.shape[0]):
        for j in range(cat.shape[1]):
            axes[1].text(j, i, symbols[int(cat[i,j])], ha="center", va="center", fontsize=5.6,
                         color="white" if cat[i,j] in {1,2} else "black")
    axes[1].set_title("Resolution-stratum interpretation")
    for label, ax in zip("ab", axes):
        add_panel_label(ax, label)
    fig.suptitle("Figure S1.12. Full candidate image-quality sensitivity audit", fontsize=13, fontweight="bold")
    return save_figure(fig, out / "Figure_S1_12_candidate_quality")


def matched_beta_matrix(alignment: pd.DataFrame, endpoints: Sequence[str], scale: str) -> tuple[np.ndarray, np.ndarray]:
    values = np.full((len(endpoints), len(PREDICTOR_ORDER)), np.nan)
    supported = np.zeros_like(values, dtype=bool)
    lookup = {(r.endpoint_id, r.predictor): r for r in alignment.itertuples(index=False)}
    for i, endpoint in enumerate(endpoints):
        for j, predictor in enumerate(PREDICTOR_ORDER):
            row = lookup[(endpoint, predictor)]
            if endpoint == "corolla_hue":
                value = row.effect_magnitude if scale == "within" else row.effect_magnitude_among
            else:
                value = row.beta_std_within if scale == "within" else row.beta_std_among
            values[i, j] = float(value)
            supported[i, j] = bool(row.within_fdr_significant_0_05 if scale == "within" else row.among_fdr_significant_0_05)
    return values, supported


def figure_s113(alignment: pd.DataFrame, contract: pd.DataFrame, out: Path) -> list[Path]:
    module_map = endpoint_module_map(contract)
    endpoints = sorted_units(alignment["endpoint_id"].unique(), module_map)
    within, ws = matched_beta_matrix(alignment, endpoints, "within")
    among, ass = matched_beta_matrix(alignment, endpoints, "among")
    limit = max(0.5, float(np.nanmax(np.abs(np.r_[within.ravel(), among.ravel()]))))
    fig, axes = plt.subplots(1, 2, figsize=(11.8, 8.0), constrained_layout=True)
    for label, ax, matrix, support, title in [
        ("a", axes[0], within, ws, "Within taxa"),
        ("b", axes[1], among, ass, "Among taxa"),
    ]:
        im = ax.imshow(matrix, cmap="RdBu_r", norm=TwoSlopeNorm(vmin=-limit, vcenter=0, vmax=limit), aspect="auto")
        ax.set_xticks(range(4), [PREDICTOR_LABEL[p] for p in PREDICTOR_ORDER])
        ax.set_yticks(range(len(endpoints)), [endpoint_label(e, True) for e in endpoints])
        for i in range(matrix.shape[0]):
            for j in range(matrix.shape[1]):
                star = "*" if support[i,j] else ""
                ax.text(j, i, f"{matrix[i,j]:.2f}{star}", ha="center", va="center", fontsize=5.8,
                        color="white" if abs(matrix[i,j]) > 0.55*limit else "black")
        ax.set_title(title + "\n* frozen BH support")
        add_panel_label(ax, label)
    plt.colorbar(im, ax=axes, fraction=0.025, pad=0.02, label="Standardized beta; hue shown as circular magnitude")
    fig.suptitle("Figure S1.13. Full matched 68-row coefficient and support map", fontsize=13, fontweight="bold")
    return save_figure(fig, out / "Figure_S1_13_matched_coefficients")


def environment_column_order() -> list[tuple[str, str, str]]:
    return [
        ("complete18_env_min5", "within_taxon", "≥5 within"),
        ("complete18_env_min5", "among_taxon", "≥5 among"),
        ("complete18_env_min2", "within_taxon", "≥2 within"),
        ("complete18_env_min2", "among_taxon", "≥2 among"),
    ]


def figure_s115(blocks: pd.DataFrame, out: Path) -> list[Path]:
    values = np.full((len(BLOCK_ORDER), 4), np.nan)
    nominal = np.zeros_like(values, dtype=bool)
    fdr = np.zeros_like(values, dtype=bool)
    for i, block in enumerate(BLOCK_ORDER):
        for j, (scope, scale, _label) in enumerate(environment_column_order()):
            row = blocks[blocks["scope"].eq(scope) & blocks["scale"].eq(scale) & blocks["block_id"].eq(block)].iloc[0]
            values[i,j] = float(row["multivariate_r2"])
            nominal[i,j] = float(row["permutation_p"]) < 0.05
            fdr[i,j] = bool(row["fdr_supported_0_05"])
    fig, ax = plt.subplots(figsize=(8.8, 5.9), constrained_layout=True)
    im = ax.imshow(values, cmap="YlGnBu", vmin=0, vmax=max(0.08, float(np.nanmax(values))), aspect="auto")
    ax.set_yticks(range(len(BLOCK_ORDER)), [BLOCK_LABEL[b] for b in BLOCK_ORDER])
    ax.set_xticks(range(4), [x[2] for x in environment_column_order()])
    for i in range(values.shape[0]):
        for j in range(values.shape[1]):
            suffix = "**" if fdr[i,j] else ("†" if nominal[i,j] else "")
            ax.text(j, i, f"{values[i,j]:.3f}{suffix}", ha="center", va="center", fontsize=7,
                    color="white" if values[i,j] > 0.55*np.nanmax(values) else "black")
    ax.set_title("Figure S1.15. Stand-alone environmental-block $R^2$\n† nominal P < 0.05; ** block-family BH support (none)")
    plt.colorbar(im, ax=ax, fraction=0.046, pad=0.02, label="Multivariate $R^2$")
    return save_figure(fig, out / "Figure_S1_15_environment_blocks")


def figure_s116(incremental: pd.DataFrame, out: Path) -> list[Path]:
    row_defs = [
        ("omnibus", "all_process_extension", "All process variables"),
        ("block_specific", "radiative_atmospheric_drying", "Radiation + VPD"),
        ("block_specific", "mechanical_exposure", "Wind"),
        ("block_specific", "growing_season_water_input", "Growing-season precipitation"),
        ("block_specific", "climatic_productivity", "Potential NPP"),
    ]
    values = np.full((len(row_defs), 4), np.nan)
    support = np.zeros_like(values, dtype=bool)
    pvals = np.full_like(values, np.nan)
    for i, (family, block, _label) in enumerate(row_defs):
        for j, (scope, scale, _label2) in enumerate(environment_column_order()):
            row = incremental[
                incremental["scope"].eq(scope)
                & incremental["scale"].eq(scale)
                & incremental["test_family"].eq(family)
                & incremental["block_id"].eq(block)
            ].iloc[0]
            values[i,j] = float(row["partial_r2"])
            support[i,j] = bool(row["supported_0_05"])
            pvals[i,j] = float(row["permutation_p"])
    fig, ax = plt.subplots(figsize=(9.0, 5.6), constrained_layout=True)
    im = ax.imshow(values, cmap="YlGnBu", vmin=0, vmax=max(0.22, float(np.nanmax(values))), aspect="auto")
    ax.set_yticks(range(len(row_defs)), [x[2] for x in row_defs])
    ax.set_xticks(range(4), [x[2] for x in environment_column_order()])
    for i in range(values.shape[0]):
        for j in range(values.shape[1]):
            star = "*" if support[i,j] else ""
            ax.text(j, i, f"{values[i,j]:.3f}{star}\nP={pvals[i,j]:.3g}", ha="center", va="center", fontsize=6.3,
                    color="white" if values[i,j] > 0.55*np.nanmax(values) else "black")
    ax.set_title("Figure S1.16. Full nested increments beyond the four-variable core\n* predeclared support rule passed")
    plt.colorbar(im, ax=ax, fraction=0.046, pad=0.02, label="Partial $R^2$")
    return save_figure(fig, out / "Figure_S1_16_incremental_tests")


def figure_s117(geometry: pd.DataFrame, out: Path) -> list[Path]:
    fig, ax = plt.subplots(figsize=(9.7, 7.2), constrained_layout=True)
    rows: list[pd.Series] = []
    for block in BLOCK_ORDER:
        for scope in ["complete18_env_min5", "complete18_env_min2"]:
            rows.append(geometry[geometry["block_id"].eq(block) & geometry["scope"].eq(scope)].iloc[0])
    y = np.arange(len(rows))[::-1]
    labels = []
    for yi, row in zip(y, rows):
        estimate = float(row["coefficient_matrix_cosine_within_vs_among"])
        low = float(row["bootstrap_ci95_low"])
        high = float(row["bootstrap_ci95_high"])
        scope = str(row["scope"])
        marker = "o" if scope.endswith("min5") else "s"
        colour = "#0072B2" if scope.endswith("min5") else "#D55E00"
        ax.errorbar(estimate, yi, xerr=[[estimate-low],[high-estimate]], fmt=marker, color=colour, capsize=3, lw=1)
        labels.append(f"{BLOCK_LABEL[str(row['block_id'])]} — {'≥5' if scope.endswith('min5') else '≥2'}")
    ax.axvline(0, color="black", ls="--", lw=0.9)
    ax.set_yticks(y, labels)
    ax.set_xlabel("Cosine similarity of within- and among-taxon coefficient matrices")
    ax.set_title("Figure S1.17. Cross-scale coefficient geometry is inconclusive\nAll bootstrap intervals include zero")
    ax.set_xlim(-0.65, 0.65)
    return save_figure(fig, out / "Figure_S1_17_cross_scale_cosines")
