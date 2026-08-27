"""Main Figures 4–6 for the frozen GEB submission release."""
from __future__ import annotations

from geb_submission_release_common import *

def plot_loading_biplot(
    ax: plt.Axes,
    scores: pd.DataFrame,
    loadings: pd.DataFrame,
    xpc: int,
    ypc: int,
    explained: Sequence[float],
    module_map: dict[str, str],
) -> None:
    xcol, ycol = f"PC{xpc}", f"PC{ypc}"
    ax.scatter(scores[xcol], scores[ycol], s=13, alpha=0.58, color="#666666", linewidths=0)
    xr = float(scores[xcol].max() - scores[xcol].min())
    yr = float(scores[ycol].max() - scores[ycol].min())
    lx = loadings[f"PC{xpc}_loading"].to_numpy(float)
    ly = loadings[f"PC{ypc}_loading"].to_numpy(float)
    nonzero_x = np.nanmax(np.abs(lx)) or 1.0
    nonzero_y = np.nanmax(np.abs(ly)) or 1.0
    scale = min(0.30 * xr / nonzero_x, 0.30 * yr / nonzero_y)
    label_magnitude = np.sqrt(lx**2 + ly**2)
    label_count = min(len(loadings), 10)
    labelled = set(loadings.iloc[np.argsort(label_magnitude)[-label_count:]]["endpoint_id"])
    for row in loadings.itertuples(index=False):
        endpoint = row.endpoint_id
        x = float(getattr(row, f"PC{xpc}_loading")) * scale
        y = float(getattr(row, f"PC{ypc}_loading")) * scale
        module = module_map.get(endpoint, "shape")
        colour = MODULE_COLOUR.get(module, "#333333")
        ax.annotate(
            "",
            xy=(x, y),
            xytext=(0, 0),
            arrowprops={"arrowstyle": "-|>", "lw": 0.8, "color": colour, "alpha": 0.9},
        )
        if endpoint in labelled:
            ax.text(x * 1.04, y * 1.04, endpoint_label(endpoint, True), color=colour, fontsize=5.8,
                    ha="left" if x >= 0 else "right", va="bottom" if y >= 0 else "top")
    ax.axhline(0, color="#BBBBBB", lw=0.6)
    ax.axvline(0, color="#BBBBBB", lw=0.6)
    ax.set_xlabel(f"PC{xpc} ({100 * explained[xpc-1]:.1f}%)")
    ax.set_ylabel(f"PC{ypc} ({100 * explained[ypc-1]:.1f}%)")


def loading_heatmap(
    ax: plt.Axes,
    loadings: pd.DataFrame,
    module_map: dict[str, str],
    title: str,
) -> None:
    temp = loadings.copy()
    temp["module"] = temp["endpoint_id"].map(module_map)
    temp["module_rank"] = temp["module"].map(module_rank)
    temp = temp.sort_values(["module_rank", "endpoint_id"])
    values = temp[["PC1_loading", "PC2_loading", "PC3_loading"]].to_numpy(float)
    limit = max(0.35, float(np.nanmax(np.abs(values))))
    im = ax.imshow(values, aspect="auto", cmap="RdBu_r", norm=TwoSlopeNorm(vmin=-limit, vcenter=0, vmax=limit))
    ax.set_xticks(range(3), ["PC1", "PC2", "PC3"])
    ax.set_yticks(range(len(temp)), [endpoint_label(v, True) for v in temp["endpoint_id"]])
    for i in range(values.shape[0]):
        for j in range(values.shape[1]):
            ax.text(j, i, f"{values[i, j]:.2f}", ha="center", va="center", fontsize=5.6,
                    color="white" if abs(values[i, j]) > 0.55 * limit else "black")
    last_module = None
    for i, module in enumerate(temp["module"]):
        if last_module is not None and module != last_module:
            ax.axhline(i - 0.5, color="black", lw=0.7)
        last_module = module
    ax.set_title(title)
    plt.colorbar(im, ax=ax, fraction=0.046, pad=0.02, label="Loading")


def figure_4(
    scores: pd.DataFrame,
    loadings: pd.DataFrame,
    explained: Sequence[float],
    module_map: dict[str, str],
    out: Path,
) -> list[Path]:
    fig = plt.figure(figsize=(12.0, 9.2), constrained_layout=True)
    gs = fig.add_gridspec(2, 2, width_ratios=[1.05, 1.0], height_ratios=[1.0, 1.0])
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[1, 0])
    ax_d = fig.add_subplot(gs[1, 1])
    plot_loading_biplot(ax_a, scores, loadings, 1, 2, explained, module_map)
    ax_a.set_title("Expanded taxon morphospace: PC1–PC2")
    plot_loading_biplot(ax_b, scores, loadings, 1, 3, explained, module_map)
    ax_b.set_title("Expanded taxon morphospace: PC1–PC3")
    x = np.arange(1, 4)
    bars = ax_c.bar(x, np.asarray(explained[:3]) * 100, color="#777777", edgecolor="black", linewidth=0.6)
    cumulative = np.cumsum(explained[:3]) * 100
    ax_c.plot(x, cumulative, marker="o", color="#0072B2", label="Cumulative")
    for bar, value in zip(bars, np.asarray(explained[:3]) * 100):
        ax_c.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.7, f"{value:.1f}%", ha="center", fontsize=8)
    ax_c.set_xticks(x, ["PC1", "PC2", "PC3"])
    ax_c.set_ylabel("Explained variance (%)")
    ax_c.set_ylim(0, max(48, float(cumulative[-1]) + 5))
    ax_c.legend(frameon=False, loc="upper left")
    ax_c.set_title(f"No dominant axis: PC1–PC3 = {cumulative[-1]:.1f}%")
    loading_heatmap(ax_d, loadings, module_map, "Endpoint loadings on the first three axes")
    for label, ax in zip("abcd", [ax_a, ax_b, ax_c, ax_d]):
        add_panel_label(ax, label)
    handles = [
        Line2D([0], [0], color=MODULE_COLOUR[m], lw=2, label=MODULE_LABEL[m])
        for m in MODULE_ORDER
    ]
    fig.legend(handles=handles, loc="lower center", ncol=5, frameon=False, bbox_to_anchor=(0.5, -0.015))
    fig.suptitle("Figure 4. Multidimensional taxon morphospace of thistle capitula", fontsize=13, fontweight="bold")
    return save_figure(fig, out / "Figure_4_morphospace")


def plot_unit_matrix(
    ax: plt.Axes,
    matrix: np.ndarray,
    units: Sequence[str],
    module_map: dict[str, str],
    title: str,
) -> None:
    masked = np.ma.masked_invalid(matrix)
    cmap = plt.get_cmap("cividis").copy()
    cmap.set_bad("#F2F2F2")
    im = ax.imshow(masked, vmin=0, vmax=max(0.65, float(np.nanmax(matrix))), cmap=cmap, aspect="equal")
    labels = [endpoint_label(v, True) for v in units]
    ax.set_xticks(range(len(units)), labels, rotation=90)
    ax.set_yticks(range(len(units)), labels)
    modules = [module_map.get(v, "") for v in units]
    for i in range(1, len(units)):
        if modules[i] != modules[i - 1]:
            ax.axhline(i - 0.5, color="white", lw=1.3)
            ax.axvline(i - 0.5, color="white", lw=1.3)
    ax.set_title(title)
    plt.colorbar(im, ax=ax, fraction=0.046, pad=0.02, label="Association strength")


def pairwise_alignment(
    pair_table: pd.DataFrame,
    scope: str,
    module_map: dict[str, str],
) -> pd.DataFrame:
    subset = pair_table[pair_table["scope"].eq(scope)].copy()
    w = subset[subset["scale"].eq("within_taxon")][["left", "right", "value"]].rename(columns={"value": "within"})
    a = subset[subset["scale"].eq("among_taxon")][["left", "right", "value"]].rename(columns={"value": "among"})
    merged = w.merge(a, on=["left", "right"], how="inner", validate="one_to_one")
    merged["same_module"] = [module_map.get(l) == module_map.get(r) for l, r in zip(merged["left"], merged["right"])]
    return merged


def plot_contrast_forest(ax: plt.Axes, summary: pd.DataFrame) -> None:
    rows: list[tuple[str, float, float, float, str, str]] = []
    for scope, threshold in [("complete18_min5", "≥5 observations/taxon"), ("complete18_min2", "≥2 observations/taxon")]:
        row = summary[summary["scope"].eq(scope)].iloc[0]
        rows.extend(
            [
                (f"{threshold}: within", row["within_module_contrast"], row["within_module_contrast_bootstrap_ci95_low"], row["within_module_contrast_bootstrap_ci95_high"], "within", scope),
                (f"{threshold}: among", row["among_module_contrast"], row["among_module_contrast_bootstrap_ci95_low"], row["among_module_contrast_bootstrap_ci95_high"], "among", scope),
            ]
        )
    y = np.arange(len(rows))[::-1]
    for yi, (label, estimate, low, high, scale, _scope) in zip(y, rows):
        marker = "o" if scale == "within" else "s"
        colour = "#0072B2" if scale == "within" else "#D55E00"
        ax.errorbar(estimate, yi, xerr=[[estimate - low], [high - estimate]], fmt=marker, color=colour,
                    capsize=3, markersize=5, lw=1.1)
    ax.axvline(0, color="#777777", lw=0.8, ls="--")
    ax.set_yticks(y, [r[0] for r in rows])
    ax.set_xlabel("Within-module minus between-module strength")
    ax.set_title("Registered-module organization and sensitivity")


def figure_5(
    pair_table: pd.DataFrame,
    summary: pd.DataFrame,
    contract: pd.DataFrame,
    out: Path,
    scope: str = "complete18_min5",
    number: str = "Figure_5_multilevel_organization",
    title_prefix: str = "Figure 5",
) -> list[Path]:
    module_map = endpoint_module_map(contract)
    subset = pair_table[pair_table["scope"].eq(scope)]
    units = sorted_units(pd.concat([subset["left"], subset["right"]]), module_map)
    within = symmetric_matrix(subset[subset["scale"].eq("within_taxon")], units)
    among = symmetric_matrix(subset[subset["scale"].eq("among_taxon")], units)
    alignment = pairwise_alignment(pair_table, scope, module_map)
    scope_row = summary[summary["scope"].eq(scope)].iloc[0]

    fig = plt.figure(figsize=(12.2, 9.7), constrained_layout=True)
    gs = fig.add_gridspec(2, 2, height_ratios=[1.08, 0.75])
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[1, 0])
    ax_d = fig.add_subplot(gs[1, 1])
    plot_unit_matrix(ax_a, within, units, module_map, "Within-taxon association strength")
    plot_unit_matrix(ax_b, among, units, module_map, "Among-taxon association strength")
    plot_contrast_forest(ax_c, summary)
    for same, group in alignment.groupby("same_module", sort=False):
        ax_d.scatter(
            group["within"],
            group["among"],
            s=28 if same else 20,
            marker="o" if same else "x",
            color="#0072B2" if same else "#777777",
            alpha=0.78,
            label="Same registered module" if same else "Between modules",
        )
    lim = max(float(alignment[["within", "among"]].max().max()) * 1.05, 0.55)
    ax_d.plot([0, lim], [0, lim], color="#AAAAAA", ls="--", lw=0.8)
    ax_d.set_xlim(0, lim)
    ax_d.set_ylim(0, lim)
    ax_d.set_xlabel("Within-taxon association strength")
    ax_d.set_ylabel("Among-taxon association strength")
    ax_d.set_title(
        "Partial cross-scale correspondence\n"
        f"Spearman = {scope_row['cross_scale_strength_matrix_spearman']:.3f} "
        f"[{scope_row['cross_scale_spearman_bootstrap_ci95_low']:.3f}, "
        f"{scope_row['cross_scale_spearman_bootstrap_ci95_high']:.3f}]"
    )
    ax_d.legend(frameon=False, loc="lower right")
    for label, ax in zip("abcd", [ax_a, ax_b, ax_c, ax_d]):
        add_panel_label(ax, label)
    fig.suptitle(
        f"{title_prefix}. Partial organization of the same 17-unit capitulum phenotype within and among taxa",
        fontsize=13,
        fontweight="bold",
    )
    return save_figure(fig, out / number)


def alignment_matrix(alignment: pd.DataFrame, contract: pd.DataFrame) -> tuple[list[str], np.ndarray]:
    module_map = endpoint_module_map(contract)
    endpoints = sorted_units(alignment["endpoint_id"].unique(), module_map)
    code = {"neither": 0, "within_only": 1, "among_only": 2, "both_scales": 3}
    matrix = np.zeros((len(endpoints), len(PREDICTOR_ORDER)), dtype=int)
    lookup = {(r.endpoint_id, r.predictor): code[r.cross_scale_class] for r in alignment.itertuples(index=False)}
    for i, endpoint in enumerate(endpoints):
        for j, predictor in enumerate(PREDICTOR_ORDER):
            matrix[i, j] = lookup[(endpoint, predictor)]
    return endpoints, matrix


def incremental_focus_matrix(incremental: pd.DataFrame) -> tuple[np.ndarray, np.ndarray]:
    columns = [
        ("complete18_env_min5", "within_taxon"),
        ("complete18_env_min5", "among_taxon"),
        ("complete18_env_min2", "within_taxon"),
        ("complete18_env_min2", "among_taxon"),
    ]
    rows = [("omnibus", "all_process_extension"), ("block_specific", "growing_season_water_input")]
    values = np.full((2, 4), np.nan)
    supported = np.zeros((2, 4), dtype=bool)
    for i, (family, block) in enumerate(rows):
        for j, (scope, scale) in enumerate(columns):
            hit = incremental[
                incremental["scope"].eq(scope)
                & incremental["scale"].eq(scale)
                & incremental["test_family"].eq(family)
                & incremental["block_id"].eq(block)
            ]
            if len(hit) != 1:
                raise ValueError(f"Expected one incremental row for {scope}, {scale}, {family}, {block}")
            row = hit.iloc[0]
            values[i, j] = float(row["partial_r2"])
            supported[i, j] = bool(row["supported_0_05"])
    return values, supported


def figure_6(
    atlas: pd.DataFrame,
    alignment: pd.DataFrame,
    incremental: pd.DataFrame,
    contract: pd.DataFrame,
    out: Path,
) -> list[Path]:
    fig = plt.figure(figsize=(12.2, 8.2), constrained_layout=True)
    gs = fig.add_gridspec(2, 2, height_ratios=[0.80, 1.2], width_ratios=[0.88, 1.12])
    ax_a = fig.add_subplot(gs[:, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[1, 1])

    counts = atlas["evidence_grade"].value_counts().reindex(EVIDENCE_GRADE_LABEL.keys()).dropna().astype(int)
    labels = [EVIDENCE_GRADE_LABEL[g] for g in counts.index]
    y = np.arange(len(counts))[::-1]
    ax_a.barh(y, counts.to_numpy(), color=[EVIDENCE_GRADE_COLOUR[g] for g in counts.index], edgecolor="black", lw=0.5)
    ax_a.set_yticks(y, labels)
    ax_a.set_xlabel("Atlas rows")
    ax_a.set_title("Evidence hierarchy is explicit, not collapsed")
    for yi, value in zip(y, counts.to_numpy()):
        ax_a.text(value + 0.12, yi, str(value), va="center", fontsize=8)
    ax_a.set_xlim(0, max(counts) + 1.4)
    candidate = atlas[atlas["evidence_grade"].str.startswith("C_")]
    candidate_text = "\n".join(
        f"• {endpoint_label(r.endpoint_id)} × {r.predictor.replace('chelsa_', '').upper()}"
        for r in candidate.itertuples(index=False)
    )
    ax_a.text(
        0.02,
        -0.09,
        "Exploratory candidate rows:\n" + candidate_text,
        transform=ax_a.transAxes,
        va="top",
        fontsize=7.5,
    )

    endpoints, support = alignment_matrix(alignment, contract)
    support_cmap = ListedColormap(["#F2F2F2", "#0072B2", "#D55E00", "#009E73"])
    norm = BoundaryNorm([-0.5, 0.5, 1.5, 2.5, 3.5], support_cmap.N)
    ax_b.imshow(support, cmap=support_cmap, norm=norm, aspect="auto")
    ax_b.set_xticks(range(4), [PREDICTOR_LABEL[p] for p in PREDICTOR_ORDER])
    ax_b.set_yticks(range(len(endpoints)), [endpoint_label(e, True) for e in endpoints])
    symbol = {0: "", 1: "W", 2: "A", 3: "W+A"}
    for i in range(support.shape[0]):
        for j in range(support.shape[1]):
            ax_b.text(j, i, symbol[int(support[i, j])], ha="center", va="center", fontsize=6.2,
                      color="white" if support[i, j] else "black", fontweight="bold")
    ax_b.set_title("Matched endpoint–climate support differs by scale")
    handles = [
        Line2D([0], [0], marker="s", linestyle="", markerfacecolor=support_cmap(i), markeredgecolor="black",
               label=label, markersize=7)
        for i, label in enumerate(["Neither", "Within only", "Among only", "Both scales"])
    ]
    ax_b.legend(handles=handles, frameon=False, ncol=4, loc="upper center", bbox_to_anchor=(0.5, -0.18))

    values, supported = incremental_focus_matrix(incremental)
    vmax = max(0.22, float(np.nanmax(values)))
    im = ax_c.imshow(values, cmap="YlGnBu", vmin=0, vmax=vmax, aspect="auto")
    ax_c.set_yticks([0, 1], ["All process variables", "Growing-season precipitation"])
    ax_c.set_xticks(
        range(4),
        ["≥5\nwithin", "≥5\namong", "≥2\nwithin", "≥2\namong"],
    )
    for i in range(values.shape[0]):
        for j in range(values.shape[1]):
            star = "*" if supported[i, j] else ""
            ax_c.text(j, i, f"{values[i, j]:.3f}{star}", ha="center", va="center", fontsize=8,
                      color="white" if values[i, j] > 0.55 * vmax else "black", fontweight="bold" if star else "normal")
    ax_c.set_title("Increment beyond BIO1/BIO4/BIO12/BIO15\n* predeclared support rule passed")
    plt.colorbar(im, ax=ax_c, fraction=0.046, pad=0.02, label="Partial $R^2$")

    for label, ax in zip("abc", [ax_a, ax_b, ax_c]):
        add_panel_label(ax, label)
    fig.suptitle("Figure 6. Evidence grade, biological scale and environmental sufficiency", fontsize=13, fontweight="bold")
    return save_figure(fig, out / "Figure_6_environmental_evidence")
