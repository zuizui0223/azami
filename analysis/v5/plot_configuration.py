"""Plot V5 comparisons on their original estimand scales."""
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

def main():
    root = Path(__file__).parent / "results/configuration"
    analogue = pd.read_csv(root / "v5_configuration_convergence.csv").query("dimensions == 3")
    overlap = pd.read_csv(root / "configuration_overlap.csv").query("dimensions == 3 and within_block")
    centroids = pd.read_csv(root / "matched_centroids.csv")
    fig, axes = plt.subplots(1, 3, figsize=(12, 3.8), layout="constrained")
    colors = ["#23736e", "#ab571e"]
    for i, row in enumerate(analogue.itertuples()):
        axes[0].plot([row.observed, row.independent_marginals], [i, i], color=".6")
        axes[0].scatter(row.observed, i, color=colors[0], label="Observed" if i == 0 else None)
        axes[0].scatter(row.independent_marginals, i, facecolor="white", edgecolor=colors[0],
                        label="Independent margins" if i == 0 else None)
    axes[0].set(yticks=[0, 1], yticklabels=["1-degree cells", "5-degree cells"], xlabel="Reduction in mean trait distance",
                title="Resemblance across taxa", ylim=(-.8, 1.45), xlim=(0,.067))
    for i, row in enumerate(overlap.itertuples()):
        axes[1].plot([row.observed_D_sim, row.observed_D_new], [i, i], color=".6")
        axes[1].scatter(row.observed_D_sim, i, facecolor="white", edgecolor=colors[0],
                        label="Native-like climate" if i == 0 else None)
        axes[1].scatter(row.observed_D_new, i, color=colors[0], label="Novel climate" if i == 0 else None)
    axes[1].set(yticks=[0, 1], yticklabels=["C. arvense","C. vulgare"], xlabel="Overlap with native heads (D)",
                title="Introduced head distributions", ylim=(-.8,1.45), xlim=(.69,.88))
    for i, (species, rows) in enumerate(centroids.groupby("species")):
        rows = rows.set_index("trait").loc[["angle", "chroma", "outline"]]
        axes[2].scatter(rows.difference_native_sd, np.arange(3)+(i-.5)*.13, color=colors[i],
                        label=species.replace("Cirsium", "C."))
    axes[2].axvline(0, color=".6", lw=1)
    axes[2].set(yticks=[0,1,2], yticklabels=["Angle","Chroma","Outline"], ylim=(-1,2.45),
                title="Direction of centroid shifts", xlabel="Novel minus native-like (native SD)")
    for ax in axes:
        ax.spines[["top","right","left"]].set_visible(False)
        ax.tick_params(axis="y", length=0)
        ax.legend(frameon=False, loc="lower center", fontsize=8)
    fig.savefig(root / "climate_and_heads.png", dpi=180)
    plt.close(fig)

if __name__ == "__main__":
    main()

