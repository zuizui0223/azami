"""Render saved linear contrasts and exposure support, without model refitting."""
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from analysis.v4.plot_colour_response import ENGLISH

ENDPOINTS = {
    "corolla_lab_lightness": "Lightness (L*)",
    "corolla_lab_chroma": "Chroma (C*)",
    "corolla_hue_sin": "HSV hue sine",
    "corolla_hue_cos": "HSV hue cosine",
    "orientation_image_vertical_angle": "Image-vertical angle (degrees)",
    "bract_projection_p95": "Relative contour projection (SI)",
    "bract_spread_fraction": "Contour spread fraction (SI)",
}


def plot(result, output):
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    for model in result["models"].values():
        for endpoint, component in model.get("components", {}).items():
            if not component.get("gradients"):
                continue
            variables = result["predictors"]
            fig, axes = plt.subplots(len(variables), 2,
                                     figsize=(12, 2.25*len(variables)), squeeze=False)
            for i, variable in enumerate(variables):
                for j, level in enumerate(("between", "within")):
                    term = variable + "__" + level
                    curve = component["gradients"][term]
                    support = result["environment_diagnostics"][endpoint]["exposure_support"][term]
                    ax = axes[i, j]
                    color = "#245C91" if j == 0 else "#BC591D"
                    ax.axhline(0, color="grey", lw=.6)
                    ax.axvline(curve["reference"], color="grey", lw=.6, ls=":")
                    ax.fill_between(curve["x"], curve["lower"], curve["upper"],
                                    color=color, alpha=.18)
                    ax.plot(curve["x"], curve["change"], color=color)
                    ax.set_xlabel(ENGLISH.get(variable, variable) +
                                  (" species mean" if j == 0 else " deviation from species mean") +
                                  "\n[stored environmental units]")
                    ax.set_ylabel("Change from median")
                    if j == 1:
                        ax.text(.02, .96,
                                f"Species spanning shown range: {support['species_spanning_q05_q95']}"
                                f"/{support['species_n']}",
                                transform=ax.transAxes, va="top", fontsize=7)
                    ax.spines[["top", "right"]].set_visible(False)
            axes[0, 0].set_title("Between species")
            axes[0, 1].set_title("Within species")
            fig.suptitle(f"{ENDPOINTS.get(endpoint, endpoint)} | n={component['n']:,}, "
                         f"species={component['n_species']}", fontsize=13, y=.998)
            fig.text(.02, .006,
                     "Linear conditional contrasts; other environmental terms held fixed. Bands: pointwise 95% slope intervals.\n"
                     "Shown ranges: marginal 5th–95th percentiles, not joint environmental support or individual predictions.\n"
                     "Colour components must be interpreted with the joint test; hue coordinates are not pigment concentrations.",
                     fontsize=8)
            fig.tight_layout(rect=(0, .05, 1, .98))
            fig.savefig(output / (endpoint + "_gradients.png"), dpi=150)
            fig.savefig(output / (endpoint + "_gradients.pdf"))
            plt.close(fig)


if __name__ == "__main__":
    plot(json.loads(Path(sys.argv[1]).read_text()), Path(sys.argv[2]))
