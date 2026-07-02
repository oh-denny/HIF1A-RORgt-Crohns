#!/usr/bin/env python3

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
import pandas as pd


PS_CMAP = LinearSegmentedColormap.from_list(
    "ps_rose_white_green", ["#a52391", "#ffffff", "#239b56"], N=256
)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Regenerate only the HIF1A KO perturbation-score panel."
    )
    parser.add_argument(
        "--result-dir",
        default="results/celloracle_hif1a_ko_1500_PS_markov_20260701_133734",
    )
    return parser.parse_args()


def savefig(fig, base):
    fig.savefig(base.with_suffix(".png"), dpi=400, bbox_inches="tight")
    fig.savefig(base.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def main():
    args = parse_args()
    result_dir = Path(args.result_dir)
    tables = result_dir / "tables"
    figures = result_dir / "figures"

    grid = pd.read_csv(
        tables / "HIF1A_KO_perturbation_score_grid.csv",
        dtype={"seurat_clusters": str},
    )
    cells = pd.read_csv(
        tables / "HIF1A_KO_markov_density_by_cell.csv",
        dtype={"seurat_clusters": str},
    )
    usable = grid[~grid["mass_filtered"]].copy()
    centers = cells.groupby("seurat_clusters")[["UMAP1", "UMAP2"]].median()

    fig, ax = plt.subplots(figsize=(6.2, 5.4))
    ax.scatter(
        cells["UMAP1"],
        cells["UMAP2"],
        s=5,
        color="#e5e7eb",
        alpha=0.65,
        linewidths=0,
        rasterized=True,
    )
    vmax = float(np.nanquantile(np.abs(usable["perturbation_score"]), 0.98))
    points = ax.scatter(
        usable["UMAP1"],
        usable["UMAP2"],
        c=usable["perturbation_score"],
        s=18,
        marker="s",
        cmap=PS_CMAP,
        vmin=-vmax,
        vmax=vmax,
        linewidths=0,
        rasterized=True,
    )

    rose = "#a52391"
    center4 = centers.loc["4"]
    center5 = centers.loc["5"]
    label_anchor = (center4 + center5) / 2
    label_position = (label_anchor["UMAP1"] - 0.9, label_anchor["UMAP2"] + 1.65)
    ax.annotate(
        "Progression inhibited\nclusters 4 and 5",
        xy=(center5["UMAP1"], center5["UMAP2"]),
        xytext=label_position,
        textcoords="data",
        color=rose,
        fontsize=8.5,
        weight="bold",
        ha="center",
        arrowprops={
            "arrowstyle": "->",
            "color": rose,
            "lw": 1.2,
            "connectionstyle": "arc3,rad=0.12",
        },
    )
    ax.annotate(
        "",
        xy=(center4["UMAP1"], center4["UMAP2"]),
        xytext=label_position,
        textcoords="data",
        arrowprops={
            "arrowstyle": "->",
            "color": rose,
            "lw": 1.2,
            "connectionstyle": "arc3,rad=-0.16",
        },
    )

    center3 = centers.loc["3"]
    ax.annotate(
        "Progression promoted\ncluster 3",
        xy=(center3["UMAP1"], center3["UMAP2"]),
        xytext=(-42, -42),
        textcoords="offset points",
        color="#15803d",
        fontsize=8.5,
        weight="bold",
        arrowprops={
            "arrowstyle": "->",
            "color": "#15803d",
            "lw": 1.2,
            "connectionstyle": "arc3,rad=0.12",
        },
    )

    colorbar = fig.colorbar(points, ax=ax, fraction=0.045, pad=0.025)
    colorbar.set_label(
        "Perturbation score\n(negative: inhibition; positive: promotion)",
        fontsize=8.5,
    )
    colorbar.ax.tick_params(labelsize=7.5, length=2)

    ax.set_title(
        "HIF1A KO perturbation score along the observed trajectory",
        fontsize=11,
    )
    ax.set_aspect("equal", adjustable="datalim")
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)

    output = figures / "HIF1A_KO_perturbation_score_PS_clusters4_5"
    savefig(fig, output)
    print(output.with_suffix(".png").resolve())
    print(output.with_suffix(".pdf").resolve())


if __name__ == "__main__":
    main()
