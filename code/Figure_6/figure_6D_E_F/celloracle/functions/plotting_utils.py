"""Publication-oriented plotting helpers shared by CellOracle pipelines."""

from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
import seaborn as sns

from .embedding_utils import arrow_scale, robust_limits
from .perturbation_utils import (
    celloracle_impact_summary,
    cluster_ps_summary,
    most_enriched_cluster,
    strongest_ps_clusters,
)


CLUSTER_COLORS = {
    "0": "#A61C4B",
    "1": "#0F6C84",
    "2": "#2C6B6B",
    "3": "#AEB8B7",
    "4": "#6C7A3A",
    "5": "#D8C7AF",
    "6": "#162E93",
}
PS_CMAP = LinearSegmentedColormap.from_list("ps_purple_white_green", ["#8e1b8e", "#ffffff", "#239b56"], N=256)
PS_CLUSTERS45_CMAP = LinearSegmentedColormap.from_list("ps_rose_white_green", ["#a52391", "#ffffff", "#239b56"], N=256)
MARKOV_CMAP = LinearSegmentedColormap.from_list("markov_blue", ["#deedf7", "#6baed6", "#08519c"], N=256)


def save_figure(fig, base, dpi=400):
    base = Path(base)
    fig.savefig(base.with_suffix(".png"), dpi=dpi, bbox_inches="tight", facecolor="white")
    fig.savefig(base.with_suffix(".pdf"), bbox_inches="tight", facecolor="white")
    plt.close(fig)


def cluster_palette(values):
    levels = sorted({str(value) for value in values}, key=lambda value: (len(value), value))
    missing = [level for level in levels if level not in CLUSTER_COLORS]
    fallback = sns.color_palette("tab10", n_colors=max(10, len(missing)))
    return {level: CLUSTER_COLORS.get(level, fallback[missing.index(level)] if level in missing else None) for level in levels}


def clean_axes(ax, x="FA1", y="FA2"):
    ax.set_xlabel(x)
    ax.set_ylabel(y)
    ax.set_aspect("equal", adjustable="datalim")
    ax.grid(False)
    sns.despine(ax=ax)


def apply_limits(ax, frame, x, y):
    ax.set_xlim(*robust_limits(frame[x]))
    ax.set_ylim(*robust_limits(frame[y]))


def clean_embedding_axis(ax):
    ax.set_aspect("equal", adjustable="datalim")
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)


def plot_propagation_magnitude(
    propagation,
    cluster_column,
    output_base,
    target_cluster=None,
    title="HIF1A KO propagation magnitude by cluster",
    value_column="mean_l2_norm",
    y_label="Mean L2 norm",
):
    """Plot unnormalized CellOracle propagation magnitude with the fixed cluster palette."""
    fig, ax = plt.subplots(figsize=(7.4, 4.6), facecolor="white")
    ax.set_facecolor("white")
    available = set(propagation[cluster_column].astype(str))
    order = [cluster for cluster in CLUSTER_COLORS if cluster in available]
    extras = sorted(available.difference(order), key=lambda value: (len(value), value))
    for cluster in order + extras:
        subset = propagation[propagation[cluster_column].astype(str) == cluster].sort_values("n_propagation")
        ax.plot(
            subset["n_propagation"], subset[value_column],
            color=CLUSTER_COLORS.get(cluster, "#4B5563"),
            linewidth=2.1 if cluster == str(target_cluster) else 1.7,
            marker="o", markersize=3.5, markeredgewidth=0,
            label=f"cluster {cluster}",
        )
    ax.set_xlabel("CellOracle propagation step", fontsize=10)
    ax.set_ylabel(y_label, fontsize=10)
    ax.set_title(title, fontsize=11)
    ax.tick_params(axis="both", labelsize=9)
    ax.grid(False)
    ax.legend(
        frameon=False, title="Cluster", title_fontsize=9, fontsize=8,
        ncol=1, bbox_to_anchor=(1.02, 1), loc="upper left",
    )
    sns.despine(ax=ax)
    fig.tight_layout()
    save_figure(fig, output_base)


def plot_clusters(frame, cluster_column, output_base):
    palette = cluster_palette(frame[cluster_column])
    fig, ax = plt.subplots(figsize=(7, 6))
    for cluster, subset in frame.groupby(cluster_column):
        ax.scatter(subset["FA1"], subset["FA2"], s=5, color=palette[str(cluster)], alpha=0.8, linewidths=0, label=str(cluster))
    for cluster, row in frame.groupby(cluster_column)[["FA1", "FA2"]].median().iterrows():
        ax.text(row["FA1"], row["FA2"], str(cluster), fontsize=10, weight="bold", ha="center", va="center")
    ax.set(xlabel="FA1", ylabel="FA2", title="HIF1A KO - FA clusters")
    ax.legend(title=cluster_column, bbox_to_anchor=(1.02, 1), loc="upper left", frameon=False, markerscale=3)
    sns.despine(ax=ax)
    save_figure(fig, output_base)


def plot_cell_arrows(vectors, cluster_column, output_base, sampled=None):
    draw = vectors if sampled is None or sampled >= len(vectors) else vectors.sample(sampled, random_state=123)
    colors = vectors[cluster_column].astype(str).map(cluster_palette(vectors[cluster_column]))
    factor = arrow_scale(vectors[["FA1", "FA2"]], vectors[["shift_FA1", "shift_FA2"]], fraction=0.055)
    fig, ax = plt.subplots(figsize=(7, 6))
    ax.scatter(vectors["FA1"], vectors["FA2"], s=4, color=colors, alpha=0.25, linewidths=0)
    ax.quiver(draw["FA1"], draw["FA2"], draw["shift_FA1"] * factor, draw["shift_FA2"] * factor,
              angles="xy", scale_units="xy", scale=1, width=0.0025, color="#111827", alpha=0.75)
    ax.set(xlabel="FA1", ylabel="FA2", title="HIF1A KO cell shift vectors")
    sns.despine(ax=ax)
    save_figure(fig, output_base)


def plot_grid_flow(embedding, grid, output_base):
    keep = grid["mass"] >= 0.01
    flow = grid.loc[keep, ["flow_FA1", "flow_FA2"]].to_numpy()
    factor = arrow_scale(embedding[["FA1", "FA2"]], flow, fraction=0.08)
    fig, ax = plt.subplots(figsize=(7, 6))
    ax.scatter(embedding["FA1"], embedding["FA2"], c="lightgray", s=5, alpha=0.25, linewidths=0)
    ax.quiver(grid.loc[keep, "grid_FA1"], grid.loc[keep, "grid_FA2"],
              grid.loc[keep, "flow_FA1"] * factor, grid.loc[keep, "flow_FA2"] * factor,
              angles="xy", scale_units="xy", scale=1, width=0.004, color="#111827")
    ax.set(xlabel="FA1", ylabel="FA2", title="HIF1A KO flow on FA grid")
    sns.despine(ax=ax)
    save_figure(fig, output_base)


def plot_delta_map(embedding, output_base):
    fig, ax = plt.subplots(figsize=(7, 6))
    points = ax.scatter(embedding["FA1"], embedding["FA2"], c=embedding["delta_l2"], s=7, cmap="magma", linewidths=0)
    fig.colorbar(points, ax=ax, label="delta L2")
    ax.set(xlabel="FA1", ylabel="FA2", title="HIF1A KO effect magnitude")
    sns.despine(ax=ax)
    save_figure(fig, output_base)


def plot_l2_distributions(embedding, cluster_column, group_column, figures_dir):
    fig, ax = plt.subplots(figsize=(7, 4))
    sns.violinplot(data=embedding, x=cluster_column, y="delta_l2", inner="box", color="#93c5fd", ax=ax)
    ax.set(xlabel="seurat_clusters", ylabel="delta L2", title="HIF1A KO delta L2 by cluster")
    sns.despine(ax=ax)
    save_figure(fig, figures_dir / "HIF1A_KO_delta_l2_by_cluster")
    fig, ax = plt.subplots(figsize=(8, 4))
    sns.boxplot(data=embedding, x=cluster_column, y="delta_l2", hue=group_column, fliersize=0.5, ax=ax)
    ax.set(xlabel="seurat_clusters", ylabel="delta L2", title="HIF1A KO delta L2 by group and cluster")
    ax.legend(title=group_column, frameon=False)
    sns.despine(ax=ax)
    save_figure(fig, figures_dir / "HIF1A_KO_delta_l2_by_group_cluster")


def plot_perturbation_score_map(score, output_base):
    fig, ax = plt.subplots(figsize=(7, 6))
    vmax = np.nanpercentile(np.abs(score["perturbation_score"]), 98)
    points = ax.scatter(score["FA1"], score["FA2"], c=score["perturbation_score"], s=7, cmap="coolwarm", vmin=-vmax, vmax=vmax, linewidths=0)
    fig.colorbar(points, ax=ax, label="perturbation score")
    ax.set(xlabel="FA1", ylabel="FA2", title="HIF1A KO perturbation score")
    sns.despine(ax=ax)
    save_figure(fig, output_base)


def delta_map(frame, base, title, x="FA1", y="FA2", highlight=None, cluster_column="seurat_clusters"):
    fig, ax = plt.subplots(figsize=(7.2, 6.2))
    ax.scatter(frame[x], frame[y], s=5, color="#d8dadd", alpha=0.55, linewidths=0)
    draw = frame if highlight is None else frame[frame[cluster_column] == highlight]
    points = ax.scatter(draw[x], draw[y], c=draw["delta_l2"], s=9 if highlight is None else 13,
                        cmap="magma", vmin=0, vmax=np.nanquantile(frame["delta_l2"], 0.98), linewidths=0, rasterized=True)
    if highlight is not None:
        center = draw[[x, y]].median()
        ax.text(center[x], center[y], f"cluster {highlight}", fontsize=9, weight="bold", ha="center", va="center",
                bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.7, "pad": 2})
    fig.colorbar(points, ax=ax, label="Magnitude of predicted transcriptomic change (delta L2)")
    ax.set_title(title)
    clean_axes(ax, x, y)
    apply_limits(ax, frame, x, y)
    save_figure(fig, base)


def arrow_plot(frame, sampled, base, title, x="FA1", y="FA2", u="shift_FA1", v="shift_FA2"):
    fig, ax = plt.subplots(figsize=(7.2, 6.2))
    ax.scatter(frame[x], frame[y], s=5, color="#d8dadd", alpha=0.55, linewidths=0)
    factor = arrow_scale(frame[[x, y]], frame[[u, v]], fraction=0.025)
    ax.quiver(sampled[x], sampled[y], sampled[u] * factor, sampled[v] * factor, angles="xy", scale_units="xy", scale=1,
              width=0.0022, headwidth=3.2, headlength=4.0, color="#111111", alpha=0.62)
    ax.set_title(title)
    clean_axes(ax, x, y)
    apply_limits(ax, frame, x, y)
    save_figure(fig, base)


def flow_plot(cells, grid, base, title, x="FA1", y="FA2", gx="grid_FA1", gy="grid_FA2", gu="flow_FA1", gv="flow_FA2"):
    kept = grid[grid["n_cells"] >= 5].copy()
    fig, ax = plt.subplots(figsize=(7.2, 6.2))
    ax.scatter(cells[x], cells[y], s=5, color="#d8dadd", alpha=0.55, linewidths=0)
    if len(kept):
        factor = arrow_scale(cells[[x, y]], kept[[gu, gv]], fraction=0.035)
        ax.quiver(kept[gx], kept[gy], kept[gu] * factor, kept[gv] * factor, angles="xy", scale_units="xy", scale=1,
                  width=0.0025, headwidth=3.2, color="#111111", alpha=0.62)
    ax.set_title(title)
    clean_axes(ax, x, y)
    apply_limits(ax, cells, x, y)
    save_figure(fig, base)
    return kept


def distribution_plot(embedding, base, cluster_column="seurat_clusters", group_column="group"):
    order = sorted(embedding[cluster_column].unique(), key=int)
    fig, ax = plt.subplots(figsize=(9, 4.8))
    sns.boxplot(data=embedding, x=cluster_column, y="delta_l2", hue=group_column, order=order,
                palette={"CD": "#c43c39", "Healthy": "#2878a8"}, showfliers=False, width=0.68, linewidth=0.9, ax=ax)
    ax.set(xlabel="Seurat cluster", ylabel="Predicted change magnitude (delta L2)",
           title="Predicted HIF1A knockout effect by disease group and cell cluster")
    ax.legend(title="Group", frameon=False)
    sns.despine(ax=ax)
    save_figure(fig, base)


def cluster_centers(oracle, cluster_column="seurat_clusters"):
    import pandas as pd
    embedding = pd.DataFrame({"UMAP1": oracle.embedding[:, 0], "UMAP2": oracle.embedding[:, 1],
                              cluster_column: oracle.adata.obs[cluster_column].astype(str).to_numpy()})
    return embedding.groupby(cluster_column)[["UMAP1", "UMAP2"]].median()


def plot_ps_panel(ax, oracle, grid, title, group=None, annotate=True, cluster_column="seurat_clusters", group_column="group"):
    embedding = np.asarray(oracle.embedding)
    obs = oracle.adata.obs
    background = np.ones(len(obs), dtype=bool) if group is None else obs[group_column].astype(str).to_numpy() == group
    ax.scatter(embedding[background, 0], embedding[background, 1], s=5, color="#e5e7eb", alpha=0.65, linewidths=0, rasterized=True)
    usable = grid[~grid["mass_filtered"]].copy()
    if group is not None:
        usable = usable[usable[group_column] == group]
    vmax = float(np.nanquantile(np.abs(grid.loc[~grid["mass_filtered"], "perturbation_score"]), 0.98))
    points = ax.scatter(usable["UMAP1"], usable["UMAP2"], c=usable["perturbation_score"], s=18, marker="s",
                        cmap=PS_CMAP, vmin=-vmax, vmax=vmax, linewidths=0, rasterized=True)
    ax.set_title(title, fontsize=11)
    clean_embedding_axis(ax)
    if annotate and group is None:
        negative, positive = strongest_ps_clusters(cluster_ps_summary(grid, cluster_column), cluster_column)
        centers = cluster_centers(oracle, cluster_column)
        for cluster, offset, label, color in zip(
            [negative, positive], [(-45, 35), (42, -38)],
            [f"Progression inhibited\ncluster {negative}", f"Progression promoted\ncluster {positive}"],
            ["#8e1b8e", "#15803d"],
        ):
            center = centers.loc[cluster]
            ax.annotate(label, xy=(center["UMAP1"], center["UMAP2"]), xytext=offset, textcoords="offset points",
                        color=color, fontsize=8.5, weight="bold",
                        arrowprops={"arrowstyle": "->", "color": color, "lw": 1.2, "connectionstyle": "arc3,rad=0.15"})
    return points


def plot_markov_panel(ax, density, results, title, density_column, annotate=True, cluster_column="seurat_clusters"):
    ax.scatter(density["UMAP1"], density["UMAP2"], s=5, color="#e5e7eb", alpha=0.55, linewidths=0, rasterized=True)
    ordered = density.iloc[np.argsort(density[density_column].to_numpy())]
    points = ax.scatter(ordered["UMAP1"], ordered["UMAP2"], c=ordered[density_column], s=8,
                        cmap=MARKOV_CMAP, vmin=0, vmax=1, linewidths=0, rasterized=True)
    ax.set_title(title, fontsize=11)
    clean_embedding_axis(ax)
    if annotate:
        cluster, _ = most_enriched_cluster(results)
        center = density[density[cluster_column] == cluster][["UMAP1", "UMAP2"]].median()
        ax.annotate(f"Predicted accumulation\ncluster {cluster}", xy=(center["UMAP1"], center["UMAP2"]),
                    xytext=(-105, 45), textcoords="offset points", color="#08519c", fontsize=8.5, weight="bold",
                    arrowprops={"arrowstyle": "->", "color": "#08519c", "lw": 1.2, "connectionstyle": "arc3,rad=-0.15"})
    return points


def add_colorbar(fig, points, ax, label, ticks=None, ticklabels=None):
    colorbar = fig.colorbar(points, ax=ax, fraction=0.045, pad=0.025)
    colorbar.set_label(label, fontsize=8.5)
    colorbar.ax.tick_params(labelsize=7.5, length=2)
    if ticks is not None:
        colorbar.set_ticks(ticks)
    if ticklabels is not None:
        colorbar.set_ticklabels(ticklabels)


def plot_single_panels(oracle, grid, density, results, figures):
    fig, ax = plt.subplots(figsize=(6.2, 5.4))
    points = plot_ps_panel(ax, oracle, grid, "HIF1A KO perturbation score along the observed trajectory")
    add_colorbar(fig, points, ax, "Perturbation score\n(negative: inhibition; positive: promotion)")
    save_figure(fig, figures / "HIF1A_KO_perturbation_score_PS")
    fig, ax = plt.subplots(figsize=(6.2, 5.4))
    points = plot_markov_panel(ax, density, results, "Predicted cell-state density after HIF1A KO Markov simulation", "markov_density_all")
    add_colorbar(fig, points, ax, "Simulated endpoint density", [0, 1], ["Low", "High"])
    save_figure(fig, figures / "HIF1A_KO_markov_simulation_density")


def plot_combined(oracle, grid, density, results, figures, markov_steps):
    fig, axes = plt.subplots(1, 2, figsize=(11.8, 5.1))
    ps_points = plot_ps_panel(axes[0], oracle, grid, "Perturbation score")
    markov_points = plot_markov_panel(axes[1], density, results, f"Markov simulation ({markov_steps} steps)", "markov_density_all")
    add_colorbar(fig, ps_points, axes[0], "PS")
    add_colorbar(fig, markov_points, axes[1], "Endpoint density", [0, 1], ["Low", "High"])
    fig.suptitle("Predicted effect of HIF1A knockout on cell-state progression", fontsize=14, y=1.01)
    fig.tight_layout()
    save_figure(fig, figures / "HIF1A_KO_PS_and_markov_combined")


def plot_groups(oracle, grid, density, results, figures, markov_steps):
    fig, axes = plt.subplots(2, 2, figsize=(10.8, 9.2))
    for row, group in enumerate(("CD", "Healthy")):
        ps_points = plot_ps_panel(axes[row, 0], oracle, grid, f"{group}: perturbation score", group=group, annotate=False)
        group_results = results[results["start_group"] == group]
        markov_points = plot_markov_panel(axes[row, 1], density, group_results,
                                          f"{group}: Markov endpoint density ({markov_steps} steps)",
                                          f"markov_density_start_{group}", annotate=False)
        add_colorbar(fig, ps_points, axes[row, 0], "PS")
        add_colorbar(fig, markov_points, axes[row, 1], "Endpoint density", [0, 1], ["Low", "High"])
    fig.suptitle("HIF1A knockout simulation stratified by disease group", fontsize=14, y=1.0)
    fig.tight_layout()
    save_figure(fig, figures / "HIF1A_KO_PS_and_markov_by_group")


def plot_ps_clusters4_5(grid, cells, output_base, cluster_column="seurat_clusters"):
    usable = grid[~grid["mass_filtered"]].copy()
    centers = cells.groupby(cluster_column)[["UMAP1", "UMAP2"]].median()
    fig, ax = plt.subplots(figsize=(6.2, 5.4))
    ax.scatter(cells["UMAP1"], cells["UMAP2"], s=5, color="#e5e7eb", alpha=0.65, linewidths=0, rasterized=True)
    vmax = float(np.nanquantile(np.abs(usable["perturbation_score"]), 0.98))
    points = ax.scatter(usable["UMAP1"], usable["UMAP2"], c=usable["perturbation_score"], s=18, marker="s",
                        cmap=PS_CLUSTERS45_CMAP, vmin=-vmax, vmax=vmax, linewidths=0, rasterized=True)
    rose = "#a52391"
    center4, center5 = centers.loc["4"], centers.loc["5"]
    anchor = (center4 + center5) / 2
    position = (anchor["UMAP1"] - 0.9, anchor["UMAP2"] + 1.65)
    ax.annotate("Progression inhibited\nclusters 4 and 5", xy=(center5["UMAP1"], center5["UMAP2"]), xytext=position,
                textcoords="data", color=rose, fontsize=8.5, weight="bold", ha="center",
                arrowprops={"arrowstyle": "->", "color": rose, "lw": 1.2, "connectionstyle": "arc3,rad=0.12"})
    ax.annotate("", xy=(center4["UMAP1"], center4["UMAP2"]), xytext=position, textcoords="data",
                arrowprops={"arrowstyle": "->", "color": rose, "lw": 1.2, "connectionstyle": "arc3,rad=-0.16"})
    center3 = centers.loc["3"]
    ax.annotate("Progression promoted\ncluster 3", xy=(center3["UMAP1"], center3["UMAP2"]), xytext=(-42, -42),
                textcoords="offset points", color="#15803d", fontsize=8.5, weight="bold",
                arrowprops={"arrowstyle": "->", "color": "#15803d", "lw": 1.2, "connectionstyle": "arc3,rad=0.12"})
    add_colorbar(fig, points, ax, "Perturbation score\n(negative: inhibition; positive: promotion)")
    ax.set_title("HIF1A KO perturbation score along the observed trajectory", fontsize=11)
    clean_embedding_axis(ax)
    save_figure(fig, output_base)


def save_primary_knockout_figures(oracle, target_gene, out_figures):
    import scanpy as sc

    sc.pl.embedding(oracle.adata, basis="umap", color=[oracle.cluster_column_name, target_gene], show=False, frameon=False, ncols=2)
    plt.tight_layout()
    plt.savefig(out_figures / "umap_cluster_and_target.png", dpi=300)
    plt.close()
    fig, ax = plt.subplots(figsize=(6, 5))
    oracle.plot_quiver(ax=ax, scale=30, s=6)
    ax.set_title(f"{target_gene} KO simulated cell shift")
    fig.tight_layout(); fig.savefig(out_figures / "hif1a_ko_quiver_cells.png", dpi=300); plt.close(fig)
    mass_filter = oracle.total_p_mass >= 0.01
    fig, ax = plt.subplots(figsize=(6, 5))
    ax.scatter(oracle.embedding[:, 0], oracle.embedding[:, 1], c="lightgray", s=5, alpha=0.25, linewidths=0)
    ax.quiver(oracle.flow_grid[mass_filter, 0], oracle.flow_grid[mass_filter, 1], oracle.flow[mass_filter, 0], oracle.flow[mass_filter, 1],
              angles="xy", scale_units="xy", scale=1, width=0.004, color="#1f2937")
    ax.set_title(f"{target_gene} KO simulated flow on grid"); ax.axis("off")
    plt.tight_layout(); fig.savefig(out_figures / "hif1a_ko_grid_arrows.png", dpi=300); plt.close(fig)
    try:
        propagation = celloracle_impact_summary(
            oracle, target_gene, max_propagation=5, norm_order=1
        )
        plot_propagation_magnitude(
            propagation,
            "cluster",
            out_figures / "n_propagation_impact",
            target_cluster="4",
            title=f"{target_gene} KO propagation magnitude by cluster",
            value_column="mean_delta_x_length",
            y_label="Mean delta X length",
        )
    except Exception as exc:
        (out_figures / "n_propagation_impact_warning.txt").write_text(str(exc))
