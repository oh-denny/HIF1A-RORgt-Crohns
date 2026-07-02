#!/usr/bin/env python3

import argparse
import json
import os
import warnings
from datetime import datetime
from pathlib import Path

PROJECT_DIR = Path(__file__).resolve().parents[2]
for name, relative in {
    "HOME": "results/celloracle_home",
    "XDG_CONFIG_HOME": "results/xdg_config",
    "XDG_CACHE_HOME": "results/xdg_cache",
    "MPLCONFIGDIR": "results/matplotlib_cache",
    "NUMBA_CACHE_DIR": "results/celloracle_numba_cache",
}.items():
    os.environ[name] = str(PROJECT_DIR / relative)
    Path(os.environ[name]).mkdir(parents=True, exist_ok=True)

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist
from sklearn.neighbors import KernelDensity, NearestNeighbors

try:
    import gimmemotifs.motif as gimmemotifs_motif

    if not hasattr(gimmemotifs_motif, "default_motifs"):
        gimmemotifs_motif.default_motifs = lambda *args, **kwargs: []
except Exception:
    pass

import celloracle as co
from celloracle.applications import Gradient_calculator, Oracle_development_module


SEED = 123
CLUSTER_COLUMN = "seurat_clusters"
GROUP_COLUMN = "group"
PS_CMAP = LinearSegmentedColormap.from_list(
    "ps_purple_white_green", ["#8e1b8e", "#ffffff", "#239b56"], N=256
)
MARKOV_CMAP = LinearSegmentedColormap.from_list(
    "markov_blue", ["#deedf7", "#6baed6", "#08519c"], N=256
)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Generate CellOracle perturbation-score and Markov panels for HIF1A KO."
    )
    parser.add_argument(
        "--oracle",
        default="results/celloracle_hif1a_ko_1500/objects/oracle_HIF1A_KO.celloracle.oracle",
    )
    parser.add_argument(
        "--out-prefix",
        default="results/celloracle_hif1a_ko_1500_PS_markov",
    )
    parser.add_argument("--pseudotime-column", default="Pseudotime")
    parser.add_argument("--n-grid", type=int, default=40)
    parser.add_argument("--grid-neighbors", type=int, default=200)
    parser.add_argument("--min-mass", type=float, default=0.01)
    parser.add_argument("--markov-steps", type=int, default=200)
    parser.add_argument("--markov-duplications", type=int, default=5)
    parser.add_argument("--n-jobs", type=int, default=1)
    return parser.parse_args()


def make_output_dir(prefix):
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    out = Path(f"{prefix}_{stamp}")
    if out.exists():
        raise FileExistsError(f"Refusing to overwrite existing output: {out}")
    for child in ("figures", "tables", "objects", "logs"):
        (out / child).mkdir(parents=True, exist_ok=False)
    return out


def savefig(fig, base):
    fig.savefig(Path(base).with_suffix(".png"), dpi=400, bbox_inches="tight")
    fig.savefig(Path(base).with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def clean_embedding_axis(ax):
    ax.set_aspect("equal", adjustable="datalim")
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)


def load_oracle(path):
    oracle = co.load_hdf5(str(path))
    required_obs = {CLUSTER_COLUMN, GROUP_COLUMN}
    missing_obs = required_obs.difference(oracle.adata.obs.columns)
    if missing_obs:
        raise KeyError(f"Missing required cell metadata: {sorted(missing_obs)}")
    if oracle.embedding_name != "X_umap":
        if "X_umap" not in oracle.adata.obsm:
            raise KeyError("X_umap is not available in the Oracle object.")
        raise ValueError(
            "The input Oracle must contain the HIF1A KO simulation calculated on X_umap."
        )
    required_attrs = (
        "delta_embedding",
        "delta_embedding_random",
        "transition_prob",
        "transition_prob_random",
    )
    missing_attrs = [name for name in required_attrs if not hasattr(oracle, name)]
    if missing_attrs:
        raise AttributeError(f"Oracle is missing simulation attributes: {missing_attrs}")
    return oracle


def calculate_celloracle_ps(oracle, args):
    gradient = Gradient_calculator(
        oracle_object=oracle,
        pseudotime_key=args.pseudotime_column,
        name="HIF1A_KO",
    )
    gradient.calculate_p_mass(
        smooth=0.8,
        n_grid=args.n_grid,
        n_neighbors=args.grid_neighbors,
        n_jobs=args.n_jobs,
    )
    gradient.calculate_mass_filter(min_mass=args.min_mass)
    gradient.transfer_data_into_grid(args={"method": "knn", "n_knn": 30})
    gradient.calculate_gradient(scale_factor="l2_norm_mean", normalization="sqrt")

    development = Oracle_development_module(
        gradient_object=gradient,
        oracle_object=oracle,
        name="HIF1A_KO",
    )
    development.calculate_inner_product()

    grid = pd.DataFrame(
        {
            "UMAP1": development.gridpoints_coordinates[:, 0],
            "UMAP2": development.gridpoints_coordinates[:, 1],
            "perturbation_score": development.inner_product,
            "perturbation_score_randomized": development.inner_product_random,
            "pseudotime": development.pseudotime_on_grid,
            "simulation_flow_UMAP1": development.flow[:, 0],
            "simulation_flow_UMAP2": development.flow[:, 1],
            "development_flow_UMAP1": development.ref_flow[:, 0],
            "development_flow_UMAP2": development.ref_flow[:, 1],
            "mass_filtered": development.mass_filter_simulation,
        }
    )
    return gradient, development, grid


def assign_grid_metadata(grid, oracle):
    embedding = np.asarray(oracle.embedding)
    nearest = NearestNeighbors(n_neighbors=1).fit(embedding)
    _, indices = nearest.kneighbors(grid[["UMAP1", "UMAP2"]].to_numpy())
    indices = indices[:, 0]
    obs = oracle.adata.obs
    grid[CLUSTER_COLUMN] = obs.iloc[indices][CLUSTER_COLUMN].astype(str).to_numpy()
    grid[GROUP_COLUMN] = obs.iloc[indices][GROUP_COLUMN].astype(str).to_numpy()
    return grid


def run_markov(oracle, args):
    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore",
            message="Functions for Markov simulation are deprecated.*",
            category=DeprecationWarning,
        )
        oracle.run_markov_chain_simulation(
            n_steps=args.markov_steps,
            n_duplication=args.markov_duplications,
            seed=SEED,
            calculate_randomized=True,
        )

    trajectories = oracle.markvov_transition_id
    random_trajectories = oracle.markvov_transition_random_id
    start = trajectories.iloc[:, 0].to_numpy(dtype=int)
    end = trajectories.iloc[:, -1].to_numpy(dtype=int)
    end_random = random_trajectories.iloc[:, -1].to_numpy(dtype=int)
    obs = oracle.adata.obs.reset_index(names="cell")
    result = pd.DataFrame(
        {
            "start_cell_index": start,
            "end_cell_index": end,
            "end_randomized_cell_index": end_random,
            "start_cell": obs.iloc[start]["cell"].to_numpy(),
            "end_cell": obs.iloc[end]["cell"].to_numpy(),
            "end_randomized_cell": obs.iloc[end_random]["cell"].to_numpy(),
            "start_cluster": obs.iloc[start][CLUSTER_COLUMN].astype(str).to_numpy(),
            "end_cluster": obs.iloc[end][CLUSTER_COLUMN].astype(str).to_numpy(),
            "end_randomized_cluster": obs.iloc[end_random][CLUSTER_COLUMN].astype(str).to_numpy(),
            "start_group": obs.iloc[start][GROUP_COLUMN].astype(str).to_numpy(),
            "end_group": obs.iloc[end][GROUP_COLUMN].astype(str).to_numpy(),
        }
    )
    return result, trajectories.to_numpy(dtype=np.int32), random_trajectories.to_numpy(
        dtype=np.int32
    )


def kde_bandwidth(embedding):
    distances, _ = NearestNeighbors(n_neighbors=min(31, len(embedding))).fit(
        embedding
    ).kneighbors()
    bandwidth = float(np.median(distances[:, -1]))
    return max(bandwidth, np.finfo(float).eps)


def normalized_density(embedding, endpoint_indices, bandwidth):
    endpoint_xy = embedding[np.asarray(endpoint_indices, dtype=int)]
    model = KernelDensity(kernel="gaussian", bandwidth=bandwidth).fit(endpoint_xy)
    log_density = model.score_samples(embedding)
    density = np.exp(log_density - np.nanmax(log_density))
    low, high = np.nanquantile(density, [0.01, 0.99])
    return np.clip((density - low) / max(high - low, 1e-12), 0, 1)


def make_markov_density_table(oracle, markov_results):
    embedding = np.asarray(oracle.embedding)
    bandwidth = kde_bandwidth(embedding)
    obs = oracle.adata.obs.reset_index(names="cell")
    out = pd.DataFrame(
        {
            "cell": obs["cell"],
            "UMAP1": embedding[:, 0],
            "UMAP2": embedding[:, 1],
            CLUSTER_COLUMN: obs[CLUSTER_COLUMN].astype(str),
            GROUP_COLUMN: obs[GROUP_COLUMN].astype(str),
        }
    )
    out["markov_density_all"] = normalized_density(
        embedding, markov_results["end_cell_index"], bandwidth
    )
    out["markov_density_randomized"] = normalized_density(
        embedding, markov_results["end_randomized_cell_index"], bandwidth
    )
    for group in ("CD", "Healthy"):
        group_end = markov_results.loc[
            markov_results["start_group"] == group, "end_cell_index"
        ]
        out[f"markov_density_start_{group}"] = normalized_density(
            embedding, group_end, bandwidth
        )
    return out, bandwidth


def cluster_centers(oracle):
    embedding = pd.DataFrame(
        {
            "UMAP1": oracle.embedding[:, 0],
            "UMAP2": oracle.embedding[:, 1],
            CLUSTER_COLUMN: oracle.adata.obs[CLUSTER_COLUMN].astype(str).to_numpy(),
        }
    )
    return embedding.groupby(CLUSTER_COLUMN)[["UMAP1", "UMAP2"]].median()


def cluster_ps_summary(grid):
    usable = grid[~grid["mass_filtered"]].copy()
    return (
        usable.groupby(CLUSTER_COLUMN)["perturbation_score"]
        .agg(n_grid_points="size", mean_PS="mean", median_PS="median")
        .reset_index()
        .sort_values("mean_PS")
    )


def transition_tables(markov_results):
    counts = pd.crosstab(
        markov_results["start_cluster"],
        markov_results["end_cluster"],
        rownames=["start_cluster"],
        colnames=["end_cluster"],
    )
    proportions = counts.div(counts.sum(axis=1), axis=0)
    group_cluster = (
        markov_results.groupby(["start_group", "end_cluster"])
        .size()
        .rename("n_trajectories")
        .reset_index()
    )
    group_cluster["proportion_within_start_group"] = group_cluster.groupby(
        "start_group"
    )["n_trajectories"].transform(lambda values: values / values.sum())
    return counts, proportions, group_cluster


def strongest_ps_clusters(ps_summary):
    negative = ps_summary.iloc[0]
    positive = ps_summary.iloc[-1]
    return str(negative[CLUSTER_COLUMN]), str(positive[CLUSTER_COLUMN])


def most_enriched_cluster(markov_results):
    start = markov_results["start_cluster"].value_counts(normalize=True)
    end = markov_results["end_cluster"].value_counts(normalize=True)
    enrichment = (end / start).replace([np.inf, -np.inf], np.nan).dropna()
    return str(enrichment.idxmax()), enrichment


def plot_ps_panel(ax, oracle, grid, title, group=None, annotate=True):
    embedding = np.asarray(oracle.embedding)
    obs = oracle.adata.obs
    background = np.ones(len(obs), dtype=bool)
    if group is not None:
        background = obs[GROUP_COLUMN].astype(str).to_numpy() == group
    ax.scatter(
        embedding[background, 0],
        embedding[background, 1],
        s=5,
        color="#e5e7eb",
        alpha=0.65,
        linewidths=0,
        rasterized=True,
    )
    usable = grid[~grid["mass_filtered"]].copy()
    if group is not None:
        usable = usable[usable[GROUP_COLUMN] == group]
    vmax = float(np.nanquantile(np.abs(grid.loc[~grid["mass_filtered"], "perturbation_score"]), 0.98))
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
    ax.set_title(title, fontsize=11)
    clean_embedding_axis(ax)

    if annotate and group is None:
        summary = cluster_ps_summary(grid)
        negative, positive = strongest_ps_clusters(summary)
        centers = cluster_centers(oracle)
        offsets = [(-45, 35), (42, -38)]
        labels = [
            f"Progression inhibited\ncluster {negative}",
            f"Progression promoted\ncluster {positive}",
        ]
        colors = ["#8e1b8e", "#15803d"]
        for cluster, offset, label, color in zip(
            [negative, positive], offsets, labels, colors
        ):
            center = centers.loc[cluster]
            ax.annotate(
                label,
                xy=(center["UMAP1"], center["UMAP2"]),
                xytext=offset,
                textcoords="offset points",
                color=color,
                fontsize=8.5,
                weight="bold",
                arrowprops={
                    "arrowstyle": "->",
                    "color": color,
                    "lw": 1.2,
                    "connectionstyle": "arc3,rad=0.15",
                },
            )
    return points


def plot_markov_panel(ax, density_table, markov_results, title, density_column, annotate=True):
    ax.scatter(
        density_table["UMAP1"],
        density_table["UMAP2"],
        s=5,
        color="#e5e7eb",
        alpha=0.55,
        linewidths=0,
        rasterized=True,
    )
    order = np.argsort(density_table[density_column].to_numpy())
    ordered = density_table.iloc[order]
    points = ax.scatter(
        ordered["UMAP1"],
        ordered["UMAP2"],
        c=ordered[density_column],
        s=8,
        cmap=MARKOV_CMAP,
        vmin=0,
        vmax=1,
        linewidths=0,
        rasterized=True,
    )
    ax.set_title(title, fontsize=11)
    clean_embedding_axis(ax)
    if annotate:
        cluster, _ = most_enriched_cluster(markov_results)
        center = (
            density_table[density_table[CLUSTER_COLUMN] == cluster][["UMAP1", "UMAP2"]]
            .median()
        )
        ax.annotate(
            f"Predicted accumulation\ncluster {cluster}",
            xy=(center["UMAP1"], center["UMAP2"]),
            xytext=(-105, 45),
            textcoords="offset points",
            color="#08519c",
            fontsize=8.5,
            weight="bold",
            arrowprops={
                "arrowstyle": "->",
                "color": "#08519c",
                "lw": 1.2,
                "connectionstyle": "arc3,rad=-0.15",
            },
        )
    return points


def add_colorbar(fig, points, ax, label, ticks=None, ticklabels=None):
    colorbar = fig.colorbar(points, ax=ax, fraction=0.045, pad=0.025)
    colorbar.set_label(label, fontsize=8.5)
    colorbar.ax.tick_params(labelsize=7.5, length=2)
    if ticks is not None:
        colorbar.set_ticks(ticks)
    if ticklabels is not None:
        colorbar.set_ticklabels(ticklabels)


def plot_single_panels(oracle, grid, density, markov_results, figures):
    fig, ax = plt.subplots(figsize=(6.2, 5.4))
    points = plot_ps_panel(
        ax,
        oracle,
        grid,
        "HIF1A KO perturbation score along the observed trajectory",
    )
    add_colorbar(fig, points, ax, "Perturbation score\n(negative: inhibition; positive: promotion)")
    savefig(fig, figures / "HIF1A_KO_perturbation_score_PS")

    fig, ax = plt.subplots(figsize=(6.2, 5.4))
    points = plot_markov_panel(
        ax,
        density,
        markov_results,
        "Predicted cell-state density after HIF1A KO Markov simulation",
        "markov_density_all",
    )
    add_colorbar(fig, points, ax, "Simulated endpoint density", [0, 1], ["Low", "High"])
    savefig(fig, figures / "HIF1A_KO_markov_simulation_density")


def plot_combined(oracle, grid, density, markov_results, figures, markov_steps):
    fig, axes = plt.subplots(1, 2, figsize=(11.8, 5.1))
    ps_points = plot_ps_panel(
        axes[0],
        oracle,
        grid,
        "Perturbation score",
    )
    markov_points = plot_markov_panel(
        axes[1],
        density,
        markov_results,
        f"Markov simulation ({markov_steps} steps)",
        "markov_density_all",
    )
    add_colorbar(fig, ps_points, axes[0], "PS")
    add_colorbar(
        fig,
        markov_points,
        axes[1],
        "Endpoint density",
        [0, 1],
        ["Low", "High"],
    )
    fig.suptitle(
        "Predicted effect of HIF1A knockout on cell-state progression",
        fontsize=14,
        y=1.01,
    )
    fig.tight_layout()
    savefig(fig, figures / "HIF1A_KO_PS_and_markov_combined")


def plot_groups(oracle, grid, density, markov_results, figures, markov_steps):
    fig, axes = plt.subplots(2, 2, figsize=(10.8, 9.2))
    for row, group in enumerate(("CD", "Healthy")):
        ps_points = plot_ps_panel(
            axes[row, 0],
            oracle,
            grid,
            f"{group}: perturbation score",
            group=group,
            annotate=False,
        )
        group_results = markov_results[markov_results["start_group"] == group]
        markov_points = plot_markov_panel(
            axes[row, 1],
            density,
            group_results,
            f"{group}: Markov endpoint density ({markov_steps} steps)",
            f"markov_density_start_{group}",
            annotate=False,
        )
        add_colorbar(fig, ps_points, axes[row, 0], "PS")
        add_colorbar(
            fig,
            markov_points,
            axes[row, 1],
            "Endpoint density",
            [0, 1],
            ["Low", "High"],
        )
    fig.suptitle(
        "HIF1A knockout simulation stratified by disease group",
        fontsize=14,
        y=1.0,
    )
    fig.tight_layout()
    savefig(fig, figures / "HIF1A_KO_PS_and_markov_by_group")


def main():
    args = parse_args()
    out = make_output_dir(args.out_prefix)
    figures = out / "figures"
    tables = out / "tables"
    objects = out / "objects"
    logs = out / "logs"

    oracle = load_oracle(args.oracle)
    gradient, development, ps_grid = calculate_celloracle_ps(oracle, args)
    ps_grid = assign_grid_metadata(ps_grid, oracle)
    ps_summary = cluster_ps_summary(ps_grid)

    markov_results, trajectories, random_trajectories = run_markov(oracle, args)
    density, bandwidth = make_markov_density_table(oracle, markov_results)
    transition_counts, transition_proportions, group_cluster = transition_tables(
        markov_results
    )

    ps_grid.to_csv(tables / "HIF1A_KO_perturbation_score_grid.csv", index=False)
    ps_summary.to_csv(tables / "HIF1A_KO_perturbation_score_by_cluster.csv", index=False)
    markov_results.to_csv(tables / "HIF1A_KO_markov_start_end_cells.csv", index=False)
    density.to_csv(tables / "HIF1A_KO_markov_density_by_cell.csv", index=False)
    transition_counts.to_csv(tables / "HIF1A_KO_markov_transition_counts.csv")
    transition_proportions.to_csv(
        tables / "HIF1A_KO_markov_transition_proportions.csv"
    )
    group_cluster.to_csv(
        tables / "HIF1A_KO_markov_endpoint_cluster_by_start_group.csv", index=False
    )
    np.savez_compressed(
        objects / "HIF1A_KO_markov_trajectories.npz",
        trajectories=trajectories,
        randomized_trajectories=random_trajectories,
        selected_cell_indices=np.asarray(oracle.ixs_markvov_simulation),
    )
    gradient.to_hdf5(str(objects / "HIF1A_KO_pseudotime_gradient.celloracle.gradient"))

    plot_single_panels(oracle, ps_grid, density, markov_results, figures)
    plot_combined(
        oracle, ps_grid, density, markov_results, figures, args.markov_steps
    )
    plot_groups(
        oracle, ps_grid, density, markov_results, figures, args.markov_steps
    )

    negative_cluster, positive_cluster = strongest_ps_clusters(ps_summary)
    enriched_cluster, enrichment = most_enriched_cluster(markov_results)
    summary = {
        "output_dir": str(out.resolve()),
        "source_oracle": str(Path(args.oracle).resolve()),
        "embedding": oracle.embedding_name,
        "n_cells": int(oracle.adata.n_obs),
        "n_markov_start_cells": int(len(oracle.ixs_markvov_simulation)),
        "n_markov_trajectories": int(len(markov_results)),
        "markov_steps": args.markov_steps,
        "markov_duplications": args.markov_duplications,
        "kde_bandwidth": bandwidth,
        "cluster_with_most_negative_mean_PS": negative_cluster,
        "cluster_with_most_positive_mean_PS": positive_cluster,
        "cluster_with_largest_markov_endpoint_enrichment": enriched_cluster,
        "largest_markov_endpoint_enrichment_ratio": float(
            enrichment.loc[enriched_cluster]
        ),
        "groups": sorted(oracle.adata.obs[GROUP_COLUMN].astype(str).unique()),
        "figures": sorted(path.name for path in figures.glob("*")),
        "tables": sorted(path.name for path in tables.glob("*")),
        "note": (
            "CellOracle marks Markov simulation as deprecated in favor of "
            "perturbation-score analysis; both are exported here to reproduce "
            "the requested paper-style comparison."
        ),
    }
    (logs / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    (logs / "run_parameters.json").write_text(json.dumps(vars(args), indent=2) + "\n")
    (logs / "methodology.txt").write_text(
        "Perturbation score (PS) was calculated with CellOracle's "
        "Gradient_calculator and Oracle_development_module. The score is the "
        "inner product between the HIF1A KO simulation flow and the pseudotime "
        "gradient on a 40 x 40 UMAP grid. Negative values indicate opposition "
        "to pseudotime progression; positive values indicate alignment.\n\n"
        "Markov trajectories were generated with Oracle.run_markov_chain_simulation "
        "from the existing HIF1A KO transition probabilities. The plotted endpoint "
        "density is a Gaussian KDE evaluated at observed cells and scaled from 0 "
        "to 1 after clipping at the 1st and 99th percentiles. Group-specific panels "
        "condition trajectories by the group of their starting cell.\n"
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
