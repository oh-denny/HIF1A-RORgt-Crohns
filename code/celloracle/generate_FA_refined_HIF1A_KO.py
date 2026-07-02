#!/usr/bin/env python3

import argparse
import json
import os
from datetime import datetime
from pathlib import Path

PROJECT_DIR = Path(__file__).resolve().parents[2]
for name, relative in {
    "MPLCONFIGDIR": "results/matplotlib_cache",
    "XDG_CACHE_HOME": "results/xdg_cache",
}.items():
    os.environ.setdefault(name, str(PROJECT_DIR / relative))
    Path(os.environ[name]).mkdir(parents=True, exist_ok=True)

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from fa2_modified import ForceAtlas2
from scipy import sparse
from sklearn.neighbors import NearestNeighbors
from sklearn.preprocessing import RobustScaler


CLUSTER_COL = "seurat_clusters"
GROUP_COL = "group"
GROUPS = ("CD", "Healthy")
SEED = 123


def parse_args():
    parser = argparse.ArgumentParser(
        description="Refined, visualization-only analysis of the CellOracle HIF1A KO results."
    )
    parser.add_argument(
        "--source-dir",
        default="results/celloracle_hif1a_ko_1500_FA_paperstyle_20260626_154852",
    )
    parser.add_argument(
        "--out-prefix",
        default="results/celloracle_hif1a_ko_1500_FA_refined",
    )
    parser.add_argument("--target-cluster", default="4")
    parser.add_argument("--related-clusters", type=int, default=3)
    parser.add_argument("--knn", type=int, default=30)
    parser.add_argument("--arrow-fraction", type=float, default=0.075)
    parser.add_argument("--fa-iterations", type=int, default=800)
    return parser.parse_args()


def make_output_dir(prefix):
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    out = Path(f"{prefix}_{stamp}")
    if out.exists():
        raise FileExistsError(f"Refusing to overwrite existing output: {out}")
    for child in ("figures", "tables", "logs"):
        (out / child).mkdir(parents=True, exist_ok=False)
    return out


def read_inputs(source):
    tables = Path(source) / "tables"
    names = {
        "embedding": "fa_embedding_cells.csv",
        "vectors": "cell_shift_vectors_FA.csv",
        "grid": "grid_shift_vectors_FA.csv",
        "propagation": "propagation_magnitude_by_cluster.csv",
        "score": "perturbation_score_FA.csv",
    }
    data = {}
    for key, name in names.items():
        path = tables / name
        if not path.exists():
            raise FileNotFoundError(path)
        data[key] = pd.read_csv(
            path,
            dtype={CLUSTER_COL: str, GROUP_COL: str},
        )
    return data


def savefig(fig, base):
    fig.savefig(Path(base).with_suffix(".png"), dpi=400, bbox_inches="tight")
    fig.savefig(Path(base).with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def clean_axes(ax, x="FA1", y="FA2"):
    ax.set_xlabel(x)
    ax.set_ylabel(y)
    ax.set_aspect("equal", adjustable="datalim")
    ax.grid(False)
    sns.despine(ax=ax)


def robust_limits(values, padding=0.04):
    lo, hi = np.nanquantile(values, [0.002, 0.998])
    span = max(hi - lo, np.finfo(float).eps)
    return lo - span * padding, hi + span * padding


def apply_limits(ax, df, x, y):
    ax.set_xlim(*robust_limits(df[x]))
    ax.set_ylim(*robust_limits(df[y]))


def arrow_scale(xy, uv, fraction=0.035):
    diag = np.linalg.norm(np.ptp(np.asarray(xy), axis=0))
    mag = np.linalg.norm(np.asarray(uv), axis=1)
    positive = mag[np.isfinite(mag) & (mag > 0)]
    ref = np.quantile(positive, 0.95) if len(positive) else 1.0
    return diag * fraction / ref


def stratified_sample(df, fraction):
    sampled = []
    for _, sub in df.groupby(CLUSTER_COL, sort=True):
        n = min(len(sub), max(1, int(round(len(sub) * fraction))))
        sampled.append(sub.sample(n=n, random_state=SEED))
    return pd.concat(sampled, ignore_index=True)


def cluster_connectivity(embedding, target, knn):
    coords = RobustScaler().fit_transform(embedding[["FA1", "FA2"]])
    labels = embedding[CLUSTER_COL].to_numpy()
    n_neighbors = min(knn + 1, len(embedding))
    indices = NearestNeighbors(n_neighbors=n_neighbors).fit(coords).kneighbors(
        return_distance=False
    )
    rows = []
    target_mask = labels == target
    target_center = np.median(coords[target_mask], axis=0)
    for cluster in sorted(set(labels), key=lambda value: int(value)):
        if cluster == target:
            continue
        cluster_mask = labels == cluster
        outgoing = sum(
            np.count_nonzero(labels[neighbors[1:]] == cluster)
            for neighbors in indices[target_mask]
        )
        incoming = sum(
            np.count_nonzero(labels[neighbors[1:]] == target)
            for neighbors in indices[cluster_mask]
        )
        rows.append(
            {
                "target_cluster": target,
                "related_cluster": cluster,
                "knn_edges_both_directions": int(outgoing + incoming),
                "edges_from_target": int(outgoing),
                "edges_to_target": int(incoming),
                "robust_centroid_distance": float(
                    np.linalg.norm(
                        target_center - np.median(coords[cluster_mask], axis=0)
                    )
                ),
            }
        )
    return pd.DataFrame(rows).sort_values(
        ["knn_edges_both_directions", "robust_centroid_distance"],
        ascending=[False, True],
    )


def knn_graph(coords, knn):
    n_neighbors = min(knn + 1, len(coords))
    distances, indices = NearestNeighbors(n_neighbors=n_neighbors).fit(coords).kneighbors()
    distances = distances[:, 1:]
    indices = indices[:, 1:]
    sigma = np.median(distances[distances > 0])
    weights = np.exp(-np.square(distances / max(sigma, np.finfo(float).eps)))
    rows = np.repeat(np.arange(len(coords)), indices.shape[1])
    graph = sparse.csr_matrix(
        (weights.ravel(), (rows, indices.ravel())),
        shape=(len(coords), len(coords)),
    )
    return graph.maximum(graph.T)


def refined_forceatlas(embedding, included, knn, iterations):
    subset = embedding[embedding[CLUSTER_COL].isin(included)].copy().reset_index(drop=True)
    initial = RobustScaler().fit_transform(subset[["FA1", "FA2"]])
    graph = knn_graph(initial, knn)
    fa = ForceAtlas2(
        outboundAttractionDistribution=False,
        linLogMode=False,
        adjustSizes=False,
        edgeWeightInfluence=1.0,
        jitterTolerance=1.0,
        barnesHutOptimize=True,
        barnesHutTheta=1.2,
        scalingRatio=2.0,
        strongGravityMode=True,
        gravity=1.0,
        verbose=False,
    )
    refined = np.asarray(fa.forceatlas2(graph, pos=initial, iterations=iterations))
    refined = RobustScaler().fit_transform(refined)
    subset["FA1_refined"] = refined[:, 0]
    subset["FA2_refined"] = refined[:, 1]
    return subset, graph


def project_vectors_to_refined(vectors, refined, knn=20):
    merged = refined[
        ["cell", "FA1_refined", "FA2_refined"]
    ].merge(vectors, on="cell", how="left", validate="one_to_one")
    old_xy = merged[["FA1", "FA2"]].to_numpy()
    new_xy = merged[["FA1_refined", "FA2_refined"]].to_numpy()
    old_scaled = RobustScaler().fit_transform(old_xy)
    neighbors = NearestNeighbors(n_neighbors=min(knn, len(merged))).fit(old_scaled)
    indices = neighbors.kneighbors(return_distance=False)
    projected = np.zeros((len(merged), 2))
    old_shift = merged[["shift_FA1", "shift_FA2"]].to_numpy()
    for i, local in enumerate(indices):
        design = np.column_stack([old_xy[local] - old_xy[i], np.ones(len(local))])
        response = new_xy[local] - new_xy[i]
        coefficients = np.linalg.lstsq(design, response, rcond=None)[0]
        projected[i] = old_shift[i] @ coefficients[:2, :]
    merged["shift_FA1_refined"] = projected[:, 0]
    merged["shift_FA2_refined"] = projected[:, 1]
    return merged


def estimate_grid(vectors, x, y, u, v, steps=32, neighbors=80):
    xy = vectors[[x, y]].to_numpy()
    uv = vectors[[u, v]].to_numpy()
    x_values = np.linspace(*robust_limits(xy[:, 0], padding=0), steps)
    y_values = np.linspace(*robust_limits(xy[:, 1], padding=0), steps)
    gx, gy = np.meshgrid(x_values, y_values)
    points = np.column_stack([gx.ravel(), gy.ravel()])
    model = NearestNeighbors(n_neighbors=min(neighbors, len(xy))).fit(xy)
    distances, indices = model.kneighbors(points)
    bandwidth = np.median(distances[:, -1])
    weights = np.exp(-0.5 * np.square(distances / max(bandwidth, 1e-12)))
    flow = (uv[indices] * weights[:, :, None]).sum(axis=1)
    flow /= np.maximum(weights.sum(axis=1)[:, None], 1e-12)
    radius = np.median(np.diff(x_values)) * 1.5
    counts = np.array(
        [
            len(found)
            for found in NearestNeighbors(radius=radius)
            .fit(xy)
            .radius_neighbors(points, return_distance=False)
        ]
    )
    return pd.DataFrame(
        {
            x.replace("FA", "grid_FA"): points[:, 0],
            y.replace("FA", "grid_FA"): points[:, 1],
            u.replace("shift_", "flow_"): flow[:, 0],
            v.replace("shift_", "flow_"): flow[:, 1],
            "n_cells": counts,
        }
    )


def delta_map(df, base, title, x="FA1", y="FA2", highlight=None):
    fig, ax = plt.subplots(figsize=(7.2, 6.2))
    ax.scatter(df[x], df[y], s=5, color="#d8dadd", alpha=0.55, linewidths=0)
    draw = df if highlight is None else df[df[CLUSTER_COL] == highlight]
    vmax = np.nanquantile(df["delta_l2"], 0.98)
    points = ax.scatter(
        draw[x],
        draw[y],
        c=draw["delta_l2"],
        s=9 if highlight is None else 13,
        cmap="magma",
        vmin=0,
        vmax=vmax,
        linewidths=0,
        rasterized=True,
    )
    if highlight is not None:
        center = draw[[x, y]].median()
        ax.text(
            center[x],
            center[y],
            f"cluster {highlight}",
            fontsize=9,
            weight="bold",
            ha="center",
            va="center",
            bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.7, "pad": 2},
        )
    fig.colorbar(points, ax=ax, label="Magnitude of predicted transcriptomic change (delta L2)")
    ax.set_title(title)
    clean_axes(ax, x, y)
    apply_limits(ax, df, x, y)
    savefig(fig, base)


def arrow_plot(df, sampled, base, title, x="FA1", y="FA2", u="shift_FA1", v="shift_FA2"):
    fig, ax = plt.subplots(figsize=(7.2, 6.2))
    ax.scatter(df[x], df[y], s=5, color="#d8dadd", alpha=0.55, linewidths=0)
    factor = arrow_scale(df[[x, y]].to_numpy(), df[[u, v]].to_numpy(), fraction=0.025)
    ax.quiver(
        sampled[x],
        sampled[y],
        sampled[u] * factor,
        sampled[v] * factor,
        angles="xy",
        scale_units="xy",
        scale=1,
        width=0.0022,
        headwidth=3.2,
        headlength=4.0,
        color="#111111",
        alpha=0.62,
    )
    ax.set_title(title)
    clean_axes(ax, x, y)
    apply_limits(ax, df, x, y)
    savefig(fig, base)


def flow_plot(
    cells,
    grid,
    base,
    title,
    x="FA1",
    y="FA2",
    gx="grid_FA1",
    gy="grid_FA2",
    gu="flow_FA1",
    gv="flow_FA2",
):
    kept = grid[grid["n_cells"] >= 5].copy()
    fig, ax = plt.subplots(figsize=(7.2, 6.2))
    ax.scatter(cells[x], cells[y], s=5, color="#d8dadd", alpha=0.55, linewidths=0)
    if len(kept):
        factor = arrow_scale(cells[[x, y]].to_numpy(), kept[[gu, gv]].to_numpy(), fraction=0.035)
        ax.quiver(
            kept[gx],
            kept[gy],
            kept[gu] * factor,
            kept[gv] * factor,
            angles="xy",
            scale_units="xy",
            scale=1,
            width=0.0025,
            headwidth=3.2,
            color="#111111",
            alpha=0.62,
        )
    ax.set_title(title)
    clean_axes(ax, x, y)
    apply_limits(ax, cells, x, y)
    savefig(fig, base)
    return kept


def distribution_plot(embedding, base):
    order = sorted(embedding[CLUSTER_COL].unique(), key=int)
    palette = {"CD": "#c43c39", "Healthy": "#2878a8"}
    fig, ax = plt.subplots(figsize=(9, 4.8))
    sns.boxplot(
        data=embedding,
        x=CLUSTER_COL,
        y="delta_l2",
        hue=GROUP_COL,
        order=order,
        palette=palette,
        showfliers=False,
        width=0.68,
        linewidth=0.9,
        ax=ax,
    )
    ax.set_xlabel("Seurat cluster")
    ax.set_ylabel("Predicted change magnitude (delta L2)")
    ax.set_title("Predicted HIF1A knockout effect by disease group and cell cluster")
    ax.legend(title="Group", frameon=False)
    sns.despine(ax=ax)
    savefig(fig, base)


def propagation_plot(propagation, included, target, base):
    fig, ax = plt.subplots(figsize=(8, 4.8))
    for cluster, sub in propagation.groupby(CLUSTER_COL):
        related = cluster in included
        ax.plot(
            sub["n_propagation"],
            sub["mean_l2_norm"],
            marker="o" if related else None,
            markersize=4,
            linewidth=2.2 if cluster == target else (1.5 if related else 0.8),
            color="#c43c39" if cluster == target else ("#2878a8" if related else "#b7babd"),
            alpha=1 if related else 0.6,
            label=f"cluster {cluster}",
        )
    ax.set_xlabel("CellOracle propagation step")
    ax.set_ylabel("Mean L2 norm")
    ax.set_title("Propagation of the predicted HIF1A knockout response across clusters")
    ax.legend(frameon=False, ncol=2, bbox_to_anchor=(1.02, 1), loc="upper left")
    sns.despine(ax=ax)
    savefig(fig, base)


def group_summary(embedding):
    return (
        embedding.groupby([CLUSTER_COL, GROUP_COL])["delta_l2"]
        .agg(n_cells="size", mean_delta_l2="mean", median_delta_l2="median", sd_delta_l2="std")
        .reset_index()
    )


def main():
    args = parse_args()
    source = Path(args.source_dir)
    out = make_output_dir(args.out_prefix)
    figures = out / "figures"
    tables = out / "tables"
    logs = out / "logs"
    data = read_inputs(source)

    embedding = data["embedding"].copy()
    vectors = data["vectors"].copy()
    grid = data["grid"].copy()
    propagation = data["propagation"].copy()
    score = data["score"].copy()
    if not embedding["cell"].equals(vectors["cell"]):
        vectors = embedding[["cell"]].merge(vectors, on="cell", validate="one_to_one")

    connectivity = cluster_connectivity(embedding, args.target_cluster, args.knn)
    related = connectivity.head(args.related_clusters)["related_cluster"].tolist()
    included = [args.target_cluster] + related
    refined, graph = refined_forceatlas(embedding, included, args.knn, args.fa_iterations)
    refined_vectors = project_vectors_to_refined(vectors, refined)

    sampled = stratified_sample(vectors, args.arrow_fraction)
    sampled_refined = stratified_sample(refined_vectors, args.arrow_fraction)
    grid_clean = grid[grid["n_cells"] >= 5].copy()
    cluster_group = group_summary(embedding)

    embedding.to_csv(tables / "fa_embedding_cells_used.csv", index=False)
    vectors.to_csv(tables / "cell_shift_vectors_FA_used.csv", index=False)
    grid_clean.to_csv(tables / "grid_shift_vectors_FA_ncells_ge5.csv", index=False)
    propagation.to_csv(tables / "propagation_magnitude_by_cluster_used.csv", index=False)
    score.to_csv(tables / "perturbation_score_FA_used.csv", index=False)
    connectivity.to_csv(tables / "cluster4_knn_connectivity.csv", index=False)
    cluster_group.to_csv(tables / "delta_l2_summary_by_cluster_group.csv", index=False)
    sampled.to_csv(tables / "cell_shift_vectors_FA_sampled.csv", index=False)
    refined.to_csv(tables / "cluster4_related_FA_refined_embedding.csv", index=False)
    refined_vectors.to_csv(tables / "cluster4_related_FA_refined_vectors.csv", index=False)
    sampled_refined.to_csv(tables / "cluster4_related_FA_refined_vectors_sampled.csv", index=False)
    sparse.save_npz(tables / "cluster4_related_knn_graph.npz", graph)

    delta_map(
        embedding,
        figures / "FA_delta_l2_clean",
        "Predicted transcriptomic impact of HIF1A knockout across all cells",
    )
    delta_map(
        embedding,
        figures / "FA_delta_l2_cluster4_highlight",
        "HIF1A knockout impact concentrated in Seurat cluster 4",
        highlight=args.target_cluster,
    )
    arrow_plot(
        vectors,
        sampled,
        figures / "FA_cell_shift_arrows_clean_sampled",
        "Predicted cell-state displacement after HIF1A knockout (stratified sample)",
    )
    flow_plot(
        embedding,
        grid,
        figures / "FA_flow_on_grid_clean",
        "Smoothed HIF1A knockout response field on the global FA embedding",
    )
    distribution_plot(embedding, figures / "delta_l2_by_cluster_group_clean")
    propagation_plot(
        propagation,
        included,
        args.target_cluster,
        figures / "propagation_magnitude_clean",
    )

    delta_map(
        refined,
        figures / "FA_refined_subset_delta_l2",
        f"HIF1A knockout impact in cluster 4 and connected clusters {', '.join(related)}",
        x="FA1_refined",
        y="FA2_refined",
    )
    delta_map(
        refined,
        figures / "FA_refined_subset_delta_l2_cluster4_highlight",
        "Refined view of the predicted HIF1A knockout impact in cluster 4",
        x="FA1_refined",
        y="FA2_refined",
        highlight=args.target_cluster,
    )
    arrow_plot(
        refined_vectors,
        sampled_refined,
        figures / "FA_refined_subset_cell_shift_arrows",
        "Projected HIF1A knockout displacement in the cluster 4 neighborhood",
        x="FA1_refined",
        y="FA2_refined",
        u="shift_FA1_refined",
        v="shift_FA2_refined",
    )
    refined_grid = estimate_grid(
        refined_vectors,
        "FA1_refined",
        "FA2_refined",
        "shift_FA1_refined",
        "shift_FA2_refined",
    )
    refined_grid_kept = flow_plot(
        refined,
        refined_grid,
        figures / "FA_refined_subset_flow_on_grid",
        "Smoothed HIF1A knockout response field in the cluster 4 neighborhood",
        x="FA1_refined",
        y="FA2_refined",
        gx="grid_FA1_refined",
        gy="grid_FA2_refined",
        gu="flow_FA1_refined",
        gv="flow_FA2_refined",
    )
    refined_grid_kept.to_csv(tables / "cluster4_related_FA_refined_grid_ncells_ge5.csv", index=False)

    for group in GROUPS:
        group_vectors = vectors[vectors[GROUP_COL] == group].copy()
        group_cells = embedding[embedding[GROUP_COL] == group].copy()
        group_sample = stratified_sample(group_vectors, args.arrow_fraction)
        group_grid = estimate_grid(group_vectors, "FA1", "FA2", "shift_FA1", "shift_FA2")
        kept = flow_plot(
            group_cells,
            group_grid,
            figures / f"FA_{group}_flow_on_grid_clean",
            f"Smoothed HIF1A knockout response field in {group} cells",
        )
        delta_map(
            group_cells,
            figures / f"FA_{group}_delta_l2_clean",
            f"Predicted HIF1A knockout impact in {group} cells",
        )
        arrow_plot(
            group_vectors,
            group_sample,
            figures / f"FA_{group}_cell_shift_arrows_clean_sampled",
            f"Predicted cell-state displacement after HIF1A knockout in {group}",
        )
        group_sample.to_csv(tables / f"{group}_cell_shift_vectors_FA_sampled.csv", index=False)
        kept.to_csv(tables / f"{group}_grid_shift_vectors_FA_ncells_ge5.csv", index=False)

        refined_group_vectors = refined_vectors[refined_vectors[GROUP_COL] == group].copy()
        refined_group_cells = refined[refined[GROUP_COL] == group].copy()
        refined_group_sample = stratified_sample(refined_group_vectors, args.arrow_fraction)
        refined_group_grid = estimate_grid(
            refined_group_vectors,
            "FA1_refined",
            "FA2_refined",
            "shift_FA1_refined",
            "shift_FA2_refined",
        )
        refined_kept = flow_plot(
            refined_group_cells,
            refined_group_grid,
            figures / f"FA_refined_subset_{group}_flow_on_grid",
            f"Cluster 4 neighborhood response field in {group} cells",
            x="FA1_refined",
            y="FA2_refined",
            gx="grid_FA1_refined",
            gy="grid_FA2_refined",
            gu="flow_FA1_refined",
            gv="flow_FA2_refined",
        )
        delta_map(
            refined_group_cells,
            figures / f"FA_refined_subset_{group}_delta_l2",
            f"HIF1A knockout impact in the cluster 4 neighborhood: {group}",
            x="FA1_refined",
            y="FA2_refined",
        )
        arrow_plot(
            refined_group_vectors,
            refined_group_sample,
            figures / f"FA_refined_subset_{group}_cell_shift_arrows",
            f"Projected HIF1A knockout displacement near cluster 4: {group}",
            x="FA1_refined",
            y="FA2_refined",
            u="shift_FA1_refined",
            v="shift_FA2_refined",
        )
        refined_group_sample.to_csv(
            tables / f"cluster4_related_{group}_FA_refined_vectors_sampled.csv", index=False
        )
        refined_kept.to_csv(
            tables / f"cluster4_related_{group}_FA_refined_grid_ncells_ge5.csv", index=False
        )

    means = embedding.groupby(CLUSTER_COL)["delta_l2"].mean().sort_values(ascending=False)
    cluster4_group = (
        embedding[embedding[CLUSTER_COL] == args.target_cluster]
        .groupby(GROUP_COL)["delta_l2"]
        .mean()
    )
    cd_minus_healthy = float(cluster4_group["CD"] - cluster4_group["Healthy"])
    figure_names = sorted(path.name for path in figures.glob("*.png"))
    summary = {
        "output_dir": str(out.resolve()),
        "source_dir": str(source.resolve()),
        "n_cells_global": int(len(embedding)),
        "n_cells_refined_subset": int(len(refined)),
        "refined_subset_clusters": included,
        "related_cluster_selection": "Top kNN cross-cluster edge counts to cluster 4 in robust-scaled global FA",
        "cluster_with_highest_mean_delta_l2": str(means.index[0]),
        "highest_mean_delta_l2": float(means.iloc[0]),
        "cluster4_mean_delta_l2_CD": float(cluster4_group["CD"]),
        "cluster4_mean_delta_l2_Healthy": float(cluster4_group["Healthy"]),
        "cluster4_mean_delta_l2_CD_minus_Healthy": cd_minus_healthy,
        "arrow_sample_fraction_per_cluster": args.arrow_fraction,
        "minimum_grid_cells": 5,
        "figures_png": figure_names,
    }
    (logs / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    (logs / "run_parameters.json").write_text(json.dumps(vars(args), indent=2) + "\n")
    (logs / "methodology.txt").write_text(
        "This is a visualization-only analysis. Original delta_l2 values were not recalculated.\n"
        "Related clusters were ranked using bidirectional cross-cluster edges in a kNN graph of\n"
        "robust-scaled global FA coordinates. ForceAtlas2 was rerun on the induced subset graph.\n"
        "Original CellOracle shift vectors were mapped to the refined embedding using local affine\n"
        "transformations fitted on 20 neighboring cells. Refined and group-specific grid flows are\n"
        "Gaussian-weighted local averages of the corresponding projected cell vectors.\n"
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
