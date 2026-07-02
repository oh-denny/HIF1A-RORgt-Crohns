#!/usr/bin/env python3

import argparse
import copy
import json
import os
from datetime import datetime
from pathlib import Path

PROJECT_DIR = Path(__file__).resolve().parents[2]
os.environ["HOME"] = str(PROJECT_DIR / "results" / "celloracle_home")
os.environ.setdefault("XDG_CONFIG_HOME", str(PROJECT_DIR / "results" / "xdg_config"))
os.environ.setdefault("XDG_CACHE_HOME", str(PROJECT_DIR / "results" / "xdg_cache"))
os.environ.setdefault("MPLCONFIGDIR", str(PROJECT_DIR / "results" / "matplotlib_cache"))
os.environ.setdefault("NUMBA_CACHE_DIR", str(PROJECT_DIR / "results" / "celloracle_numba_cache"))
for env_name in ["HOME", "XDG_CONFIG_HOME", "XDG_CACHE_HOME", "MPLCONFIGDIR", "NUMBA_CACHE_DIR"]:
    Path(os.environ[env_name]).mkdir(parents=True, exist_ok=True)

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
from sklearn.neighbors import NearestNeighbors

try:
    import gimmemotifs.motif as gimmemotifs_motif
    if not hasattr(gimmemotifs_motif, "default_motifs"):
        gimmemotifs_motif.default_motifs = lambda *args, **kwargs: []
except Exception:
    pass

import celloracle as co


FA_KEYS = ["X_draw_graph_fa", "X_fa", "X_forceatlas2", "X_force_directed", "FA", "FA2"]


def parse_args():
    parser = argparse.ArgumentParser(description="Generate FA paper-style HIF1A KO visualizations.")
    parser.add_argument("--oracle", default="results/celloracle_hif1a_ko_1500/objects/oracle_HIF1A_KO.celloracle.oracle")
    parser.add_argument("--links", default="results/celloracle_hif1a_ko_1500/objects/links_filtered.celloracle.links")
    parser.add_argument("--cell-metrics", default="results/celloracle_hif1a_ko_1500/tables/cell_delta_metrics.csv")
    parser.add_argument("--out-dir", default="results/celloracle_hif1a_ko_1500_FA_paperstyle")
    parser.add_argument("--target-gene", default="HIF1A")
    parser.add_argument("--cluster-column", default="seurat_clusters")
    parser.add_argument("--group-column", default="group")
    parser.add_argument("--n-jobs", type=int, default=1)
    parser.add_argument("--sample-arrows", type=int, default=700)
    parser.add_argument("--grid-steps", type=int, default=40)
    parser.add_argument("--grid-neighbors", type=int, default=200)
    return parser.parse_args()


def unique_out_dir(base_dir):
    base = Path(base_dir)
    if not base.exists():
        return base
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    return base.with_name(f"{base.name}_{stamp}")


def savefig(fig, out_base):
    out_base = Path(out_base)
    fig.savefig(out_base.with_suffix(".png"), dpi=400, bbox_inches="tight")
    fig.savefig(out_base.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def visual_arrow_scale(xy, uv, fraction=0.06):
    xy = np.asarray(xy)
    uv = np.asarray(uv)
    diag = np.linalg.norm(np.nanmax(xy, axis=0) - np.nanmin(xy, axis=0))
    mag = np.linalg.norm(uv, axis=1)
    ref = np.nanpercentile(mag[mag > 0], 95) if np.any(mag > 0) else 1.0
    if not np.isfinite(ref) or ref == 0:
        ref = 1.0
    return (diag * fraction) / ref


def load_oracle(path):
    return co.load_hdf5(str(path))


def get_or_make_fa(adata, logs_dir):
    for key in FA_KEYS:
        if key in adata.obsm.keys():
            fa = np.asarray(adata.obsm[key])
            adata.obsm["X_fa"] = fa[:, :2]
            return "recovered", key

    if "neighbors" not in adata.uns:
        if "X_pca" not in adata.obsm:
            sc.tl.pca(adata, n_comps=min(50, adata.shape[1] - 1), svd_solver="arpack")
        sc.pp.neighbors(adata, n_neighbors=30, n_pcs=min(50, adata.obsm["X_pca"].shape[1]))

    try:
        sc.tl.draw_graph(adata, layout="fa", random_state=123)
        adata.obsm["X_fa"] = np.asarray(adata.obsm["X_draw_graph_fa"])[:, :2]
        return "calculated", "X_draw_graph_fa"
    except Exception as exc:
        (logs_dir / "forceatlas2_fallback_used.txt").write_text(
            "sc.tl.draw_graph(layout='fa') failed. Used layout='fr' as force-directed fallback.\n"
            f"Error: {exc}\n"
        )
        sc.tl.draw_graph(adata, layout="fr", random_state=123)
        adata.obsm["X_fa"] = np.asarray(adata.obsm["X_draw_graph_fr"])[:, :2]
        return "calculated_fallback_fr", "X_draw_graph_fr"


def cluster_palette(values):
    levels = sorted(pd.Series(values).astype(str).unique(), key=lambda x: (len(x), x))
    colors = sns.color_palette("tab10", n_colors=max(10, len(levels)))
    return {level: colors[i % len(colors)] for i, level in enumerate(levels)}


def make_embedding_table(oracle, metrics, cluster_col, group_col):
    adata = oracle.adata
    fa = np.asarray(adata.obsm["X_fa"])
    df = pd.DataFrame({
        "cell": adata.obs_names,
        "FA1": fa[:, 0],
        "FA2": fa[:, 1],
        cluster_col: adata.obs[cluster_col].astype(str).values,
        group_col: adata.obs[group_col].astype(str).values if group_col in adata.obs else "NA",
    })
    metrics = metrics.rename(columns={"cell_id": "cell"})
    keep = ["cell", "delta_l2", "target_imputed", "target_simulated"]
    return df.merge(metrics[keep], on="cell", how="left")


def recalc_fa_shift_and_grid(oracle, target_gene, n_jobs, grid_steps, grid_neighbors):
    oracle.embedding_name = "X_fa"
    oracle.embedding = np.asarray(oracle.adata.obsm["X_fa"])
    oracle.simulate_shift(
        perturb_condition={target_gene: 0.0},
        GRN_unit="cluster",
        n_propagation=3,
        ignore_warning=True,
    )
    oracle.estimate_transition_prob(
        n_neighbors=200,
        knn_random=True,
        sampled_fraction=1.0,
        n_jobs=n_jobs,
        calculate_randomized=True,
    )
    oracle.calculate_embedding_shift(sigma_corr=0.05)
    oracle.calculate_grid_arrows(
        smooth=0.8,
        steps=(grid_steps, grid_steps),
        n_neighbors=grid_neighbors,
        n_jobs=n_jobs,
    )
    return oracle


def make_cell_vectors(oracle, emb_df, cluster_col, group_col):
    shift = np.asarray(oracle.delta_embedding)
    out = emb_df[["cell", "FA1", "FA2", "delta_l2", cluster_col, group_col]].copy()
    out["shift_FA1"] = shift[:, 0]
    out["shift_FA2"] = shift[:, 1]
    out["FA1_end"] = out["FA1"] + out["shift_FA1"]
    out["FA2_end"] = out["FA2"] + out["shift_FA2"]
    return out[["cell", "FA1", "FA2", "shift_FA1", "shift_FA2", "FA1_end", "FA2_end", "delta_l2", cluster_col, group_col]]


def make_grid_vectors(oracle):
    grid = np.asarray(oracle.flow_grid)
    flow = np.asarray(oracle.flow)
    mass = np.asarray(oracle.total_p_mass)
    fa = np.asarray(oracle.adata.obsm["X_fa"])
    step = np.median(np.diff(np.unique(np.round(grid[:, 0], 8))))
    if not np.isfinite(step) or step <= 0:
        step = 0.5
    nbrs = NearestNeighbors(radius=step * 1.5).fit(fa)
    n_cells = np.array([len(x) for x in nbrs.radius_neighbors(grid, return_distance=False)])
    return pd.DataFrame({
        "grid_FA1": grid[:, 0],
        "grid_FA2": grid[:, 1],
        "flow_FA1": flow[:, 0],
        "flow_FA2": flow[:, 1],
        "mass": mass,
        "n_cells": n_cells,
    })


def propagation_summary(oracle, target_gene, cluster_col):
    rows = []
    clusters = oracle.adata.obs[cluster_col].astype(str).values
    for n_prop in range(6):
        if n_prop == 0 and hasattr(oracle, "_Oracle__simulate_shift"):
            oracle._Oracle__simulate_shift(
                perturb_condition={target_gene: 0.0},
                GRN_unit="cluster",
                n_propagation=0,
                ignore_warning=True,
                n_min=0,
                n_max=6,
            )
        else:
            oracle.simulate_shift(
                perturb_condition={target_gene: 0.0},
                GRN_unit="cluster",
                n_propagation=n_prop,
                ignore_warning=True,
            )
        delta = oracle.adata.to_df(layer="delta_X").values
        l2 = np.linalg.norm(delta, axis=1)
        df = pd.DataFrame({cluster_col: clusters, "l2": l2})
        for cluster, sub in df.groupby(cluster_col):
            rows.append({
                cluster_col: cluster,
                "n_propagation": n_prop,
                "mean_l2_norm": sub["l2"].mean(),
                "median_l2_norm": sub["l2"].median(),
                "sd_l2_norm": sub["l2"].std(),
                "n_cells": len(sub),
            })
    return pd.DataFrame(rows)


def find_pseudotime_column(adata):
    candidates = ["pseudotime", "dpt_pseudotime", "monocle_pseudotime", "slingshot_pseudotime", "palantir_pseudotime"]
    lower_map = {c.lower(): c for c in adata.obs.columns}
    for candidate in candidates:
        if candidate.lower() in lower_map:
            return lower_map[candidate.lower()]
    for col in adata.obs.columns:
        if "pseudotime" in col.lower():
            return col
    return None


def perturbation_score(oracle, cell_vectors, cluster_col, group_col, logs_dir):
    ptime_col = find_pseudotime_column(oracle.adata)
    if ptime_col is None:
        (logs_dir / "perturbation_score_not_generated.txt").write_text(
            "Perturbation score was not generated because no pseudotime column was found in adata.obs.\n"
        )
        return None, False, None

    ptime = pd.to_numeric(oracle.adata.obs[ptime_col], errors="coerce").to_numpy()
    fa = np.asarray(oracle.adata.obsm["X_fa"])
    nn = NearestNeighbors(n_neighbors=min(31, len(fa))).fit(fa)
    _, idx = nn.kneighbors(fa)
    traj = np.zeros_like(fa)
    for i, neigh in enumerate(idx):
        forward = neigh[np.isfinite(ptime[neigh]) & np.isfinite(ptime[i]) & (ptime[neigh] > ptime[i])]
        if len(forward) == 0:
            forward = neigh[1:min(6, len(neigh))]
        else:
            forward = forward[:5]
        traj[i, :] = fa[forward].mean(axis=0) - fa[i]

    shift = cell_vectors[["shift_FA1", "shift_FA2"]].to_numpy()
    denom = np.linalg.norm(shift, axis=1) * np.linalg.norm(traj, axis=1)
    score = np.divide((shift * traj).sum(axis=1), denom, out=np.full(len(denom), np.nan), where=denom > 0)
    out = pd.DataFrame({
        "cell": oracle.adata.obs_names,
        "FA1": fa[:, 0],
        "FA2": fa[:, 1],
        "shift_FA1": shift[:, 0],
        "shift_FA2": shift[:, 1],
        "trajectory_FA1": traj[:, 0],
        "trajectory_FA2": traj[:, 1],
        "perturbation_score": score,
        cluster_col: oracle.adata.obs[cluster_col].astype(str).values,
        group_col: oracle.adata.obs[group_col].astype(str).values if group_col in oracle.adata.obs else "NA",
    })
    return out, True, ptime_col


def plot_clusters(df, cluster_col, out_base):
    palette = cluster_palette(df[cluster_col])
    fig, ax = plt.subplots(figsize=(7, 6))
    for cluster, sub in df.groupby(cluster_col):
        ax.scatter(sub["FA1"], sub["FA2"], s=5, color=palette[str(cluster)], alpha=0.8, linewidths=0, label=str(cluster))
    centroids = df.groupby(cluster_col)[["FA1", "FA2"]].median()
    for cluster, row in centroids.iterrows():
        ax.text(row["FA1"], row["FA2"], str(cluster), fontsize=10, weight="bold", ha="center", va="center")
    ax.set_xlabel("FA1")
    ax.set_ylabel("FA2")
    ax.set_title("HIF1A KO - FA clusters")
    ax.legend(title=cluster_col, bbox_to_anchor=(1.02, 1), loc="upper left", frameon=False, markerscale=3)
    sns.despine(ax=ax)
    savefig(fig, out_base)


def plot_cell_arrows(vectors, cluster_col, out_base, sampled=None):
    plot_df = vectors if sampled is None or sampled >= len(vectors) else vectors.sample(sampled, random_state=123)
    palette = cluster_palette(vectors[cluster_col])
    colors = vectors[cluster_col].astype(str).map(palette)
    xy = vectors[["FA1", "FA2"]].to_numpy()
    uv_all = vectors[["shift_FA1", "shift_FA2"]].to_numpy()
    scale_factor = visual_arrow_scale(xy, uv_all, fraction=0.055)
    fig, ax = plt.subplots(figsize=(7, 6))
    ax.scatter(vectors["FA1"], vectors["FA2"], s=4, color=colors, alpha=0.25, linewidths=0)
    ax.quiver(plot_df["FA1"], plot_df["FA2"],
              plot_df["shift_FA1"] * scale_factor,
              plot_df["shift_FA2"] * scale_factor,
              angles="xy", scale_units="xy", scale=1, width=0.0025, color="#111827", alpha=0.75)
    ax.set_xlabel("FA1")
    ax.set_ylabel("FA2")
    ax.set_title("HIF1A KO cell shift vectors")
    sns.despine(ax=ax)
    savefig(fig, out_base)


def plot_grid_flow(emb_df, grid_df, out_base):
    keep = grid_df["mass"] >= 0.01
    xy = grid_df.loc[keep, ["grid_FA1", "grid_FA2"]].to_numpy()
    uv = grid_df.loc[keep, ["flow_FA1", "flow_FA2"]].to_numpy()
    scale_factor = visual_arrow_scale(emb_df[["FA1", "FA2"]].to_numpy(), uv, fraction=0.08)
    fig, ax = plt.subplots(figsize=(7, 6))
    ax.scatter(emb_df["FA1"], emb_df["FA2"], c="lightgray", s=5, alpha=0.25, linewidths=0)
    ax.quiver(grid_df.loc[keep, "grid_FA1"], grid_df.loc[keep, "grid_FA2"],
              grid_df.loc[keep, "flow_FA1"] * scale_factor,
              grid_df.loc[keep, "flow_FA2"] * scale_factor,
              angles="xy", scale_units="xy", scale=1, width=0.004, color="#111827")
    ax.set_xlabel("FA1")
    ax.set_ylabel("FA2")
    ax.set_title("HIF1A KO flow on FA grid")
    sns.despine(ax=ax)
    savefig(fig, out_base)


def plot_delta_map(emb_df, out_base):
    fig, ax = plt.subplots(figsize=(7, 6))
    sca = ax.scatter(emb_df["FA1"], emb_df["FA2"], c=emb_df["delta_l2"], s=7, cmap="magma", linewidths=0)
    fig.colorbar(sca, ax=ax, label="delta L2")
    ax.set_xlabel("FA1")
    ax.set_ylabel("FA2")
    ax.set_title("HIF1A KO effect magnitude")
    sns.despine(ax=ax)
    savefig(fig, out_base)


def plot_l2_distributions(emb_df, cluster_col, group_col, figures_dir):
    fig, ax = plt.subplots(figsize=(7, 4))
    sns.violinplot(data=emb_df, x=cluster_col, y="delta_l2", inner="box", color="#93c5fd", ax=ax)
    ax.set_xlabel("seurat_clusters")
    ax.set_ylabel("delta L2")
    ax.set_title("HIF1A KO delta L2 by cluster")
    sns.despine(ax=ax)
    savefig(fig, figures_dir / "HIF1A_KO_delta_l2_by_cluster")

    fig, ax = plt.subplots(figsize=(8, 4))
    sns.boxplot(data=emb_df, x=cluster_col, y="delta_l2", hue=group_col, fliersize=0.5, ax=ax)
    ax.set_xlabel("seurat_clusters")
    ax.set_ylabel("delta L2")
    ax.set_title("HIF1A KO delta L2 by group and cluster")
    ax.legend(title=group_col, frameon=False)
    sns.despine(ax=ax)
    savefig(fig, figures_dir / "HIF1A_KO_delta_l2_by_group_cluster")


def plot_propagation(prop_df, cluster_col, out_base):
    fig, ax = plt.subplots(figsize=(7, 4))
    sns.lineplot(data=prop_df, x="n_propagation", y="mean_l2_norm", hue=cluster_col, marker="o", ax=ax)
    ax.set_xlabel("n_propagation")
    ax.set_ylabel("mean L2 norm")
    ax.set_title("HIF1A KO propagation magnitude by cluster")
    ax.legend(title=cluster_col, frameon=False, bbox_to_anchor=(1.02, 1), loc="upper left")
    sns.despine(ax=ax)
    savefig(fig, out_base)


def plot_perturbation_score(score_df, out_base):
    fig, ax = plt.subplots(figsize=(7, 6))
    vmax = np.nanpercentile(np.abs(score_df["perturbation_score"]), 98)
    sca = ax.scatter(score_df["FA1"], score_df["FA2"], c=score_df["perturbation_score"],
                     s=7, cmap="coolwarm", vmin=-vmax, vmax=vmax, linewidths=0)
    fig.colorbar(sca, ax=ax, label="perturbation score")
    ax.set_xlabel("FA1")
    ax.set_ylabel("FA2")
    ax.set_title("HIF1A KO perturbation score")
    sns.despine(ax=ax)
    savefig(fig, out_base)


def main():
    args = parse_args()
    out_dir = unique_out_dir(args.out_dir)
    figures_dir = out_dir / "figures"
    tables_dir = out_dir / "tables"
    objects_dir = out_dir / "objects"
    logs_dir = out_dir / "logs"
    for path in [figures_dir, tables_dir, objects_dir, logs_dir]:
        path.mkdir(parents=True, exist_ok=False)

    oracle = load_oracle(Path(args.oracle))
    _ = Path(args.links)  # Kept for explicit provenance; links are already embedded in the fitted Oracle object.
    metrics = pd.read_csv(args.cell_metrics)

    fa_status, fa_source = get_or_make_fa(oracle.adata, logs_dir)
    emb_df = make_embedding_table(oracle, metrics, args.cluster_column, args.group_column)
    emb_df.to_csv(tables_dir / "fa_embedding_cells.csv", index=False)

    prop_df = propagation_summary(oracle, args.target_gene, args.cluster_column)
    prop_df.to_csv(tables_dir / "propagation_magnitude_by_cluster.csv", index=False)

    oracle = recalc_fa_shift_and_grid(oracle, args.target_gene, args.n_jobs, args.grid_steps, args.grid_neighbors)
    oracle.to_hdf5(str(objects_dir / "oracle_HIF1A_KO_FA.celloracle.oracle"))

    cell_vectors = make_cell_vectors(oracle, emb_df, args.cluster_column, args.group_column)
    cell_vectors.to_csv(tables_dir / "cell_shift_vectors_FA.csv", index=False)

    grid_vectors = make_grid_vectors(oracle)
    grid_vectors.to_csv(tables_dir / "grid_shift_vectors_FA.csv", index=False)

    score_df, score_generated, ptime_col = perturbation_score(oracle, cell_vectors, args.cluster_column, args.group_column, logs_dir)
    if score_generated:
        score_df.to_csv(tables_dir / "perturbation_score_FA.csv", index=False)

    plot_clusters(emb_df, args.cluster_column, figures_dir / "FA_clusters_HIF1A_KO")
    plot_cell_arrows(cell_vectors, args.cluster_column, figures_dir / "FA_HIF1A_KO_cell_shift_arrows")
    plot_cell_arrows(cell_vectors, args.cluster_column, figures_dir / "FA_HIF1A_KO_cell_shift_arrows_sampled", sampled=args.sample_arrows)
    plot_grid_flow(emb_df, grid_vectors, figures_dir / "FA_HIF1A_KO_flow_on_grid")
    plot_delta_map(emb_df, figures_dir / "FA_HIF1A_KO_delta_l2_map")
    plot_l2_distributions(emb_df, args.cluster_column, args.group_column, figures_dir)
    plot_propagation(prop_df, args.cluster_column, figures_dir / "HIF1A_KO_propagation_magnitude_by_cluster")
    if score_generated:
        plot_perturbation_score(score_df, figures_dir / "FA_HIF1A_KO_perturbation_score")

    cluster_l2 = emb_df.groupby(args.cluster_column)["delta_l2"].mean().sort_values(ascending=False)
    pivot = prop_df.pivot(index=args.cluster_column, columns="n_propagation", values="mean_l2_norm")
    increase = (pivot[5] - pivot[0]).sort_values(ascending=False)

    summary = {
        "out_dir": str(out_dir),
        "tables": sorted([p.name for p in tables_dir.glob("*.csv")]),
        "figures": sorted([p.name for p in figures_dir.glob("*") if p.suffix in [".png", ".pdf"]]),
        "fa_status": fa_status,
        "fa_source": fa_source,
        "perturbation_score_generated": score_generated,
        "pseudotime_column": ptime_col,
        "n_cells": int(oracle.adata.shape[0]),
        "top_mean_delta_l2_clusters": cluster_l2.head(3).to_dict(),
        "largest_propagation_increase_cluster": increase.index[0],
        "largest_propagation_increase": float(increase.iloc[0]),
    }
    (logs_dir / "summary.json").write_text(json.dumps(summary, indent=2))
    (logs_dir / "run_parameters.json").write_text(json.dumps(vars(args), indent=2))

    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
