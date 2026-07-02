#!/usr/bin/env python3

import argparse
import json
import os
from pathlib import Path

PROJECT_DIR = Path(__file__).resolve().parents[2]
os.environ["HOME"] = str(PROJECT_DIR / "results" / "celloracle_home")
os.environ.setdefault("XDG_CONFIG_HOME", str(PROJECT_DIR / "results" / "xdg_config"))
os.environ.setdefault("XDG_CACHE_HOME", str(PROJECT_DIR / "results" / "xdg_cache"))
os.environ.setdefault("MPLCONFIGDIR", str(PROJECT_DIR / "results" / "matplotlib_cache"))
os.environ.setdefault("NUMBA_CACHE_DIR", str(PROJECT_DIR / "results" / "celloracle_numba_cache"))
for env_name in ["HOME", "XDG_CONFIG_HOME", "XDG_CACHE_HOME", "MPLCONFIGDIR", "NUMBA_CACHE_DIR"]:
    Path(os.environ[env_name]).mkdir(parents=True, exist_ok=True)

import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import io

# CellOracle 0.20 imports an old gimmemotifs symbol at package import time.
# The knockout workflow below uses prebuilt promoter base-GRNs, not motif scans.
try:
    import gimmemotifs.motif as gimmemotifs_motif
    if not hasattr(gimmemotifs_motif, "default_motifs"):
        gimmemotifs_motif.default_motifs = lambda *args, **kwargs: []
except Exception:
    pass

import celloracle as co


def parse_args():
    parser = argparse.ArgumentParser(
        description="Run CellOracle in silico knockout analysis for HIF1A."
    )
    parser.add_argument("--export-dir", default="celloracle_export_hif1a")
    parser.add_argument("--out-dir", default="results/celloracle_hif1a_ko")
    parser.add_argument("--target-gene", default="HIF1A")
    parser.add_argument("--cluster-column", default="seurat_clusters")
    parser.add_argument("--group-column", default="group")
    parser.add_argument("--base-grn", default=None,
                        help="Optional parquet/csv/pickle CellOracle TF_info_matrix.")
    parser.add_argument("--species", default="human", choices=["human", "mouse"])
    parser.add_argument("--genome", default="hg38",
                        help="human: hg19 or hg38; mouse: mm10 or mm39.")
    parser.add_argument("--n-pca", type=int, default=50)
    parser.add_argument("--k", type=int, default=None)
    parser.add_argument("--alpha", type=float, default=10)
    parser.add_argument("--bagging-number", type=int, default=20)
    parser.add_argument("--link-p", type=float, default=0.001)
    parser.add_argument("--link-threshold-number", type=int, default=2000)
    parser.add_argument("--n-propagation", type=int, default=3)
    parser.add_argument("--transition-neighbors", type=int, default=200)
    parser.add_argument("--sampled-fraction", type=float, default=1.0)
    parser.add_argument("--n-jobs", type=int, default=-1)
    parser.add_argument("--test-mode", action="store_true",
                        help="Run CellOracle get_links on one cluster only.")
    return parser.parse_args()


def resolve_gene(gene, genes):
    gene = gene.strip()
    candidates = [gene, gene.upper(), gene.capitalize()]
    for candidate in candidates:
        if candidate in genes:
            return candidate
    raise ValueError(f"Target gene {gene!r} was not found in exported genes.")


def load_export(export_dir, cluster_column):
    export_dir = Path(export_dir)
    counts = io.mmread(export_dir / "counts_genes_x_cells.mtx").tocsr()
    genes = pd.read_csv(export_dir / "genes.tsv", header=None)[0].astype(str).tolist()
    cells = pd.read_csv(export_dir / "cells.tsv", header=None)[0].astype(str).tolist()
    metadata = pd.read_csv(export_dir / "metadata.csv")
    umap = pd.read_csv(export_dir / "umap.csv")

    metadata = metadata.set_index("cell_id").loc[cells]
    umap = umap.set_index("cell_id").loc[cells]

    if cluster_column not in metadata.columns:
        raise ValueError(f"Missing cluster column in metadata: {cluster_column}")

    adata = ad.AnnData(X=counts.T.tocsr())
    adata.obs_names = cells
    adata.var_names = genes
    adata.obs = metadata
    adata.obs[cluster_column] = adata.obs[cluster_column].astype(str).astype("category")
    adata.obsm["X_umap"] = umap[["UMAP_1", "UMAP_2"]].to_numpy()
    return adata


def load_base_grn(args):
    if args.base_grn is not None:
        path = Path(args.base_grn)
        if path.suffix == ".parquet":
            return pd.read_parquet(path)
        if path.suffix in {".csv", ".txt", ".tsv"}:
            sep = "\t" if path.suffix in {".txt", ".tsv"} else ","
            return pd.read_csv(path, sep=sep)
        if path.suffix in {".pkl", ".pickle"}:
            return pd.read_pickle(path)
        raise ValueError(f"Unsupported base GRN format: {path.suffix}")

    if args.species == "human":
        version = f"{args.genome}_gimmemotifsv5_fpr2"
        return co.data.load_human_promoter_base_GRN(version=version)
    if args.species == "mouse":
        version = f"{args.genome}_gimmemotifsv5_fpr2"
        return co.data.load_mouse_promoter_base_GRN(version=version)
    raise ValueError(f"Unsupported species: {args.species}")


def save_links_tables(links, out_tables):
    raw_parts = []
    filtered_parts = []
    for cluster, df in links.links_dict.items():
        tmp = df.copy()
        tmp["cluster"] = cluster
        raw_parts.append(tmp)
    if hasattr(links, "filtered_links"):
        for cluster, df in links.filtered_links.items():
            tmp = df.copy()
            tmp["cluster"] = cluster
            filtered_parts.append(tmp)

    if raw_parts:
        pd.concat(raw_parts, axis=0).to_csv(out_tables / "links_raw.csv", index=False)
    if filtered_parts:
        pd.concat(filtered_parts, axis=0).to_csv(out_tables / "links_filtered.csv", index=False)

    try:
        links.get_network_score()
        links.merged_score.to_csv(out_tables / "network_scores.csv", index=False)
    except Exception as exc:
        (out_tables / "network_scores_warning.txt").write_text(str(exc))


def save_delta_tables(oracle, target_gene, cluster_column, group_column, out_tables):
    delta = oracle.adata.to_df(layer="delta_X")
    imputed = oracle.adata.to_df(layer="imputed_count")
    simulated = oracle.adata.to_df(layer="simulated_count")

    summary = pd.DataFrame({
        "gene": delta.columns,
        "mean_delta": delta.mean(axis=0).values,
        "median_delta": delta.median(axis=0).values,
        "mean_abs_delta": delta.abs().mean(axis=0).values,
        "mean_imputed": imputed.mean(axis=0).values,
        "mean_simulated": simulated.mean(axis=0).values,
    })
    summary["abs_rank"] = summary["mean_abs_delta"].rank(ascending=False, method="dense")
    summary.sort_values("mean_abs_delta", ascending=False).to_csv(
        out_tables / "gene_delta_summary.csv", index=False
    )

    obs = oracle.adata.obs.copy()
    vector_length = np.linalg.norm(delta.values, ord=2, axis=1)
    cell_metrics = pd.DataFrame({
        "cell_id": oracle.adata.obs_names,
        "delta_l2": vector_length,
        "target_imputed": imputed[target_gene].values,
        "target_simulated": simulated[target_gene].values,
        cluster_column: obs[cluster_column].astype(str).values,
    })
    if group_column in obs.columns:
        cell_metrics[group_column] = obs[group_column].astype(str).values
    cell_metrics.to_csv(out_tables / "cell_delta_metrics.csv", index=False)

    group_cols = [cluster_column]
    if group_column in cell_metrics.columns:
        group_cols.append(group_column)
    cluster_summary = (
        cell_metrics
        .groupby(group_cols, observed=True)
        .agg(cells=("cell_id", "size"),
             mean_delta_l2=("delta_l2", "mean"),
             median_delta_l2=("delta_l2", "median"),
             mean_target_imputed=("target_imputed", "mean"))
        .reset_index()
    )
    cluster_summary.to_csv(out_tables / "cluster_delta_summary.csv", index=False)


def save_figures(oracle, target_gene, out_figures):
    sc.pl.embedding(
        oracle.adata,
        basis="umap",
        color=[oracle.cluster_column_name, target_gene],
        show=False,
        frameon=False,
        ncols=2,
    )
    plt.tight_layout()
    plt.savefig(out_figures / "umap_cluster_and_target.png", dpi=300)
    plt.close()

    fig, ax = plt.subplots(figsize=(6, 5))
    oracle.plot_quiver(ax=ax, scale=30, s=6)
    ax.set_title(f"{target_gene} KO simulated cell shift")
    fig.tight_layout()
    fig.savefig(out_figures / "hif1a_ko_quiver_cells.png", dpi=300)
    plt.close(fig)

    min_mass = 0.01
    mass_filter = oracle.total_p_mass >= min_mass
    fig, ax = plt.subplots(figsize=(6, 5))
    ax.scatter(
        oracle.embedding[:, 0],
        oracle.embedding[:, 1],
        c="lightgray",
        s=5,
        alpha=0.25,
        linewidths=0,
    )
    ax.quiver(
        oracle.flow_grid[mass_filter, 0],
        oracle.flow_grid[mass_filter, 1],
        oracle.flow[mass_filter, 0],
        oracle.flow[mass_filter, 1],
        angles="xy",
        scale_units="xy",
        scale=1,
        width=0.004,
        color="#1f2937",
    )
    ax.set_title(f"{target_gene} KO simulated flow on grid")
    ax.axis("off")
    plt.tight_layout()
    fig.savefig(out_figures / "hif1a_ko_grid_arrows.png", dpi=300)
    plt.close(fig)

    try:
        fig = oracle.estimate_impact_of_perturbations_under_various_ns(
            perturb_condition={target_gene: 0.0},
            n_prop_max=5,
            GRN_unit="cluster",
            figsize=[7, 3],
        )
        fig.savefig(out_figures / "n_propagation_impact.png", dpi=300)
        plt.close(fig)
    except Exception as exc:
        (out_figures / "n_propagation_impact_warning.txt").write_text(str(exc))


def main():
    args = parse_args()
    out_dir = Path(args.out_dir)
    out_objects = out_dir / "objects"
    out_tables = out_dir / "tables"
    out_figures = out_dir / "figures"
    out_logs = out_dir / "logs"
    for path in [out_objects, out_tables, out_figures, out_logs]:
        path.mkdir(parents=True, exist_ok=True)

    adata = load_export(args.export_dir, args.cluster_column)
    target_gene = resolve_gene(args.target_gene, adata.var_names)

    adata.write_h5ad(out_objects / "input_for_celloracle.h5ad")

    base_grn = load_base_grn(args)
    base_grn.to_csv(out_tables / "base_grn_used.csv", index=False)

    oracle = co.Oracle()
    oracle.import_anndata_as_raw_count(
        adata=adata,
        cluster_column_name=args.cluster_column,
        embedding_name="X_umap",
    )
    oracle.import_TF_data(TF_info_matrix=base_grn)
    oracle.to_hdf5(str(out_objects / "oracle_imported.celloracle.oracle"))

    oracle.perform_PCA(n_components=args.n_pca)
    k = args.k if args.k is not None else max(20, int(oracle.adata.shape[0] * 0.025))
    oracle.knn_imputation(
        n_pca_dims=args.n_pca,
        k=k,
        balanced=True,
        b_sight=k * 8,
        b_maxl=k * 4,
        n_jobs=args.n_jobs,
    )
    oracle.to_hdf5(str(out_objects / "oracle_preprocessed.celloracle.oracle"))

    links = oracle.get_links(
        cluster_name_for_GRN_unit=args.cluster_column,
        alpha=args.alpha,
        bagging_number=args.bagging_number,
        verbose_level=1,
        test_mode=args.test_mode,
        n_jobs=args.n_jobs,
    )
    links.filter_links(
        p=args.link_p,
        weight="coef_abs",
        threshold_number=args.link_threshold_number,
    )
    links.to_hdf5(str(out_objects / "links_filtered.celloracle.links"))
    save_links_tables(links, out_tables)

    oracle.get_cluster_specific_TFdict_from_Links(links_object=links)
    oracle.fit_GRN_for_simulation(
        GRN_unit="cluster",
        alpha=args.alpha,
        use_cluster_specific_TFdict=True,
        verbose_level=1,
    )
    oracle.to_hdf5(str(out_objects / "oracle_fitted_grn.celloracle.oracle"))

    oracle.simulate_shift(
        perturb_condition={target_gene: 0.0},
        GRN_unit="cluster",
        n_propagation=args.n_propagation,
        ignore_warning=True,
    )
    oracle.estimate_transition_prob(
        n_neighbors=args.transition_neighbors,
        knn_random=True,
        sampled_fraction=args.sampled_fraction,
        n_jobs=args.n_jobs,
        calculate_randomized=True,
    )
    oracle.calculate_embedding_shift(sigma_corr=0.05)
    oracle.calculate_grid_arrows(smooth=0.8, steps=(40, 40), n_neighbors=200, n_jobs=args.n_jobs)
    oracle.calculate_mass_filter(min_mass=0.01)
    oracle.to_hdf5(str(out_objects / f"oracle_{target_gene}_KO.celloracle.oracle"))

    save_delta_tables(oracle, target_gene, args.cluster_column, args.group_column, out_tables)
    save_figures(oracle, target_gene, out_figures)

    run_info = {
        "celloracle_version": co.__version__,
        "target_gene": target_gene,
        "n_cells": int(oracle.adata.shape[0]),
        "n_genes": int(oracle.adata.shape[1]),
        "cluster_column": args.cluster_column,
        "k": int(k),
        "alpha": args.alpha,
        "bagging_number": args.bagging_number,
        "n_propagation": args.n_propagation,
    }
    (out_logs / "run_info.json").write_text(json.dumps(run_info, indent=2))
    print(json.dumps(run_info, indent=2))


if __name__ == "__main__":
    main()
