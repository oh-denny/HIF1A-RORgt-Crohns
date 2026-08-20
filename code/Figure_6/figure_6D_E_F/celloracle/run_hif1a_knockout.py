#!/usr/bin/env python3

import argparse
import json
from pathlib import Path

from functions.io_utils import configure_runtime_environment

PROJECT_DIR = Path(__file__).resolve().parents[4]
configure_runtime_environment(PROJECT_DIR, include_celloracle_cache=True)

import matplotlib

matplotlib.use("Agg")

from functions.celloracle_utils import (
    load_base_grn, load_export, patch_gimmemotifs_compatibility, resolve_gene,
    save_delta_tables, save_links_tables,
)
from functions.plotting_utils import save_primary_knockout_figures

patch_gimmemotifs_compatibility()
import celloracle as co


def parse_args():
    parser = argparse.ArgumentParser(description="Run CellOracle in silico knockout analysis for HIF1A.")
    parser.add_argument("--export-dir", default="celloracle_export_hif1a")
    parser.add_argument("--out-dir", default="results/celloracle_hif1a_ko")
    parser.add_argument("--target-gene", default="HIF1A")
    parser.add_argument("--cluster-column", default="seurat_clusters")
    parser.add_argument("--group-column", default="group")
    parser.add_argument("--base-grn", default=None, help="Optional parquet/csv/pickle CellOracle TF_info_matrix.")
    parser.add_argument("--species", default="human", choices=["human", "mouse"])
    parser.add_argument("--genome", default="hg38", help="human: hg19 or hg38; mouse: mm10 or mm39.")
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
    parser.add_argument("--test-mode", action="store_true", help="Run CellOracle get_links on one cluster only.")
    return parser.parse_args()


def main():
    args = parse_args()
    out_dir = Path(args.out_dir)
    out_objects, out_tables = out_dir / "objects", out_dir / "tables"
    out_figures, out_logs = out_dir / "figures", out_dir / "logs"
    for path in (out_objects, out_tables, out_figures, out_logs):
        path.mkdir(parents=True, exist_ok=True)

    adata = load_export(args.export_dir, args.cluster_column)
    target_gene = resolve_gene(args.target_gene, adata.var_names)
    adata.write_h5ad(out_objects / "input_for_celloracle.h5ad")
    base_grn = load_base_grn(args)
    base_grn.to_csv(out_tables / "base_grn_used.csv", index=False)

    oracle = co.Oracle()
    oracle.import_anndata_as_raw_count(adata=adata, cluster_column_name=args.cluster_column, embedding_name="X_umap")
    oracle.import_TF_data(TF_info_matrix=base_grn)
    oracle.to_hdf5(str(out_objects / "oracle_imported.celloracle.oracle"))
    oracle.perform_PCA(n_components=args.n_pca)
    k = args.k if args.k is not None else max(20, int(oracle.adata.shape[0] * 0.025))
    oracle.knn_imputation(
        n_pca_dims=args.n_pca, k=k, balanced=True, b_sight=k * 8,
        b_maxl=k * 4, n_jobs=args.n_jobs,
    )
    oracle.to_hdf5(str(out_objects / "oracle_preprocessed.celloracle.oracle"))

    links = oracle.get_links(
        cluster_name_for_GRN_unit=args.cluster_column, alpha=args.alpha,
        bagging_number=args.bagging_number, verbose_level=1,
        test_mode=args.test_mode, n_jobs=args.n_jobs,
    )
    links.filter_links(p=args.link_p, weight="coef_abs", threshold_number=args.link_threshold_number)
    links.to_hdf5(str(out_objects / "links_filtered.celloracle.links"))
    save_links_tables(links, out_tables)
    oracle.get_cluster_specific_TFdict_from_Links(links_object=links)
    oracle.fit_GRN_for_simulation(GRN_unit="cluster", alpha=args.alpha, use_cluster_specific_TFdict=True, verbose_level=1)
    oracle.to_hdf5(str(out_objects / "oracle_fitted_grn.celloracle.oracle"))

    oracle.simulate_shift(
        perturb_condition={target_gene: 0.0}, GRN_unit="cluster",
        n_propagation=args.n_propagation, ignore_warning=True,
    )
    oracle.estimate_transition_prob(
        n_neighbors=args.transition_neighbors, knn_random=True,
        sampled_fraction=args.sampled_fraction, n_jobs=args.n_jobs,
        calculate_randomized=True,
    )
    oracle.calculate_embedding_shift(sigma_corr=0.05)
    oracle.calculate_grid_arrows(smooth=0.8, steps=(40, 40), n_neighbors=200, n_jobs=args.n_jobs)
    oracle.calculate_mass_filter(min_mass=0.01)
    oracle.to_hdf5(str(out_objects / f"oracle_{target_gene}_KO.celloracle.oracle"))

    save_delta_tables(oracle, target_gene, args.cluster_column, args.group_column, out_tables)
    save_primary_knockout_figures(oracle, target_gene, out_figures)
    run_info = {
        "celloracle_version": co.__version__, "target_gene": target_gene,
        "n_cells": int(oracle.adata.shape[0]), "n_genes": int(oracle.adata.shape[1]),
        "cluster_column": args.cluster_column, "k": int(k), "alpha": args.alpha,
        "bagging_number": args.bagging_number, "n_propagation": args.n_propagation,
    }
    (out_logs / "run_info.json").write_text(json.dumps(run_info, indent=2))
    print(json.dumps(run_info, indent=2))


if __name__ == "__main__":
    main()
