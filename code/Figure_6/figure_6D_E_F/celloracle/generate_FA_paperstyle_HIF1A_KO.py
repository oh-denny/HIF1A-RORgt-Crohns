#!/usr/bin/env python3

import argparse
import json
from pathlib import Path

from functions.io_utils import configure_runtime_environment, unique_output_dir

PROJECT_DIR = Path(__file__).resolve().parents[4]
configure_runtime_environment(PROJECT_DIR, include_celloracle_cache=True)

import matplotlib

matplotlib.use("Agg")
import pandas as pd

from functions.celloracle_utils import load_oracle
from functions.embedding_utils import get_or_make_fa, make_cell_vectors, make_embedding_table, make_grid_vectors
from functions.perturbation_utils import perturbation_score, propagation_summary, recalculate_fa_shift_and_grid
from functions.plotting_utils import (
    plot_cell_arrows, plot_clusters, plot_delta_map, plot_grid_flow,
    plot_l2_distributions, plot_perturbation_score_map, plot_propagation_magnitude,
)


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


def main():
    args = parse_args()
    out_dir = unique_output_dir(args.out_dir)
    figures_dir, tables_dir = out_dir / "figures", out_dir / "tables"
    objects_dir, logs_dir = out_dir / "objects", out_dir / "logs"
    for path in (figures_dir, tables_dir, objects_dir, logs_dir):
        path.mkdir(parents=True, exist_ok=False)

    oracle = load_oracle(Path(args.oracle))
    _ = Path(args.links)  # Provenance only; links are embedded in the fitted Oracle.
    metrics = pd.read_csv(args.cell_metrics)
    fa_status, fa_source = get_or_make_fa(oracle.adata, logs_dir)
    embedding = make_embedding_table(oracle, metrics, args.cluster_column, args.group_column)
    embedding.to_csv(tables_dir / "fa_embedding_cells.csv", index=False)

    propagation = propagation_summary(oracle, args.target_gene, args.cluster_column)
    propagation.to_csv(tables_dir / "propagation_magnitude_by_cluster.csv", index=False)
    oracle = recalculate_fa_shift_and_grid(oracle, args.target_gene, args.n_jobs, args.grid_steps, args.grid_neighbors)
    oracle.to_hdf5(str(objects_dir / "oracle_HIF1A_KO_FA.celloracle.oracle"))

    cell_vectors = make_cell_vectors(oracle, embedding, args.cluster_column, args.group_column)
    cell_vectors.to_csv(tables_dir / "cell_shift_vectors_FA.csv", index=False)
    grid_vectors = make_grid_vectors(oracle)
    grid_vectors.to_csv(tables_dir / "grid_shift_vectors_FA.csv", index=False)
    score, score_generated, pseudotime_column = perturbation_score(
        oracle, cell_vectors, args.cluster_column, args.group_column, logs_dir
    )
    if score_generated:
        score.to_csv(tables_dir / "perturbation_score_FA.csv", index=False)

    plot_clusters(embedding, args.cluster_column, figures_dir / "FA_clusters_HIF1A_KO")
    plot_cell_arrows(cell_vectors, args.cluster_column, figures_dir / "FA_HIF1A_KO_cell_shift_arrows")
    plot_cell_arrows(cell_vectors, args.cluster_column, figures_dir / "FA_HIF1A_KO_cell_shift_arrows_sampled", sampled=args.sample_arrows)
    plot_grid_flow(embedding, grid_vectors, figures_dir / "FA_HIF1A_KO_flow_on_grid")
    plot_delta_map(embedding, figures_dir / "FA_HIF1A_KO_delta_l2_map")
    plot_l2_distributions(embedding, args.cluster_column, args.group_column, figures_dir)
    plot_propagation_magnitude(
        propagation, args.cluster_column, figures_dir / "HIF1A_KO_propagation_magnitude_by_cluster", target_cluster="4"
    )
    if score_generated:
        plot_perturbation_score_map(score, figures_dir / "FA_HIF1A_KO_perturbation_score")

    cluster_l2 = embedding.groupby(args.cluster_column)["delta_l2"].mean().sort_values(ascending=False)
    pivot = propagation.pivot(index=args.cluster_column, columns="n_propagation", values="mean_l2_norm")
    increase = (pivot[5] - pivot[0]).sort_values(ascending=False)
    summary = {
        "out_dir": str(out_dir),
        "tables": sorted(path.name for path in tables_dir.glob("*.csv")),
        "figures": sorted(path.name for path in figures_dir.glob("*") if path.suffix in [".png", ".pdf"]),
        "fa_status": fa_status, "fa_source": fa_source,
        "perturbation_score_generated": score_generated, "pseudotime_column": pseudotime_column,
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
