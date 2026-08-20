#!/usr/bin/env python3

import argparse
import json
from pathlib import Path

from functions.io_utils import configure_runtime_environment, make_timestamped_output_dir

PROJECT_DIR = Path(__file__).resolve().parents[4]
configure_runtime_environment(PROJECT_DIR, include_celloracle_cache=False)

import matplotlib

matplotlib.use("Agg")
from scipy import sparse

from functions.embedding_utils import (
    cluster_connectivity, estimate_grid, project_vectors_to_refined,
    refined_forceatlas, stratified_sample,
)
from functions.io_utils import read_refined_inputs
from functions.perturbation_utils import group_summary
from functions.plotting_utils import (
    arrow_plot, delta_map, distribution_plot, flow_plot, plot_propagation_magnitude,
)


CLUSTER_COL = "seurat_clusters"
GROUP_COL = "group"
GROUPS = ("CD", "Healthy")


def parse_args():
    parser = argparse.ArgumentParser(description="Refined, visualization-only analysis of the CellOracle HIF1A KO results.")
    parser.add_argument("--source-dir", default="results/celloracle_hif1a_ko_1500_FA_paperstyle_20260626_154852")
    parser.add_argument("--out-prefix", default="results/celloracle_hif1a_ko_1500_FA_refined")
    parser.add_argument("--target-cluster", default="4")
    parser.add_argument("--related-clusters", type=int, default=3)
    parser.add_argument("--knn", type=int, default=30)
    parser.add_argument("--arrow-fraction", type=float, default=0.075)
    parser.add_argument("--fa-iterations", type=int, default=800)
    return parser.parse_args()


def main():
    args = parse_args()
    source = Path(args.source_dir)
    out = make_timestamped_output_dir(args.out_prefix, ("figures", "tables", "logs"))
    figures, tables, logs = out / "figures", out / "tables", out / "logs"
    data = read_refined_inputs(source, CLUSTER_COL, GROUP_COL)
    embedding, vectors = data["embedding"].copy(), data["vectors"].copy()
    grid, propagation, score = data["grid"].copy(), data["propagation"].copy(), data["score"].copy()
    if not embedding["cell"].equals(vectors["cell"]):
        vectors = embedding[["cell"]].merge(vectors, on="cell", validate="one_to_one")

    connectivity = cluster_connectivity(embedding, args.target_cluster, args.knn, CLUSTER_COL)
    related = connectivity.head(args.related_clusters)["related_cluster"].tolist()
    included = [args.target_cluster] + related
    refined, graph = refined_forceatlas(embedding, included, args.knn, args.fa_iterations, CLUSTER_COL)
    refined_vectors = project_vectors_to_refined(vectors, refined)
    sampled = stratified_sample(vectors, args.arrow_fraction, CLUSTER_COL)
    sampled_refined = stratified_sample(refined_vectors, args.arrow_fraction, CLUSTER_COL)
    grid_clean = grid[grid["n_cells"] >= 5].copy()
    cluster_group = group_summary(embedding, CLUSTER_COL, GROUP_COL)

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

    delta_map(embedding, figures / "FA_delta_l2_clean", "Predicted transcriptomic impact of HIF1A knockout across all cells")
    delta_map(embedding, figures / "FA_delta_l2_cluster4_highlight", "HIF1A knockout impact concentrated in Seurat cluster 4", highlight=args.target_cluster)
    arrow_plot(vectors, sampled, figures / "FA_cell_shift_arrows_clean_sampled", "Predicted cell-state displacement after HIF1A knockout (stratified sample)")
    flow_plot(embedding, grid, figures / "FA_flow_on_grid_clean", "Smoothed HIF1A knockout response field on the global FA embedding")
    distribution_plot(embedding, figures / "delta_l2_by_cluster_group_clean", CLUSTER_COL, GROUP_COL)
    plot_propagation_magnitude(
        propagation, CLUSTER_COL, figures / "propagation_magnitude_clean", target_cluster=args.target_cluster,
        title="Propagation of the predicted HIF1A knockout response across clusters",
    )

    delta_map(refined, figures / "FA_refined_subset_delta_l2", f"HIF1A knockout impact in cluster 4 and connected clusters {', '.join(related)}", x="FA1_refined", y="FA2_refined")
    delta_map(refined, figures / "FA_refined_subset_delta_l2_cluster4_highlight", "Refined view of the predicted HIF1A knockout impact in cluster 4", x="FA1_refined", y="FA2_refined", highlight=args.target_cluster)
    arrow_plot(refined_vectors, sampled_refined, figures / "FA_refined_subset_cell_shift_arrows", "Projected HIF1A knockout displacement in the cluster 4 neighborhood", x="FA1_refined", y="FA2_refined", u="shift_FA1_refined", v="shift_FA2_refined")
    refined_grid = estimate_grid(refined_vectors, "FA1_refined", "FA2_refined", "shift_FA1_refined", "shift_FA2_refined")
    refined_grid_kept = flow_plot(
        refined, refined_grid, figures / "FA_refined_subset_flow_on_grid",
        "Smoothed HIF1A knockout response field in the cluster 4 neighborhood",
        x="FA1_refined", y="FA2_refined", gx="grid_FA1_refined", gy="grid_FA2_refined",
        gu="flow_FA1_refined", gv="flow_FA2_refined",
    )
    refined_grid_kept.to_csv(tables / "cluster4_related_FA_refined_grid_ncells_ge5.csv", index=False)

    for group in GROUPS:
        group_vectors = vectors[vectors[GROUP_COL] == group].copy()
        group_cells = embedding[embedding[GROUP_COL] == group].copy()
        group_sample = stratified_sample(group_vectors, args.arrow_fraction, CLUSTER_COL)
        group_grid = estimate_grid(group_vectors, "FA1", "FA2", "shift_FA1", "shift_FA2")
        kept = flow_plot(group_cells, group_grid, figures / f"FA_{group}_flow_on_grid_clean", f"Smoothed HIF1A knockout response field in {group} cells")
        delta_map(group_cells, figures / f"FA_{group}_delta_l2_clean", f"Predicted HIF1A knockout impact in {group} cells")
        arrow_plot(group_vectors, group_sample, figures / f"FA_{group}_cell_shift_arrows_clean_sampled", f"Predicted cell-state displacement after HIF1A knockout in {group}")
        group_sample.to_csv(tables / f"{group}_cell_shift_vectors_FA_sampled.csv", index=False)
        kept.to_csv(tables / f"{group}_grid_shift_vectors_FA_ncells_ge5.csv", index=False)

        refined_group_vectors = refined_vectors[refined_vectors[GROUP_COL] == group].copy()
        refined_group_cells = refined[refined[GROUP_COL] == group].copy()
        refined_group_sample = stratified_sample(refined_group_vectors, args.arrow_fraction, CLUSTER_COL)
        refined_group_grid = estimate_grid(refined_group_vectors, "FA1_refined", "FA2_refined", "shift_FA1_refined", "shift_FA2_refined")
        refined_kept = flow_plot(
            refined_group_cells, refined_group_grid, figures / f"FA_refined_subset_{group}_flow_on_grid",
            f"Cluster 4 neighborhood response field in {group} cells",
            x="FA1_refined", y="FA2_refined", gx="grid_FA1_refined", gy="grid_FA2_refined",
            gu="flow_FA1_refined", gv="flow_FA2_refined",
        )
        delta_map(refined_group_cells, figures / f"FA_refined_subset_{group}_delta_l2", f"HIF1A knockout impact in the cluster 4 neighborhood: {group}", x="FA1_refined", y="FA2_refined")
        arrow_plot(refined_group_vectors, refined_group_sample, figures / f"FA_refined_subset_{group}_cell_shift_arrows", f"Projected HIF1A knockout displacement near cluster 4: {group}", x="FA1_refined", y="FA2_refined", u="shift_FA1_refined", v="shift_FA2_refined")
        refined_group_sample.to_csv(tables / f"cluster4_related_{group}_FA_refined_vectors_sampled.csv", index=False)
        refined_kept.to_csv(tables / f"cluster4_related_{group}_FA_refined_grid_ncells_ge5.csv", index=False)

    means = embedding.groupby(CLUSTER_COL)["delta_l2"].mean().sort_values(ascending=False)
    cluster4_group = embedding[embedding[CLUSTER_COL] == args.target_cluster].groupby(GROUP_COL)["delta_l2"].mean()
    summary = {
        "output_dir": str(out.resolve()), "source_dir": str(source.resolve()),
        "n_cells_global": int(len(embedding)), "n_cells_refined_subset": int(len(refined)),
        "refined_subset_clusters": included,
        "related_cluster_selection": "Top kNN cross-cluster edge counts to cluster 4 in robust-scaled global FA",
        "cluster_with_highest_mean_delta_l2": str(means.index[0]), "highest_mean_delta_l2": float(means.iloc[0]),
        "cluster4_mean_delta_l2_CD": float(cluster4_group["CD"]),
        "cluster4_mean_delta_l2_Healthy": float(cluster4_group["Healthy"]),
        "cluster4_mean_delta_l2_CD_minus_Healthy": float(cluster4_group["CD"] - cluster4_group["Healthy"]),
        "arrow_sample_fraction_per_cluster": args.arrow_fraction, "minimum_grid_cells": 5,
        "figures_png": sorted(path.name for path in figures.glob("*.png")),
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
