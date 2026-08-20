#!/usr/bin/env python3

import argparse
import json
from pathlib import Path

from functions.io_utils import configure_runtime_environment, make_timestamped_output_dir

PROJECT_DIR = Path(__file__).resolve().parents[4]
configure_runtime_environment(PROJECT_DIR, include_celloracle_cache=True)

import matplotlib

matplotlib.use("Agg")
import numpy as np

from functions.celloracle_utils import load_validated_markov_oracle
from functions.perturbation_utils import (
    assign_grid_metadata, calculate_celloracle_ps, cluster_ps_summary,
    make_markov_density_table, most_enriched_cluster, run_markov,
    strongest_ps_clusters, transition_tables,
)
from functions.plotting_utils import plot_combined, plot_groups, plot_single_panels


CLUSTER_COLUMN = "seurat_clusters"
GROUP_COLUMN = "group"


def parse_args():
    parser = argparse.ArgumentParser(description="Generate CellOracle perturbation-score and Markov panels for HIF1A KO.")
    parser.add_argument("--oracle", default="results/celloracle_hif1a_ko_1500/objects/oracle_HIF1A_KO.celloracle.oracle")
    parser.add_argument("--out-prefix", default="results/celloracle_hif1a_ko_1500_PS_markov")
    parser.add_argument("--pseudotime-column", default="Pseudotime")
    parser.add_argument("--n-grid", type=int, default=40)
    parser.add_argument("--grid-neighbors", type=int, default=200)
    parser.add_argument("--min-mass", type=float, default=0.01)
    parser.add_argument("--markov-steps", type=int, default=200)
    parser.add_argument("--markov-duplications", type=int, default=5)
    parser.add_argument("--n-jobs", type=int, default=1)
    return parser.parse_args()


def main():
    args = parse_args()
    out = make_timestamped_output_dir(args.out_prefix, ("figures", "tables", "objects", "logs"))
    figures, tables, objects, logs = out / "figures", out / "tables", out / "objects", out / "logs"

    oracle = load_validated_markov_oracle(args.oracle, CLUSTER_COLUMN, GROUP_COLUMN)
    gradient, development, ps_grid = calculate_celloracle_ps(oracle, args)
    ps_grid = assign_grid_metadata(ps_grid, oracle, CLUSTER_COLUMN, GROUP_COLUMN)
    ps_summary = cluster_ps_summary(ps_grid, CLUSTER_COLUMN)
    markov_results, trajectories, random_trajectories = run_markov(oracle, args, CLUSTER_COLUMN, GROUP_COLUMN)
    density, bandwidth = make_markov_density_table(oracle, markov_results, CLUSTER_COLUMN, GROUP_COLUMN)
    transition_counts, transition_proportions, group_cluster = transition_tables(markov_results)

    ps_grid.to_csv(tables / "HIF1A_KO_perturbation_score_grid.csv", index=False)
    ps_summary.to_csv(tables / "HIF1A_KO_perturbation_score_by_cluster.csv", index=False)
    markov_results.to_csv(tables / "HIF1A_KO_markov_start_end_cells.csv", index=False)
    density.to_csv(tables / "HIF1A_KO_markov_density_by_cell.csv", index=False)
    transition_counts.to_csv(tables / "HIF1A_KO_markov_transition_counts.csv")
    transition_proportions.to_csv(tables / "HIF1A_KO_markov_transition_proportions.csv")
    group_cluster.to_csv(tables / "HIF1A_KO_markov_endpoint_cluster_by_start_group.csv", index=False)
    np.savez_compressed(
        objects / "HIF1A_KO_markov_trajectories.npz",
        trajectories=trajectories, randomized_trajectories=random_trajectories,
        selected_cell_indices=np.asarray(oracle.ixs_markvov_simulation),
    )
    gradient.to_hdf5(str(objects / "HIF1A_KO_pseudotime_gradient.celloracle.gradient"))

    plot_single_panels(oracle, ps_grid, density, markov_results, figures)
    plot_combined(oracle, ps_grid, density, markov_results, figures, args.markov_steps)
    plot_groups(oracle, ps_grid, density, markov_results, figures, args.markov_steps)

    negative_cluster, positive_cluster = strongest_ps_clusters(ps_summary, CLUSTER_COLUMN)
    enriched_cluster, enrichment = most_enriched_cluster(markov_results)
    summary = {
        "output_dir": str(out.resolve()), "source_oracle": str(Path(args.oracle).resolve()),
        "embedding": oracle.embedding_name, "n_cells": int(oracle.adata.n_obs),
        "n_markov_start_cells": int(len(oracle.ixs_markvov_simulation)),
        "n_markov_trajectories": int(len(markov_results)),
        "markov_steps": args.markov_steps, "markov_duplications": args.markov_duplications,
        "kde_bandwidth": bandwidth,
        "cluster_with_most_negative_mean_PS": negative_cluster,
        "cluster_with_most_positive_mean_PS": positive_cluster,
        "cluster_with_largest_markov_endpoint_enrichment": enriched_cluster,
        "largest_markov_endpoint_enrichment_ratio": float(enrichment.loc[enriched_cluster]),
        "groups": sorted(oracle.adata.obs[GROUP_COLUMN].astype(str).unique()),
        "figures": sorted(path.name for path in figures.glob("*")),
        "tables": sorted(path.name for path in tables.glob("*")),
        "note": "CellOracle marks Markov simulation as deprecated in favor of perturbation-score analysis; both are exported here to reproduce the requested paper-style comparison.",
    }
    (logs / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    (logs / "run_parameters.json").write_text(json.dumps(vars(args), indent=2) + "\n")
    (logs / "methodology.txt").write_text(
        "Perturbation score (PS) was calculated with CellOracle's Gradient_calculator and "
        "Oracle_development_module. The score is the inner product between the HIF1A KO "
        "simulation flow and the pseudotime gradient on a 40 x 40 UMAP grid. Negative values "
        "indicate opposition to pseudotime progression; positive values indicate alignment.\n\n"
        "Markov trajectories were generated with Oracle.run_markov_chain_simulation from the "
        "existing HIF1A KO transition probabilities. The plotted endpoint density is a Gaussian "
        "KDE evaluated at observed cells and scaled from 0 to 1 after clipping at the 1st and "
        "99th percentiles. Group-specific panels condition trajectories by the group of their starting cell.\n"
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
