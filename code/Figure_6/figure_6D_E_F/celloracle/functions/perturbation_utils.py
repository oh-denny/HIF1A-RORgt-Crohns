"""CellOracle perturbation, propagation, score, and Markov calculations."""

import warnings

import numpy as np
import pandas as pd
from sklearn.neighbors import KernelDensity, NearestNeighbors


def recalculate_fa_shift_and_grid(oracle, target_gene, n_jobs, grid_steps, grid_neighbors):
    oracle.embedding_name = "X_fa"
    oracle.embedding = np.asarray(oracle.adata.obsm["X_fa"])
    oracle.simulate_shift(
        perturb_condition={target_gene: 0.0}, GRN_unit="cluster", n_propagation=3, ignore_warning=True
    )
    oracle.estimate_transition_prob(
        n_neighbors=200, knn_random=True, sampled_fraction=1.0,
        n_jobs=n_jobs, calculate_randomized=True,
    )
    oracle.calculate_embedding_shift(sigma_corr=0.05)
    oracle.calculate_grid_arrows(
        smooth=0.8, steps=(grid_steps, grid_steps), n_neighbors=grid_neighbors, n_jobs=n_jobs
    )
    return oracle


def propagation_summary(oracle, target_gene, cluster_column, max_propagation=5):
    """Calculate real per-cluster delta-X magnitudes without curve normalization."""
    rows = []
    clusters = oracle.adata.obs[cluster_column].astype(str).values
    for n_propagation in range(max_propagation + 1):
        if n_propagation == 0 and hasattr(oracle, "_Oracle__simulate_shift"):
            oracle._Oracle__simulate_shift(
                perturb_condition={target_gene: 0.0}, GRN_unit="cluster",
                n_propagation=0, ignore_warning=True, n_min=0, n_max=6,
            )
        else:
            oracle.simulate_shift(
                perturb_condition={target_gene: 0.0}, GRN_unit="cluster",
                n_propagation=n_propagation, ignore_warning=True,
            )
        l2 = np.linalg.norm(oracle.adata.to_df(layer="delta_X").values, axis=1)
        frame = pd.DataFrame({cluster_column: clusters, "l2": l2})
        for cluster, subset in frame.groupby(cluster_column):
            rows.append(
                {
                    cluster_column: cluster,
                    "n_propagation": n_propagation,
                    "mean_l2_norm": subset["l2"].mean(),
                    "median_l2_norm": subset["l2"].median(),
                    "sd_l2_norm": subset["l2"].std(),
                    "n_cells": len(subset),
                }
            )
    return pd.DataFrame(rows)


def celloracle_impact_summary(oracle, target_gene, max_propagation=5, norm_order=1):
    """Reproduce CellOracle's built-in propagation-impact values as a tidy table."""
    lengths = []
    for n_propagation in range(max_propagation + 1):
        oracle._Oracle__simulate_shift(
            perturb_condition={target_gene: 0.0},
            GRN_unit=None,
            n_propagation=n_propagation,
            ignore_warning=False,
            use_randomized_GRN=False,
            n_min=0,
            n_max=max_propagation + 1,
        )
        delta = oracle.adata.to_df(layer="delta_X")
        lengths.append(np.linalg.norm(delta, ord=norm_order, axis=1))

    values = pd.DataFrame(lengths).transpose()
    values.columns = range(max_propagation + 1)
    values["cluster"] = oracle.adata.obs[oracle.cluster_column_name].astype(str).values
    summary = values.groupby("cluster").mean().reset_index()
    return summary.melt(
        id_vars="cluster",
        var_name="n_propagation",
        value_name="mean_delta_x_length",
    )


def find_pseudotime_column(adata):
    candidates = ["pseudotime", "dpt_pseudotime", "monocle_pseudotime", "slingshot_pseudotime", "palantir_pseudotime"]
    lower_map = {column.lower(): column for column in adata.obs.columns}
    for candidate in candidates:
        if candidate.lower() in lower_map:
            return lower_map[candidate.lower()]
    for column in adata.obs.columns:
        if "pseudotime" in column.lower():
            return column
    return None


def perturbation_score(oracle, cell_vectors, cluster_column, group_column, logs_dir):
    pseudotime_column = find_pseudotime_column(oracle.adata)
    if pseudotime_column is None:
        (logs_dir / "perturbation_score_not_generated.txt").write_text(
            "Perturbation score was not generated because no pseudotime column was found in adata.obs.\n"
        )
        return None, False, None
    pseudotime = pd.to_numeric(oracle.adata.obs[pseudotime_column], errors="coerce").to_numpy()
    fa = np.asarray(oracle.adata.obsm["X_fa"])
    _, indices = NearestNeighbors(n_neighbors=min(31, len(fa))).fit(fa).kneighbors(fa)
    trajectory = np.zeros_like(fa)
    for i, neighbors in enumerate(indices):
        forward = neighbors[
            np.isfinite(pseudotime[neighbors]) & np.isfinite(pseudotime[i]) & (pseudotime[neighbors] > pseudotime[i])
        ]
        forward = neighbors[1:min(6, len(neighbors))] if len(forward) == 0 else forward[:5]
        trajectory[i, :] = fa[forward].mean(axis=0) - fa[i]
    shift = cell_vectors[["shift_FA1", "shift_FA2"]].to_numpy()
    denominator = np.linalg.norm(shift, axis=1) * np.linalg.norm(trajectory, axis=1)
    score = np.divide(
        (shift * trajectory).sum(axis=1), denominator,
        out=np.full(len(denominator), np.nan), where=denominator > 0,
    )
    output = pd.DataFrame(
        {
            "cell": oracle.adata.obs_names,
            "FA1": fa[:, 0], "FA2": fa[:, 1],
            "shift_FA1": shift[:, 0], "shift_FA2": shift[:, 1],
            "trajectory_FA1": trajectory[:, 0], "trajectory_FA2": trajectory[:, 1],
            "perturbation_score": score,
            cluster_column: oracle.adata.obs[cluster_column].astype(str).values,
            group_column: oracle.adata.obs[group_column].astype(str).values if group_column in oracle.adata.obs else "NA",
        }
    )
    return output, True, pseudotime_column


def calculate_celloracle_ps(oracle, args):
    from celloracle.applications import Gradient_calculator, Oracle_development_module

    gradient = Gradient_calculator(oracle_object=oracle, pseudotime_key=args.pseudotime_column, name="HIF1A_KO")
    gradient.calculate_p_mass(smooth=0.8, n_grid=args.n_grid, n_neighbors=args.grid_neighbors, n_jobs=args.n_jobs)
    gradient.calculate_mass_filter(min_mass=args.min_mass)
    gradient.transfer_data_into_grid(args={"method": "knn", "n_knn": 30})
    gradient.calculate_gradient(scale_factor="l2_norm_mean", normalization="sqrt")
    development = Oracle_development_module(gradient_object=gradient, oracle_object=oracle, name="HIF1A_KO")
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


def assign_grid_metadata(grid, oracle, cluster_column="seurat_clusters", group_column="group"):
    indices = NearestNeighbors(n_neighbors=1).fit(np.asarray(oracle.embedding)).kneighbors(
        grid[["UMAP1", "UMAP2"]].to_numpy(), return_distance=False
    )[:, 0]
    obs = oracle.adata.obs
    grid[cluster_column] = obs.iloc[indices][cluster_column].astype(str).to_numpy()
    grid[group_column] = obs.iloc[indices][group_column].astype(str).to_numpy()
    return grid


def run_markov(oracle, args, cluster_column="seurat_clusters", group_column="group", seed=123):
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message="Functions for Markov simulation are deprecated.*", category=DeprecationWarning)
        oracle.run_markov_chain_simulation(
            n_steps=args.markov_steps, n_duplication=args.markov_duplications,
            seed=seed, calculate_randomized=True,
        )
    trajectories = oracle.markvov_transition_id
    randomized = oracle.markvov_transition_random_id
    start = trajectories.iloc[:, 0].to_numpy(dtype=int)
    end = trajectories.iloc[:, -1].to_numpy(dtype=int)
    end_random = randomized.iloc[:, -1].to_numpy(dtype=int)
    obs = oracle.adata.obs.reset_index(names="cell")
    result = pd.DataFrame(
        {
            "start_cell_index": start, "end_cell_index": end,
            "end_randomized_cell_index": end_random,
            "start_cell": obs.iloc[start]["cell"].to_numpy(),
            "end_cell": obs.iloc[end]["cell"].to_numpy(),
            "end_randomized_cell": obs.iloc[end_random]["cell"].to_numpy(),
            "start_cluster": obs.iloc[start][cluster_column].astype(str).to_numpy(),
            "end_cluster": obs.iloc[end][cluster_column].astype(str).to_numpy(),
            "end_randomized_cluster": obs.iloc[end_random][cluster_column].astype(str).to_numpy(),
            "start_group": obs.iloc[start][group_column].astype(str).to_numpy(),
            "end_group": obs.iloc[end][group_column].astype(str).to_numpy(),
        }
    )
    return result, trajectories.to_numpy(dtype=np.int32), randomized.to_numpy(dtype=np.int32)


def kde_bandwidth(embedding):
    distances, _ = NearestNeighbors(n_neighbors=min(31, len(embedding))).fit(embedding).kneighbors()
    return max(float(np.median(distances[:, -1])), np.finfo(float).eps)


def normalized_density(embedding, endpoint_indices, bandwidth):
    endpoints = embedding[np.asarray(endpoint_indices, dtype=int)]
    log_density = KernelDensity(kernel="gaussian", bandwidth=bandwidth).fit(endpoints).score_samples(embedding)
    density = np.exp(log_density - np.nanmax(log_density))
    low, high = np.nanquantile(density, [0.01, 0.99])
    return np.clip((density - low) / max(high - low, 1e-12), 0, 1)


def make_markov_density_table(oracle, results, cluster_column="seurat_clusters", group_column="group"):
    embedding = np.asarray(oracle.embedding)
    bandwidth = kde_bandwidth(embedding)
    obs = oracle.adata.obs.reset_index(names="cell")
    output = pd.DataFrame(
        {
            "cell": obs["cell"], "UMAP1": embedding[:, 0], "UMAP2": embedding[:, 1],
            cluster_column: obs[cluster_column].astype(str), group_column: obs[group_column].astype(str),
        }
    )
    output["markov_density_all"] = normalized_density(embedding, results["end_cell_index"], bandwidth)
    output["markov_density_randomized"] = normalized_density(embedding, results["end_randomized_cell_index"], bandwidth)
    for group in ("CD", "Healthy"):
        group_end = results.loc[results["start_group"] == group, "end_cell_index"]
        output[f"markov_density_start_{group}"] = normalized_density(embedding, group_end, bandwidth)
    return output, bandwidth


def cluster_ps_summary(grid, cluster_column="seurat_clusters"):
    return (
        grid[~grid["mass_filtered"]].groupby(cluster_column)["perturbation_score"]
        .agg(n_grid_points="size", mean_PS="mean", median_PS="median")
        .reset_index().sort_values("mean_PS")
    )


def transition_tables(results):
    counts = pd.crosstab(results["start_cluster"], results["end_cluster"], rownames=["start_cluster"], colnames=["end_cluster"])
    proportions = counts.div(counts.sum(axis=1), axis=0)
    group_cluster = results.groupby(["start_group", "end_cluster"]).size().rename("n_trajectories").reset_index()
    group_cluster["proportion_within_start_group"] = group_cluster.groupby("start_group")["n_trajectories"].transform(lambda values: values / values.sum())
    return counts, proportions, group_cluster


def strongest_ps_clusters(summary, cluster_column="seurat_clusters"):
    return str(summary.iloc[0][cluster_column]), str(summary.iloc[-1][cluster_column])


def most_enriched_cluster(results):
    start = results["start_cluster"].value_counts(normalize=True)
    end = results["end_cluster"].value_counts(normalize=True)
    enrichment = (end / start).replace([np.inf, -np.inf], np.nan).dropna()
    return str(enrichment.idxmax()), enrichment
