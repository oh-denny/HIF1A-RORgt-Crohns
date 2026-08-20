"""Embedding, vector-projection, and grid-flow helpers."""

import numpy as np
import pandas as pd
from scipy import sparse
from sklearn.neighbors import NearestNeighbors
from sklearn.preprocessing import RobustScaler


FA_KEYS = ["X_draw_graph_fa", "X_fa", "X_forceatlas2", "X_force_directed", "FA", "FA2"]


def get_or_make_fa(adata, logs_dir):
    """Recover an existing ForceAtlas embedding or calculate a deterministic one."""
    import scanpy as sc

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


def make_embedding_table(oracle, metrics, cluster_column, group_column):
    adata = oracle.adata
    fa = np.asarray(adata.obsm["X_fa"])
    table = pd.DataFrame(
        {
            "cell": adata.obs_names,
            "FA1": fa[:, 0],
            "FA2": fa[:, 1],
            cluster_column: adata.obs[cluster_column].astype(str).values,
            group_column: adata.obs[group_column].astype(str).values if group_column in adata.obs else "NA",
        }
    )
    metrics = metrics.rename(columns={"cell_id": "cell"})
    return table.merge(
        metrics[["cell", "delta_l2", "target_imputed", "target_simulated"]],
        on="cell",
        how="left",
    )


def make_cell_vectors(oracle, embedding, cluster_column, group_column):
    shift = np.asarray(oracle.delta_embedding)
    out = embedding[["cell", "FA1", "FA2", "delta_l2", cluster_column, group_column]].copy()
    out["shift_FA1"] = shift[:, 0]
    out["shift_FA2"] = shift[:, 1]
    out["FA1_end"] = out["FA1"] + out["shift_FA1"]
    out["FA2_end"] = out["FA2"] + out["shift_FA2"]
    return out[
        ["cell", "FA1", "FA2", "shift_FA1", "shift_FA2", "FA1_end", "FA2_end", "delta_l2", cluster_column, group_column]
    ]


def make_grid_vectors(oracle):
    grid = np.asarray(oracle.flow_grid)
    flow = np.asarray(oracle.flow)
    mass = np.asarray(oracle.total_p_mass)
    fa = np.asarray(oracle.adata.obsm["X_fa"])
    step = np.median(np.diff(np.unique(np.round(grid[:, 0], 8))))
    if not np.isfinite(step) or step <= 0:
        step = 0.5
    neighbors = NearestNeighbors(radius=step * 1.5).fit(fa)
    n_cells = np.array([len(x) for x in neighbors.radius_neighbors(grid, return_distance=False)])
    return pd.DataFrame(
        {
            "grid_FA1": grid[:, 0],
            "grid_FA2": grid[:, 1],
            "flow_FA1": flow[:, 0],
            "flow_FA2": flow[:, 1],
            "mass": mass,
            "n_cells": n_cells,
        }
    )


def robust_limits(values, padding=0.04):
    lo, hi = np.nanquantile(values, [0.002, 0.998])
    span = max(hi - lo, np.finfo(float).eps)
    return lo - span * padding, hi + span * padding


def arrow_scale(xy, uv, fraction=0.035):
    diag = np.linalg.norm(np.ptp(np.asarray(xy), axis=0))
    mag = np.linalg.norm(np.asarray(uv), axis=1)
    positive = mag[np.isfinite(mag) & (mag > 0)]
    ref = np.quantile(positive, 0.95) if len(positive) else 1.0
    return diag * fraction / ref


def stratified_sample(frame, fraction, cluster_column="seurat_clusters", seed=123):
    sampled = []
    for _, subset in frame.groupby(cluster_column, sort=True):
        n = min(len(subset), max(1, int(round(len(subset) * fraction))))
        sampled.append(subset.sample(n=n, random_state=seed))
    return pd.concat(sampled, ignore_index=True)


def cluster_connectivity(embedding, target, knn, cluster_column="seurat_clusters"):
    coords = RobustScaler().fit_transform(embedding[["FA1", "FA2"]])
    labels = embedding[cluster_column].to_numpy()
    indices = NearestNeighbors(n_neighbors=min(knn + 1, len(embedding))).fit(coords).kneighbors(return_distance=False)
    rows = []
    target_mask = labels == target
    target_center = np.median(coords[target_mask], axis=0)
    for cluster in sorted(set(labels), key=lambda value: int(value)):
        if cluster == target:
            continue
        cluster_mask = labels == cluster
        outgoing = sum(np.count_nonzero(labels[neighbors[1:]] == cluster) for neighbors in indices[target_mask])
        incoming = sum(np.count_nonzero(labels[neighbors[1:]] == target) for neighbors in indices[cluster_mask])
        rows.append(
            {
                "target_cluster": target,
                "related_cluster": cluster,
                "knn_edges_both_directions": int(outgoing + incoming),
                "edges_from_target": int(outgoing),
                "edges_to_target": int(incoming),
                "robust_centroid_distance": float(np.linalg.norm(target_center - np.median(coords[cluster_mask], axis=0))),
            }
        )
    return pd.DataFrame(rows).sort_values(
        ["knn_edges_both_directions", "robust_centroid_distance"], ascending=[False, True]
    )


def knn_graph(coords, knn):
    distances, indices = NearestNeighbors(n_neighbors=min(knn + 1, len(coords))).fit(coords).kneighbors()
    distances = distances[:, 1:]
    indices = indices[:, 1:]
    sigma = np.median(distances[distances > 0])
    weights = np.exp(-np.square(distances / max(sigma, np.finfo(float).eps)))
    rows = np.repeat(np.arange(len(coords)), indices.shape[1])
    graph = sparse.csr_matrix((weights.ravel(), (rows, indices.ravel())), shape=(len(coords), len(coords)))
    return graph.maximum(graph.T)


def refined_forceatlas(embedding, included, knn, iterations, cluster_column="seurat_clusters"):
    from fa2_modified import ForceAtlas2

    subset = embedding[embedding[cluster_column].isin(included)].copy().reset_index(drop=True)
    initial = RobustScaler().fit_transform(subset[["FA1", "FA2"]])
    graph = knn_graph(initial, knn)
    forceatlas = ForceAtlas2(
        outboundAttractionDistribution=False, linLogMode=False, adjustSizes=False,
        edgeWeightInfluence=1.0, jitterTolerance=1.0, barnesHutOptimize=True,
        barnesHutTheta=1.2, scalingRatio=2.0, strongGravityMode=True,
        gravity=1.0, verbose=False,
    )
    refined = np.asarray(forceatlas.forceatlas2(graph, pos=initial, iterations=iterations))
    refined = RobustScaler().fit_transform(refined)
    subset["FA1_refined"] = refined[:, 0]
    subset["FA2_refined"] = refined[:, 1]
    return subset, graph


def project_vectors_to_refined(vectors, refined, knn=20):
    merged = refined[["cell", "FA1_refined", "FA2_refined"]].merge(
        vectors, on="cell", how="left", validate="one_to_one"
    )
    old_xy = merged[["FA1", "FA2"]].to_numpy()
    new_xy = merged[["FA1_refined", "FA2_refined"]].to_numpy()
    old_scaled = RobustScaler().fit_transform(old_xy)
    indices = NearestNeighbors(n_neighbors=min(knn, len(merged))).fit(old_scaled).kneighbors(return_distance=False)
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
    distances, indices = NearestNeighbors(n_neighbors=min(neighbors, len(xy))).fit(xy).kneighbors(points)
    bandwidth = np.median(distances[:, -1])
    weights = np.exp(-0.5 * np.square(distances / max(bandwidth, 1e-12)))
    flow = (uv[indices] * weights[:, :, None]).sum(axis=1)
    flow /= np.maximum(weights.sum(axis=1)[:, None], 1e-12)
    radius = np.median(np.diff(x_values)) * 1.5
    counts = np.array(
        [len(found) for found in NearestNeighbors(radius=radius).fit(xy).radius_neighbors(points, return_distance=False)]
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
