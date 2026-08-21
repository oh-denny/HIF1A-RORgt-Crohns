"""Embedding, vector-projection, and grid-flow helpers."""

import numpy as np
import pandas as pd
from sklearn.neighbors import NearestNeighbors


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


def arrow_scale(xy, uv, fraction=0.035):
    diag = np.linalg.norm(np.ptp(np.asarray(xy), axis=0))
    mag = np.linalg.norm(np.asarray(uv), axis=1)
    positive = mag[np.isfinite(mag) & (mag > 0)]
    ref = np.quantile(positive, 0.95) if len(positive) else 1.0
    return diag * fraction / ref
