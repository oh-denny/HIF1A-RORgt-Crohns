"""CellOracle object loading, validation, and export helpers."""

import os
from pathlib import Path

import numpy as np
import pandas as pd


def patch_gimmemotifs_compatibility():
    """Apply import-time compatibility shims required by CellOracle 0.20."""
    cache_root = os.environ.get("XDG_CACHE_HOME")
    if cache_root:
        try:
            import appdirs

            def pipeline_cache_dir(appname=None, appauthor=None, version=None, opinion=True):
                path = Path(cache_root)
                if appname:
                    path /= appname
                if version:
                    path /= version
                return str(path)

            appdirs.user_cache_dir = pipeline_cache_dir
        except Exception:
            pass
    try:
        import gimmemotifs.motif as gimmemotifs_motif

        if not hasattr(gimmemotifs_motif, "default_motifs"):
            gimmemotifs_motif.default_motifs = lambda *args, **kwargs: []
    except Exception:
        pass


def load_oracle(path):
    patch_gimmemotifs_compatibility()
    import celloracle as co

    return co.load_hdf5(str(path))


def load_validated_markov_oracle(path, cluster_column="seurat_clusters", group_column="group"):
    oracle = load_oracle(path)
    missing_obs = {cluster_column, group_column}.difference(oracle.adata.obs.columns)
    if missing_obs:
        raise KeyError(f"Missing required cell metadata: {sorted(missing_obs)}")
    if oracle.embedding_name != "X_umap":
        if "X_umap" not in oracle.adata.obsm:
            raise KeyError("X_umap is not available in the Oracle object.")
        raise ValueError("The input Oracle must contain the HIF1A KO simulation calculated on X_umap.")
    required = ("delta_embedding", "delta_embedding_random", "transition_prob", "transition_prob_random")
    missing = [name for name in required if not hasattr(oracle, name)]
    if missing:
        raise AttributeError(f"Oracle is missing simulation attributes: {missing}")
    return oracle


def resolve_gene(gene, genes):
    gene = gene.strip()
    for candidate in (gene, gene.upper(), gene.capitalize()):
        if candidate in genes:
            return candidate
    raise ValueError(f"Target gene {gene!r} was not found in exported genes.")


def load_export(export_dir, cluster_column):
    import anndata as ad
    from scipy import io

    export_dir = Path(export_dir)
    counts = io.mmread(export_dir / "counts_genes_x_cells.mtx").tocsr()
    genes = pd.read_csv(export_dir / "genes.tsv", header=None)[0].astype(str).tolist()
    cells = pd.read_csv(export_dir / "cells.tsv", header=None)[0].astype(str).tolist()
    metadata = pd.read_csv(export_dir / "metadata.csv").set_index("cell_id").loc[cells]
    umap = pd.read_csv(export_dir / "umap.csv").set_index("cell_id").loc[cells]
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
            return pd.read_csv(path, sep="\t" if path.suffix in {".txt", ".tsv"} else ",")
        if path.suffix in {".pkl", ".pickle"}:
            return pd.read_pickle(path)
        raise ValueError(f"Unsupported base GRN format: {path.suffix}")
    patch_gimmemotifs_compatibility()
    import celloracle as co

    version = f"{args.genome}_gimmemotifsv5_fpr2"
    if args.species == "human":
        return co.data.load_human_promoter_base_GRN(version=version)
    if args.species == "mouse":
        return co.data.load_mouse_promoter_base_GRN(version=version)
    raise ValueError(f"Unsupported species: {args.species}")


def save_links_tables(links, out_tables):
    raw_parts = []
    filtered_parts = []
    for cluster, frame in links.links_dict.items():
        part = frame.copy()
        part["cluster"] = cluster
        raw_parts.append(part)
    if hasattr(links, "filtered_links"):
        for cluster, frame in links.filtered_links.items():
            part = frame.copy()
            part["cluster"] = cluster
            filtered_parts.append(part)
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
    summary = pd.DataFrame(
        {
            "gene": delta.columns,
            "mean_delta": delta.mean(axis=0).values,
            "median_delta": delta.median(axis=0).values,
            "mean_abs_delta": delta.abs().mean(axis=0).values,
            "mean_imputed": imputed.mean(axis=0).values,
            "mean_simulated": simulated.mean(axis=0).values,
        }
    )
    summary["abs_rank"] = summary["mean_abs_delta"].rank(ascending=False, method="dense")
    summary.sort_values("mean_abs_delta", ascending=False).to_csv(out_tables / "gene_delta_summary.csv", index=False)
    obs = oracle.adata.obs.copy()
    cell_metrics = pd.DataFrame(
        {
            "cell_id": oracle.adata.obs_names,
            "delta_l2": np.linalg.norm(delta.values, ord=2, axis=1),
            "target_imputed": imputed[target_gene].values,
            "target_simulated": simulated[target_gene].values,
            cluster_column: obs[cluster_column].astype(str).values,
        }
    )
    if group_column in obs.columns:
        cell_metrics[group_column] = obs[group_column].astype(str).values
    cell_metrics.to_csv(out_tables / "cell_delta_metrics.csv", index=False)
    group_columns = [cluster_column]
    if group_column in cell_metrics.columns:
        group_columns.append(group_column)
    cluster_summary = (
        cell_metrics.groupby(group_columns, observed=True)
        .agg(
            cells=("cell_id", "size"),
            mean_delta_l2=("delta_l2", "mean"),
            median_delta_l2=("delta_l2", "median"),
            mean_target_imputed=("target_imputed", "mean"),
        )
        .reset_index()
    )
    cluster_summary.to_csv(out_tables / "cluster_delta_summary.csv", index=False)
