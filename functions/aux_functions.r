

library(Seurat)
library(Matrix)
library(dplyr)
library(tibble)
library(DESeq2)
library(ggplot2)
library(fgsea)
library(msigdbr)
library(GSVA)
library(ComplexHeatmap)
library(circlize)
library(grid)




run_pb_deseq_cluster <- function(
    cluster_use,
    counts_mat,
    meta_all,
    min_cells = 9,
    outdir = "pseudobulk_clusters"
) {
  
  message("Rodando pseudobulk para cluster ", cluster_use)
  
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  
  cells_use <- meta_all %>%
    filter(cluster == cluster_use) %>%
    pull(cell)
  
  message("Número de células no cluster ", cluster_use, ": ", length(cells_use))
  
  counts_cl <- counts_mat[, cells_use, drop = FALSE]
  
  meta_cl <- meta_all %>%
    filter(cell %in% cells_use) %>%
    arrange(match(cell, colnames(counts_cl)))
  
  stopifnot(all(meta_cl$cell == colnames(counts_cl)))
  
  meta_cl <- meta_cl %>%
    mutate(
      pb_id = paste(sample_id, group, sep = "__")
    )
  
  n_cells_cl <- meta_cl %>%
    as.data.frame() %>%
    group_by(pb_id, sample_id, group) %>%
    summarise(
      n_cells = n(),
      .groups = "drop"
    ) %>%
    arrange(group, sample_id)
  
  print(n_cells_cl)
  print(table(n_cells_cl$group))
  
  mm <- Matrix::sparse.model.matrix(~ 0 + pb_id, data = meta_cl)
  colnames(mm) <- gsub("^pb_id", "", colnames(mm))
  
  pb_counts_cl <- counts_cl %*% mm
  
  pb_meta_cl <- n_cells_cl %>%
    filter(pb_id %in% colnames(pb_counts_cl)) %>%
    arrange(match(pb_id, colnames(pb_counts_cl)))
  
  pb_counts_cl <- pb_counts_cl[, pb_meta_cl$pb_id, drop = FALSE]
  
  stopifnot(all(colnames(pb_counts_cl) == pb_meta_cl$pb_id))
  
  pb_meta_cl_filt <- pb_meta_cl %>%
    filter(n_cells >= min_cells)
  
  pb_counts_cl_filt <- pb_counts_cl[, pb_meta_cl_filt$pb_id, drop = FALSE]
  
  message("Após filtro min_cells = ", min_cells)
  print(table(pb_meta_cl_filt$group))
  
  if (length(unique(pb_meta_cl_filt$group)) < 2) {
    stop("Cluster ", cluster_use, " não tem os dois grupos após filtro.")
  }
  
  if (any(table(pb_meta_cl_filt$group) < 2)) {
    warning("Cluster ", cluster_use, " tem menos de 2 amostras em algum grupo. Resultado será exploratório.")
  }
  
  counts_use <- pb_counts_cl_filt
  
  meta_use <- pb_meta_cl_filt %>%
    as.data.frame()
  
  rownames(meta_use) <- meta_use$pb_id
  meta_use$group <- factor(meta_use$group, levels = c("Healthy", "CD"))
  
  counts_use <- counts_use[, rownames(meta_use), drop = FALSE]
  
  stopifnot(all(colnames(counts_use) == rownames(meta_use)))
  
  keep_genes <- rowSums(counts_use >= 10) >= 2
  counts_use <- counts_use[keep_genes, ]
  
  dds <- DESeqDataSetFromMatrix(
    countData = round(as.matrix(counts_use)),
    colData = meta_use,
    design = ~ group
  )
  
  dds <- DESeq(dds)
  
  res <- results(
    dds,
    contrast = c("group", "CD", "Healthy")
  )
  
  de <- as.data.frame(res) %>%
    rownames_to_column("gene") %>%
    arrange(padj)
  
  write.csv(
    de,
    file = file.path(outdir, paste0("DESeq2_cluster", cluster_use, "_CD_vs_Healthy.csv")),
    row.names = FALSE
  )
  
  write.csv(
    pb_meta_cl_filt,
    file = file.path(outdir, paste0("metadata_cluster", cluster_use, ".csv")),
    row.names = FALSE
  )
  
  return(list(
    cluster = cluster_use,
    de = de,
    dds = dds,
    meta_use = meta_use,
    counts_use = counts_use,
    pb_meta = pb_meta_cl_filt,
    pb_counts = pb_counts_cl_filt
  ))
}






run_hallmark_fgsea <- function(de, cluster_use, outdir = "pseudobulk_clusters") {
  
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  
  ranks <- de %>%
    dplyr::filter(!is.na(stat)) %>%
    dplyr::arrange(desc(stat)) %>%
    dplyr::select(gene, stat) %>%
    tibble::deframe()
  
  # remover genes duplicados, se houver
  ranks <- ranks[!duplicated(names(ranks))]
  
  msig_hallmark <- msigdbr(
    species = "Homo sapiens",
    category = "H"
  )
  
  hallmark_sets <- split(
    msig_hallmark$gene_symbol,
    msig_hallmark$gs_name
  )
  
  fgsea_res <- fgsea(
    pathways = hallmark_sets,
    stats = ranks,
    minSize = 10,
    maxSize = 500,
    eps = 0
  ) %>%
    as.data.frame() %>%
    dplyr::arrange(padj)
  
  # criar versão salvável do resultado
  fgsea_to_save <- fgsea_res %>%
    dplyr::mutate(
      leadingEdge = sapply(leadingEdge, function(x) paste(x, collapse = ";"))
    )
  
  write.csv(
    fgsea_to_save,
    file = file.path(outdir, paste0("Hallmark_fgsea_cluster", cluster_use, "_CD_vs_Healthy.csv")),
    row.names = FALSE
  )
  
  return(list(
    ranks = ranks,
    hallmark_sets = hallmark_sets,
    fgsea = fgsea_res,
    fgsea_table = fgsea_to_save
  ))
}


fgsea_c5 <- run_hallmark_fgsea(
  de = res_c5$de,
  cluster_use = "5"
)

fgsea_c4 <- run_hallmark_fgsea(
  de = res_c4$de,
  cluster_use = "4"
)

fgsea_c5$fgsea %>%
  dplyr::select(pathway, NES, pval, padj, size) %>%
  head(20)

fgsea_c4$fgsea %>%
  dplyr::select(pathway, NES, pval, padj, size) %>%
  head(20)



plot_enrichment_cluster <- function(fgsea_obj, pathway_name, cluster_use) {
  
  p <- plotEnrichment(
    pathway = fgsea_obj$hallmark_sets[[pathway_name]],
    stats = fgsea_obj$ranks
  ) +
    labs(
      title = pathway_name,
      subtitle = paste0("Human cluster ", cluster_use, " (CD vs Healthy)"),
      x = "Genes ranked by DESeq2 statistic",
      y = "Running enrichment score"
    ) +
    theme_classic(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5)
    )
  
  return(p)
}











run_pseudobulk_ssgsea_heatmap = function(
  seurat_obj,
  clusters_use = c("4", "5", "6"),
  assay = "RNA",
  cluster_colname = "seurat_clusters",
  sample_colname = "sample",
  group_colname = "group",
  group_levels = c("Healthy", "CD"),
  min_cells = 9,
  min_genes_per_set = 10,
  pathways_use = c(
    "HALLMARK_HYPOXIA",
    "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
    "HALLMARK_INTERFERON_GAMMA_RESPONSE",
    "HALLMARK_IL2_STAT5_SIGNALING",
    "HALLMARK_IL6_JAK_STAT3_SIGNALING",
    "HALLMARK_INFLAMMATORY_RESPONSE",
    "HALLMARK_GLYCOLYSIS",
    "HALLMARK_MTORC1_SIGNALING",
    "HALLMARK_APOPTOSIS",
    "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
    "HALLMARK_CHOLESTEROL_HOMEOSTASIS"
  ),
  group_colors = c(
    "CD" = "#E9969B",
    "Healthy" = "#2D6F6B"
  ),
  cluster_colors = c(
    "4" = "#B2182B",
    "5" = "#2166AC",
    "6" = "#4D9221"
  ),
  zscore_limit = 2
) {

  DefaultAssay(seurat_obj) = assay

  seurat_pb = JoinLayers(
    object = seurat_obj,
    assay = assay
  )

  counts_mat = GetAssayData(
    seurat_pb,
    assay = assay,
    layer = "counts"
  )

  meta_all = seurat_pb@meta.data %>%
    as.data.frame() %>%
    tibble::rownames_to_column("cell") %>%
    dplyr::mutate(
      cluster = as.character(.data[[cluster_colname]]),
      sample_id = as.character(.data[[sample_colname]]),
      group = factor(
        .data[[group_colname]],
        levels = group_levels
      )
    )

  meta_selected = meta_all %>%
    dplyr::filter(cluster %in% clusters_use) %>%
    dplyr::mutate(
      pb_id = paste0(
        "C",
        cluster,
        "__",
        sample_id,
        "__",
        group
      )
    )

  counts_selected = counts_mat[, meta_selected$cell, drop = FALSE]

  meta_selected = meta_selected %>%
    dplyr::arrange(
      match(cell, colnames(counts_selected))
    )

  stopifnot(
    all(meta_selected$cell == colnames(counts_selected))
  )

  design_pb = Matrix::sparse.model.matrix(
    ~ 0 + pb_id,
    data = meta_selected
  )

  colnames(design_pb) = gsub(
    "^pb_id",
    "",
    colnames(design_pb)
  )

  pb_counts = counts_selected %*% design_pb

  pb_meta = meta_selected %>%
    dplyr::group_by(
      pb_id,
      cluster,
      sample_id,
      group
    ) %>%
    dplyr::summarise(
      n_cells = dplyr::n(),
      .groups = "drop"
    ) %>%
    dplyr::arrange(
      match(pb_id, colnames(pb_counts))
    )

  pb_counts = pb_counts[
    ,
    pb_meta$pb_id,
    drop = FALSE
  ]

  stopifnot(
    all(colnames(pb_counts) == pb_meta$pb_id)
  )

  pb_meta_filt = pb_meta %>%
    dplyr::filter(
      n_cells >= min_cells
    )

  pb_counts_filt = pb_counts[
    ,
    pb_meta_filt$pb_id,
    drop = FALSE
  ]
  dge = edgeR::DGEList(
    counts = pb_counts_filt
  )

  keep_genes = edgeR::filterByExpr(
    dge,
    group = paste(
      pb_meta_filt$cluster,
      pb_meta_filt$group
    )
  )

  dge = dge[
    keep_genes,
    ,
    keep.lib.sizes = FALSE
  ]

  dge = edgeR::calcNormFactors(dge)

  logCPM = edgeR::cpm(
    dge,
    log = TRUE,
    prior.count = 1
  )
  msig_hallmark = msigdbr::msigdbr(
    species = "Homo sapiens",
    category = "H"
  )

  hallmark_sets = split(
    msig_hallmark$gene_symbol,
    msig_hallmark$gs_name
  )

  hallmark_sets_use = hallmark_sets[pathways_use]

  hallmark_sets_use = lapply(
    hallmark_sets_use,
    function(gene_set) {
      intersect(
        gene_set,
        rownames(logCPM)
      )
    }
  )

  hallmark_sets_use = hallmark_sets_use[
    lengths(hallmark_sets_use) >= min_genes_per_set
  ]

  ssgsea_param = GSVA::ssgseaParam(
    exprData = as.matrix(logCPM),
    geneSets = hallmark_sets_use,
    alpha = 0.25,
    normalize = FALSE
  )

  ssgsea_scores = GSVA::gsva(
    ssgsea_param
  )

  meta_hm = pb_meta_filt %>%
    as.data.frame() %>%
    dplyr::filter(
      pb_id %in% colnames(ssgsea_scores)
    ) %>%
    dplyr::mutate(
      cluster = factor(
        as.character(cluster),
        levels = clusters_use
      ),
      group = factor(
        as.character(group),
        levels = c("CD", "Healthy")
      ),
      sample_id = as.character(sample_id),
      cluster_group = paste0(
        "C",
        cluster,
        "_",
        group
      )
    ) %>%
    dplyr::arrange(
      cluster,
      group,
      sample_id
    )

  ssgsea_samples = ssgsea_scores[
    ,
    meta_hm$pb_id,
    drop = FALSE
  ]

  stopifnot(
    all(colnames(ssgsea_samples) == meta_hm$pb_id)
  )

  mat_samples_z = t(
    scale(
      t(ssgsea_samples)
    )
  )

  mat_samples_z[is.na(mat_samples_z)] = 0

  mat_samples_z[
    mat_samples_z > zscore_limit
  ] = zscore_limit

  mat_samples_z[
    mat_samples_z < -zscore_limit
  ] = -zscore_limit

  rownames(mat_samples_z) = rownames(mat_samples_z) %>%
    gsub(
      "^HALLMARK_",
      "",
      .
    ) %>%
    gsub(
      "_",
      " ",
      .
    )

  col_fun = circlize::colorRamp2(
    c(
      -zscore_limit,
      0,
      zscore_limit
    ),
    c(
      "#2166AC",
      "white",
      "#B2182B"
    )
  )

  ha_samples = ComplexHeatmap::HeatmapAnnotation(
    Cluster = meta_hm$cluster,
    Group = meta_hm$group,
    `n cells` = ComplexHeatmap::anno_barplot(
      meta_hm$n_cells,
      gp = grid::gpar(
        fill = "grey70",
        col = NA
      ),
      height = grid::unit(
        0.8,
        "cm"
      )
    ),
    col = list(
      Cluster = cluster_colors,
      Group = group_colors
    ),
    annotation_name_side = "left",
    annotation_name_gp = grid::gpar(
      fontsize = 9
    )
  )

  column_split_df = data.frame(
    Cluster = meta_hm$cluster,
    Group = meta_hm$group
  )

  heatmap = ComplexHeatmap::Heatmap(
    mat_samples_z,
    name = "ssGSEA\nz-score",
    col = col_fun,
    top_annotation = ha_samples,
    column_split = column_split_df,
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    cluster_column_slices = FALSE,
    show_column_names = FALSE,
    show_row_names = TRUE,
    row_names_gp = grid::gpar(
      fontsize = 12
    ),
    border = TRUE,
    row_title = NULL,
    column_title = "Hallmark pathway activity across samples",
    column_title_gp = grid::gpar(
      fontsize = 12,
      fontface = "bold"
    ),
    heatmap_legend_param = list(
      title = "ssGSEA\nz-score",
      title_gp = grid::gpar(
        fontsize = 9,
        fontface = "bold"
      ),
      labels_gp = grid::gpar(
        fontsize = 8
      )
    )
  )

  return(
    list(
      pseudobulk_counts = pb_counts,
      pseudobulk_metadata = pb_meta,
      pseudobulk_counts_filtered = pb_counts_filt,
      pseudobulk_metadata_filtered = pb_meta_filt,
      dge = dge,
      logCPM = logCPM,
      hallmark_gene_sets = hallmark_sets_use,
      ssgsea_scores = ssgsea_scores,
      heatmap_metadata = meta_hm,
      heatmap_zscore_matrix = mat_samples_z,
      heatmap = heatmap
    )
  )
}
