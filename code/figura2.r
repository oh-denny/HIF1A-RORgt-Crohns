


library(Seurat)
library(Matrix)
library(dplyr)
library(tibble)
library(DESeq2)
library(ggplot2)
library(purrr)

# cluster RORγt+ Treg
cluster_ror <- "4"

meta <- SC_obj_sub@meta.data %>%
  as.data.frame() %>%
  rownames_to_column("cell") %>%
  mutate(
    cluster = as.character(seurat_clusters),
    sample_id = sample,
    group = factor(group, levels = c("Healthy", "CD"))
  )

cells_cluster4 <- meta %>%
  filter(cluster == cluster_ror) %>%
  pull(cell)

length(cells_cluster4)
table(meta$group[meta$cell %in% cells_cluster4])
table(meta$sample_id[meta$cell %in% cells_cluster4], meta$group[meta$cell %in% cells_cluster4])


DefaultAssay(SC_obj_sub) <- "RNA"

SC_obj_pb <- JoinLayers(
  object = SC_obj_sub,
  assay = "RNA"
)

Layers(SC_obj_pb[["RNA"]])

counts_mat <- GetAssayData(
  SC_obj_pb,
  assay = "RNA",
  layer = "counts"
)


meta <- SC_obj_pb@meta.data %>%
  as.data.frame() %>%
  tibble::rownames_to_column("cell") %>%
  mutate(
    cluster = as.character(seurat_clusters),
    sample_id = sample,
    group = factor(group, levels = c("Healthy", "CD"))
  )


cluster_ror <- "4"

cells_cluster4 <- meta %>%
  filter(cluster == cluster_ror) %>%
  pull(cell)

counts_c4 <- counts_mat[, cells_cluster4]

meta_c4 <- meta %>%
  filter(cell %in% cells_cluster4) %>%
  arrange(match(cell, colnames(counts_c4)))





meta_c4 <- meta_c4 %>%
  dplyr::mutate(
    pb_id = paste(sample_id, group, sep = "__")
  )
# checar número de células por amostra
n_cells_c4 <- meta_c4 %>%
  as.data.frame() %>%
  dplyr::group_by(pb_id, sample_id, group) %>%
  dplyr::summarise(
    n_cells = dplyr::n(),
    .groups = "drop"
  ) %>%
  dplyr::arrange(group, sample_id)

n_cells_c4
n_cells_c4
table(n_cells_c4$group)



mm <- Matrix::sparse.model.matrix(~ 0 + pb_id, data = meta_c4)

colnames(mm) <- gsub("^pb_id", "", colnames(mm))

pb_counts_c4 <- counts_c4 %*% mm

pb_meta_c4 <- n_cells_c4 %>%
  dplyr::filter(pb_id %in% colnames(pb_counts_c4)) %>%
  dplyr::arrange(match(pb_id, colnames(pb_counts_c4)))

pb_counts_c4 <- pb_counts_c4[, pb_meta_c4$pb_id]

stopifnot(all(colnames(pb_counts_c4) == pb_meta_c4$pb_id))

dim(pb_counts_c4)
pb_meta_c4

min_cells <- 9

pb_meta_c4_filt <- pb_meta_c4 %>%
  dplyr::filter(n_cells >= min_cells)

pb_counts_c4_filt <- pb_counts_c4[, pb_meta_c4_filt$pb_id]

table(pb_meta_c4_filt$group)
pb_meta_c4_filt



library(DESeq2)

counts_use <- pb_counts_c4_filt

meta_use <- pb_meta_c4_filt %>%
  as.data.frame()

rownames(meta_use) <- meta_use$pb_id

meta_use$group <- factor(meta_use$group, levels = c("Healthy", "CD"))

# garantir mesma ordem
counts_use <- counts_use[, rownames(meta_use)]

stopifnot(all(colnames(counts_use) == rownames(meta_use)))

# filtrar genes pouco expressos
keep_genes <- rowSums(counts_use >= 10) >= 2
counts_use <- counts_use[keep_genes, ]

dds_c4 <- DESeqDataSetFromMatrix(
  countData = round(as.matrix(counts_use)),
  colData = meta_use,
  design = ~ group
)

dds_c4 <- DESeq(dds_c4)

res_c4 <- results(
  dds_c4,
  contrast = c("group", "CD", "Healthy")
)

de_c4 <- as.data.frame(res_c4) %>%
  tibble::rownames_to_column("gene") %>%
  dplyr::arrange(padj)

head(de_c4, 30)



genes_check <- c(
  "FOXP3", "RORC", "IL2RA", "CTLA4", "ICOS",
  "PDCD1", "LAG3", "HAVCR2",
  "AHR", "HIF1A", "RORA",
  "CCR6", "CCR9", "ITGAE", "ITGB7",
  "CCR7", "CXCR3", "CCR5", "CXCR4",
  "S1PR1", "SELL", "CD44"
)

de_c4 %>%
  dplyr::filter(gene %in% genes_check) %>%
  dplyr::arrange(desc(log2FoldChange))



homing_genes <- c(
  "CCR6", "CCR9", "ITGAE", "ITGB7",
  "CCR7", "CXCR3", "CCR5", "CXCR4",
  "S1PR1", "SELL", "CD44", "SELPLG",
  "ITGAL", "ITGB2"
)
de_c4_homing <- de_c4 %>%
  dplyr::filter(gene %in% homing_genes) %>%
  dplyr::mutate(
    sig = dplyr::case_when(
      padj < 0.05 & log2FoldChange > 0 ~ "Up in CD",
      padj < 0.05 & log2FoldChange < 0 ~ "Down in CD",
      TRUE ~ "NS"
    )
  )

de_c4_homing


p_c4_homing <- de_c4_homing %>%
  dplyr::mutate(
    gene = factor(gene, levels = gene[order(log2FoldChange)])
  ) %>%
  ggplot(aes(x = gene, y = log2FoldChange, color = sig)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4) +
  geom_point(size = 3) +
  coord_flip() +
  theme_classic(base_size = 13) +
  labs(
    x = NULL,
    y = "log2FC CD vs Healthy",
    title = "Cluster 4 RORγt⁺ Treg-like homing genes"
  ) +
  theme(
    axis.text.y = element_text(face = "italic"),
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.title = element_blank()
  )

p_c4_homing


de_c4_plot <- de_c4 %>%
  dplyr::mutate(
    neglog10padj = -log10(padj),
    sig = dplyr::case_when(
      padj < 0.05 & log2FoldChange > 0.25 ~ "Up in CD",
      padj < 0.05 & log2FoldChange < -0.25 ~ "Down in CD",
      TRUE ~ "NS"
    )
  )

p_volcano_c4 <- ggplot(
  de_c4_plot,
  aes(x = log2FoldChange, y = -log10(pvalue), color = sig)
) +
  geom_point(size = 1.3, alpha = 0.8) +
  geom_vline(xintercept = c(-0.25, 0.25), linetype = "dashed", linewidth = 0.3) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", linewidth = 0.3) +
  theme_classic(base_size = 13) +
  labs(
    x = "log2FC CD vs Healthy",
    y = "-log10 adjusted p-value",
    title = "Pseudobulk DE in cluster 4"
  ) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.title = element_blank()
  )

p_volcano_c4





library(DESeq2)
library(ggplot2)

vsd_c4 <- vst(dds_c4, blind = FALSE)

pcaData <- plotPCA(vsd_c4, intgroup = c("group"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(PC1, PC2, color = group, label = name)) +
  geom_point(size = 4) +
  ggrepel::geom_text_repel(size = 3) +
  theme_classic(base_size = 13) +
  labs(
    x = paste0("PC1: ", percentVar[1], "% variance"),
    y = paste0("PC2: ", percentVar[2], "% variance"),
    title = "Pseudobulk PCA - cluster 4"
  )





genes_fig2 <- c(
  "CXCR4", "ITGAE", "CCR6", "CCR7",
  "IL2RA", "FOXP3", "RORA", "HIF1A",
  "SOCS1", "DDIT4", "AREG", "LGALS3", "GZMB"
)

de_fig2 <- de_c4 %>%
  dplyr::filter(gene %in% genes_fig2) %>%
  dplyr::mutate(
    sig = dplyr::case_when(
      padj < 0.05 & log2FoldChange > 0 ~ "FDR < 0.05, up in CD",
      padj < 0.05 & log2FoldChange < 0 ~ "FDR < 0.05, down in CD",
      pvalue < 0.05 ~ "nominal p < 0.05",
      TRUE ~ "NS"
    ),
    gene = factor(gene, levels = gene[order(log2FoldChange)])
  )

p_fig2_lollipop <- ggplot(
  de_fig2,
  aes(x = gene, y = log2FoldChange)
) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4) +
  geom_segment(aes(x = gene, xend = gene, y = 0, yend = log2FoldChange), linewidth = 0.5) +
  geom_point(size = 3) +
  coord_flip() +
  theme_classic(base_size = 20) +
  labs(
    x = NULL,
    y = "log2FC CD vs Healthy",
    title = "Human RORC⁺ Treg-like cluster 4"
  ) +
  theme(
    axis.text.y = element_text(),
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.title = element_blank()
  )

p_fig2_lollipop












#### enriquecimento
library(dplyr)
library(msigdbr)
library(fgsea)
library(ggplot2)
library(tibble)

# ranking usando estatística Wald do DESeq2
ranks_c4 <- de_c4 %>%
  dplyr::filter(!is.na(stat)) %>%
  dplyr::arrange(desc(stat)) %>%
  dplyr::select(gene, stat) %>%
  tibble::deframe()

# Hallmark humano
msig_hallmark <- msigdbr(
  species = "Homo sapiens",
  category = "H"
)

hallmark_sets <- split(
  msig_hallmark$gene_symbol,
  msig_hallmark$gs_name
)

fgsea_hallmark_c4 <- fgsea(
  pathways = hallmark_sets,
  stats = ranks_c4,
  minSize = 10,
  maxSize = 500,
  eps = 0
) %>%
  as.data.frame() %>%
  arrange(padj)

fgsea_hallmark_c4 %>%
  select(pathway, NES, pval, padj, size) %>%
  head(20)





library(fgsea)
library(ggplot2)

# enrichment plot da via de hipóxia
p_hypoxia <- plotEnrichment(
  pathway = hallmark_sets[["HALLMARK_HYPOXIA"]],
  stats = ranks_c4
) +
  labs(
    title = "HALLMARK_HYPOXIA",
    subtitle = "Human RORC+ Treg-like cluster 4 (CD vs Healthy)",
    x = "Genes ranked by DESeq2 statistic",
    y = "Running enrichment score"
  ) +
  theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5)
  )

p_hypoxia



p_tnfa <- plotEnrichment(
  pathway = hallmark_sets[["HALLMARK_TNFA_SIGNALING_VIA_NFKB"]],
  stats = ranks_c4
) +
  labs(
    title = "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
    subtitle = "Human RORC+ Treg-like cluster 4 (CD vs Healthy)",
    x = "Genes ranked by DESeq2 statistic",
    y = "Running enrichment score"
  ) +
  theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5)
  )

p_tnfa




sig_hallmarks <- fgsea_hallmark_c4 %>%
  dplyr::filter(padj < 0.05) %>%
  dplyr::arrange(desc(NES)) %>%
  dplyr::pull(pathway)

sig_hallmarks



library(DESeq2)

vsd_c4 <- vst(dds_c4, blind = FALSE)
mat_vst <- assay(vsd_c4)

dim(mat_vst)
head(colnames(mat_vst))


library(GSVA)

# garantir que é matriz numérica genes x samples
mat_vst <- as.matrix(mat_vst)

# manter só genes presentes nos sets
sig_hallmark_sets <- hallmark_sets[sig_hallmarks]

sig_hallmark_sets <- lapply(sig_hallmark_sets, function(x) {
  intersect(x, rownames(mat_vst))
})

# remover vias que ficaram com poucos genes
sig_hallmark_sets <- sig_hallmark_sets[lengths(sig_hallmark_sets) >= 10]

# novo formato da API GSVA
param_ssgsea <- ssgseaParam(
  exprData = mat_vst,
  geneSets = sig_hallmark_sets,
  alpha = 0.25,
  normalize = TRUE
)

ssgsea_res <- gsva(param_ssgsea)

dim(ssgsea_res)
ssgsea_res[1:5, 1:5]





library(ComplexHeatmap)
library(circlize)
library(grid)
library(dplyr)

# usar apenas vias significativas, ordenadas pelo NES do fgsea
sig_paths <- fgsea_hallmark_c4 %>%
  dplyr::filter(padj < 0.05) %>%
  dplyr::arrange(desc(NES)) %>%
  dplyr::pull(pathway)

sig_paths <- intersect(sig_paths, rownames(ssgsea_res))

mat_hm <- ssgsea_res[sig_paths, , drop = FALSE]


# meta_use precisa ter rownames = pb_id, igual às colunas do ssgsea_res
meta_hm <- meta_use[colnames(mat_hm), , drop = FALSE]

meta_hm$group <- factor(meta_hm$group, levels = c("CD", "Healthy"))

# ordenar CD primeiro, Healthy depois
ord <- order(meta_hm$group)

mat_hm <- mat_hm[, ord, drop = FALSE]
meta_hm <- meta_hm[ord, , drop = FALSE]


mat_hm_z <- t(scale(t(mat_hm)))
mat_hm_z[is.na(mat_hm_z)] <- 0



pathway_labels <- rownames(mat_hm_z) %>%
  gsub("^HALLMARK_", "", .) %>%
  gsub("_", " ", .)


group_col <- c(
  "CD" = "#E9969B",
  "Healthy" = "#2D6F6B"
)

col_fun <- circlize::colorRamp2(
  c(-2, 0, 2),
  c("#2166AC", "white", "#B2182B")
)

ha <- HeatmapAnnotation(
  Group = meta_hm$group,
  `n cells` = anno_barplot(
    meta_hm$n_cells,
    gp = gpar(fill = "grey70", col = NA),
    height = unit(1, "cm")
  ),
  col = list(Group = group_col),
  annotation_name_side = "left",
  annotation_name_gp = gpar(fontsize = 9),
  show_annotation_name = TRUE
)

ht <- Heatmap(
  mat_hm_z,
  name = "ssGSEA\nz-score",
  col = col_fun,
  top_annotation = ha,
  column_split = meta_hm$group,
  cluster_columns = TRUE,
  cluster_column_slices = FALSE,
  cluster_rows = TRUE,
  row_labels = pathway_labels,
  row_names_gp = gpar(fontsize = 9),
  show_column_names = FALSE,
  border = TRUE,
  row_title = NULL,
  column_title = NULL,
  heatmap_legend_param = list(
    title = "ssGSEA\nz-score",
    title_gp = gpar(fontsize = 9, fontface = "bold"),
    labels_gp = gpar(fontsize = 8)
  )
)

draw(ht)





















###### ultimo



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

DefaultAssay(SC_obj_sub) <- "RNA"

# Junta layers do Seurat v5
SC_obj_pb <- JoinLayers(
  object = SC_obj_sub,
  assay = "RNA"
)

counts_mat <- GetAssayData(
  SC_obj_pb,
  assay = "RNA",
  layer = "counts"
)

meta_all <- SC_obj_pb@meta.data %>%
  as.data.frame() %>%
  rownames_to_column("cell") %>%
  mutate(
    cluster = as.character(seurat_clusters),
    sample_id = sample,
    group = factor(group, levels = c("Healthy", "CD"))
  )

stopifnot(all(colnames(counts_mat) %in% meta_all$cell))



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


res_c5 <- run_pb_deseq_cluster(
  cluster_use = "5",
  counts_mat = counts_mat,
  meta_all = meta_all,
  min_cells = 9,
  outdir = "pseudobulk_clusters"
)

res_c6 <- run_pb_deseq_cluster(
  cluster_use = "6",
  counts_mat = counts_mat,
  meta_all = meta_all,
  min_cells = 9,
  outdir = "pseudobulk_clusters"
)



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

fgsea_c6 <- run_hallmark_fgsea(
  de = res_c6$de,
  cluster_use = "6"
)

fgsea_c5$fgsea %>%
  dplyr::select(pathway, NES, pval, padj, size) %>%
  head(20)

fgsea_c6$fgsea %>%
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



p_hypoxia_c5 <- plot_enrichment_cluster(fgsea_c5, "HALLMARK_HYPOXIA", "5")
p_tnfa_c5    <- plot_enrichment_cluster(fgsea_c5, "HALLMARK_TNFA_SIGNALING_VIA_NFKB", "5")

p_hypoxia_c6 <- plot_enrichment_cluster(fgsea_c6, "HALLMARK_HYPOXIA", "6")
p_tnfa_c6    <- plot_enrichment_cluster(fgsea_c6, "HALLMARK_TNFA_SIGNALING_VIA_NFKB", "6")

p_hypoxia_c5
p_tnfa_c5
p_hypoxia_c6
p_tnfa_c6




run_ssgsea_heatmap_cluster <- function(
    dds,
    meta_use,
    fgsea_obj,
    cluster_use,
    outdir = "pseudobulk_clusters"
) {
  
  sig_paths <- fgsea_obj$fgsea %>%
    filter(padj < 0.05) %>%
    arrange(desc(NES)) %>%
    pull(pathway)
  
  if (length(sig_paths) < 2) {
    warning("Poucas vias significativas para cluster ", cluster_use)
  }
  
  vsd <- vst(dds, blind = FALSE)
  mat_vst <- assay(vsd)
  mat_vst <- as.matrix(mat_vst)
  
  sig_hallmark_sets <- fgsea_obj$hallmark_sets[sig_paths]
  
  sig_hallmark_sets <- lapply(sig_hallmark_sets, function(x) {
    intersect(x, rownames(mat_vst))
  })
  
  sig_hallmark_sets <- sig_hallmark_sets[lengths(sig_hallmark_sets) >= 10]
  
  param_ssgsea <- ssgseaParam(
    exprData = mat_vst,
    geneSets = sig_hallmark_sets,
    alpha = 0.25,
    normalize = TRUE
  )
  
  ssgsea_res <- gsva(param_ssgsea)
  
  meta_hm <- meta_use[colnames(ssgsea_res), , drop = FALSE]
  meta_hm$group <- factor(meta_hm$group, levels = c("CD", "Healthy"))
  
  ord <- order(meta_hm$group)
  
  mat_hm <- ssgsea_res[, ord, drop = FALSE]
  meta_hm <- meta_hm[ord, , drop = FALSE]
  
  mat_hm_z <- t(scale(t(mat_hm)))
  mat_hm_z[is.na(mat_hm_z)] <- 0
  
  pathway_labels <- rownames(mat_hm_z) %>%
    gsub("^HALLMARK_", "", .) %>%
    gsub("_", " ", .)
  
  group_col <- c(
    "CD" = "#E9969B",
    "Healthy" = "#2D6F6B"
  )
  
  col_fun <- circlize::colorRamp2(
    c(-2, 0, 2),
    c("#2166AC", "white", "#B2182B")
  )
  
  ha <- HeatmapAnnotation(
    Group = meta_hm$group,
    `n cells` = anno_barplot(
      meta_hm$n_cells,
      gp = gpar(fill = "grey70", col = NA),
      height = unit(1, "cm")
    ),
    col = list(Group = group_col),
    annotation_name_side = "left",
    annotation_name_gp = gpar(fontsize = 9),
    show_annotation_name = TRUE
  )
  
  ht <- Heatmap(
    mat_hm_z,
    name = "ssGSEA\nz-score",
    col = col_fun,
    top_annotation = ha,
    column_split = meta_hm$group,
    cluster_columns = TRUE,
    cluster_column_slices = FALSE,
    cluster_rows = TRUE,
    row_labels = pathway_labels,
    row_names_gp = gpar(fontsize = 9),
    show_column_names = FALSE,
    border = TRUE,
    row_title = NULL,
    column_title = paste0("Cluster ", cluster_use),
    column_title_gp = gpar(fontsize = 12, fontface = "bold"),
    heatmap_legend_param = list(
      title = "ssGSEA\nz-score",
      title_gp = gpar(fontsize = 9, fontface = "bold"),
      labels_gp = gpar(fontsize = 8)
    )
  )
  
  pdf(
    file.path(outdir, paste0("ssGSEA_Hallmark_cluster", cluster_use, "_ComplexHeatmap.pdf")),
    width = 6,
    height = 5
  )
  draw(ht)
  dev.off()
  
  return(list(
    ssgsea = ssgsea_res,
    heatmap = ht,
    mat_z = mat_hm_z,
    meta_hm = meta_hm
  ))
}

ssgsea_c5 <- run_ssgsea_heatmap_cluster(
  dds = res_c5$dds,
  meta_use = res_c5$meta_use,
  fgsea_obj = fgsea_c5,
  cluster_use = "5"
)

draw(ssgsea_c5$heatmap)

ssgsea_c6 <- run_ssgsea_heatmap_cluster(
  dds = res_c6$dds,
  meta_use = res_c6$meta_use,
  fgsea_obj = fgsea_c6,
  cluster_use = "6"
)

draw(ssgsea_c6$heatmap)








library(Seurat)
library(Matrix)
library(dplyr)
library(tibble)
library(edgeR)
library(msigdbr)
library(GSVA)
library(ComplexHeatmap)
library(circlize)
library(grid)

DefaultAssay(SC_obj_sub) <- "RNA"

SC_obj_pb <- JoinLayers(
  object = SC_obj_sub,
  assay = "RNA"
)

counts_mat <- GetAssayData(
  SC_obj_pb,
  assay = "RNA",
  layer = "counts"
)

meta_all <- SC_obj_pb@meta.data %>%
  as.data.frame() %>%
  rownames_to_column("cell") %>%
  mutate(
    cluster = as.character(seurat_clusters),
    sample_id = sample,
    group = factor(group, levels = c("Healthy", "CD"))
  )

clusters_use <- c("4", "5", "6")

meta_456 <- meta_all %>%
  filter(cluster %in% clusters_use) %>%
  mutate(
    pb_id = paste0("C", cluster, "__", sample_id, "__", group)
  )

counts_456 <- counts_mat[, meta_456$cell]

meta_456 <- meta_456 %>%
  arrange(match(cell, colnames(counts_456)))

stopifnot(all(meta_456$cell == colnames(counts_456)))


mm <- Matrix::sparse.model.matrix(~ 0 + pb_id, data = meta_456)
colnames(mm) <- gsub("^pb_id", "", colnames(mm))

pb_counts_456 <- counts_456 %*% mm

pb_meta_456 <- meta_456 %>%
  group_by(pb_id, cluster, sample_id, group) %>%
  summarise(
    n_cells = n(),
    .groups = "drop"
  ) %>%
  arrange(match(pb_id, colnames(pb_counts_456)))

pb_counts_456 <- pb_counts_456[, pb_meta_456$pb_id]

stopifnot(all(colnames(pb_counts_456) == pb_meta_456$pb_id))

pb_meta_456
table(pb_meta_456$cluster, pb_meta_456$group)

min_cells <- 9

pb_meta_456_filt <- pb_meta_456 %>%
  filter(n_cells >= min_cells)

pb_counts_456_filt <- pb_counts_456[, pb_meta_456_filt$pb_id]

table(pb_meta_456_filt$cluster, pb_meta_456_filt$group)


dge <- edgeR::DGEList(counts = pb_counts_456_filt)

keep_genes <- edgeR::filterByExpr(
  dge,
  group = paste(pb_meta_456_filt$cluster, pb_meta_456_filt$group)
)

dge <- dge[keep_genes, , keep.lib.sizes = FALSE]
dge <- edgeR::calcNormFactors(dge)

logCPM_456 <- edgeR::cpm(
  dge,
  log = TRUE,
  prior.count = 1
)

dim(logCPM_456)



msig_hallmark <- msigdbr(
  species = "Homo sapiens",
  category = "H"
)

hallmark_sets <- split(
  msig_hallmark$gene_symbol,
  msig_hallmark$gs_name
)

# escolha vias de interesse
pathways_use <- c(
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
)

hallmark_sets_use <- hallmark_sets[pathways_use]

hallmark_sets_use <- lapply(hallmark_sets_use, function(x) {
  intersect(x, rownames(logCPM_456))
})

hallmark_sets_use <- hallmark_sets_use[lengths(hallmark_sets_use) >= 10]

param_ssgsea_456 <- ssgseaParam(
  exprData = as.matrix(logCPM_456),
  geneSets = hallmark_sets_use,
  alpha = 0.25,
  normalize = FALSE
)

ssgsea_456 <- gsva(param_ssgsea_456)

dim(ssgsea_456)

meta_hm <- pb_meta_456_filt %>%
  mutate(
    cluster_group = paste0("C", cluster, "_", group)
  )

rownames(meta_hm) <- meta_hm$pb_id

meta_hm <- pb_meta_456_filt %>%
  as.data.frame() %>%
  dplyr::filter(pb_id %in% colnames(ssgsea_456)) %>%
  dplyr::mutate(
    cluster = as.character(cluster),
    group = as.character(group),
    cluster_group = paste0("C", cluster, "_", group)
  )

# ordenar metadata exatamente como as colunas do ssGSEA
meta_hm <- meta_hm %>%
  dplyr::arrange(match(pb_id, colnames(ssgsea_456)))

ssgsea_456 <- ssgsea_456[, meta_hm$pb_id, drop = FALSE]

stopifnot(all(meta_hm$pb_id == colnames(ssgsea_456)))

table(meta_hm$cluster_group)

cluster_group_levels <- c(
  "C4_CD", "C4_Healthy",
  "C5_CD", "C5_Healthy",
  "C6_CD", "C6_Healthy"
)

mat_avg <- sapply(cluster_group_levels, function(cg) {
  
  cols <- meta_hm$pb_id[meta_hm$cluster_group == cg]
  cols <- intersect(cols, colnames(ssgsea_456))
  
  if (length(cols) == 0) {
    warning("Nenhuma amostra encontrada para: ", cg)
    return(rep(NA_real_, nrow(ssgsea_456)))
  }
  
  rowMeans(ssgsea_456[, cols, drop = FALSE], na.rm = TRUE)
})

mat_avg <- as.matrix(mat_avg)
rownames(mat_avg) <- rownames(ssgsea_456)
colnames(mat_avg) <- cluster_group_levels


keep_cols <- colSums(is.na(mat_avg)) < nrow(mat_avg)

mat_avg <- mat_avg[, keep_cols, drop = FALSE]

mat_avg

mat_avg_z <- t(scale(t(mat_avg)))
mat_avg_z[is.na(mat_avg_z)] <- 0

rownames(mat_avg_z) <- rownames(mat_avg_z) %>%
  gsub("^HALLMARK_", "", .) %>%
  gsub("_", " ", .)


library(ComplexHeatmap)
library(circlize)
library(grid)

annotation_df <- data.frame(
  Cluster = gsub("^C([0-9]+)_.*", "\\1", colnames(mat_avg_z)),
  Group = gsub("^C[0-9]+_", "", colnames(mat_avg_z))
)

rownames(annotation_df) <- colnames(mat_avg_z)

ha <- HeatmapAnnotation(
  Cluster = annotation_df$Cluster,
  Group = annotation_df$Group,
  col = list(
    Cluster = c(
      "4" = "#B2182B",
      "5" = "#2166AC",
      "6" = "#4D9221"
    ),
    Group = c(
      "CD" = "#E9969B",
      "Healthy" = "#2D6F6B"
    )
  ),
  annotation_name_side = "left"
)

col_fun <- circlize::colorRamp2(
  c(-2, 0, 2),
  c("#2166AC", "white", "#B2182B")
)

ht_combined <- Heatmap(
  mat_avg_z,
  name = "ssGSEA\nz-score",
  col = col_fun,
  top_annotation = ha,
  column_split = annotation_df$Cluster,
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  cluster_column_slices = FALSE,
  show_column_names = TRUE,
  show_row_names = TRUE,
  row_names_gp = gpar(fontsize = 9),
  column_names_gp = gpar(fontsize = 9),
  column_names_rot = 45,
  border = TRUE,
  row_title = NULL,
  column_title = "Hallmark pathway activity across clusters",
  column_title_gp = gpar(fontsize = 12, fontface = "bold")
)

draw(ht_combined)






library(dplyr)
library(ComplexHeatmap)
library(circlize)
library(grid)

# metadata alinhada com as colunas do ssGSEA
meta_hm <- pb_meta_456_filt %>%
  as.data.frame() %>%
  dplyr::filter(pb_id %in% colnames(ssgsea_456)) %>%
  dplyr::mutate(
    cluster = factor(as.character(cluster), levels = c("4", "5", "6")),
    group = factor(as.character(group), levels = c("CD", "Healthy")),
    sample_id = as.character(sample_id),
    cluster_group = paste0("C", cluster, "_", group)
  ) %>%
  dplyr::arrange(cluster, group, sample_id)

# alinhar matriz com metadata
ssgsea_samples <- ssgsea_456[, meta_hm$pb_id, drop = FALSE]


mat_samples_z <- t(scale(t(ssgsea_samples)))
mat_samples_z[is.na(mat_samples_z)] <- 0

# limitar extremos para melhorar visualização
mat_samples_z[mat_samples_z > 2] <- 2
mat_samples_z[mat_samples_z < -2] <- -2

# limpar nomes das vias
rownames(mat_samples_z) <- rownames(mat_samples_z) %>%
  gsub("^HALLMARK_", "", .) %>%
  gsub("_", " ", .)

group_col <- c(
  "CD" = "#E9969B",
  "Healthy" = "#2D6F6B"
)

cluster_col <- c(
  "4" = "#B2182B",
  "5" = "#2166AC",
  "6" = "#4D9221"
)

col_fun <- circlize::colorRamp2(
  c(-2, 0, 2),
  c("#2166AC", "white", "#B2182B")
)

ha_samples <- HeatmapAnnotation(
  Cluster = meta_hm$cluster,
  Group = meta_hm$group,
  `n cells` = anno_barplot(
    meta_hm$n_cells,
    gp = gpar(fill = "grey70", col = NA),
    height = unit(0.8, "cm")
  ),
  col = list(
    Cluster = cluster_col,
    Group = group_col
  ),
  annotation_name_side = "left",
  annotation_name_gp = gpar(fontsize = 9)
)

column_split_df <- data.frame(
  Cluster = meta_hm$cluster,
  Group = meta_hm$group
)

ht_samples <- Heatmap(
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
  row_names_gp = gpar(fontsize = 12),
  border = TRUE,
  row_title = NULL,
  column_title = "Hallmark pathway activity across samples",
  column_title_gp = gpar(fontsize = 12, fontface = "bold"),
  heatmap_legend_param = list(
    title = "ssGSEA\nz-score",
    title_gp = gpar(fontsize = 9, fontface = "bold"),
    labels_gp = gpar(fontsize = 8)
  )
)

draw(ht_samples)
