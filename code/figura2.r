


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


























