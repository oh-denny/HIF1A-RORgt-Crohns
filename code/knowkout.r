#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
})

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || all(is.na(x))) y else x
}

script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1] %||% "code/knowkout.r")
project_dir <- normalizePath(file.path(dirname(script_file), ".."), mustWork = FALSE)
if (!dir.exists(file.path(project_dir, "data"))) {
  project_dir <- normalizePath(".", mustWork = TRUE)
}
setwd(project_dir)

target_gene <- Sys.getenv("CELLO_TARGET_GENE", unset = "HIF1A")
target_gene <- toupper(target_gene)
n_features <- as.integer(Sys.getenv("CELLO_N_FEATURES", unset = "3000"))

seurat_rds <- Sys.getenv(
  "CELLO_SEURAT_RDS",
  unset = "data/seu_obj/SC_obj_sub_filtered.rds"
)
monocle_rds <- Sys.getenv(
  "CELLO_MONOCLE_RDS",
  unset = "data/seu_obj/monocle2_cds_ordered.rds"
)
outdir <- Sys.getenv(
  "CELLO_EXPORT_DIR",
  unset = "celloracle_export_hif1a"
)

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

message("Loading Seurat object: ", seurat_rds)
SC_obj_sub <- readRDS(seurat_rds)

DefaultAssay(SC_obj_sub) <- "RNA"
if ("Assay5" %in% class(SC_obj_sub[["RNA"]]) && length(Layers(SC_obj_sub[["RNA"]])) > 1) {
  SC_obj_sub <- JoinLayers(SC_obj_sub, assay = "RNA")
}

counts_mat <- tryCatch(
  GetAssayData(SC_obj_sub, assay = "RNA", layer = "counts"),
  error = function(e) GetAssayData(SC_obj_sub, assay = "RNA", slot = "counts")
)
counts_mat <- as(counts_mat, "dgCMatrix")

if (!target_gene %in% rownames(counts_mat)) {
  stop("Target gene not found in RNA counts: ", target_gene)
}

message("Selecting genes for CellOracle.")
DefaultAssay(SC_obj_sub) <- "SCT"
variable_genes <- VariableFeatures(SC_obj_sub)
variable_genes <- intersect(variable_genes, rownames(counts_mat))

if (length(variable_genes) < n_features) {
  detected_cells <- Matrix::rowSums(counts_mat > 0)
  detected_cells <- sort(detected_cells, decreasing = TRUE)
  variable_genes <- unique(c(variable_genes, names(detected_cells)))
}

genes_keep <- unique(c(
  target_gene,
  "FOXP3", "IL2RA", "CTLA4", "IKZF2", "ICOS", "PDCD1", "LAG3", "HAVCR2",
  "RORC", "IL23R", "CCR6", "MAF", "AHR", "RORA", "IL10", "IL17A",
  "IL17F", "IFNG", "TNF",
  variable_genes
))
genes_keep <- intersect(genes_keep, rownames(counts_mat))
genes_keep <- genes_keep[seq_len(min(length(genes_keep), n_features))]

if (!target_gene %in% genes_keep) {
  genes_keep <- unique(c(target_gene, genes_keep))
}

counts_mat <- counts_mat[genes_keep, , drop = FALSE]

if (file.exists(monocle_rds)) {
  message("Adding Monocle pseudotime metadata: ", monocle_rds)
  cds <- readRDS(monocle_rds)
  pdata <- Biobase::pData(cds)
  common_cells <- intersect(colnames(SC_obj_sub), rownames(pdata))
  SC_obj_sub$Pseudotime <- NA_real_
  SC_obj_sub$State <- NA
  SC_obj_sub@meta.data[common_cells, "Pseudotime"] <- pdata[common_cells, "Pseudotime"]
  SC_obj_sub@meta.data[common_cells, "State"] <- pdata[common_cells, "State"]
}

meta <- SC_obj_sub@meta.data
meta$cell_id <- rownames(meta)
meta <- meta[colnames(counts_mat), , drop = FALSE]
meta$cell_id <- rownames(meta)

if (!"seurat_clusters" %in% colnames(meta)) {
  stop("metadata column 'seurat_clusters' is required for CellOracle GRN units.")
}

umap_mat <- Embeddings(SC_obj_sub, reduction = "umap")
umap_mat <- umap_mat[colnames(counts_mat), , drop = FALSE]
umap_df <- data.frame(
  cell_id = rownames(umap_mat),
  UMAP_1 = umap_mat[, 1],
  UMAP_2 = umap_mat[, 2],
  row.names = NULL
)

write.csv(meta, file = file.path(outdir, "metadata.csv"), row.names = FALSE)
write.csv(umap_df, file = file.path(outdir, "umap.csv"), row.names = FALSE)
writeLines(rownames(counts_mat), con = file.path(outdir, "genes.tsv"))
writeLines(colnames(counts_mat), con = file.path(outdir, "cells.tsv"))
Matrix::writeMM(counts_mat, file = file.path(outdir, "counts_genes_x_cells.mtx"))

summary_df <- data.frame(
  item = c("cells", "genes", "target_gene", "clusters", "groups"),
  value = c(
    ncol(counts_mat),
    nrow(counts_mat),
    target_gene,
    paste(names(table(meta$seurat_clusters)), table(meta$seurat_clusters), sep = ":", collapse = ";"),
    if ("group" %in% colnames(meta)) paste(names(table(meta$group)), table(meta$group), sep = ":", collapse = ";") else NA
  )
)
write.csv(summary_df, file = file.path(outdir, "export_summary.csv"), row.names = FALSE)

message("CellOracle export written to: ", normalizePath(outdir))
message("Cells: ", ncol(counts_mat), " | Genes: ", nrow(counts_mat), " | Target: ", target_gene)






library(tidyverse)

cluster_delta = read_csv(
  "./results/celloracle_hif1a_ko_1500/tables/cluster_delta_summary.csv"
)

gene_delta = read_csv(
  "./results/celloracle_hif1a_ko_1500/tables/gene_delta_summary.csv"
)

gene_delta %>%
  arrange(desc(abs(mean_delta))) %>%
  head(30)

cluster_delta

FeaturePlot(SC_obj_sub, features = "HIF1A")





library(tidyverse)
library(grid)




library(Seurat)
library(tidyverse)




library(tidyverse)
library(ggrepel)
library(viridis)
library(grid)

# caminho da rodada final recomendada
res_dir <- "/Users/denny/Desktop/niels/prof_niels/results/celloracle_hif1a_ko_1500_FA_paperstyle_20260626_154852"

tab_dir <- file.path(res_dir, "tables")
fig_dir <- file.path(res_dir, "figures_R")

dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# carregar tabelas
fa_cells <- read_csv(file.path(tab_dir, "fa_embedding_cells.csv"))
cell_shift <- read_csv(file.path(tab_dir, "cell_shift_vectors_FA.csv"))
grid_shift <- read_csv(file.path(tab_dir, "grid_shift_vectors_FA.csv"))
propagation <- read_csv(file.path(tab_dir, "propagation_magnitude_by_cluster.csv"))
pert_score <- read_csv(file.path(tab_dir, "perturbation_score_FA.csv"))

# conferir colunas
names(fa_cells)
names(cell_shift)
names(grid_shift)
names(propagation)
names(pert_score)




fa_cells <- fa_cells %>%
  mutate(seurat_clusters = as.factor(seurat_clusters))

cluster_centers <- fa_cells %>%
  group_by(seurat_clusters) %>%
  summarise(
    FA1 = median(FA1, na.rm = TRUE),
    FA2 = median(FA2, na.rm = TRUE),
    .groups = "drop"
  )

p_fa_clusters <- ggplot(fa_cells, aes(x = FA1, y = FA2)) +
  geom_point(aes(color = seurat_clusters), size = 0.7, alpha = 0.75) +
  geom_text_repel(
    data = cluster_centers,
    aes(label = seurat_clusters),
    color = "black",
    size = 4,
    fontface = "bold",
    max.overlaps = Inf
  ) +
  theme_classic(base_size = 14) +
  labs(
    title = "HIF1A KO - FA clusters",
    x = "FA1",
    y = "FA2",
    color = "seurat_clusters"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 18),
    legend.position = "right"
  )

p_fa_clusters

ggsave(file.path(fig_dir, "FA_clusters_HIF1A_KO_R.png"), p_fa_clusters, width = 8, height = 6, dpi = 300)
ggsave(file.path(fig_dir, "FA_clusters_HIF1A_KO_R.pdf"), p_fa_clusters, width = 8, height = 6)






cell_shift <- cell_shift %>%
  mutate(seurat_clusters = as.factor(seurat_clusters))

arrow_scale <- 1

p_arrows_all <- ggplot() +
  geom_point(
    data = cell_shift,
    aes(x = FA1, y = FA2, color = seurat_clusters),
    size = 0.45,
    alpha = 0.25
  ) +
  geom_segment(
    data = cell_shift,
    aes(
      x = FA1,
      y = FA2,
      xend = FA1 + shift_FA1 * arrow_scale,
      yend = FA2 + shift_FA2 * arrow_scale
    ),
    arrow = arrow(length = unit(0.04, "inches")),
    linewidth = 0.25,
    alpha = 0.65,
    color = "#111827"
  ) +
  theme_classic(base_size = 14) +
  labs(
    title = "HIF1A KO cell shift vectors",
    x = "FA1",
    y = "FA2"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 18),
    legend.position = "none"
  )

p_arrows_all

ggsave(file.path(fig_dir, "FA_HIF1A_KO_cell_shift_arrows_R.png"), p_arrows_all, width = 8, height = 6, dpi = 300)
ggsave(file.path(fig_dir, "FA_HIF1A_KO_cell_shift_arrows_R.pdf"), p_arrows_all, width = 8, height = 6)





set.seed(123)

cell_shift_sampled <- cell_shift %>%
  group_by(seurat_clusters) %>%
  slice_sample(prop = 0.15) %>%
  ungroup()

arrow_scale <- 0.7

p_arrows_sampled <- ggplot() +
  geom_point(
    data = cell_shift,
    aes(x = FA1, y = FA2, color = seurat_clusters),
    size = 0.45,
    alpha = 0.18
  ) +
  geom_segment(
    data = cell_shift_sampled,
    aes(
      x = FA1,
      y = FA2,
      xend = FA1 + shift_FA1 * arrow_scale,
      yend = FA2 + shift_FA2 * arrow_scale
    ),
    arrow = arrow(length = unit(0.045, "inches")),
    linewidth = 0.35,
    alpha = 0.75,
    color = "#111827"
  ) +
  theme_classic(base_size = 14) +
  labs(
    title = "HIF1A KO cell shift vectors",
    x = "FA1",
    y = "FA2"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 18),
    legend.position = "none"
  )

p_arrows_sampled

ggsave(file.path(fig_dir, "FA_HIF1A_KO_cell_shift_arrows_sampled_R.png"), p_arrows_sampled, width = 8, height = 6, dpi = 300)
ggsave(file.path(fig_dir, "FA_HIF1A_KO_cell_shift_arrows_sampled_R.pdf"), p_arrows_sampled, width = 8, height = 6)






p_delta_map <- ggplot(fa_cells, aes(x = FA1, y = FA2)) +
  geom_point(aes(color = delta_l2), size = 0.75, alpha = 0.9) +
  scale_color_viridis_c(option = "magma", name = "delta L2") +
  theme_classic(base_size = 14) +
  labs(
    title = "HIF1A KO effect magnitude",
    x = "FA1",
    y = "FA2"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 18),
    legend.position = "right"
  )

p_delta_map

ggsave(file.path(fig_dir, "FA_HIF1A_KO_delta_l2_map_R.png"), p_delta_map, width = 8, height = 6, dpi = 300)
ggsave(file.path(fig_dir, "FA_HIF1A_KO_delta_l2_map_R.pdf"), p_delta_map, width = 8, height = 6)






propagation <- propagation %>%
  mutate(seurat_clusters = as.factor(seurat_clusters))

p_propagation <- ggplot(
  propagation,
  aes(
    x = n_propagation,
    y = mean_l2_norm,
    color = seurat_clusters,
    group = seurat_clusters
  )
) +
  geom_line(linewidth = 1) +
  geom_point(size = 2.5) +
  theme_classic(base_size = 14) +
  labs(
    title = "HIF1A KO propagation magnitude by cluster",
    x = "n_propagation",
    y = "mean L2 norm",
    color = "seurat_clusters"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 18),
    legend.position = "right"
  )

p_propagation

ggsave(file.path(fig_dir, "HIF1A_KO_propagation_magnitude_by_cluster_R.png"), p_propagation, width = 8, height = 6, dpi = 300)
ggsave(file.path(fig_dir, "HIF1A_KO_propagation_magnitude_by_cluster_R.pdf"), p_propagation, width = 8, height = 6)