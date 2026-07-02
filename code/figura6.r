
library(Seurat)
library(dplyr)
library(tibble)
library(ComplexHeatmap)
library(circlize)
library(grid)

DefaultAssay(SC_obj_pb) <- "SCT"

# genes por categoria
genes_treg <- c(
  "FOXP3", "RORC", "IL2RA", "CTLA4", "ICOS", "PDCD1",
  "LAG3", "HAVCR2"
)

genes_hif_axis <- c(
  "HIF1A", "RORA", "AHR"
)

genes_tissue_inflam <- c(
  "CXCR4", "ITGAE",
  "DDIT4", "SOCS1", "AREG", "LGALS3", "GZMB"
)

genes_mito <- c(
  "PPARGC1A",  # PGC-1alpha
  "TFAM",
  "CPT1A",
  "BNIP3",
  "OPA1",
  "MFN1",
  "MFN2",
  "DNM1L",     # DRP1
  "FIS1",
  "MFF",
  "MIEF1",
  "MIEF2"
)

genes_B <- c(
  genes_treg,
  genes_hif_axis,
  genes_tissue_inflam,
  genes_mito
)

# manter apenas genes presentes no objeto
genes_B_present <- genes_B[genes_B %in% rownames(SC_obj_pb)]

# ver quais genes não foram encontrados
setdiff(genes_B, genes_B_present)

genes_B <- genes_B_present
genes_B


df_B <- FetchData(
  SC_obj_pb,
  vars = c(genes_B, "seurat_clusters")
) %>%
  tibble::rownames_to_column("cell") %>%
  dplyr::mutate(
    cluster = as.character(seurat_clusters)
  )

avg_B <- df_B %>%
  dplyr::select(cluster, all_of(genes_B)) %>%
  dplyr::group_by(cluster) %>%
  dplyr::summarise(
    dplyr::across(all_of(genes_B), ~ mean(.x, na.rm = TRUE)),
    .groups = "drop"
  )

avg_B <- avg_B %>%
  dplyr::mutate(cluster_num = as.numeric(cluster)) %>%
  dplyr::arrange(cluster_num) %>%
  dplyr::select(-cluster_num)

mat_B <- avg_B %>%
  tibble::column_to_rownames("cluster") %>%
  as.matrix() %>%
  t()

# z-score por gene
mat_B_z <- t(scale(t(mat_B)))
mat_B_z[is.na(mat_B_z)] <- 0

# limitar escala
mat_B_z[mat_B_z > 2] <- 2
mat_B_z[mat_B_z < -2] <- -2

colnames(mat_B_z) <- paste0("Cluster ", colnames(mat_B_z))


cluster_ids <- gsub("Cluster ", "", colnames(mat_B_z))


ha_B <- HeatmapAnnotation(
  `Cluster identity` = cluster_ids,
  col = list(
    `Cluster identity` = c(
      "4" = "#6C7A3A",
      "5" = "#D8C7AF",
      "6" = "#162E93",
      "1" = "#0F6C84",
      "2" = "#2C6B6B",
      "3" = "#AEB8B7",
      "0" = "#A61C4B"
    )
  ),
  annotation_name_side = "left",
  annotation_name_gp = gpar(fontsize = 9)
)


gene_category <- dplyr::case_when(
  rownames(mat_B_z) %in% genes_treg ~ "Treg / suppressive",
  rownames(mat_B_z) %in% genes_hif_axis ~ "HIF / ROR / AHR axis",
  rownames(mat_B_z) %in% genes_tissue_inflam ~ "Tissue adaptation / inflammation",
  rownames(mat_B_z) %in% genes_mito ~ "Mitochondrial metabolism",
  TRUE ~ "Other"
)

gene_category <- factor(
  gene_category,
  levels = c(
    "Treg / suppressive",
    "HIF / ROR / AHR axis",
    "Tissue adaptation / inflammation",
    "Mitochondrial metabolism"
  )
)

col_fun_B <- circlize::colorRamp2(
  c(-2, 0, 2),
  c("#2166AC", "white", "#B2182B")
)

ht_B <- Heatmap(
  mat_B_z,
  name = "Average\nexpression\nz-score",
  col = col_fun_B,
  top_annotation = ha_B,
  row_split = gene_category,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_column_names = TRUE,
  show_row_names = TRUE,
  row_names_gp = gpar(fontsize = 9, fontface = "italic"),
  column_names_gp = gpar(fontsize = 9),
  column_names_rot = 45,
  border = TRUE,
  row_title_gp = gpar(fontsize = 9, fontface = "bold"),
  row_title_rot = 0,
  row_title = NULL,
  column_title = "Human CD4+ T cell clusters",
  column_title_gp = gpar(fontsize = 11, fontface = "bold"),
  heatmap_legend_param = list(
    title = "Average\nexpression\nz-score",
    title_gp = gpar(fontsize = 9, fontface = "bold"),
    labels_gp = gpar(fontsize = 8)
  )
)

draw(ht_B)










library(Seurat)
library(msigdbr)
library(dplyr)
library(ggplot2)
library(tidyr)

DefaultAssay(SC_obj_pb) <- "SCT"

# garantir que as células batem entre Seurat e Monocle
cells_use <- intersect(colnames(SC_obj_sub), colnames(cds))

SC_tmp <- subset(SC_obj_sub, cells = cells_use)
cds_tmp <- cds[, cells_use]

# extrair pseudotime do Monocle2
df_traj <- pData(cds_tmp) %>%
  as.data.frame() %>%
  tibble::rownames_to_column("cell") %>%
  dplyr::select(
    cell,
    Pseudotime,
    State,
    seurat_clusters,
    group
  ) %>%
  dplyr::mutate(
    seurat_clusters = as.character(seurat_clusters),
    group = factor(group, levels = c("Healthy", "CD"))
  )

head(df_traj)



coords <- reducedDimS(cds_tmp) %>%
  t() %>%
  as.data.frame()

colnames(coords)[1:2] <- c("Component1", "Component2")

coords <- coords %>%
  tibble::rownames_to_column("cell")

df_traj <- df_traj %>%
  dplyr::left_join(coords, by = "cell")

head(df_traj)


hallmarks_interest <- c(
  "HALLMARK_HYPOXIA",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "HALLMARK_IL2_STAT5_SIGNALING"
)

msig_h <- msigdbr(
  species = "Homo sapiens",
  category = "H"
)

hallmark_list <- msig_h %>%
  dplyr::filter(gs_name %in% hallmarks_interest) %>%
  split(.$gs_name) %>%
  lapply(function(x) unique(x$gene_symbol))

# manter só genes presentes no objeto
hallmark_list <- lapply(hallmark_list, function(x) {
  intersect(x, rownames(SC_tmp))
})

# remover vias com poucos genes presentes
hallmark_list <- hallmark_list[lengths(hallmark_list) >= 10]

names(hallmark_list)
lengths(hallmark_list)


library(UCell)

SC_tmp <- UCell::AddModuleScore_UCell(
  SC_tmp,
  features = hallmark_list
)

ucell_cols <- paste0(names(hallmark_list), "_UCell")

ucell_cols
ucell_cols %in% colnames(SC_tmp@meta.data)

# colunas criadas pelo AddModuleScore
score_cols_raw <- paste0("HM", seq_along(hallmark_list))

score_cols_raw
score_cols_raw %in% colnames(SC_tmp@meta.data)


df_scores <- SC_tmp@meta.data %>%
  as.data.frame() %>%
  tibble::rownames_to_column("cell") %>%
  dplyr::select(cell, dplyr::all_of(score_cols_raw))

# renomear sem usar match()
colnames(df_scores) <- c("cell", names(hallmark_list))

head(df_scores)


df_traj_scores <- df_traj %>%
  dplyr::left_join(df_scores, by = "cell")

head(df_traj_scores)



ggplot(df_traj_scores, aes(x = Component1, y = Component2, color = HALLMARK_HYPOXIA)) +
  geom_point(size = 0.5, alpha = 0.8) +
  scale_color_viridis_c(option = "magma", name = "Hypoxia") +
  theme_classic(base_size = 14) +
  labs(
    x = "Component 1",
    y = "Component 2"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16)
  )




df_long_scores <- df_traj_scores %>%
  dplyr::select(
    cell,
    group,
    seurat_clusters,
    State,
    Pseudotime,
    dplyr::all_of(names(hallmark_list))
  ) %>%
  tidyr::pivot_longer(
    cols = dplyr::all_of(names(hallmark_list)),
    names_to = "pathway",
    values_to = "score"
  ) %>%
  dplyr::mutate(
    pathway = gsub("^HALLMARK_", "", pathway),
    pathway = gsub("_", " ", pathway),
    pathway = stringr::str_to_title(tolower(pathway))
  )

head(df_long_scores)


p_E <- ggplot(
  df_long_scores,
  aes(x = Pseudotime, y = score, color = pathway)
) +
  geom_smooth(
    se = FALSE,
    method = "loess",
    linewidth = 1.1,
    span = 0.7
  ) +
  theme_classic(base_size = 14) +
  labs(
    x = "Pseudotime",
    y = "Module score",
    title = "Inflammatory and hypoxia programs along pseudotime"
  ) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.title = element_blank()
  )

p_E

ggsave(
  "Figure6E_hallmark_scores_along_pseudotime.pdf",
  p_E,
  width = 6,
  height = 4
)


p_E_group <- ggplot(
  df_long_scores,
  aes(x = Pseudotime, y = score, color = group)
) +
  geom_smooth(
    se = FALSE,
    method = "loess",
    linewidth = 1,
    span = 0.7
  ) +
  facet_wrap(~ pathway, scales = "free_y", ncol = 2) +
  scale_color_manual(
    values = c(
      "Healthy" = "#2D6F6B",
      "CD" = "#E9969B"
    )
  ) +
  theme_classic(base_size = 13) +
  labs(
    x = "Pseudotime",
    y = "Module score",
    title = "Pathway activity along pseudotime"
  ) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.title = element_blank(),
    strip.text = element_text(face = "bold")
  )

p_E_group

ggsave(
  "Figure6E_hallmark_scores_along_pseudotime_by_group.pdf",
  p_E_group,
  width = 7,
  height = 5
)


ucell_cols <- paste0(names(hallmark_list), "_UCell")

ucell_cols
ucell_cols %in% colnames(SC_tmp@meta.data)


df_ucell <- SC_tmp@meta.data %>%
  as.data.frame() %>%
  tibble::rownames_to_column("cell") %>%
  dplyr::select(cell, sample, group, seurat_clusters, dplyr::all_of(ucell_cols))

# renomear
colnames(df_ucell) <- c(
  "cell", "sample", "group", "seurat_clusters",
  names(hallmark_list)
)

df_traj_ucell <- df_traj %>%
  dplyr::left_join(
    df_ucell %>% dplyr::select(cell, sample, dplyr::all_of(names(hallmark_list))),
    by = "cell"
  )


ggplot(df_traj_ucell, aes(x = Component1, y = Component2, color = HALLMARK_TNFA_SIGNALING_VIA_NFKB)) +
  geom_point(size = 0.5, alpha = 0.8) +
  scale_color_viridis_c(option = "magma", name = "TNFa") +
  theme_classic(base_size = 14) +
  labs(
    x = "Component 1",
    y = "Component 2"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16)
  )


df_long_ucell <- df_traj_ucell %>%
  tidyr::pivot_longer(
    cols = dplyr::all_of(names(hallmark_list)),
    names_to = "pathway",
    values_to = "score"
  ) %>%
  dplyr::mutate(
    pathway = gsub("^HALLMARK_", "", pathway),
    pathway = gsub("_", " ", pathway),
    pathway = stringr::str_to_title(tolower(pathway))
  )








