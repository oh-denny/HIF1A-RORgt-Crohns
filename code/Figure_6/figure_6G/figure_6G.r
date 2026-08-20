library(Seurat)
library(dplyr)
library(tibble)
library(ComplexHeatmap)
library(circlize)
library(grid)

DefaultAssay(SC_obj_pb) = "SCT"

genes_treg = c(
  "FOXP3", "RORC", "IL2RA", "CTLA4", "ICOS", "PDCD1",
  "LAG3", "HAVCR2"
)

genes_hif_axis = c(
  "HIF1A", "RORA", "AHR"
)

genes_tissue_inflam = c(
  "CXCR4", "ITGAE",
  "DDIT4", "SOCS1", "AREG", "LGALS3", "GZMB"
)

genes_mito = c(
  "PPARGC1A", 
  "TFAM",
  "CPT1A",
  "BNIP3",
  "OPA1",
  "MFN1",
  "MFN2",
  "DNM1L",
  "FIS1",
  "MFF",
  "MIEF1",
  "MIEF2"
)

genes_B = c(
  genes_treg,
  genes_hif_axis,
  genes_tissue_inflam,
  genes_mito
)

genes_B_present <- genes_B[genes_B %in% rownames(SC_obj_pb)]
setdiff(genes_B, genes_B_present)

genes_B = genes_B_present
genes_B

df_B = FetchData(
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

avg_B = avg_B %>%
  dplyr::mutate(cluster_num = as.numeric(cluster)) %>%
  dplyr::arrange(cluster_num) %>%
  dplyr::select(-cluster_num)

mat_B = avg_B %>%
  tibble::column_to_rownames("cluster") %>%
  as.matrix() %>%
  t()

mat_B_z = t(scale(t(mat_B)))
mat_B_z[is.na(mat_B_z)] = 0
mat_B_z[mat_B_z > 2] = 2
mat_B_z[mat_B_z < -2] = -2

colnames(mat_B_z) = paste0("Cluster ", colnames(mat_B_z))
cluster_ids = gsub("Cluster ", "", colnames(mat_B_z))

ha_B = HeatmapAnnotation(
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

gene_category = dplyr::case_when(
  rownames(mat_B_z) %in% genes_treg ~ "Treg / suppressive",
  rownames(mat_B_z) %in% genes_hif_axis ~ "HIF / ROR / AHR axis",
  rownames(mat_B_z) %in% genes_tissue_inflam ~ "Tissue adaptation / inflammation",
  rownames(mat_B_z) %in% genes_mito ~ "Mitochondrial metabolism",
  TRUE ~ "Other"
)

gene_category = factor(
  gene_category,
  levels = c(
    "Treg / suppressive",
    "HIF / ROR / AHR axis",
    "Tissue adaptation / inflammation",
    "Mitochondrial metabolism"
  )
)

col_fun_B = circlize::colorRamp2(
  c(-2, 0, 2),
  c("#2166AC", "white", "#B2182B")
)

ht_B = Heatmap(
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

tiff("figures/figure6B.tiff", width = 5, height = 6, units = "in", res = 300)
draw(ht_B)
dev.off()
