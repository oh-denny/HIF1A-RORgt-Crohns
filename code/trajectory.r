library(Seurat)
library(tidyverse)
library(magrittr)
library(monocle)

saveRDS(SC_obj_sub, "data/seu_obj/SC_obj_sub_filtered.rds")
SC_obj_sub = readRDS("data/seu_obj/SC_obj_sub_filtered.rds")

source("functions/patch/patch.r")


treg_core = c("FOXP3", "IL2RA", "CTLA4", "IKZF2", "ICOS", "PDCD1", "LAG3", "HAVCR2", "TIGIT")
rorgt_program = c("RORC", "IL23R", "CCR6", "MAF", "AHR")
inflam_genes = c("HIF1A", "IL10", "IL17A", "IL17F", "IFNG", "TNF")

DefaultAssay(SC_obj_sub) = "SCT"
Idents(SC_obj_sub) = "seurat_clusters"

markers = FindAllMarkers(
  SC_obj_sub,
  assay = "SCT",
  slot = "data",
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25,
  recorrect_umi = FALSE
)

ordering_genes_markers = markers |>
  dplyr::filter(p_val_adj < 0.05) |>
  dplyr::pull(gene) |>
  unique()

ordering_genes_markers = ordering_genes_markers[
  ordering_genes_markers %in% rownames(cds)
]

head(markers)
table(markers$cluster)

ordering_genes_markers = markers |>
  dplyr::filter(p_val_adj < 0.05) |>
  dplyr::group_by(cluster) |>
  dplyr::slice_max(order_by = avg_log2FC, n = 50) |>
  dplyr::pull(gene) |>
  unique()

ordering_genes_markers = ordering_genes_markers[
  ordering_genes_markers %in% rownames(cds)
]


cds = as.CellDataSet(SC_obj_sub)

patch_dplyr_for_monocle2()

cds = estimateSizeFactors(cds)
cds = estimateDispersions(cds)

disp_table = dispersionTable(cds)

cds = setOrderingFilter(cds, ordering_genes_markers)
cds = reduceDimension(cds, max_components = 2, method = "DDRTree")
cds = orderCells(cds)

patch_dplyr_select_rename_monocle2()
patch_plot_cell_trajectory_monocle2()
patch_plot_cell_trajectory_monocle2_v2()
patch_plot_cell_trajectory_all_select()
patch_plot_cell_trajectory_remaining_select()

monocle::plot_cell_trajectory(cds, color_by = "seurat_clusters")
monocle::plot_cell_trajectory(cds, color_by = "Pseudotime")
monocle::plot_cell_trajectory(cds, color_by = "seurat_clusters")



table(
  pData(cds)$State,
  pData(cds)$seurat_clusters
)


plot_genes_in_pseudotime(
    cds[c("HIF1A"),],
    color_by="State"
)
plot_genes_in_pseudotime(
  cds["FOXP3",],
  color_by = "State"
)









ordering_genes = subset(
  disp_table,
  mean_expression >= 0.1 &
    dispersion_empirical >= 1 * dispersion_fit
)

ordering_genes = ordering_genes$gene_id

cds = setOrderingFilter(cds, ordering_genes)

plot_ordering_genes(cds)

cds = reduceDimension(
  cds,
  max_components = 2,
  method = "DDRTree"
)

patch_igraph_for_monocle2()
patch_project2MST_monocle2_v2()

cds <- orderCells(cds)






grep("select_", deparse(body(monocle:::plot_cell_trajectory)), value = TRUE)
monocle::plot_cell_trajectory(cds, color_by = "Pseudotime")
monocle::plot_cell_trajectory(cds, color_by = "seurat_clusters")


monocle::plot_genes_in_pseudotime(
  cds[traj_genes, ],
  color_by = "seurat_clusters"
)



DefaultAssay(SC_obj_sub) <- "RNA"

DefaultAssay(SC_obj_sub) <- "RNA"

SC_obj_sub[["RNA"]] <- JoinLayers(SC_obj_sub[["RNA"]])

SC_obj_sub <- NormalizeData(SC_obj_sub)

Idents(SC_obj_sub) <- "seurat_clusters"

SC_obj_sub = PrepSCTFindMarkers(SC_obj_sub, assay = "SCT")

markers <- FindAllMarkers(
  SC_obj_sub,
  assay = "SCT",
  slot = "data",
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)


