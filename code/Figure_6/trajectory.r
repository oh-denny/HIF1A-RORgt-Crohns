library(Seurat)
library(monocle)
library(dplyr)
library(ggplot2)

source("functions/patch/patch.r")
patch_monocle2()

seurat_rds = "data/seu_obj/SC_obj_sub_filtered.rds"

if (!exists("SC_obj_sub")) {
  SC_obj_sub = readRDS(seurat_rds)
}

DefaultAssay(SC_obj_sub) = "SCT"
Idents(SC_obj_sub) = "seurat_clusters"

cds = as.CellDataSet(SC_obj_sub)
patch_monocle2()

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
  filter(p_val_adj < 0.05) |>
  group_by(cluster) |>
  slice_max(order_by = avg_log2FC, n = 50, with_ties = FALSE) |>
  pull(gene) |>
  unique()

ordering_genes_markers = intersect(ordering_genes_markers, rownames(cds))


cds = estimateSizeFactors(cds)
cds = estimateDispersions(cds)
cds = setOrderingFilter(cds, ordering_genes_markers)

cds = reduceDimension(
  cds,
  max_components = 2,
  method = "DDRTree"
)

patch_monocle2()
cds = monocle::orderCells(cds)

saveRDS(cds, "data/seu_obj/monocle2_cds_ordered.rds")
