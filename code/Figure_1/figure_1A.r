library(Seurat)
library(ggplot2)
library(grid)

seurat_list = readRDS("data/seu_obj/seurat_list.rds")

seurat_list = lapply(seurat_list, function(obj) {
  obj[["percent.mt"]] = PercentageFeatureSet(obj, pattern = "^MT-")
  obj
})

VlnPlot(
  seurat_list[[1]],
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  ncol = 3
)

seurat_list_filtered = lapply(seurat_list, function(obj) {
  subset(
    obj,
    subset =
      nFeature_RNA > 500 &
      nFeature_RNA < 6000 &
      percent.mt < 15
  )
})

seurat_list_sct = lapply(seurat_list_filtered, function(obj) {
  SCTransform(
    obj,
    assay = "RNA",
    new.assay.name = "SCT",
    vars.to.regress = "percent.mt",
    vst.flavor = "v2",
    verbose = TRUE
  )
})


SC_obj = merge(
  x = seurat_list_sct[[1]],
  y = seurat_list_sct[-1],
  add.cell.ids = names(seurat_list_sct)
)

DefaultAssay(SC_obj) = "SCT"
saveRDS(SC_obj, "data/seu_obj/SC_obj_merged_SCT.rds")

SC_obj = readRDS("data/seu_obj/SC_obj_merged_SCT.rds")


table(SC_obj$group)
table(SC_obj$sample)

features = SelectIntegrationFeatures(
  object.list = seurat_list_sct,
  nfeatures = 3000
)

features = SelectIntegrationFeatures(SC_obj, nfeatures = 3000)
VariableFeatures(SC_obj) = features
VariableFeatures(SC_obj) = rownames(SC_obj[["SCT"]]@scale.data)

SC_obj = RunPCA(SC_obj, assay = "SCT", verbose = FALSE, features = features)
ElbowPlot(SC_obj)

SC_obj = FindNeighbors(SC_obj, dims = 1:15)
SC_obj = FindClusters(SC_obj, resolution = 0.5)
SC_obj = RunUMAP(SC_obj, dims = 1:15)

saveRDS(SC_obj, "data/seu_obj/SC_obj_clustered.rds")
SC_obj = readRDS("data/seu_obj/SC_obj_clustered.rds")



#### FOXP3+ subclusters
SC_obj_sub = subset(
  SC_obj,
  idents = c("1"),
)

SC_obj_sub$seurat_clusters = droplevels(SC_obj_sub$seurat_clusters)

DefaultAssay(SC_obj_sub) = "SCT"

SC_obj_sub = RunPCA(
  SC_obj_sub,
  assay = "SCT",
  features = VariableFeatures(SC_obj_sub),
  verbose = FALSE
)
ElbowPlot(SC_obj_sub)

SC_obj_sub = FindNeighbors(SC_obj_sub, dims = 1:10)
SC_obj_sub = FindClusters(SC_obj_sub,resolution = 0.5)
SC_obj_sub = RunUMAP(SC_obj_sub, dims = 1:10)

cols = c(
  "#A61C4B",
  "#0F6C84",
  "#2C6B6B",
  "#AEB8B7",
  "#6C7A3A",
  "#D8C7AF",
  "#162E93"
)

DimPlot(SC_obj_sub, group.by = "group")

tiff("figures/umap_foxp3_subclusters.tiff", width = 5, height = 4, units = "in", res = 300)
DimPlot(SC_obj_sub, group.by = "SCT_snn_res.0.5",
  label = FALSE,
  cols = cols,
  pt.size = 1.5) +
  NoAxes() +
  NoLegend() +
  theme_void() + 
  annotation_custom(
    grob = linesGrob(
      x = unit(c(0.08, 0.08), "npc"),
      y = unit(c(0.08, 0.28), "npc"),
      gp = gpar(lwd = 3)
    )
  ) +
  annotation_custom(
    grob = linesGrob(
      x = unit(c(0.08, 0.28), "npc"),
      y = unit(c(0.08, 0.08), "npc"),
      gp = gpar(lwd = 3)
    )
  ) +
  annotation_custom(
    grob = textGrob(
      "UMAP2",
      x = unit(0.035, "npc"),
      y = unit(0.18, "npc"),
      rot = 90,
      gp = gpar(fontsize = 15)
    )
  ) +
  annotation_custom(
    grob = textGrob(
      "UMAP1",
      x = unit(0.18, "npc"),
      y = unit(0.035, "npc"),
      gp = gpar(fontsize = 15)
    )
  ) +
   annotation_custom(
    segmentsGrob(
      x0 = unit(0.08, "npc"), y0 = unit(0.08, "npc"),
      x1 = unit(0.08, "npc"), y1 = unit(0.28, "npc"),
      arrow = arrow(length = unit(0.18, "cm")),
      gp = gpar(lwd = 3)
    )
  ) +
  annotation_custom(
    segmentsGrob(
      x0 = unit(0.08, "npc"), y0 = unit(0.08, "npc"),
      x1 = unit(0.28, "npc"), y1 = unit(0.08, "npc"),
      arrow = arrow(length = unit(0.18, "cm")),
      gp = gpar(lwd = 3)
    )
  )
dev.off()
