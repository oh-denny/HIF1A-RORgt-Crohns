
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


