

source("functions/aux_functions.r")

ssgsea_results = run_pseudobulk_ssgsea_heatmap(
  seurat_obj = SC_obj_sub,
  clusters_use = c("4", "5", "6"),
  min_cells = 9
)

tiff("figures/GSEA_1f.tiff", width = 10, height = 7, units = "in", res = 300)
ComplexHeatmap::draw(
  ssgsea_results$heatmap,
  heatmap_legend_side = "bottom",
  annotation_legend_side = "bottom",
  padding = grid::unit(
    c(2, 20, 2, 50),
    "mm"
  )
)
dev.off()
