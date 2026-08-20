

source("functions/aux_functions.r")

ssgsea_results = run_pseudobulk_ssgsea_heatmap(
  seurat_obj = SC_obj_sub,
  clusters_use = c("4", "5", "6"),
  min_cells = 9
)

tiff("figures/GSEA_1f.tiff", width = 9, height = 5, units = "in", res = 300)
ComplexHeatmap::draw(
  ssgsea_results$heatmap
)
dev.off()
