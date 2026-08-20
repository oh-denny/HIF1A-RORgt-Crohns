
library(viridis)

cds = readRDS("data/seu_obj/monocle2_cds_ordered.rds")

table(
  pData(cds)$State,
  pData(cds)$seurat_clusters
)

pData(cds)$seurat_clusters = factor(
  pData(cds)$seurat_clusters,
  levels = names(cols)
)

p_pseudotime = monocle::plot_cell_trajectory(
  cds,
  color_by = "Pseudotime"
) +
  scale_color_viridis_c(option = "magma", direction = 1) +
  theme_classic(base_size = 14) +
  labs(color = "Pseudotime")

tiff("figures/pseudotime.tiff", width = 5, height = 4, units = "in", res = 300)
p_pseudotime
dev.off()

p_state = monocle::plot_cell_trajectory(
  cds,
  color_by = "State"
) +
  #scale_color_manual(values = c("Healthy" = "#1F5C5C", "CD" = "#E18487")) +
  theme_classic(base_size = 14) +
  labs(color = "State")

tiff("figures/state.tiff", width = 5, height = 4, units = "in", res = 300)
p_state
dev.off()


