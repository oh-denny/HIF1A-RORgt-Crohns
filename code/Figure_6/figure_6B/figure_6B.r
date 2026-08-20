

cds = readRDS("data/seu_obj/monocle2_cds_ordered.rds")

table(
  pData(cds)$State,
  pData(cds)$seurat_clusters
)

cols = c(
  "0" = "#A61C4B",
  "1" = "#0F6C84",
  "2" = "#2C6B6B",
  "3" = "#AEB8B7",
  "4" = "#6C7A3A",
  "5" = "#D8C7AF",
  "6" = "#162E93"
)

pData(cds)$seurat_clusters <- factor(
  pData(cds)$seurat_clusters,
  levels = names(cols)
)


p_cluster = monocle::plot_cell_trajectory(
  cds,
  color_by = "seurat_clusters"
) +
  scale_color_manual(values = cols) +
  theme_classic(base_size = 14) +
  labs(color = "Cluster")

tiff("figures/cluster.tiff", width = 5, height = 4, units = "in", res = 300)
p_cluster
dev.off()

p_group = monocle::plot_cell_trajectory(
  cds,
  color_by = "group"
) +
  scale_color_manual(values = c("Healthy" = "#1F5C5C", "CD" = "#E18487")) +
  theme_classic(base_size = 14) +
  labs(color = "group")

tiff("figures/group.tiff", width = 5, height = 4, units = "in", res = 300)
p_group
dev.off()
