

library(Seurat)

SC_obj_sub$group = factor(SC_obj_sub$group, levels = c("CD", "Healthy"))
genes = c("IL2RA", "IKZF2", "ICOS", "CTLA4", "PDCD1", "FOXP3", "LAG3", "HAVCR2", "AHR", "HIF1A", "RORA")


tiff("figures/figure1c.tiff", width = 8, height = 5, units = "in", res = 300)            
VlnPlot(
  SC_obj_sub,
  features = c("IL2RA", "IKZF2", "ICOS", "CTLA4", "PDCD1", "FOXP3", "LAG3", "HAVCR2", "AHR", "HIF1A", "RORA"),
  stack = T,
  split.by = "group"
) +
    scale_fill_manual(values = c("Healthy" = "#1F5C5C", "CD" = "#E18487"))
dev.off()
