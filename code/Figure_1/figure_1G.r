
source("functions/aux_functions.r")

DefaultAssay(SC_obj_sub) <- "RNA"

SC_obj_pb <- JoinLayers(
  object = SC_obj_sub,
  assay = "RNA"
)
counts_mat <- GetAssayData(
  SC_obj_pb,
  assay = "RNA",
  layer = "counts"
)

meta_all <- SC_obj_pb@meta.data %>%
  as.data.frame() %>%
  rownames_to_column("cell") %>%
  mutate(
    cluster = as.character(seurat_clusters),
    sample_id = sample,
    group = factor(group, levels = c("Healthy", "CD"))
  )

res_c5 <- run_pb_deseq_cluster(
  cluster_use = "5",
  counts_mat = counts_mat,
  meta_all = meta_all,
  min_cells = 9,
  outdir = "pseudobulk_clusters"
)

res_c4 <- run_pb_deseq_cluster(
  cluster_use = "4",
  counts_mat = counts_mat,
  meta_all = meta_all,
  min_cells = 9,
  outdir = "pseudobulk_clusters"
)

p_hypoxia_c5 <- plot_enrichment_cluster(fgsea_c5, "HALLMARK_HYPOXIA", "5")
p_tnfa_c5    <- plot_enrichment_cluster(fgsea_c5, "HALLMARK_TNFA_SIGNALING_VIA_NFKB", "5")

p_hypoxia_c4 <- plot_enrichment_cluster(fgsea_c4, "HALLMARK_HYPOXIA", "4")
p_tnfa_c4    <- plot_enrichment_cluster(fgsea_c4, "HALLMARK_TNFA_SIGNALING_VIA_NFKB", "4")

tiff("figures/hypoxia_c5.tiff", width = 6, height = 4, units = "in", res = 300)
p_hypoxia_c5
dev.off()

tiff("figures/tnfa_c5.tiff", width = 6, height = 4, units = "in", res = 300)
p_tnfa_c5
dev.off()

tiff("figures/hypoxia_c4.tiff", width = 6, height = 4, units = "in", res = 300)
p_hypoxia_c4
dev.off()

tiff("figures/tnfa_c4.tiff", width = 6, height = 4, units = "in", res = 300)
p_tnfa_c4
dev.off()

