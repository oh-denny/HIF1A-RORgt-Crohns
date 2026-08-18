
   
library(Seurat)
library(dplyr)
library(ggplot2)

feature = "RORgt_Treg_Score1" 
          
cols = c(
  "0" = "#A61C4B",
  "1" = "#0F6C84",
  "2" = "#2C6B6B",
  "3" = "#AEB8B7",
  "4" = "#6C7A3A",
  "5" = "#D8C7AF",
  "6" = "#162E93"
)

rorgt_genes = c("RORC", "IL23R", "CCR6", "MAF", "AHR")
rorgt_genes = rorgt_genes[rorgt_genes %in% rownames(SC_obj_sub)]
rorgt_genes %in% rownames(SC_obj_sub)

SC_obj_sub <- AddModuleScore(
  object = SC_obj_sub,
  features = list(rorgt_genes),
  name = "RORgt_Treg_Score"
)


umap_df = as.data.frame(Embeddings(SC_obj_sub, reduction = "umap"))
colnames(umap_df)[1:2] = c("UMAP_1", "UMAP_2")

meta_df = FetchData(
  SC_obj_sub,
  vars = c(feature, "group", "seurat_clusters")
)

df_plot = cbind(umap_df[rownames(meta_df), ], meta_df)

colnames(df_plot)[colnames(df_plot) == feature] <- "RORgt_score"

df_plot = df_plot %>%
  mutate(
    seurat_clusters = factor(seurat_clusters, levels = names(cols)),
    group = factor(group, levels = c("Healthy", "CD"))
  )


df_plot = df_plot %>%
  group_by(group) %>%
  mutate(
    RORgt_cutoff = quantile(RORgt_score, 0.95, na.rm = TRUE),
    RORgt_high = RORgt_score >= RORgt_cutoff
  ) %>%
  ungroup()


p_rorgt_simple = ggplot(df_plot, aes(x = UMAP_1, y = UMAP_2)) +
  stat_density_2d(
    data = df_plot %>% filter(RORgt_high),
    aes(fill = after_stat(level)),
    geom = "polygon",
    alpha = 0.35,
    color = NA,
    bins = 6
  ) +
  geom_point(
    color = "black",
    size = 0.25,
    alpha = 0.45
  ) +
  facet_wrap(~ group, nrow = 1) +
  scale_fill_gradient(
    low = "#D8C7AF",
    high = "#A61C4B",
    name = "RORγt-high\ndensity"
  ) +
  coord_equal() +
  labs(
    title = "Density of RORγt-high FOXP3+ cells",
    x = "UMAP 1",
    y = "UMAP 2"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    strip.text = element_text(face = "bold", size = 14),
    legend.position = "right"
  )


tiff("figures/figure1b.tiff", width = 7, height = 5.5, units = "in", res = 300)
p_rorgt_simple
dev.off()
