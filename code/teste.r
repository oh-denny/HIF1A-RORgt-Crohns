



library(Seurat)
library(dplyr)
library(ggplot2)

feature <- "RORgt_Treg_Score1" 
# feature <- "RORC"  # humano, caso seja gene
# feature <- "Rorc"  # mouse, caso seja gene

cols <- c(
  "0" = "#A61C4B",
  "1" = "#0F6C84",
  "2" = "#2C6B6B",
  "3" = "#AEB8B7",
  "4" = "#6C7A3A",
  "5" = "#D8C7AF",
  "6" = "#162E93"
)

umap_df <- as.data.frame(Embeddings(SC_obj_sub, reduction = "umap"))
colnames(umap_df)[1:2] <- c("UMAP_1", "UMAP_2")

meta_df <- FetchData(
  SC_obj_sub,
  vars = c(feature, "group", "seurat_clusters")
)

df_plot <- cbind(umap_df[rownames(meta_df), ], meta_df)

colnames(df_plot)[colnames(df_plot) == feature] <- "RORgt_score"

df_plot <- df_plot %>%
  mutate(
    seurat_clusters = factor(seurat_clusters, levels = names(cols)),
    group = factor(group, levels = c("Healthy", "CD"))
  )


df_plot <- df_plot %>%
  group_by(group) %>%
  mutate(
    RORgt_cutoff = quantile(RORgt_score, 0.95, na.rm = TRUE),
    RORgt_high = RORgt_score >= RORgt_cutoff
  ) %>%
  ungroup()


df_density <- df_plot %>%
  filter(RORgt_high) %>%
  group_by(group, seurat_clusters) %>%
  filter(n() >= 20) %>%
  ungroup()

p_rorgt_region <- ggplot(df_plot, aes(x = UMAP_1, y = UMAP_2)) +
  stat_density_2d(
    data = df_density,
    aes(
      fill = seurat_clusters,
      group = interaction(group, seurat_clusters)
    ),
    geom = "polygon",
    alpha = 0.25,
    color = NA,
    bins = 4,
    contour_var = "ndensity"
  ) +
  geom_point(
    color = "black",
    size = 0.25,
    alpha = 0.45
  ) +
  #facet_wrap(~ group, nrow = 1) +
  scale_fill_manual(values = cols, drop = FALSE) +
  coord_equal() +
  labs(
    title = "Spatial distribution of RORγt-high FOXP3+ cells",
    fill = "Cluster",
    x = "UMAP 1",
    y = "UMAP 2"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    strip.text = element_text(face = "bold", size = 14),
    legend.position = "bottom"
  )

p_rorgt_region






library(ggforce)

df_hull <- df_plot %>%
  filter(RORgt_high) %>%
  group_by(group, seurat_clusters) %>%
  filter(n() >= 10) %>%
  ungroup()

p_rorgt_hull <- ggplot(df_plot, aes(x = UMAP_1, y = UMAP_2)) +
  geom_mark_hull(
    data = df_hull,
    aes(
      fill = seurat_clusters,
      color = seurat_clusters,
      group = interaction(group, seurat_clusters)
    ),
    alpha = 0.18,
    expand = unit(2, "mm"),
    radius = unit(2, "mm"),
    show.legend = TRUE
  ) +
  geom_point(
    color = "black",
    size = 0.25,
    alpha = 0.50
  ) +
  #facet_wrap(~ group, nrow = 1) +
  scale_fill_manual(values = cols, drop = FALSE) +
  scale_color_manual(values = cols, drop = FALSE) +
  coord_equal() +
  labs(
    title = "RORγt-high regions across FOXP3+ subclusters",
    fill = "Cluster",
    color = "Cluster",
    x = "UMAP 1",
    y = "UMAP 2"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    strip.text = element_text(face = "bold", size = 14),
    legend.position = "bottom"
  )

p_rorgt_hull



p_rorgt_simple <- ggplot(df_plot, aes(x = UMAP_1, y = UMAP_2)) +
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

p_rorgt_simple








p_rorgt_simple <- ggplot(df_plot, aes(x = UMAP_1, y = UMAP_2)) +
  stat_density_2d(
    data = df_plot %>% filter(RORgt_high),
    aes(fill = after_stat(level)),
    geom = "polygon",
    contour = TRUE,
    contour_var = "ndensity",
    bins = 10,
    alpha = 0.75,
    color = NA,
    h = c(0.5, 0.5)
  ) +
  geom_point(
    color = "black",
    size = 0.22,
    alpha = 0.38
  ) +
  facet_wrap(~ group, nrow = 1) +
  scale_fill_gradientn(
    colors = c("#F1E4D4", "#E0A08F", "#C94F6B", "#A61C4B"),
    limits = c(0, 1),
    oob = scales::squish,
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

p_rorgt_simple
