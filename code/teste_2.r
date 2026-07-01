


library(tidyverse)
library(grid)
library(ggrepel)
library(viridis)

res_dir <- "/Users/denny/Desktop/niels/prof_niels/results/celloracle_hif1a_ko_1500_FA_paperstyle_20260626_154852"
tab_dir <- file.path(res_dir, "tables")
fig_dir <- file.path(res_dir, "figures_R_improved")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

fa_cells <- read_csv(file.path(tab_dir, "fa_embedding_cells.csv"))
cell_shift <- read_csv(file.path(tab_dir, "cell_shift_vectors_FA.csv"))
grid_shift <- read_csv(file.path(tab_dir, "grid_shift_vectors_FA.csv"))

fa_cells <- fa_cells %>%
  mutate(seurat_clusters = as.factor(seurat_clusters))

cell_shift <- cell_shift %>%
  mutate(seurat_clusters = as.factor(seurat_clusters))

set.seed(123)

cell_shift_sampled <- cell_shift %>%
  group_by(seurat_clusters) %>%
  slice_sample(prop = 0.08) %>%
  ungroup()

arrow_scale <- 0.35

p_arrows_clean <- ggplot() +
  geom_point(
    data = fa_cells,
    aes(x = FA1, y = FA2),
    color = "grey85",
    size = 0.4,
    alpha = 0.5
  ) +
  geom_segment(
    data = cell_shift_sampled,
    aes(
      x = FA1,
      y = FA2,
      xend = FA1 + shift_FA1 * arrow_scale,
      yend = FA2 + shift_FA2 * arrow_scale
    ),
    arrow = arrow(length = unit(0.035, "inches")),
    linewidth = 0.25,
    alpha = 0.65,
    color = "black"
  ) +
  coord_equal() +
  theme_classic(base_size = 14) +
  labs(
    title = "HIF1A KO simulated shift",
    x = "FA1",
    y = "FA2"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 18)
  )

p_arrows_clean

ggsave(file.path(fig_dir, "FA_HIF1A_KO_cell_shift_arrows_clean.png"),
       p_arrows_clean, width = 7, height = 6, dpi = 300)



p_delta_highlight <- ggplot(fa_cells, aes(x = FA1, y = FA2)) +
  geom_point(
    aes(color = delta_l2),
    size = 0.65,
    alpha = 0.9
  ) +
  geom_point(
    data = fa_cells %>% filter(seurat_clusters == "4"),
    aes(x = FA1, y = FA2),
    shape = 21,
    color = "black",
    fill = NA,
    size = 1.1,
    stroke = 0.25,
    alpha = 0.55
  ) +
  scale_color_viridis_c(option = "magma", name = "delta L2") +
  coord_equal() +
  theme_classic(base_size = 14) +
  labs(
    title = "HIF1A KO effect magnitude",
    subtitle = "Cluster 4 highlighted",
    x = "FA1",
    y = "FA2"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 18),
    plot.subtitle = element_text(hjust = 0.5, size = 12)
  )

p_delta_highlight

ggsave(file.path(fig_dir, "FA_HIF1A_KO_delta_l2_cluster4_highlight.png"),
       p_delta_highlight, width = 7, height = 6, dpi = 300)