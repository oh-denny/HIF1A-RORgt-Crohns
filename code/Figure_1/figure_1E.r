

library(Seurat)
library(dplyr)
library(ggplot2)

DefaultAssay(SC_obj_sub) = "SCT"
Idents(SC_obj_sub) = "SCT_snn_res.0.5"

df_hif = FetchData(
  SC_obj_sub,
  vars = c("HIF1A", "SCT_snn_res.0.5", "group")
)

df_hif = df_hif %>%
  mutate(
    HIF1A_pos = HIF1A > 0
  )

prop_hif_group = df_hif %>%
  group_by(SCT_snn_res.0.5, group) %>%
  summarise(
    total_cells = n(),
    HIF1A_pos_cells = sum(HIF1A_pos),
    percent_HIF1A_pos = 100 * HIF1A_pos_cells / total_cells,
    .groups = "drop"
  )

prop_hif


donut_df = df_hif %>%
  filter(HIF1A_pos) %>%
  group_by(SCT_snn_res.0.5, group) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(SCT_snn_res.0.5) %>%
  mutate(percent = 100 * n / sum(n))

tiff("figures/figure_1e.tiff", width = 6, height = 7, units = "in", res = 300)
ggplot(donut_df, aes(x = 2, y = percent, fill = group)) +
  geom_col(color = "white", width = 1) +
  coord_polar(theta = "y") +
  xlim(0.5, 2.5) +
  facet_wrap(~ SCT_snn_res.0.5, nrow = 2) +
  scale_fill_manual(values = c("Healthy" = "#1F5C5C", "CD" = "#E18487")) +
  theme_void() +
  labs(
    fill = "Group",
    title = "Distribution of HIF1A+ cells by group in each subcluster"
  ) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    strip.text = element_text(face = "bold", size = 12)
  )
dev.off()
