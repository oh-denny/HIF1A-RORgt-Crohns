


library(Seurat)
library(dplyr)
library(ggplot2)
library(ggalluvial)

DefaultAssay(SC_obj_sub) = "SCT"
Idents(SC_obj_sub) = "SCT_snn_res.0.5"

df_hif = FetchData(
  SC_obj_sub,
  vars = c("HIF1A", "SCT_snn_res.0.5", "group", "sample")
)

df_hif = df_hif %>%
  mutate(
    HIF1A_pos = HIF1A > 0,
    HIF1A_status = ifelse(HIF1A_pos, "HIF1A+", "HIF1A-"),
    SCT_snn_res.0.5 = factor(SCT_snn_res.0.5, levels = sort(unique(SCT_snn_res.0.5))),
    group = factor(group, levels = c("Healthy", "CD")),
    HIF1A_status = factor(HIF1A_status, levels = c("HIF1A-", "HIF1A+"))
  )

alluvial_df = df_hif %>%
  count(group, SCT_snn_res.0.5, HIF1A_status, name = "n")

alluvial = ggplot(
  alluvial_df,
  aes(
    axis1 = group,
    axis2 = SCT_snn_res.0.5,
    axis3 = HIF1A_status,
    y = n
  )
) +
  geom_alluvium(
    aes(fill = group),
    alpha = 0.75,
    width = 1/12,
    color = "gray80"
  ) +
  geom_stratum(
    width = 1/12,
    color = "gray40",
    fill = "gray95"
  ) +
   geom_stratum(
    aes(
      fill = after_stat(
        ifelse(x == 2, as.character(stratum), NA)
      )
    ),
    width = 1/12,
    color = "gray40"
  ) +
  scale_x_discrete(
    limits = c("Group", "FOXP3+ subcluster", "HIF1A status"),
    expand = c(0.08, 0.08)
  ) +
  scale_fill_manual(
    values = c(
      "Healthy" = "#1F5C5C",
      "CD" = "#E18487",

      # cores dos subclusters
      "0" = "#A61C4B",
      "1" = "#0F6C84",
      "2" = "#2C6B6B",
      "3" = "#AEB8B7",
      "4" = "#6C7A3A",
      "5" = "#D8C7AF",
      "6" = "#162E93"
    ),
    na.value = "gray95"
  ) +
  labs(
    title = "Flow of FOXP3+ cells by group, subcluster and HIF1A expression",
    y = "Number of cells",
    fill = "Group"
  ) +
  theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.title.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )


tiff("figures/figure_1d.tiff", width = 6, height = 7, units = "in", res = 300)
alluvial
dev.off()
