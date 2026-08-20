

min_nonzero <- min(
  c(rorgt_summary$Healthy, rorgt_summary$CD)[
    c(rorgt_summary$Healthy, rorgt_summary$CD) > 0
  ]
)

pseudocount <- min_nonzero / 2

pseudocount

rorgt_summary <- rorgt_summary %>%
  mutate(
    log2FC_plot = log2(
      (CD + pseudocount) /
      (Healthy + pseudocount)
    )
  )

rorgt_summary


cols <- c(
  "0" = "#A61C4B",
  "1" = "#0F6C84",
  "2" = "#2C6B6B",
  "3" = "#AEB8B7",
  "4" = "#6C7A3A",
  "5" = "#D8C7AF",
  "6" = "#162E93"
)

tiff("figures/supp_1.tiff", width = 5, height = 5.5, units = "in", res = 300)
ggplot(
  rorgt_summary,
  aes(
    x = seurat_clusters,
    y = log2FC_plot,
    fill = seurat_clusters
  )
) +
  geom_col(
    width = 0.7,
    color = "black",
    linewidth = 0.3
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.5
  ) +
  scale_fill_manual(values = cols) +
  labs(
    x = expression(FOXP3^"+"*" subcluster"),
    y = expression(log[2]*"FC of ROR"*gamma*"t-high cell abundance (CD / Healthy)")
  ) +
  guides(fill = "none") +
  theme_classic(base_size = 18)
dev.off()




rorgt_sample <- df_plot %>%
  group_by(sample, group) %>%
  summarise(
    mean_RORgt_score = mean(RORgt_score, na.rm = TRUE),
    median_RORgt_score = median(RORgt_score, na.rm = TRUE),
    .groups = "drop"
  )


tiff("figures/supp2.tiff", width = 5, height = 5.5, units = "in", res = 300)
ggplot(
  rorgt_sample,
  aes(x = group, y = mean_RORgt_score)
) +
  geom_boxplot(
    width = 0.5,
    outlier.shape = NA
  ) +
  geom_jitter(
    width = 0.08,
    size = 3
  ) +
  labs(
    x = NULL,
    y = "Mean RORγt Treg score"
  ) +
  theme_classic(base_size = 18)
dev.off()



tiff("figures/supp3.tiff", width = 5, height = 5.5, units = "in", res = 300)
ggplot(
  rorgt_donor,
  aes(x = group, y = pct_RORgt_high)
) +
  geom_boxplot(
    width = 0.45,
    outlier.shape = NA
  ) +
  geom_jitter(
    width = 0.08,
    size = 3
  ) +
  labs(
    x = NULL,
    y = expression("ROR"*gamma*"t-high FOXP3+"*" cells (%)")
  ) +
  theme_classic(base_size = 18)
dev.off()
