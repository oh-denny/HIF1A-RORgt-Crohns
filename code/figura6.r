












df_long_scores <- df_traj_scores %>%
  dplyr::select(
    cell,
    group,
    seurat_clusters,
    State,
    Pseudotime,
    dplyr::all_of(names(hallmark_list))
  ) %>%
  tidyr::pivot_longer(
    cols = dplyr::all_of(names(hallmark_list)),
    names_to = "pathway",
    values_to = "score"
  ) %>%
  dplyr::mutate(
    pathway = gsub("^HALLMARK_", "", pathway),
    pathway = gsub("_", " ", pathway),
    pathway = stringr::str_to_title(tolower(pathway))
  )

head(df_long_scores)


p_E <- ggplot(
  df_long_scores,
  aes(x = Pseudotime, y = score, color = pathway)
) +
  geom_smooth(
    se = FALSE,
    method = "loess",
    linewidth = 1.1,
    span = 0.7
  ) +
  theme_classic(base_size = 14) +
  labs(
    x = "Pseudotime",
    y = "Module score",
    title = "Inflammatory and hypoxia programs along pseudotime"
  ) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.title = element_blank()
  )

p_E

ggsave(
  "Figure6E_hallmark_scores_along_pseudotime.pdf",
  p_E,
  width = 6,
  height = 4
)


p_E_group <- ggplot(
  df_long_scores,
  aes(x = Pseudotime, y = score, color = group)
) +
  geom_smooth(
    se = FALSE,
    method = "loess",
    linewidth = 1,
    span = 0.7
  ) +
  facet_wrap(~ pathway, scales = "free_y", ncol = 2) +
  scale_color_manual(
    values = c(
      "Healthy" = "#2D6F6B",
      "CD" = "#E9969B"
    )
  ) +
  theme_classic(base_size = 13) +
  labs(
    x = "Pseudotime",
    y = "Module score",
    title = "Pathway activity along pseudotime"
  ) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.title = element_blank(),
    strip.text = element_text(face = "bold")
  )

p_E_group

ggsave(
  "Figure6E_hallmark_scores_along_pseudotime_by_group.pdf",
  p_E_group,
  width = 7,
  height = 5
)


ucell_cols <- paste0(names(hallmark_list), "_UCell")

ucell_cols
ucell_cols %in% colnames(SC_tmp@meta.data)


df_ucell <- SC_tmp@meta.data %>%
  as.data.frame() %>%
  tibble::rownames_to_column("cell") %>%
  dplyr::select(cell, sample, group, seurat_clusters, dplyr::all_of(ucell_cols))

# renomear
colnames(df_ucell) <- c(
  "cell", "sample", "group", "seurat_clusters",
  names(hallmark_list)
)

df_traj_ucell <- df_traj %>%
  dplyr::left_join(
    df_ucell %>% dplyr::select(cell, sample, dplyr::all_of(names(hallmark_list))),
    by = "cell"
  )


ggplot(df_traj_ucell, aes(x = Component1, y = Component2, color = HALLMARK_TNFA_SIGNALING_VIA_NFKB)) +
  geom_point(size = 0.5, alpha = 0.8) +
  scale_color_viridis_c(option = "magma", name = "TNFa") +
  theme_classic(base_size = 14) +
  labs(
    x = "Component 1",
    y = "Component 2"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16)
  )


df_long_ucell <- df_traj_ucell %>%
  tidyr::pivot_longer(
    cols = dplyr::all_of(names(hallmark_list)),
    names_to = "pathway",
    values_to = "score"
  ) %>%
  dplyr::mutate(
    pathway = gsub("^HALLMARK_", "", pathway),
    pathway = gsub("_", " ", pathway),
    pathway = stringr::str_to_title(tolower(pathway))
  )








