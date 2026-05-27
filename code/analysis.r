
library(Seurat)
library(ggplot2)
library(grid)

seurat_list = readRDS("data/seu_obj/seurat_list.rds")

seurat_list = lapply(seurat_list, function(obj) {
  obj[["percent.mt"]] = PercentageFeatureSet(obj, pattern = "^MT-")
  obj
})

VlnPlot(
  seurat_list[[1]],
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  ncol = 3
)

seurat_list_filtered = lapply(seurat_list, function(obj) {
  subset(
    obj,
    subset =
      nFeature_RNA > 500 &
      nFeature_RNA < 6000 &
      percent.mt < 15
  )
})

seurat_list_sct = lapply(seurat_list_filtered, function(obj) {
  SCTransform(
    obj,
    assay = "RNA",
    new.assay.name = "SCT",
    vars.to.regress = "percent.mt",
    vst.flavor = "v2",
    verbose = TRUE
  )
})


SC_obj = merge(
  x = seurat_list_sct[[1]],
  y = seurat_list_sct[-1],
  add.cell.ids = names(seurat_list_sct)
)

DefaultAssay(SC_obj) = "SCT"
saveRDS(SC_obj, "data/seu_obj/SC_obj_merged_SCT.rds")

SC_obj = readRDS("data/seu_obj/SC_obj_merged_SCT.rds")


table(SC_obj$group)
table(SC_obj$sample)

features = SelectIntegrationFeatures(
  object.list = seurat_list_sct,
  nfeatures = 3000
)

features = SelectIntegrationFeatures(SC_obj, nfeatures = 3000)
VariableFeatures(SC_obj) = features
VariableFeatures(SC_obj) = rownames(SC_obj[["SCT"]]@scale.data)

SC_obj = RunPCA(SC_obj, assay = "SCT", verbose = FALSE, features = features)
ElbowPlot(SC_obj)

SC_obj = FindNeighbors(SC_obj, dims = 1:15)
SC_obj = FindClusters(SC_obj, resolution = 0.5)
SC_obj = RunUMAP(SC_obj, dims = 1:15)

saveRDS(SC_obj, "data/seu_obj/SC_obj_clustered.rds")
SC_obj = readRDS("data/seu_obj/SC_obj_clustered.rds")


httpgd::hgd()

DimPlot(SC_obj, group.by = "group")
DimPlot(SC_obj, group.by = "seurat_clusters", label = TRUE)




### MY TEST
library(Seurat)
library(ggplot2)
library(ggforce)
library(dplyr)

Idents(SC_obj) <- "seurat_clusters"

umap_df <- as.data.frame(Embeddings(SC_obj, "umap"))
colnames(umap_df)[1:2] <- c("UMAP_1", "UMAP_2")
umap_df$cluster <- Idents(SC_obj)

cluster1_df = umap_df |> 
  filter(cluster == "1")

nrow(cluster1_df)

p <- FeaturePlot(
  SC_obj,
  features = "FOXP3",
  cols = c("lightgrey", "#408A71"),
  order = TRUE,
  min.cutoff = 0,
  max.cutoff = "q95"
) +
  NoAxes() +
  NoLegend() +
  theme_void() +
  geom_density_2d(
    data = cluster1_df,
    aes(x = UMAP_1, y = UMAP_2),
    inherit.aes = FALSE,
    color = "black",
    linewidth = .5,
    linetype = "dashed",
    bins = 5
  ) +
  annotate(
    "text",
    x = -8.8, y = 4,
    label = "FOXP3+ T-Reg\nCluster 1",
    size = 4
  )

p + 
  annotation_custom(
    grob = linesGrob(
      x = unit(c(0.08, 0.08), "npc"),
      y = unit(c(0.08, 0.28), "npc"),
      gp = gpar(lwd = 3)
    )
  ) +
  annotation_custom(
    grob = linesGrob(
      x = unit(c(0.08, 0.28), "npc"),
      y = unit(c(0.08, 0.08), "npc"),
      gp = gpar(lwd = 3)
    )
  ) +
  annotation_custom(
    grob = textGrob(
      "UMAP2",
      x = unit(0.035, "npc"),
      y = unit(0.18, "npc"),
      rot = 90,
      gp = gpar(fontsize = 15)
    )
  ) +
  annotation_custom(
    grob = textGrob(
      "UMAP1",
      x = unit(0.18, "npc"),
      y = unit(0.035, "npc"),
      gp = gpar(fontsize = 15)
    )
  ) +
   annotation_custom(
    segmentsGrob(
      x0 = unit(0.08, "npc"), y0 = unit(0.08, "npc"),
      x1 = unit(0.08, "npc"), y1 = unit(0.28, "npc"),
      arrow = arrow(length = unit(0.18, "cm")),
      gp = gpar(lwd = 3)
    )
  ) +
  annotation_custom(
    segmentsGrob(
      x0 = unit(0.08, "npc"), y0 = unit(0.08, "npc"),
      x1 = unit(0.28, "npc"), y1 = unit(0.08, "npc"),
      arrow = arrow(length = unit(0.18, "cm")),
      gp = gpar(lwd = 3)
    )
  )







#### reclsutering without clusters 10 and 12

SC_obj_sub = subset(
  SC_obj,
  idents = c("1"),
)

SC_obj_sub$seurat_clusters = droplevels(SC_obj_sub$seurat_clusters)

DefaultAssay(SC_obj_sub) = "SCT"

SC_obj_sub = RunPCA(
  SC_obj_sub,
  assay = "SCT",
  features = VariableFeatures(SC_obj_sub),
  verbose = FALSE
)
ElbowPlot(SC_obj_sub)

SC_obj_sub = FindNeighbors(SC_obj_sub, dims = 1:10)
SC_obj_sub = FindClusters(SC_obj_sub,resolution = 0.5)
SC_obj_sub = RunUMAP(SC_obj_sub, dims = 1:10)

cols = c(
  "#A61C4B",
  "#0F6C84",
  "#2C6B6B",
  "#AEB8B7",
  "#6C7A3A",
  "#D8C7AF",
  "#162E93"
)

DimPlot(SC_obj_sub, group.by = "group")
DimPlot(SC_obj_sub, group.by = "seurat_clusters",
  label = FALSE,
  cols = cols) +
  NoAxes() +
  NoLegend() +
  theme_void() + 
  annotation_custom(
    grob = linesGrob(
      x = unit(c(0.08, 0.08), "npc"),
      y = unit(c(0.08, 0.28), "npc"),
      gp = gpar(lwd = 3)
    )
  ) +
  annotation_custom(
    grob = linesGrob(
      x = unit(c(0.08, 0.28), "npc"),
      y = unit(c(0.08, 0.08), "npc"),
      gp = gpar(lwd = 3)
    )
  ) +
  annotation_custom(
    grob = textGrob(
      "UMAP2",
      x = unit(0.035, "npc"),
      y = unit(0.18, "npc"),
      rot = 90,
      gp = gpar(fontsize = 15)
    )
  ) +
  annotation_custom(
    grob = textGrob(
      "UMAP1",
      x = unit(0.18, "npc"),
      y = unit(0.035, "npc"),
      gp = gpar(fontsize = 15)
    )
  ) +
   annotation_custom(
    segmentsGrob(
      x0 = unit(0.08, "npc"), y0 = unit(0.08, "npc"),
      x1 = unit(0.08, "npc"), y1 = unit(0.28, "npc"),
      arrow = arrow(length = unit(0.18, "cm")),
      gp = gpar(lwd = 3)
    )
  ) +
  annotation_custom(
    segmentsGrob(
      x0 = unit(0.08, "npc"), y0 = unit(0.08, "npc"),
      x1 = unit(0.28, "npc"), y1 = unit(0.08, "npc"),
      arrow = arrow(length = unit(0.18, "cm")),
      gp = gpar(lwd = 3)
    )
  )







FeaturePlot(SC_obj_sub, features = c("FOXP3", "HIF1A"), 
            cols = c("lightgrey", "#000B58"),
            order = TRUE,
            max.cutoff = "q95")

head(SC_obj_sub)

df = table(SC_obj_sub$group) |> data.frame()
df$Var1 = factor(df$Var1, levels = c("Healthy", "CD"))

ggplot(df, aes(x = Var1, y = Freq, fill = Var1)) +
  geom_bar(stat = "identity", width = 0.4) +
  theme_test(18) +
  labs(x = "Group", y = "Number of Cells", title = "Cell Count per Group in Cluster 1") +
  theme(legend.position = "none")+
  scale_fill_manual(values = c("Healthy" = "#E18487", "CD" = "#1F5C5C"))





### HIF1A+ T-Reg cluster 10
library(Seurat)
library(dplyr)
library(ggplot2)

DefaultAssay(SC_obj_sub) <- "SCT"
Idents(SC_obj_sub) <- "seurat_clusters"



df_hif <- FetchData(
  SC_obj_sub,
  vars = c("HIF1A", "seurat_clusters", "group")
)

df_hif <- df_hif %>%
  mutate(
    HIF1A_pos = HIF1A > 0
  )

prop_hif_group <- df_hif %>%
  group_by(seurat_clusters, group) %>%
  summarise(
    total_cells = n(),
    HIF1A_pos_cells = sum(HIF1A_pos),
    percent_HIF1A_pos = 100 * HIF1A_pos_cells / total_cells,
    .groups = "drop"
  )

prop_hif


donut_df = df_hif %>%
  filter(HIF1A_pos) %>%
  group_by(seurat_clusters, group) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(seurat_clusters) %>%
  mutate(percent = 100 * n / sum(n))

ggplot(donut_df, aes(x = 2, y = percent, fill = group)) +
  geom_col(color = "white", width = 1) +
  coord_polar(theta = "y") +
  xlim(0.5, 2.5) +
  facet_wrap(~ seurat_clusters, nrow = 1) +
  scale_fill_manual(values = c("Healthy" = "#E18487", "CD" = "#1F5C5C")) +
  theme_void() +
  labs(
    fill = "Group",
    title = "Distribution of HIF1A+ cells by group in each subcluster"
  ) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    strip.text = element_text(face = "bold", size = 12)
  )
