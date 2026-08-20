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