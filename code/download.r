

setwd("prof_niels")

library(GEOquery)
library(Seurat)
library(tidyverse)


dir.create("GSE209832", showWarnings = FALSE)
setwd("GSE209832")

getGEOSuppFiles("GSE209832", makeDirectory = TRUE)

untar("data/GSE209832/GSE209832_RAW.tar", exdir = "data/out/H5_files")

h5_files = list.files(
  "data/out/H5_files",
  pattern = "filtered_feature_bc_matrix.h5$",
  full.names = TRUE
)


sample_names = basename(h5_files) |>
  gsub("_GEX_.*", "", x = _)

sample_names

metadata = tibble(
  sample = c(
    "GSM6401749", "GSM6401750", "GSM6401751",
    "GSM6401752", "GSM6401753", "GSM6401754",
    "GSM6401755", "GSM6401756", "GSM6401757",
    "GSM6401758", "GSM6401759"
  ),
  group = c(
    "CD", "Healthy", "CD",
    "CD", "CD", "Healthy",
    "Healthy", "CD", "Healthy",
    "Healthy", "Healthy"
  )
)


seurat_list = lapply(h5_files, function(file) {
  gsm = basename(file) |> gsub("_GEX_.*", "", x = _)
  counts = Read10X_h5(file)
  obj = CreateSeuratObject(
    counts = counts,
    project = gsm,
    min.cells = 3,
    min.features = 200
  )
  obj$sample = gsm
  obj$group = metadata$group[match(gsm, metadata$sample)]
  return(obj)
})

names(seurat_list) = sample_names
saveRDS(seurat_list, "data/seu_obj/seurat_list.rds")