
setwd("prof_niels")



seurat_list <- readRDS("data/seu_obj/GSE216910_scATAC_seurat_list.rds")

sapply(seurat_list, ncol)
sapply(seurat_list, nrow)
