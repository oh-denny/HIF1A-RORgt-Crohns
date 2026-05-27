
setwd("prof_niels")



source("functions/atac_functions/functions.r")


dir.create("data/out/scATAC_files", recursive = TRUE, showWarnings = FALSE)
dir.create("data/seu_obj", recursive = TRUE, showWarnings = FALSE)

untar(
  tarfile = "data/sc_atac_mice/GSE216910_RAW.tar",
  exdir = "data/out/scATAC_files"
)

base_dir = "data/out/scATAC_files"

mtx_files = list.files(
  base_dir,
  pattern = "\\.mtx\\.gz$",
  full.names = TRUE
)

sample_names = basename(mtx_files) |>
  gsub("\\.mtx\\.gz$", "", x = _)

metadata = tibble(
  sample = c(
    "GSM6697676_Spleen_Colon_Treg",
    "GSM6697677_Gata3_WT_KO_counts",
    "GSM6697678_B6_Cast_F1_Treg_only",
    "GSM6697678_B6_Cast_F1_Treg_Tconv",
    "GSM6697679_FoxP3_WT_KO_Treg_Tconv"
  ),
  experiment = c(
    "Spleen_Colon_Treg",
    "Gata3_WT_KO",
    "B6_Cast_F1_Treg_only",
    "B6_Cast_F1_Treg_Tconv",
    "FoxP3_WT_KO_Treg_Tconv"
  )
)


seurat_list = lapply(sample_names, function(smp) {

  message("\n==============================")
  message("Processando: ", smp)
  message("==============================")

  prefix <- smp |>
    gsub("_counts$", "", x = _)

  short_prefix <- prefix |>
    gsub("_Treg_Tconv$", "", x = _)

  mtx_file <- file.path(base_dir, paste0(smp, ".mtx.gz"))

  peak_file <- list.files(
    base_dir,
    pattern = paste0("^", prefix, "_peaks\\.(tsv|txt)\\.gz$"),
    full.names = TRUE
  )[1]

  barcode_file <- list.files(
    base_dir,
    pattern = paste0("^", prefix, "_barcodes\\.(tsv|txt)\\.gz$"),
    full.names = TRUE
  )[1]

  frag_file <- list.files(
    base_dir,
    pattern = paste0("^", prefix, ".*fragments\\.(tsv|bed)\\.gz$"),
    full.names = TRUE
  )

  if (length(frag_file) == 0 && short_prefix != prefix) {
    frag_file <- list.files(
      base_dir,
      pattern = paste0("^", short_prefix, ".*fragments\\.(tsv|bed)\\.gz$"),
      full.names = TRUE
    )
  }

  meta_file <- list.files(
    base_dir,
    pattern = paste0("^", prefix, "_Metadata\\.txt\\.gz$"),
    full.names = TRUE
  )[1]

  if ((is.na(meta_file) || !file.exists(meta_file)) && short_prefix != prefix) {
    meta_file <- list.files(
      base_dir,
      pattern = paste0("^", short_prefix, "_Metadata\\.txt\\.gz$"),
      full.names = TRUE
    )[1]
  }

  message("MTX: ", mtx_file)
  message("PEAK: ", peak_file)
  message("BARCODE: ", barcode_file)
  message("FRAG: ", paste(frag_file, collapse = "; "))
  message("META: ", meta_file)

  if (is.na(peak_file) || !file.exists(peak_file)) {
    stop("Peak file não encontrado para: ", smp)
  }

  if (is.na(barcode_file) || !file.exists(barcode_file)) {
    stop("Barcode file não encontrado para: ", smp)
  }

  counts <- readMM(mtx_file)

  peak_names <- clean_peaks(peak_file)

  if (length(peak_names) != nrow(counts)) {
    stop(
      "Número de peaks não bate com linhas da matriz em ", smp,
      "\nPeaks: ", length(peak_names),
      "\nMatriz: ", nrow(counts)
    )
  }

  barcodes <- read_tsv(
    barcode_file,
    col_names = FALSE,
    show_col_types = FALSE
  )$X1

  if (length(barcodes) != ncol(counts)) {
    stop(
      "Número de barcodes não bate com colunas da matriz em ", smp,
      "\nBarcodes: ", length(barcodes),
      "\nMatriz: ", ncol(counts)
    )
  }

  rownames(counts) <- peak_names
  colnames(counts) <- barcodes

  chrom_assay <- create_chromatin_assay(counts, frag_file, smp)

  obj <- CreateSeuratObject(
    counts = chrom_assay,
    assay = "peaks",
    project = smp
  )

  obj$sample <- smp
  obj$experiment <- metadata$experiment[match(smp, metadata$sample)]
  obj$fragment_files <- ifelse(
    length(frag_file) > 0,
    paste(basename(frag_file), collapse = ";"),
    NA_character_
  )

  if (!is.na(meta_file) && file.exists(meta_file)) {

    cell_meta <- read_cell_metadata(meta_file)

    common_cells <- intersect(colnames(obj), rownames(cell_meta))

    message("Células em comum com metadata: ", length(common_cells))

    if (length(common_cells) > 0) {

      obj <- subset(obj, cells = common_cells)

      obj <- AddMetaData(
        obj,
        metadata = cell_meta[common_cells, , drop = FALSE]
      )

    } else {

      message("Metadata não bateu com os barcodes para: ", smp)
      message("Mantendo objeto sem metadata adicional.")
    }
  }

  return(obj)
})

names(seurat_list) <- sample_names

saveRDS(
  seurat_list,
  "data/seu_obj/GSE216910_scATAC_seurat_list.rds"
)
