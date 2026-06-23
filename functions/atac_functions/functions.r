


library(Signac)
library(Seurat)
library(Matrix)
library(tidyverse)
library(GenomicRanges)

clean_peaks <- function(peak_file) {

  peak_names <- readLines(peak_file)

  peak_names <- gsub("^b6-", "", peak_names)
  peak_names <- gsub("^cast-", "", peak_names)
  peak_names <- gsub(":", "-", peak_names)
  peak_names <- gsub("\t", "-", peak_names)

  peak_names <- peak_names[!grepl(
    "^chr\\s|^seqnames|^chrom|^start|^end",
    peak_names,
    ignore.case = TRUE
  )]

  peak_names <- peak_names[grepl("^chr.+-[0-9]+-[0-9]+$", peak_names)]

  return(peak_names)
}

read_cell_metadata <- function(meta_file) {
  cell_meta <- read.delim(
    gzfile(meta_file),
    header = FALSE,
    sep = "\t",
    stringsAsFactors = FALSE,
    fill = TRUE,
    check.names = FALSE
  )

  header <- trimws(as.character(cell_meta[1, ]))
  header <- header[!is.na(header) & !header %in% c("", "NA")]

  cell_meta <- cell_meta[-1, , drop = FALSE]

  if (ncol(cell_meta) == length(header) + 1) {
    colnames(cell_meta) <- c("barcode", make.names(header, unique = TRUE))
  } else {
    colnames(cell_meta) <- c(
      "barcode",
      paste0("metadata_", seq_len(ncol(cell_meta) - 1))
    )
  }

  rownames(cell_meta) <- cell_meta$barcode
  cell_meta$barcode <- NULL

  return(cell_meta)
}

create_chromatin_assay <- function(counts, frag_file, smp) {
  if (length(frag_file) == 1) {
    return(CreateChromatinAssay(
      counts = counts,
      sep = c("-", "-"),
      fragments = frag_file,
      min.cells = 10,
      min.features = 200,
      validate.fragments = FALSE
    ))
  }

  if (length(frag_file) > 1) {
    message("Mais de um fragment file para: ", smp)
    message("Arquivos: ", paste(basename(frag_file), collapse = ", "))
    message("Criando ChromatinAssay sem fragments para evitar erro do Signac.")
  } else {
    message("Sem fragment file para: ", smp)
    message("Criando ChromatinAssay sem fragments.")
  }

  CreateChromatinAssay(
    counts = counts,
    sep = c("-", "-"),
    min.cells = 10,
    min.features = 200
  )
}