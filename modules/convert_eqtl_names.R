#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("Error: No tissue name provided. Usage: Rscript convert_eqtl_names.R <TISSUE>", call. = FALSE)
}

tissue <- args[1]
cat("--- Running gene name conversion for tissue:", tissue, "---\n")

suppressPackageStartupMessages({
  library(AnnotationHub)
  library(AnnotationDbi)
})

convert_eqtl_genes <- function(file_path, is_boxplot_file = FALSE) {
  if (!file.exists(file_path)) {
    warning("File not found, skipping: ", file_path)
    return(NULL)
  }
  
  cat("Processing:", basename(file_path), "\n")
  
  # ── Smart Header Detection ──────────────────────────────────────────────────
  first_line <- readLines(file_path, n = 1)
  has_header <- !grepl("ENSG[0-9]+", first_line)
  
  data_table <- read.table(file_path, sep="\t", header=has_header, stringsAsFactors=FALSE, check.names=FALSE)
  
  # ── Identify Gene Column ────────────────────────────────────────────────────
  if (has_header) {
    if ("geneid" %in% colnames(data_table)) { id_col <- "geneid" }
    else if ("gene" %in% colnames(data_table)) { id_col <- "gene" }
    else if ("gene_id" %in% colnames(data_table)) { id_col <- "gene_id" }
    else { id_col <- 1 }
    ensembl_ids <- data_table[[id_col]]
  } else {
    id_col <- 1
    ensembl_ids <- data_table[, 1]
  }
  
  # Clean Ensembl IDs for database matching
  ensembl_ids_clean <- gsub("\\..*$", "", ensembl_ids)
  
  # Safety check: Verify valid keys exist
  if (!any(grepl("^ENSG", ensembl_ids_clean))) {
    warning("No valid Ensembl IDs found in column ", id_col, " of ", basename(file_path), ". Skipping conversion.")
    return(NULL)
  }
  
  # ── Annotation Query ────────────────────────────────────────────────────────
  ah <- AnnotationHub()
  human_org <- query(ah, c("OrgDb", "Homo sapiens"))
  orgdb <- ah[[names(human_org)[1]]]
  
  gene_symbols <- mapIds(
    orgdb, 
    keys = ensembl_ids_clean, 
    keytype = "ENSEMBL", 
    column = "SYMBOL", 
    multiVals = "first"
  )
  
  # Fallback to the original Ensembl ID (with version) if no mapping symbol exists
  final_symbols <- ifelse(is.na(gene_symbols), ensembl_ids, gene_symbols)
  
  # ── Overwrite Column In Place ───────────────────────────────────────────────
  if (is.numeric(id_col)) {
    data_table[, id_col] <- final_symbols
  } else {
    data_table[[id_col]] <- final_symbols
  }
  
  # Force standard headers on headerless boxplot files so downstream scripts parse cleanly
  if (!has_header && is_boxplot_file) {
    colnames(data_table) <- c("geneid", "snpid", if(ncol(data_table) > 2) colnames(data_table)[3:ncol(data_table)] else NULL)
    has_header <- TRUE
  }
  
  # Save back cleanly
  write.table(data_table, file_path, sep="\t", row.names=FALSE, col.names=has_header, quote=FALSE)
  cat("Successfully updated:", basename(file_path), "\n\n")
}

# ==============================================================================
# EXECUTION
# ==============================================================================
base_dir <- paste0("~/donglab/data/target_ALS/", tissue, "/eQTL/results/")
sub_dir  <- paste0(base_dir, tissue, "_eQTL/")

# 1. Standard summary tables
convert_eqtl_genes(paste0(base_dir, tissue, "_eQTL.full_annotated.txt"), is_boxplot_file=FALSE)
convert_eqtl_genes(paste0(base_dir, tissue, "_eQTL.FDR0.05.txt"), is_boxplot_file=FALSE)
convert_eqtl_genes(paste0(base_dir, tissue, "_eQTL.lead_snps.txt"), is_boxplot_file=FALSE)

# 2. Meta-analysis path handling
meta_path <- paste0(sub_dir, "target_ALS_Combined_Meta_eQTL.txt")
if (!file.exists(meta_path)) { meta_path <- paste0(base_dir, "target_ALS_Combined_Meta_eQTL.txt") }
convert_eqtl_genes(meta_path, is_boxplot_file=FALSE)

# 3. Boxplot input pairs file (Processed uniformly inline)
convert_eqtl_genes(paste0(base_dir, tissue, "_eQTL.top_for_boxplot.txt"), is_boxplot_file=TRUE)
