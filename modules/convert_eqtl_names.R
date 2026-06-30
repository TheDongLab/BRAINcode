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
  # If the first line contains "ENSG" or a decimal version of it, it's data (headerless)
  has_header <- !grepl("ENSG[0-9]+", first_line)
  
  data_table <- read.table(file_path, sep="\t", header=has_header, stringsAsFactors=FALSE, check.names=FALSE)
  
  # ── Identify Gene Column ────────────────────────────────────────────────────
  if (has_header) {
    if ("geneid" %in% colnames(data_table)) { id_col <- "geneid" }
    else if ("gene" %in% colnames(data_table)) { id_col <- "gene" }
    else if ("gene_id" %in% colnames(data_table)) { id_col <- "gene_id" }
    else { id_col <- 1 } # Fallback to first column
    ensembl_ids <- data_table[[id_col]]
  } else {
    id_col <- 1 # No header means column 1 is the Ensembl ID
    ensembl_ids <- data_table[, 1]
  }
  
  # Clean Ensembl IDs
  ensembl_ids_clean <- gsub("\\..*$", "", ensembl_ids)
  
  # Safety check: Ensure we actually have Ensembl keys to query
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
  
  final_symbols <- ifelse(is.na(gene_symbols), ensembl_ids, gene_symbols)
  final_symbols <- make.unique(final_symbols)
  
  # ── Reconstruct Table ───────────────────────────────────────────────────────
  if (is_boxplot_file) {
    # If the file didn't have headers, assign standard ones for the new boxplot script
    if (!has_header) {
      colnames(data_table) <- c("geneid", "snpid", if(ncol(data_table) > 2) colnames(data_table)[3:ncol(data_table)] else NULL)
      id_col <- "geneid"
    }
    
    # Inject the symbol column right next to geneid
    data_table$gene_symbol <- final_symbols
    remaining_cols <- setdiff(colnames(data_table), c(id_col, "gene_symbol"))
    data_table <- data_table[, c(id_col, "gene_symbol", remaining_cols), drop = FALSE]
    header_to_write <- TRUE # Force header to TRUE so the boxplot script's fread(header=TRUE) works!
    
  } else {
    # Standard reporting summaries just overwrite the column inline
    if (is.numeric(id_col)) {
      data_table[, id_col] <- final_symbols
    } else {
      data_table[[id_col]] <- final_symbols
    }
    header_to_write <- has_header
  }
  
  # Save back
  write.table(data_table, file_path, sep="\t", row.names=FALSE, col.names=header_to_write, quote=FALSE)
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

# 3. Boxplot input pairs file (Now with smart auto-header upgrading)
convert_eqtl_genes(paste0(base_dir, tissue, "_eQTL.top_for_boxplot.txt"), is_boxplot_file=TRUE)
