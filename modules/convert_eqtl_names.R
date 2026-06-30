#!/usr/bin/env Rscript

# Capture arguments passed from the bash shell
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

# Reusable function adapted to append gene symbols instead of overwriting raw IDs
convert_eqtl_genes <- function(file_path, is_boxplot_file = FALSE) {
  if (!file.exists(file_path)) {
    warning("File not found, skipping: ", file_path)
    return(NULL)
  }
  
  cat("Converting names inside:", basename(file_path), "\n")
  data_table <- read.table(file_path, sep="\t", header=TRUE, stringsAsFactors=FALSE, check.names=FALSE)
  
  # Identify the Ensembl ID column
  if ("geneid" %in% colnames(data_table)) {
    ensembl_ids <- data_table$geneid
    id_col <- "geneid"
  } else if ("gene" %in% colnames(data_table)) {
    ensembl_ids <- data_table$gene
    id_col <- "gene"
  } else if ("gene_id" %in% colnames(data_table)) {
    ensembl_ids <- data_table$gene_id
    id_col <- "gene_id"
  } else {
    ensembl_ids <- rownames(data_table)
    id_col <- NULL
  }
  
  # Strip version decimals (e.g., ENSG00000001460.19 -> ENSG00000001460)
  ensembl_ids_clean <- gsub("\\..*$", "", ensembl_ids)
  
  # Query AnnotationHub for Human OrgDb
  ah <- AnnotationHub()
  human_org <- query(ah, c("OrgDb", "Homo sapiens"))
  orgdb <- ah[[names(human_org)[1]]]
  
  # Map names safely
  gene_symbols <- mapIds(
    orgdb, 
    keys = ensembl_ids_clean, 
    keytype = "ENSEMBL", 
    column = "SYMBOL", 
    multiVals = "first"
  )
  
  # Fall back to original ID if symbol is missing
  final_symbols <- ifelse(is.na(gene_symbols), ensembl_ids, gene_symbols)
  final_symbols <- make.unique(final_symbols)
  
  if (is_boxplot_file) {
    # CRITICAL FIX FOR NEW BOXPLOT SCRIPT:
    # Keep the raw 'geneid' intact for matrix matching, and append the 'gene_symbol' column
    data_table$gene_symbol <- final_symbols
    
    # Reorder columns to ensure layout matches: geneid, gene_symbol, snpid, etc.
    remaining_cols <- setdiff(colnames(data_table), c(id_col, "gene_symbol"))
    data_table <- data_table[, c(id_col, "gene_symbol", remaining_cols), drop = FALSE]
    
  } else {
    # Standard results files can continue overwriting the ID for clean reporting
    if (!is.null(id_col)) {
      data_table[[id_col]] <- final_symbols
    } else {
      rownames(data_table) <- final_symbols
    }
  }
  
  # Save the table right back to its original location
  write.table(data_table, file_path, sep="\t", row.names=FALSE, quote=FALSE)
  cat("Successfully updated:", basename(file_path), "\n\n")
}

# ==============================================================================
# EXECUTION
# ==============================================================================
# Base directory points directly to results
base_dir <- paste0("~/donglab/data/target_ALS/", tissue, "/eQTL/results/")
# Subdirectory where Matrix eQTL / postprocess dumps files before cleanup
sub_dir  <- paste0(base_dir, tissue, "_eQTL/")

# 1. Overwrite standard summary tables in place (check both paths safely)
convert_eqtl_genes(paste0(base_dir, tissue, "_eQTL.full_annotated.txt"), is_boxplot_file=FALSE)
convert_eqtl_genes(paste0(base_dir, tissue, "_eQTL.FDR0.05.txt"), is_boxplot_file=FALSE)
convert_eqtl_genes(paste0(base_dir, tissue, "_eQTL.lead_snps.txt"), is_boxplot_file=FALSE)

# Target the meta-analysis file inside the active subdirectory
meta_path <- paste0(sub_dir, "target_ALS_Combined_Meta_eQTL.txt")
if (!file.exists(meta_path)) {
  meta_path <- paste0(base_dir, "target_ALS_Combined_Meta_eQTL.txt") # Fallback if already cleaned
}
convert_eqtl_genes(meta_path, is_boxplot_file=FALSE)

# 2. Add symbols to the boxplot pairs file while PRESERVING raw Ensembl IDs
convert_eqtl_genes(paste0(base_dir, tissue, "_eQTL.top_for_boxplot.txt"), is_boxplot_file=TRUE)
