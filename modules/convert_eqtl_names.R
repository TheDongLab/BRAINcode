#!/usr/bin/env Rscript

# Capture arguments passed from the bash shell
args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
  stop("Error: No tissue name provided. Usage: Rscript convert_eqtl_names.R <TISSUE>", call. = FALSE)
}

# The first argument passed will be our tissue variable
tissue <- args[1]
cat("--- Running gene name conversion for tissue:", tissue, "---\n")

library(AnnotationHub)
library(AnnotationDbi)

# Reusable function adapted to accept full exact paths
convert_eqtl_genes <- function(file_path, output_path) {
  if (!file.exists(file_path)) {
    warning("File not found, skipping: ", file_path)
    return(NULL)
  }
  
  cat("Processing:", basename(file_path), "\n")
  data_table <- read.table(file_path, sep="\t", header=TRUE, stringsAsFactors=FALSE, check.names=FALSE)
  
  # Identify Ensembl ID column (added 'geneid' to the check)
  if ("geneid" %in% colnames(data_table)) {
    ensembl_ids <- data_table$geneid
  } else if ("gene" %in% colnames(data_table)) {
    ensembl_ids <- data_table$gene
  } else if ("gene_id" %in% colnames(data_table)) {
    ensembl_ids <- data_table$gene_id
  } else {
    ensembl_ids <- rownames(data_table)
  }
  
  # Strip version decimals (e.g., ENSG00000001460.19 -> ENSG00000001460)
  ensembl_ids_clean <- gsub("\\..*$", "", ensembl_ids)
  
  # Query AnnotationHub
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
  
  # Handle NAs and uniqueness
  final_names <- ifelse(is.na(gene_symbols), ensembl_ids, gene_symbols)
  final_names <- make.unique(final_names)
  
  # Inject column at the front
  data_table$gene_symbol <- final_names
  col_order <- c("gene_symbol", colnames(data_table)[colnames(data_table) != "gene_symbol"])
  data_table <- data_table[, col_order]
  
  write.table(data_table, output_path, sep="\t", row.names=FALSE, quote=FALSE)
  cat("Saved to:", basename(output_path), "\n\n")
}

# ==============================================================================
# DYNAMIC PATH EXECUTION
# ==============================================================================
# Construct the path dynamically based on the passed tissue string
base_dir <- paste0("~/donglab/data/target_ALS/", tissue, "/eQTL/results/")

# Define the targets dynamically using the tissue variable
convert_eqtl_genes(
  file_path   = paste0(base_dir, tissue, "_eQTL.FDR0.05.txt"),
  output_path = paste0(base_dir, tissue, "_eQTL.FDR0.05.with_symbols.txt")
)

convert_eqtl_genes(
  file_path   = paste0(base_dir, tissue, "_eQTL.lead_snps.txt"),
  output_path = paste0(base_dir, tissue, "_eQTL.lead_snps.with_symbols.txt")
)

# Meta-analysis file (assuming the filename layout remains consistent)
convert_eqtl_genes(
  file_path   = paste0(base_dir, "target_ALS_Combined_Meta_eQTL.txt"),
  output_path = paste0(base_dir, "target_ALS_Combined_Meta_eQTL.with_symbols.txt")
)
