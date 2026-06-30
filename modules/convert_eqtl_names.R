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
convert_eqtl_genes <- function(file_path) {
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
  
  # If a symbol is missing, fall back to the original Ensembl ID
  final_names <- ifelse(is.na(gene_symbols), ensembl_ids, gene_symbols)
  final_names <- make.unique(final_names)
  
  # OVERWRITE the identity column directly so plotting tools pick it up automatically!
  if (!is.null(id_col)) {
    data_table[[id_col]] <- final_names
  } else {
    rownames(data_table) <- final_names
  }
  
  # Save the table right back to its original location
  write.table(data_table, file_path, sep="\t", row.names=FALSE, quote=FALSE)
  cat("Successfully updated:", basename(file_path), "\n\n")
}

# ==============================================================================
# EXECUTION
# ==============================================================================
args <- commandArgs(trailingOnly = TRUE)
tissue <- args[1]
base_dir <- paste0("~/donglab/data/target_ALS/", tissue, "/eQTL/results/")

# Overwrite the files in place
convert_eqtl_genes(paste0(base_dir, tissue, "_eQTL.full_annotated.txt"))
convert_eqtl_genes(paste0(base_dir, tissue, "_eQTL.FDR0.05.txt"))
convert_eqtl_genes(paste0(base_dir, tissue, "_eQTL.lead_snps.txt"))
convert_eqtl_genes(paste0(base_dir, "target_ALS_Combined_Meta_eQTL.txt"))
