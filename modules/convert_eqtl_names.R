#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("Error: No tissue name provided. Usage: Rscript convert_eqtl_names.R <TISSUE> [standard|interaction]", call. = FALSE)
}

tissue   <- args[1]
run_type <- if (length(args) >= 2) args[2] else "standard"

cat("--- Running gene name conversion ---\n")
cat("Tissue:", tissue, "\n")
cat("Run Type:", run_type, "\n---\n")

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
  
  # Clean Ensembl IDs for database matching (Strips trailing .versions or suffixes)
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
  
  # ── Structural Output Logic ──────────────────────────────────────────────────
  if (is_boxplot_file) {
    new_table <- data.frame(
      geneid      = ensembl_ids,
      gene_symbol = final_symbols,
      snpid       = data_table[, 2],
      stringsAsFactors = FALSE
    )
    if (ncol(data_table) > 2) {
      new_table <- cbind(new_table, data_table[, 3:ncol(data_table), drop=FALSE])
    }
    data_table <- new_table
    has_header <- TRUE
  } else {
    if (is.numeric(id_col)) {
      data_table[, id_col] <- final_symbols
    } else {
      data_table[[id_col]] <- final_symbols
    }
  }
  
  write.table(data_table, file_path, sep="\t", row.names=FALSE, col.names=has_header, quote=FALSE)
  cat("Successfully updated:", basename(file_path), "\n\n")
}

# ==============================================================================
# TARGET EXECUTION PATHS
# ==============================================================================
if (run_type == "interaction") {
  base_dir  <- paste0("~/donglab/data/target_ALS/", tissue, "/eQTL/interaction_results/")
  sub_dir   <- paste0(base_dir, tissue, "_eQTL/")
  meta_name <- "target_ALS_Combined_Meta_Interaction_eQTL.txt"
} else {
  base_dir  <- paste0("~/donglab/data/target_ALS/", tissue, "/eQTL/results/")
  sub_dir   <- paste0(base_dir, tissue, "_eQTL/")
  meta_name <- "target_ALS_Combined_Meta_eQTL.txt"
}

# Run processing on exactly what is needed for this specific pipeline run
convert_eqtl_genes(paste0(base_dir, tissue, "_eQTL.cis.txt"), is_boxplot_file=FALSE)
convert_eqtl_genes(paste0(base_dir, tissue, "_eQTL.full_annotated.txt"), is_boxplot_file=FALSE)
convert_eqtl_genes(paste0(base_dir, tissue, "_eQTL.FDR0.05.txt"), is_boxplot_file=FALSE)
convert_eqtl_genes(paste0(base_dir, tissue, "_eQTL.lead_snps.txt"), is_boxplot_file=FALSE)

# Boxplot coordinate matching file
convert_eqtl_genes(paste0(base_dir, tissue, "_eQTL.top_for_boxplot.txt"), is_boxplot_file=TRUE)

# Meta-analysis file (checks sub-directory first, falls back to base)
meta_path <- paste0(sub_dir, meta_name)
if (!file.exists(meta_path)) { meta_path <- paste0(base_dir, meta_name) }
convert_eqtl_genes(meta_path, is_boxplot_file=FALSE)
