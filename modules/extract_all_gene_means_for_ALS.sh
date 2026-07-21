#!/bin/bash
#SBATCH --job-name=target_als_deseq2
#SBATCH --output=/home/zw529/donglab/data/target_ALS/target_als_deseq.out
#SBATCH --error=/home/zw529/donglab/data/target_ALS/target_als_deseq.err
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=70G

# TOP-LEVEL directory containing all tissue subfolders
export BASE_DIR="/home/zw529/donglab/data/target_ALS"
export METADATA="/home/zw529/donglab/data/target_ALS/targetALS_rnaseq_metadata.csv"

# Output Paths
export MEANS_CSV="/home/zw529/donglab/data/target_ALS/all_genes_mean_tpm_all_tissues.csv"
export DESEQ2_CSV="/home/zw529/donglab/data/target_ALS/deseq2_als_vs_non_als_all_tissues.csv"

module load R

echo "Job started at: $(date)"
echo "Searching across ALL tissue types under $BASE_DIR..."
echo "----------------------------------------"

Rscript -e '
library(DESeq2)
library(BiocParallel)

register(MulticoreParam(16))

base_dir    <- Sys.getenv("BASE_DIR")
meta_path   <- Sys.getenv("METADATA")
means_path    <- Sys.getenv("MEANS_CSV")
deseq2_path <- Sys.getenv("DESEQ2_CSV")

# 1. Search recursively across ALL tissue folders for normalization.tab files
tab_files <- list.files(base_dir, pattern = "normalization\\.tab$", recursive = TRUE, full.names = TRUE)
cat("Found a total of", length(tab_files), "normalization.tab files across all tissues.\n")

if (length(tab_files) == 0) stop("Error: No normalization.tab files found!")

# Extract sample IDs (e.g., CGND_HRA_00618)
raw_sample_ids <- basename(dirname(tab_files))
clean_expr_ids <- gsub("[_-]", ".", gsub(" ", "", raw_sample_ids))

# Extract Tissue Type from path (e.g., Target_ALS/Cerebellum/... -> Cerebellum)
extract_tissue <- function(path) {
    parts <- unlist(strsplit(path, "/"))
    idx <- which(tolower(parts) == "target_als")
    if (length(idx) > 0 && (idx[1] + 1) <= length(parts)) {
        return(parts[idx[1] + 1])
    }
    return("Unknown")
}

tissue_types <- sapply(tab_files, extract_tissue)

# Deduplicate directory entries if necessary
if (any(duplicated(clean_expr_ids))) {
    cat("Removing duplicate directory entries...\n")
    keep_idx       <- !duplicated(clean_expr_ids)
    tab_files      <- tab_files[keep_idx]
    clean_expr_ids <- clean_expr_ids[keep_idx]
    tissue_types   <- tissue_types[keep_idx]
}

# 2. Build Count and TPM Matrices
cat("Assembling matrix across all discovered tissues...\n")
first_df <- read.delim(tab_files[1], sep="\t", header=TRUE, stringsAsFactors=FALSE)
gene_ids <- first_df$gene_id

raw_counts_list <- list()
tpm_list        <- list()

for (i in seq_along(tab_files)) {
    s_id <- clean_expr_ids[i]
    df   <- read.delim(tab_files[i], sep="\t", header=TRUE, stringsAsFactors=FALSE)
    df   <- df[match(gene_ids, df$gene_id), ]
    
    raw_counts_list[[s_id]] <- df$raw_count
    tpm_list[[s_id]]        <- df$TPM
}

raw_matrix <- do.call(cbind, raw_counts_list)
tpm_matrix <- do.call(cbind, tpm_list)

rownames(raw_matrix) <- gene_ids
rownames(tpm_matrix) <- gene_ids

# Clean NAs
raw_matrix[is.na(raw_matrix)] <- 0
tpm_matrix[is.na(tpm_matrix)] <- 0

# 3. SAMPLE QC: Drop zero-read samples (e.g. GBB_12_17_CBLL)
sample_depths <- colSums(raw_matrix)
zero_samples  <- names(sample_depths[sample_depths == 0])

if (length(zero_samples) > 0) {
    cat("Removing", length(zero_samples), "zero-read sample(s):", paste(zero_samples, collapse=", "), "\n")
    valid_cols <- !(colnames(raw_matrix) %in% zero_samples)
    raw_matrix <- raw_matrix[, valid_cols, drop=FALSE]
    tpm_matrix <- tpm_matrix[, valid_cols, drop=FALSE]
}

cat("Final matrix dimensions:", nrow(raw_matrix), "genes across", ncol(raw_matrix), "samples.\n")

# 4. Process Metadata & Covariates
meta <- read.delim(meta_path, sep=",", header=TRUE, quote="\"", fill=TRUE, check.names=FALSE, stringsAsFactors=FALSE)
meta$clean_id <- gsub("[_-]", ".", gsub(" ", "", meta[["externalsampleid"]]))

# Keep metadata entries matching matrix
meta <- meta[meta$clean_id %in% colnames(raw_matrix), ]
if (any(duplicated(meta$clean_id))) {
    meta <- meta[!duplicated(meta$clean_id), ]
}

meta$subject_group <- trimws(gsub("[\r\n\t]+", " ", meta$subject_group))
als_keywords <- "ALS|MND|Amyotrophic|Motor Neuron"
status_vec   <- ifelse(grepl(als_keywords, meta$subject_group, ignore.case=TRUE), "ALS", "NonALS")
meta$als_status <- factor(status_vec, levels = c("NonALS", "ALS"))

# Attach Tissue Type mapping
sample_tissue_map <- setNames(tissue_types, clean_expr_ids)
meta$tissue_type  <- factor(make.names(sample_tissue_map[meta$clean_id]))

rownames(meta) <- meta$clean_id

# Align Matrix and Metadata
common_ids <- intersect(meta$clean_id, colnames(raw_matrix))
meta       <- meta[common_ids, ]
raw_matrix <- raw_matrix[, common_ids, drop=FALSE]
tpm_matrix <- tpm_matrix[, common_ids, drop=FALSE]

stopifnot(all(colnames(raw_matrix) == rownames(meta)))

# =========================================================================
# FILE 1: Group Mean TPM
# =========================================================================
cat("Calculating mean TPM expression values...\n")
als_samples     <- meta$als_status == "ALS"
non_als_samples <- meta$als_status == "NonALS"

mean_als     <- rowMeans(tpm_matrix[, als_samples, drop=FALSE], na.rm=TRUE)
mean_non_als <- rowMeans(tpm_matrix[, non_als_samples, drop=FALSE], na.rm=TRUE)

final_means_list <- data.frame(
    gene_id          = rownames(tpm_matrix),
    mean_TPM_ALS     = mean_als,
    mean_TPM_Non_ALS = mean_non_als,
    stringsAsFactors = FALSE
)

write.csv(final_means_list, means_path, row.names = FALSE)
cat("[File 1 Saved]: Multi-tissue TPM means written to ->", means_path, "\n\n")

# =========================================================================
# FILE 2: Multi-Tissue DESeq2 (Controlling for Tissue Type)
# =========================================================================
cat("Running DESeq2 controlling for tissue type (~ tissue_type + als_status)...\n")

count_data <- round(as.matrix(raw_matrix))
count_data <- count_data[rowSums(count_data) > 0, ]

# Multi-factor design controlling for tissue differences
dds <- DESeqDataSetFromMatrix(
    countData = count_data,
    colData   = meta,
    design    = ~ tissue_type + als_status
)

# Robust Geometric Mean Size Factors
cts <- counts(dds)
geo_means <- apply(cts, 1, function(row) {
    non_zero <- row[row > 0]
    if (length(non_zero) == 0) return(0)
    exp(mean(log(non_zero)))
})

# Estimate size factors manually and pass directly to DESeq() without sfType
dds <- estimateSizeFactors(dds, geoMeans = geo_means)
dds <- DESeq(dds, parallel = TRUE)

# Extract ALS vs NonALS results
res <- results(dds, contrast=c("als_status", "ALS", "NonALS"), parallel=TRUE)
res_df <- as.data.frame(res)
res_df$gene_id <- rownames(res_df)

col_order <- c("gene_id", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")
other_cols <- setdiff(colnames(res_df), col_order)
res_df <- res_df[, c(col_order, other_cols)]
res_df <- res_df[order(res_df$padj, na.last = TRUE), ]

write.csv(res_df, deseq2_path, row.names = FALSE)
cat("[File 2 Saved]: Multi-tissue DESeq2 results written to ->", deseq2_path, "\n\n")

cat("Top 10 Differentially Expressed Genes (across all tissues):\n")
print(head(res_df, 10), row.names = FALSE)
'

echo "----------------------------------------"
echo "Job finished at: $(date)"
