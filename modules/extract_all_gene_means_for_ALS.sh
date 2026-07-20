#!/bin/bash
#SBATCH --job-name=target_als_deseq2
#SBATCH --output=/home/zw529/donglab/data/target_ALS/target_als_deseq.out
#SBATCH --error=/home/zw529/donglab/data/target_ALS/target_als_deseq.err
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=80G

export PROCESSED_DIR="/home/zw529/donglab/data/target_ALS/Cerebellum/RNAseq/Processed"
export METADATA="/home/zw529/donglab/data/target_ALS/targetALS_rnaseq_metadata.csv"

# Output Paths
export MEANS_CSV="/home/zw529/donglab/data/target_ALS/all_genes_mean_tpm.csv"
export DESEQ2_CSV="/home/zw529/donglab/data/target_ALS/deseq2_als_vs_non_als_results.csv"

module load R

echo "Job started at: $(date)"
echo "Building matrices directly from normalization.tab files..."
echo "----------------------------------------"

Rscript -e '
library(DESeq2)
library(BiocParallel)

register(MulticoreParam(8))

processed_dir <- Sys.getenv("PROCESSED_DIR")
meta_path     <- Sys.getenv("METADATA")
means_path    <- Sys.getenv("MEANS_CSV")
deseq2_path   <- Sys.getenv("DESEQ2_CSV")

# 1. Find all normalization.tab files dynamically
tab_files <- list.files(processed_dir, pattern = "normalization\\.tab$", recursive = TRUE, full.names = TRUE)
cat("Found", length(tab_files), "sample normalization.tab files.\n")

if (length(tab_files) == 0) {
    stop("Error: No normalization.tab files found under ", processed_dir)
}

# Extract sample IDs from parent directory names (e.g., CGND_HRA_00618)
sample_ids <- basename(dirname(tab_files))

# 2. Read all files into memory and assemble raw_count and TPM matrices
cat("Assembling Raw Count and TPM matrices...\n")
first_df <- read.delim(tab_files[1], sep="\t", header=TRUE, stringsAsFactors=FALSE)
gene_ids <- first_df$gene_id

raw_counts_list <- list()
tpm_list        <- list()

for (i in seq_along(tab_files)) {
    s_id <- sample_ids[i]
    df   <- read.delim(tab_files[i], sep="\t", header=TRUE, stringsAsFactors=FALSE)
    
    # Ensure genes are in the exact same order
    df <- df[match(gene_ids, df$gene_id), ]
    
    raw_counts_list[[s_id]] <- df$raw_count
    tpm_list[[s_id]]        <- df$TPM
}

raw_matrix <- do.call(cbind, raw_counts_list)
tpm_matrix <- do.call(cbind, tpm_list)

rownames(raw_matrix) <- gene_ids
rownames(tpm_matrix) <- gene_ids

# 3. Load & Sanitize Metadata
meta <- read.delim(meta_path, sep=",", header=TRUE, quote="\"", fill=TRUE, check.names=FALSE, stringsAsFactors=FALSE)

clean_meta_ids <- gsub("[_-]", ".", gsub(" ", "", meta[["externalsampleid"]]))
clean_expr_ids <- gsub("[_-]", ".", gsub(" ", "", colnames(raw_matrix)))

colnames(raw_matrix) <- clean_expr_ids
colnames(tpm_matrix) <- clean_expr_ids

# Match metadata rows to matrix columns
meta <- meta[clean_meta_ids %in% clean_expr_ids, ]
clean_meta_ids <- gsub("[_-]", ".", gsub(" ", "", meta[["externalsampleid"]]))

raw_matrix <- raw_matrix[, clean_meta_ids, drop=FALSE]
tpm_matrix <- tpm_matrix[, clean_meta_ids, drop=FALSE]

# Categorize ALS vs Non-ALS
meta$subject_group <- trimws(gsub("[\r\n\t]+", " ", meta$subject_group))
als_keywords <- "ALS|MND|Amyotrophic|Motor Neuron"
meta$als_status <- ifelse(grepl(als_keywords, meta$subject_group, ignore.case=TRUE), "ALS", "Non-ALS")
meta$als_status <- factor(meta$als_status, levels = c("Non-ALS", "ALS"))
rownames(meta)  <- clean_meta_ids

# Define sample logical vectors
als_samples     <- meta$als_status == "ALS"
non_als_samples <- meta$als_status == "Non-ALS"

# =========================================================================
# FILE 1: Save Group TPM Means Output
# =========================================================================
cat("Calculating mean TPM expression values...\n")
mean_als     <- rowMeans(tpm_matrix[, als_samples, drop=FALSE], na.rm=TRUE)
mean_non_als <- rowMeans(tpm_matrix[, non_als_samples, drop=FALSE], na.rm=TRUE)

final_means_list <- data.frame(
    gene_id          = rownames(tpm_matrix),
    mean_TPM_ALS     = mean_als,
    mean_TPM_Non_ALS = mean_non_als,
    stringsAsFactors = FALSE
)

write.csv(final_means_list, means_path, row.names = FALSE)
cat("[File 1 Saved]: Means list written to ->", means_path, "\n\n")

# =========================================================================
# FILE 2: Run True DESeq2 using Raw Count Matrix
# =========================================================================
cat("Running standard DESeq2 analysis on true raw counts...\n")

# Ensure count matrix is integer
count_data <- round(as.matrix(raw_matrix))

dds <- DESeqDataSetFromMatrix(
    countData = count_data,
    colData   = meta,
    design    = ~ als_status
)

# Filter low-expression genes across samples
keep <- rowSums(counts(dds)) >= 10
dds  <- dds[keep, ]

# Execute DESeq2 pipeline
dds <- DESeq(dds, parallel = TRUE)

# Extract standard results table (ALS vs Non-ALS)
res <- results(dds, contrast=c("als_status", "ALS", "Non-ALS"), parallel=TRUE)
res_df <- as.data.frame(res)
res_df$gene_id <- rownames(res_df)

# Reorder columns: gene_id first
col_order <- c("gene_id", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")
other_cols <- setdiff(colnames(res_df), col_order)
res_df <- res_df[, c(col_order, other_cols)]

# Sort by adjusted p-value
res_df <- res_df[order(res_df$padj, na.last = TRUE), ]

write.csv(res_df, deseq2_path, row.names = FALSE)
cat("[File 2 Saved]: Traditional DESeq2 results written to ->", deseq2_path, "\n\n")

# Preview top hits
options(width = 150)
cat("Top 10 Differentially Expressed Genes (by padj):\n")
print(head(res_df, 10), row.names = FALSE)
'

echo "----------------------------------------"
echo "Job finished at: $(date)"
