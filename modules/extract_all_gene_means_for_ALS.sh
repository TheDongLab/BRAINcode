#!/bin/bash
#SBATCH --job-name=deseq2_and_means
#SBATCH --output=/home/zw529/donglab/data/target_ALS/deseq2_and_means_%j.out
#SBATCH --error=/home/zw529/donglab/data/target_ALS/deseq2_and_means_%j.err
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G

export METADATA="/home/zw529/donglab/data/target_ALS/targetALS_rnaseq_metadata.csv"
export EXPR_MATRIX="/home/zw529/donglab/data/target_ALS/QTL/expression_matrix.txt"

# Output File Paths
export MEANS_CSV="/home/zw529/donglab/data/target_ALS/all_genes_mean_tpm.csv"
export DESEQ2_CSV="/home/zw529/donglab/data/target_ALS/deseq2_als_vs_non_als_results.csv"

module load R

echo "Job started at: $(date)"
echo "Running DESeq2 & generating both separate output tables..."
echo "----------------------------------------"

Rscript -e '
library(DESeq2)
library(BiocParallel)

register(MulticoreParam(8))

meta_path   <- Sys.getenv("METADATA")
expr_path   <- Sys.getenv("EXPR_MATRIX")
means_path  <- Sys.getenv("MEANS_CSV")
deseq2_path <- Sys.getenv("DESEQ2_CSV")

# 1. Load Data
meta <- read.delim(meta_path, sep=",", header=TRUE, quote="\"", fill=TRUE, check.names=FALSE, stringsAsFactors=FALSE)
expr <- read.table(expr_path, header=TRUE, row.names=1, check.names=FALSE)

# 2. Sanitize and Align Sample IDs
clean_meta_ids <- gsub("[_-]", ".", gsub(" ", "", meta[["externalsampleid"]]))
clean_expr_ids <- gsub("[_-]", ".", gsub(" ", "", colnames(expr)))
colnames(expr) <- clean_expr_ids

# Filter & align metadata rows to match matrix columns exactly
meta <- meta[clean_meta_ids %in% clean_expr_ids, ]
clean_meta_ids <- gsub("[_-]", ".", gsub(" ", "", meta[["externalsampleid"]]))
expr_matched <- expr[, clean_meta_ids, drop=FALSE]

# 3. Categorize into ALS vs Non-ALS using Regex
meta$subject_group <- trimws(gsub("[\r\n\t]+", " ", meta$subject_group))
als_keywords <- "ALS|MND|Amyotrophic|Motor Neuron"
meta$als_status <- ifelse(grepl(als_keywords, meta$subject_group, ignore.case=TRUE), "ALS", "Non-ALS")

# Set Non-ALS as baseline
meta$als_status <- factor(meta$als_status, levels = c("Non-ALS", "ALS"))
rownames(meta)  <- clean_meta_ids

# Define logical indices for groups
als_samples     <- meta$als_status == "ALS"
non_als_samples <- meta$als_status == "Non-ALS"

# =========================================================================
# FILE 1: Save Means Table (Original Matrix Scale)
# =========================================================================
cat("Calculating mean expression values for ALS and Non-ALS groups...\n")
mean_als     <- rowMeans(expr_matched[, als_samples, drop=FALSE], na.rm=TRUE)
mean_non_als <- rowMeans(expr_matched[, non_als_samples, drop=FALSE], na.rm=TRUE)

final_means_list <- data.frame(
    gene_id          = rownames(expr_matched),
    mean_TPM_ALS     = mean_als,
    mean_TPM_Non_ALS = mean_non_als,
    stringsAsFactors = FALSE
)

write.csv(final_means_list, means_path, row.names = FALSE)
cat("[File 1 Saved]: Means list written to ->", means_path, "\n\n")


# =========================================================================
# FILE 2: Run DESeq2 (With Count Sanitization)
# =========================================================================
cat("Preparing expression matrix for DESeq2...\n")
raw_mat <- as.matrix(expr_matched)

# Check if data contains negative values (e.g., log-transformed or Z-scored)
if (min(raw_mat, na.rm = TRUE) < 0) {
    cat("Warning: Detected negative values in matrix. Converting from log-scale to positive counts for DESeq2...\n")
    # Shift minimum value to 0 or exponentiate if standard log2(x + 1)
    raw_mat <- raw_mat - min(raw_mat, na.rm = TRUE)
}

# Ensure positive non-zero integers required by DESeq2 model
expr_counts <- round(raw_mat)

cat("Setting up DESeqDataSet and running DESeq()...\n")
dds <- DESeqDataSetFromMatrix(
    countData = expr_counts,
    colData   = meta,
    design    = ~ als_status
)

# Filter low-count genes
keep <- rowSums(counts(dds)) >= 10
dds  <- dds[keep, ]

# Run pipeline
dds <- DESeq(dds, parallel = TRUE)

# Extract DESeq2 results table
res <- results(dds, contrast=c("als_status", "ALS", "Non-ALS"), parallel=TRUE)
res_df <- as.data.frame(res)
res_df$gene_id <- rownames(res_df)

# Reorder columns cleanly
col_order <- c("gene_id", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")
other_cols <- setdiff(colnames(res_df), col_order)
res_df <- res_df[, c(col_order, other_cols)]

# Sort by adjusted p-value
res_df <- res_df[order(res_df$padj, na.last = TRUE), ]

write.csv(res_df, deseq2_path, row.names = FALSE)
cat("[File 2 Saved]: Traditional DESeq2 results written to ->", deseq2_path, "\n\n")

# Preview top hits in output log
options(width = 150)
cat("Top 10 Differentially Expressed Genes (by padj):\n")
print(head(res_df, 10), row.names = FALSE)
'

echo "----------------------------------------"
echo "Job finished at: $(date)"
