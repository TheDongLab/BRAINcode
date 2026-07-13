#!/bin/bash
#SBATCH --job-name=all_gene_means
#SBATCH --output=/home/zw529/donglab/data/target_ALS/all_gene_means_%j.out
#SBATCH --error=/home/zw529/donglab/data/target_ALS/all_gene_means_%j.err
#SBATCH --time=00:20:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G

export METADATA="/home/zw529/donglab/data/target_ALS/targetALS_rnaseq_metadata.csv"
export EXPR_MATRIX="/home/zw529/donglab/data/target_ALS/QTL/expression_matrix.txt"
export OUTPUT_FILE="/home/zw529/donglab/data/target_ALS/all_genes_mean_tpm.csv"

module load R

echo "Job started at: $(date)"
echo "Calculating group means for all genes (Handling Complex Subgroups)..."
echo "----------------------------------------"

Rscript -e '
meta_path <- Sys.getenv("METADATA")
expr_path <- Sys.getenv("EXPR_MATRIX")
out_path  <- Sys.getenv("OUTPUT_FILE")

# 1. Load Data
meta <- read.delim(meta_path, sep=",", header=TRUE, quote="\"", fill=TRUE, check.names=FALSE, stringsAsFactors=FALSE)
expr <- read.table(expr_path, header=TRUE, row.names=1, check.names=FALSE)

# 2. Sanitize and Map IDs (Preserving exact previous script logic)
clean_meta_ids <- gsub("[_-]", ".", gsub(" ", "", meta[["externalsampleid"]]))
clean_expr_ids <- gsub("[_-]", ".", gsub(" ", "", colnames(expr)))
colnames(expr) <- clean_expr_ids

# Filter metadata to match available matrix columns
meta <- meta[clean_meta_ids %in% clean_expr_ids, ]
clean_meta_ids <- gsub("[_-]", ".", gsub(" ", "", meta[["externalsampleid"]]))
expr_matched <- expr[, clean_meta_ids, drop=FALSE]

# 3. Clean and Standardize Categories into ALS vs Non-ALS using Regex
meta$subject_group <- trimws(gsub("[\r\n\t]+", " ", meta$subject_group))

# This regex pattern dynamically flags any group containing ALS, MND, or Amyotrophic variations
als_keywords <- "ALS|MND|Amyotrophic|Motor Neuron"
meta$als_status <- ifelse(grepl(als_keywords, meta$subject_group, ignore.case=TRUE), "ALS", "Non-ALS")

# --- QC Verification Printouts ---
options(width = 150)
cat("\n=== GROUP MAPPING VERIFICATION ===\n")
print(unique(meta[, c("subject_group", "als_status")]), row.names=FALSE)
cat("==================================\n\n")

cat("Sample counts per group:\n")
print(table(meta$als_status))
cat("\n")

# 4. Separate columns and calculate Means
als_samples     <- meta$als_status == "ALS"
non_als_samples <- meta$als_status == "Non-ALS"

cat("Computing matrix-wide row means...\n")
mean_als     <- rowMeans(expr_matched[, als_samples, drop=FALSE], na.rm=TRUE)
mean_non_als <- rowMeans(expr_matched[, non_als_samples, drop=FALSE], na.rm=TRUE)

# 5. Compile and Output Results
final_list <- data.frame(
    gene_id          = rownames(expr_matched),
    mean_TPM_ALS     = mean_als,
    mean_TPM_Non_ALS = mean_non_als,
    stringsAsFactors = FALSE
)

# Save to CSV
write.csv(final_list, out_path, row.names = FALSE)
cat("Full list successfully saved to:", out_path, "\n\n")

# Preview output in the log file
cat("First 15 rows of the generated list:\n")
print(head(final_list, 15), row.names=FALSE)
'

echo "----------------------------------------"
echo "Job finished at: $(date)"
