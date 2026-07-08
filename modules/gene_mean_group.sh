#!/bin/bash
#SBATCH --job-name=gene_mean_group
#SBATCH --output=/home/zw529/donglab/data/target_ALS/gene_mean_group_%j.out
#SBATCH --error=/home/zw529/donglab/data/target_ALS/gene_mean_group_%j.err
#SBATCH --time=00:20:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G

# Check if a Gene ID argument was passed
if [ -z "$1" ]; then
    echo "Error: No Gene ID provided."
    echo "Usage: sbatch gene_mean_group.sh <GENE_ID>"
    exit 1
fi

export METADATA="/home/zw529/donglab/data/target_ALS/targetALS_rnaseq_metadata.csv"
export EXPR_MATRIX="/home/zw529/donglab/data/target_ALS/QTL/expression_matrix.txt"
export GENE_ID="$1"

module load R

echo "Job started at: $(date)"
echo "Target Gene (Pattern Match): ${GENE_ID}"
echo "----------------------------------------"

Rscript -e '
meta_path <- Sys.getenv("METADATA")
expr_path <- Sys.getenv("EXPR_MATRIX")
gene_id   <- Sys.getenv("GENE_ID")

meta <- read.delim(meta_path, sep=",", header=TRUE, quote="\"", fill=TRUE, check.names=FALSE, stringsAsFactors=FALSE)
expr <- read.table(expr_path, header=TRUE, row.names=1, check.names=FALSE)

matched_rows <- grep(gene_id, rownames(expr), value=TRUE)
if (length(matched_rows) == 0) {
    stop(paste("Error: No row matching pattern", gene_id, "found in expression matrix."))
}
actual_row_name <- matched_rows[1]
cat("Matched row name in matrix:", actual_row_name, "\n")

# Extract expression data
gene_counts <- as.numeric(expr[actual_row_name, ])
names(gene_counts) <- colnames(expr)

# --- SANITIZE IDS: Convert hyphens and underscores to periods using lowercase header ---
clean_meta_ids <- gsub("[_-]", ".", gsub(" ", "", meta[["externalsampleid"]]))
clean_expr_ids <- gsub("[_-]", ".", gsub(" ", "", names(gene_counts)))
names(gene_counts) <- clean_expr_ids

# Map data using sanitized IDs
meta$expression <- gene_counts[match(clean_meta_ids, names(gene_counts))]

# QC Check
valid_counts <- sum(!is.na(meta$expression))
cat("Total metadata rows:", nrow(meta), "\n")
cat("Successfully matched matrix samples:", valid_counts, "\n\n")

if (valid_counts == 0) {
    cat("Sample matrix column example:", head(names(gene_counts), 3), "\n")
    cat("Metadata ID example:", head(clean_meta_ids, 3), "\n")
    stop("Error: 0 rows matched after converting hyphens/underscores.")
}

# Compute and output the categorical means using lowercase group header
results <- aggregate(expression ~ subject_group, data=meta, FUN=mean, na.rm=TRUE)
print(results, row.names=FALSE)
'

echo "----------------------------------------"
echo "Job finished at: $(date)"
