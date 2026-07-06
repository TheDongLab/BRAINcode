#!/bin/bash
#SBATCH --job-name=gene_mean_group
#SBATCH --output=/home/zw529/donglab/data/target_ALS/gene_mean_group_%j.out
#SBATCH --error=/home/zw529/donglab/data/target_ALS/gene_mean_group_%j.err
#SBATCH --time=00:20:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G

# Check if a Gene ID argument was passed
if [ -z "$1" ]; then
    echo "Error: No Gene ID provided."
    echo "Usage: sbatch gene_mean_group.sh <GENE_ID>"
    exit 1
fi

# Fixed dataset paths
export METADATA="/home/zw529/donglab/data/target_ALS/collections.postmortem_tissue_core.rnaseq_metadata.csv"
export EXPR_MATRIX="/home/zw529/donglab/data/target_ALS/QTL/expression_matrix.txt"
export GENE_ID="$1"

# Load R module
module load R

echo "Job started at: $(date)"
echo "Target Gene (Pattern Match): ${GENE_ID}"
echo "----------------------------------------"

# Run the R processing logic
Rscript -e '
meta_path <- Sys.getenv("METADATA")
expr_path <- Sys.getenv("EXPR_MATRIX")
gene_id   <- Sys.getenv("GENE_ID")

# Load datasets
meta <- read.csv(meta_path, check.names=FALSE, stringsAsFactors=FALSE)
expr <- read.table(expr_path, header=TRUE, row.names=1, check.names=FALSE)

# Find the row name that contains the gene string
matched_rows <- grep(gene_id, rownames(expr), value=TRUE)

if (length(matched_rows) == 0) {
    stop(paste("Error: No row matching pattern", gene_id, "found in expression matrix."))
} else if (length(matched_rows) > 1) {
    cat("Warning: Multiple matches found:", paste(matched_rows, collapse=", "), "\n")
    cat("Using the first match:", matched_rows[1], "\n\n")
}

actual_row_name <- matched_rows[1]
cat("Matched row name in matrix:", actual_row_name, "\n\n")

# Pull row and map to metadata via sample ID
gene_counts <- as.numeric(expr[actual_row_name, ])
names(gene_counts) <- colnames(expr)
meta$expression <- gene_counts[match(meta[["Externalsampleid"]], names(gene_counts))]

# Compute and output the categorical means
results <- aggregate(expression ~ `Subject Group`, data=meta, FUN=mean, na.rm=TRUE)
print(results, row.names=FALSE)
'

echo "----------------------------------------"
echo "Job finished at: $(date)"
