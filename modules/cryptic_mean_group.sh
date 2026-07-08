#!/bin/bash
#SBATCH --job-name=cryptic_mean_group
#SBATCH --output=/home/zw529/donglab/data/target_ALS/cryptic_mean_group_%j.out
#SBATCH --error=/home/zw529/donglab/data/target_ALS/cryptic_mean_group_%j.err
#SBATCH --time=00:20:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G

if [ -z "$1" ]; then
    echo "Error: No Coordinate Pattern provided (e.g., chr19:17642)."
    echo "Usage: sbatch get_splicing_means.sh <COORD_PATTERN>"
    exit 1
fi

# Paths updated for the splicing matrix
export METADATA="/home/zw529/donglab/data/target_ALS/collections.postmortem_tissue_core.rnaseq_metadata.csv"
export SPLICE_MATRIX="/home/zw529/donglab/data/target_ALS/QTL/splicing_matrix.txt"
export COORD_PATTERN="$1"

module load R

echo "Job started at: $(date)"
echo "Target Junction Pattern: ${COORD_PATTERN}"
echo "----------------------------------------"

Rscript -e '
meta_path   <- Sys.getenv("METADATA")
splice_path <- Sys.getenv("SPLICE_MATRIX")
coord_pat   <- Sys.getenv("COORD_PATTERN")

meta <- read.csv(meta_path, check.names=FALSE, stringsAsFactors=FALSE)
# Using read.table assuming splicing_matrix.txt is tab-separated with junction_id as column 1
expr <- read.table(splice_path, header=TRUE, row.names=1, check.names=FALSE)

matched_rows <- grep(coord_pat, rownames(expr), value=TRUE)
if (length(matched_rows) == 0) {
    stop(paste("Error: No junctions matching pattern", coord_pat, "found in matrix."))
}

# Print out matches so you can see all junctions passing through the cryptic region
cat("Found", length(matched_rows), "matching junction(s):\n")
print(matched_rows)
cat("\nUsing the primary match:", matched_rows[1], "\n\n")

actual_row_name <- matched_rows[1]

# Extract values
junc_counts <- as.numeric(expr[actual_row_name, ])
names(junc_counts) <- colnames(expr)

# Sanitize sample IDs (hyphens, underscores, or periods converted to standard dots)
clean_meta_ids <- gsub("[_-]", ".", gsub(" ", "", meta[["Externalsampleid"]]))
clean_expr_ids <- gsub("[_-]", ".", gsub(" ", "", names(junc_counts)))
names(junc_counts) <- clean_expr_ids

# Map data 
meta$splicing_value <- junc_counts[match(clean_meta_ids, names(junc_counts))]

# QC Check
valid_counts <- sum(!is.na(meta$splicing_value))
cat("Total metadata rows:", nrow(meta), "\n")
cat("Successfully matched samples:", valid_counts, "\n\n")

if (valid_counts == 0) {
    stop("Error: 0 samples matched between metadata and splicing matrix.")
}

# Compute and output the categorical means
results <- aggregate(splicing_value ~ `Subject Group`, data=meta, FUN=mean, na.rm=TRUE)
print(results, row.names=FALSE)
'

echo "----------------------------------------"
echo "Job finished at: $(date)"
