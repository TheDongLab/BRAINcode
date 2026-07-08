#!/bin/bash
#SBATCH --job-name=cryptic_mean_group
#SBATCH --output=/home/zw529/donglab/data/target_ALS/cryptic_mean_group_%j.out
#SBATCH --error=/home/zw529/donglab/data/target_ALS/cryptic_mean_group_%j.err
#SBATCH --time=00:20:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=28G

if [ -z "$1" ] || [ -z "$2" ]; then
    echo "Error: Missing arguments."
    echo "Usage: sbatch cryptic_mean_group.sh <CHROMOSOME> <COORD_PATTERN> [STRAND]"
    echo "Examples:"
    echo "  sbatch cryptic_mean_group.sh chr8 79611214               (Prefix pattern, any strand)"
    echo "  sbatch cryptic_mean_group.sh chr8 79611214-79616822 '+'  (Exact match, plus strand)"
    exit 1
fi

# Updated metadata file path here
export METADATA="/home/zw529/donglab/data/target_ALS/targetALS_rnaseq_metadata.csv"
export SPLICE_MATRIX="/home/zw529/donglab/data/target_ALS/QTL/splicing_matrix.txt"
export CHR="$1"
export COORD_PAT="$2"
export STRAND="$3"

module load R

echo "Job started at: $(date)"
echo "Target: ${CHR} looking for coordinates matching: ${COORD_PAT}"
if [ ! -z "${STRAND}" ]; then echo "Strand Filter: ${STRAND}"; fi
echo "----------------------------------------"

Rscript -e '
meta_path     <- Sys.getenv("METADATA")
splice_path   <- Sys.getenv("SPLICE_MATRIX")
target_chr    <- Sys.getenv("CHR")
coord_pat     <- Sys.getenv("COORD_PAT")
target_strand <- Sys.getenv("STRAND")

if (target_strand == "") {
    target_strand <- "[+-]"
}

# Dynamic exact vs. prefix mapping with literal string escape for plus signs
if (grepl("-", coord_pat)) {
    safe_strand <- gsub("\\+", "\\\\+", target_strand)
    regex_pattern <- paste0("^", target_chr, ":", safe_strand, ":", coord_pat, "$")
} else {
    safe_strand <- gsub("\\+", "\\\\+", target_strand)
    regex_pattern <- paste0("^", target_chr, ":", safe_strand, ":", coord_pat)
}

meta <- read.csv(meta_path, check.names=FALSE, stringsAsFactors=FALSE)
expr <- read.table(splice_path, header=TRUE, row.names=1, check.names=FALSE)

matched_rows <- grep(regex_pattern, rownames(expr), value=TRUE)
if (length(matched_rows) == 0) {
    stop(paste("Error: No junctions matching regex", regex_pattern, "found in matrix."))
}

cat("Found", length(matched_rows), "matching junction(s):\n")
print(matched_rows)
cat("\nUsing the primary match:", matched_rows[1], "\n\n")

actual_row_name <- matched_rows[1]

# Extract values
junc_counts <- as.numeric(expr[actual_row_name, ])
names(junc_counts) <- colnames(expr)

# Sanitize sample IDs
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
