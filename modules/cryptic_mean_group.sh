#!/bin/bash
#SBATCH --job-name=cryptic_mean_group
#SBATCH --output=/home/zw529/donglab/data/target_ALS/cryptic_mean_group_%j.out
#SBATCH --error=/home/zw529/donglab/data/target_ALS/cryptic_mean_group_%j.err
#SBATCH --time=00:20:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=28G

if [ -z "$1" ] || [ -z "$2" ]; then
    echo "Error: Missing arguments."
    echo "Usage: sbatch cryptic_mean_group.sh <CHROMOSOME> <COORD_PATTERN> [STRAND]"
    exit 1
fi

export METADATA="/home/zw529/donglab/data/target_ALS/targetALS_rnaseq_metadata.csv"
export SPLICE_MATRIX="/home/zw529/donglab/data/target_ALS/QTL/splicing_matrix.txt"
export CHR="$1"
export COORD_PAT="$2"
export STRAND="$3"

module load R

echo "Job started at: $(date)"
echo "Target: ${CHR} looking for coordinates matching: ${COORD_PAT}"
echo "----------------------------------------"

Rscript -e '
meta_path     <- Sys.getenv("METADATA")
splice_path   <- Sys.getenv("SPLICE_MATRIX")
target_chr    <- Sys.getenv("CHR")
coord_pat     <- Sys.getenv("COORD_PAT")
target_strand <- Sys.getenv("STRAND")

if (target_strand == "") { target_strand <- "[+-]" }

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
junc_counts <- as.numeric(expr[actual_row_name, ])
names(junc_counts) <- colnames(expr)

clean_meta_ids <- gsub("[_-]", ".", gsub(" ", "", meta[["externalsampleid"]]))
clean_expr_ids <- gsub("[_-]", ".", gsub(" ", "", names(junc_counts)))
names(junc_counts) <- clean_expr_ids

meta$splicing_value <- junc_counts[match(clean_meta_ids, names(junc_counts))]

valid_counts <- sum(!is.na(meta$splicing_value))
cat("Total metadata rows:", nrow(meta), "\n")
cat("Successfully matched samples:", valid_counts, "\n\n")

# --- CLEAN AND MERGE CATEGORIES ---
# 1. Strip raw trailing newlines and whitespace
meta$subject_group <- trimws(gsub("[\r\n\t]+", " ", meta$subject_group))

# 2. Merge Control duplicates
meta$subject_group[meta$subject_group == "Non Neurological Control"] <- "Non-Neurological Control"

# 3. Merge ALS Spectrum duplicates
meta$subject_group[meta$subject_group == "ALS Spectrum MND, Other Neurological Diseases"] <- "ALS Spectrum MND, Other Neurological Disorders"

# 4. Merge Other Neuro duplicates
meta$subject_group[meta$subject_group == "Other Neurological Disorders"] <- "Other Neurological Disorders"

# Compute the categorical means
results <- aggregate(splicing_value ~ subject_group, data=meta, FUN=mean, na.rm=TRUE)

# Force R to extend its print console width so columns stay side-by-side
options(width = 150)
print(results, row.names=FALSE, right=FALSE)
'

echo "----------------------------------------"
echo "Job finished at: $(date)"
