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
    exit 1
fi

export METADATA="/home/zw529/donglab/data/target_ALS/targetALS_rnaseq_metadata.csv"
export SPLICE_MATRIX="/home/zw529/donglab/data/target_ALS/QTL/splicing_matrix.txt"
export CHR="$1"
export COORD_PAT="$2"
export STRAND="$3"

module load R

echo "Job started at: $(date)"
echo "Target: ${CHR} scanning for any matching pattern containing: ${COORD_PAT}"
echo "----------------------------------------"

Rscript -e '
meta_path     <- Sys.getenv("METADATA")
splice_path   <- Sys.getenv("SPLICE_MATRIX")
target_chr    <- Sys.getenv("CHR")
coord_pat     <- Sys.getenv("COORD_PAT")
target_strand <- Sys.getenv("STRAND")

if (target_strand == "") { target_strand <- "[+-]" }

# Build pattern to scan for the coordinate pattern anywhere inside the junction string
safe_strand <- gsub("\\+", "\\\\+", target_strand)
regex_pattern <- paste0("^", target_chr, ":", safe_strand, ":.*", coord_pat)

meta <- read.csv(meta_path, check.names=FALSE, stringsAsFactors=FALSE)
expr <- read.table(splice_path, header=TRUE, row.names=1, check.names=FALSE)

matched_rows <- grep(regex_pattern, rownames(expr), value=TRUE)
if (length(matched_rows) == 0) {
    stop(paste("Error: No junctions matching regex", regex_pattern, "found in matrix."))
}

cat("Found", length(matched_rows), "matching junction(s) total.\n")
cat("------------------------------------------------------------\n\n")

# --- CLEAN AND MERGE CATEGORIES ---
meta$subject_group <- trimws(gsub("[\r\n\t]+", " ", meta$subject_group))
meta$subject_group[meta$subject_group == "Non Neurological Control"] <- "Non-Neurological Control"
meta$subject_group[meta$subject_group == "ALS Spectrum MND, Other Neurological Diseases"] <- "ALS Spectrum MND, Other Neurological Disorders"
meta$subject_group[meta$subject_group == "Other Neurological Disorders"] <- "Other Neurological Disorders"

# Sanitize metadata sample IDs once
clean_meta_ids <- gsub("[_-]", ".", gsub(" ", "", meta[["externalsampleid"]]))

# Loop through every single matched junction and calculate values independently
for (i in 1:length(matched_rows)) {
    junc_name <- matched_rows[i]
    
    cat("============================================================\n")
    cat(" JUNCTION ", i, " OF ", length(matched_rows), ": ", junc_name, "\n", sep="")
    cat("============================================================\n")
    
    # Extract values for this specific row
    junc_counts <- as.numeric(expr[junc_name, ])
    names(junc_counts) <- colnames(expr)
    
    # Sanitize matrix expression IDs
    clean_expr_ids <- gsub("[_-]", ".", gsub(" ", "", names(junc_counts)))
    names(junc_counts) <- clean_expr_ids
    
    # Map data dynamically to a temporary loop column to avoid cross-contamination
    meta$tmp_splicing_value <- junc_counts[match(clean_meta_ids, names(junc_counts))]
    
    valid_counts <- sum(!is.na(meta$tmp_splicing_value))
    cat("Successfully matched samples:", valid_counts, "out of", nrow(meta), "\n\n")
    
    # Compute the categorical means
    results <- aggregate(tmp_splicing_value ~ subject_group, data=meta, FUN=mean, na.rm=TRUE)
    colnames(results)[2] <- "mean_splicing_value"
    
    # Print the table side-by-side cleanly
    options(width = 150)
    print(results, row.names=FALSE, right=FALSE)
    cat("\n\n")
}
'

echo "----------------------------------------"
echo "Job finished at: $(date)"
