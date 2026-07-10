#!/bin/bash
#SBATCH --job-name=circ_mean_group
#SBATCH --output=/home/zw529/donglab/data/target_ALS/circ_mean_group_%j.out
#SBATCH --error=/home/zw529/donglab/data/target_ALS/circ_mean_group_%j.err
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=28G

if [ -z "$1" ] || [ -z "$2" ]; then
    echo "Error: Missing arguments."
    echo "Usage: sbatch circRNA_mean_group.sh <CHROMOSOME> <COORD_PATTERN> [STRAND]"
    exit 1
fi

export METADATA="/home/zw529/donglab/data/target_ALS/targetALS_rnaseq_metadata.csv"
export CIRC_MATRIX="/home/zw529/donglab/data/target_ALS/QTL/circ_matrix.txt"
export TDP43_METADATA="/home/zw529/donglab/data/target_ALS/aligned_rnaseq_tdp43_by_subject.csv"
export CHR="$1"
export COORD_PAT="$2"
export STRAND="$3"

module load R

echo "Job started at: $(date)"
echo "Target: ${CHR} scanning for any matching pattern containing: ${COORD_PAT}"
echo "----------------------------------------"

Rscript -e '
meta_path     <- Sys.getenv("METADATA")
circ_path     <- Sys.getenv("CIRC_MATRIX")
tdp43_path    <- Sys.getenv("TDP43_METADATA")
target_chr    <- Sys.getenv("CHR")
coord_pat     <- Sys.getenv("COORD_PAT")
target_strand <- Sys.getenv("STRAND")

if (target_strand == "") { target_strand <- "[+-]" }

# --- FIX: Rebuilt regex for the chr:start-end:strand format ---
safe_strand <- gsub("\\+", "\\\\+", target_strand)
regex_pattern <- paste0("^", target_chr, ":.*", coord_pat, ".*:", safe_strand, "$")

meta <- read.csv(meta_path, check.names=FALSE, stringsAsFactors=FALSE)
expr <- read.table(circ_path, header=TRUE, row.names=1, check.names=FALSE)
tdp43 <- read.csv(tdp43_path, check.names=FALSE, stringsAsFactors=FALSE)

matched_rows <- grep(regex_pattern, rownames(expr), value=TRUE)
if (length(matched_rows) == 0) {
    stop(paste("Error: No circRNAs matching regex", regex_pattern, "found in matrix."))
}

cat("Found", length(matched_rows), "matching circRNA(s) total.\n")
cat("------------------------------------------------------------\n\n")

# --- CLEAN AND MERGE CATEGORIES ---
meta$subject_group <- trimws(gsub("[\r\n\t]+", " ", meta$subject_group))
meta$subject_group[meta$subject_group == "Non Neurological Control"] <- "Non-Neurological Control"
meta$subject_group[meta$subject_group == "ALS Spectrum MND, Other Neurological Diseases"] <- "ALS Spectrum MND, Other Neurological Disorders"
meta$subject_group[meta$subject_group == "Other Neurological Disorders"] <- "Other Neurological Disorders"

# Sanitize metadata sample IDs once
clean_meta_ids <- gsub("[_-]", ".", gsub(" ", "", meta[["externalsampleid"]]))

# Loop through every single matched circRNA and calculate values independently
for (i in 1:length(matched_rows)) {
    circ_name <- matched_rows[i]
    
    cat("============================================================\n")
    cat(" circRNA ", i, " OF ", length(matched_rows), ": ", circ_name, "\n", sep="")
    cat("============================================================\n")
    
    # Extract values for this specific row
    circ_counts <- as.numeric(expr[circ_name, ])
    names(circ_counts) <- colnames(expr)
    
    # Sanitize matrix expression IDs
    clean_expr_ids <- gsub("[_-]", ".", gsub(" ", "", names(circ_counts)))
    names(circ_counts) = clean_expr_ids
    
    # Map data dynamically to a temporary loop column to avoid cross-contamination
    meta$tmp_circ_value <- circ_counts[match(clean_meta_ids, names(circ_counts))]
    
    valid_counts <- sum(!is.na(meta$tmp_circ_value))
    cat("Successfully matched samples:", valid_counts, "out of", nrow(meta), "\n\n")
    
    # Compute the categorical means
    results <- aggregate(tmp_circ_value ~ subject_group, data=meta, FUN=mean, na.rm=TRUE)
    colnames(results)[2] <- "mean_circ_percentage"
    
    # Print the table side-by-side cleanly
    options(width = 150)
    print(results, row.names=FALSE, right=FALSE)
    cat("\n\n")
}

# ==============================================================================
# --- TDP43 SPECIFIC ANALYSIS ---
# ==============================================================================
cat("\n\n")
cat("============================================================\n")
cat(" RUNNING INTERSECTION WITH TDP43 COHORT                     \n")
cat("============================================================\n\n")

# Clean up whitespace and ensure missing scores are marked explicitly
tdp43$Neuronal_TDP43_Score <- trimws(gsub("[\r\n\t]+", " ", tdp43$Neuronal_TDP43_Score))
tdp43$Neuronal_TDP43_Score[tdp43$Neuronal_TDP43_Score == ""] <- "Unknown/Missing"

# --- ENFORCE EXPLICIT ORDERING ---
custom_levels <- c(
    "Frequent", 
    "Moderate", 
    "Sparse", 
    "Sparse Residual MN (with pTDP-43 Cytoplasmic Inclusions)", 
    "Absent"
)

# Catch any unexpected values using safe indices (avoids Bash symbol syntax conflicts)
extra_idx <- is.na(match(tdp43$Neuronal_TDP43_Score, custom_levels))
extra_levels <- unique(tdp43$Neuronal_TDP43_Score[extra_idx])
all_levels <- c(custom_levels, extra_levels)

tdp43$Neuronal_TDP43_Score <- factor(tdp43$Neuronal_TDP43_Score, levels = all_levels)

# Sanitize TDP43 file sample IDs (mapping via RNAseq_Sample_ID)
clean_tdp43_ids <- gsub("[_-]", ".", gsub(" ", "", tdp43[["RNAseq_Sample_ID"]]))

# Re-sanitize the global matrix expression IDs just to ensure clean scope
clean_expr_ids_tdp43 <- gsub("[_-]", ".", gsub(" ", "", colnames(expr)))

for (i in 1:length(matched_rows)) {
    circ_name <- matched_rows[i]
    
    cat("============================================================\n")
    cat(" TDP43 COHORT - circRNA ", i, " OF ", length(matched_rows), ": ", circ_name, "\n", sep="")
    cat("============================================================\n")
    
    # Extract values for this specific row
    circ_counts <- as.numeric(expr[circ_name, ])
    names(circ_counts) <- clean_expr_ids_tdp43
    
    # Map splicing values directly to the TDP43 structure using sanitized IDs
    tdp43$tmp_circ_value <- circ_counts[match(clean_tdp43_ids, names(circ_counts))]
    
    valid_counts_tdp43 <- sum(!is.na(tdp43$tmp_circ_value))
    cat("Successfully matched TDP43 cohort samples:", valid_counts_tdp43, "out of", nrow(tdp43), "\n\n")
    
    if (valid_counts_tdp43 > 0) {
        # Compute the means grouped by Neuronal_TDP43_Score
        tdp43_results <- aggregate(tmp_circ_value ~ Neuronal_TDP43_Score, data=tdp43, FUN=mean, na.rm=TRUE)
        colnames(tdp43_results)[2] <- "mean_circ_percentage"
        
        # Print the side-by-side TDP43 table cleanly
        options(width = 150)
        print(tdp43_results, row.names=FALSE, right=FALSE)
    } else {
        cat("Warning: No matching sample expression data found for the TDP43 cohort.\n")
    }
    cat("\n\n")
}
'

echo "----------------------------------------"
echo "Job finished at: $(date)"
