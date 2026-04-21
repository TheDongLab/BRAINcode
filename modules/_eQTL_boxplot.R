#!/usr/bin/env Rscript
###########################################
# _eQTL_boxplot.R - FIXED v2 (Robust Covariate Handling)
###########################################
 
suppressPackageStartupMessages({
  library(data.table)
})
 
args <- commandArgs(TRUE)
if (length(args) < 7) {
  stop("Usage: Rscript _eQTL_boxplot.R <pairs_file> <snp_file> <expr_file> <cov_file> <snp_loc> <out_dir> <tissue>")
}
 
GSfile    <- args[1]
snp_file  <- args[2]
expr_file <- args[3]
cov_file  <- args[4]
snp_loc   <- args[5]
out_dir   <- args[6]
tissue    <- args[7]
 
out_file_std <- file.path(out_dir, paste0(tissue, "_all_sig_eQTL_boxplots.pdf"))
out_file_col <- file.path(out_dir, paste0(tissue, "_all_sig_eQTL_boxplots_colored.pdf"))
 
message("# Loading genomic matrices...")
pairs    <- fread(GSfile, header=FALSE, col.names=c("geneid","snpid"))
snp_mat  <- fread(snp_file, header=TRUE)
expr_mat <- fread(expr_file, header=TRUE)
cov_mat  <- fread(cov_file, header=TRUE)
 
message(sprintf("# Covariate matrix: %d rows x %d columns", nrow(cov_mat), ncol(cov_mat)))
message(sprintf("# Covariate first column name: '%s'", names(cov_mat)[1]))
 
# --- SAMPLE ALIGNMENT ---
# Extract sample IDs from matrix column names (skip first column = row ID)
snp_samples  <- names(snp_mat)[-1]
expr_samples <- names(expr_mat)[-1]
cov_samples  <- names(cov_mat)[-1]
 
message(sprintf("# SNP matrix has %d samples", length(snp_samples)))
message(sprintf("# Expression matrix has %d samples", length(expr_samples)))
message(sprintf("# Covariate matrix has %d samples", length(cov_samples)))
 
# Find samples present in ALL three matrices (tissue-specific filtering)
common_samples <- Reduce(intersect, list(snp_samples, expr_samples, cov_samples))
message(sprintf("# Found %d samples common to all matrices (tissue-specific)", length(common_samples)))
 
if (length(common_samples) < 5) {
  stop("ERROR: Fewer than 5 common samples. Check data alignment.")
}
 
# Set keys for fast lookups
setkey(snp_mat, snpid)
setkey(expr_mat, geneid)
 
# --- COVARIATE / STATUS EXTRACTION ---
# Get the first column name (the row identifier)
cov_first_col <- names(cov_mat)[1]
message(sprintf("# Looking for 'is_als' row in column '%s'...", cov_first_col))
 
# Try to find the 'is_als' row using get() function
als_row <- NULL
tryCatch({
  als_row <- cov_mat[get(cov_first_col) == "is_als", ]
  if (nrow(als_row) > 0) {
    message(sprintf("# Found 'is_als' row"))
  }
}, error = function(e) {
  message(sprintf("# Error searching for 'is_als': %s", e$message))
})
 
# Fallback: try direct indexing if get() didn't work
if (is.null(als_row) || nrow(als_row) == 0) {
  message("# Attempting fallback method to find 'is_als'...")
  # Try accessing first column directly
  first_col_data <- cov_mat[[cov_first_col]]
  als_idx <- which(first_col_data == "is_als")
  
  if (length(als_idx) > 0) {
    als_row <- cov_mat[als_idx[1], ]
    message(sprintf("# Found 'is_als' at row %d using fallback method", als_idx[1]))
  }
}
 
# Final fallback: use all purple if not found
if (is.null(als_row) || nrow(als_row) == 0) {
  message("# WARNING: 'is_als' row not found in covariates. Using all purple points.")
  als_vals <- rep(0, length(common_samples))
} else {
  # Extract values in the order of common_samples
  tryCatch({
    als_vals <- as.numeric(als_row[, ..common_samples, with=FALSE])
    message(sprintf("# Successfully extracted ALS status for %d samples", length(als_vals)))
  }, error = function(e) {
    message(sprintf("# Error extracting ALS values: %s. Using all purple.", e$message))
    als_vals <<- rep(0, length(common_samples))
  })
}
 
p_colors <- ifelse(als_vals > 0.5, "red", "#9932CC")
p_shapes <- ifelse(als_vals > 0.5, 16, 18)
p_sizes  <- ifelse(als_vals > 0.5, 0.7, 1.1)
 
run_plotting <- function(pdf_path, use_status_colors = FALSE) {
    pdf(pdf_path, width=12, height=5)
    
    plots_made <- 0
    pairs_skipped <- 0
    
    for (i in seq_len(nrow(pairs))) {
        G <- pairs$geneid[i]
        S <- pairs$snpid[i]
        
        # Retrieve values by column name, ensuring order matches common_samples EXACTLY
        snp_row <- snp_mat[snpid == S]
        expr_row <- expr_mat[geneid == G]
        
        if (nrow(snp_row) == 0 || nrow(expr_row) == 0) {
            pairs_skipped <- pairs_skipped + 1
            next
        }
        
        # Extract values in the SAME order as common_samples
        snp_vals  <- as.numeric(snp_row[, ..common_samples, with=FALSE])
        expr_vals <- as.numeric(expr_row[, ..common_samples, with=FALSE])
        
        # Sanity check: all values should be numbers
        if (all(is.na(snp_vals)) || all(is.na(expr_vals))) {
            pairs_skipped <- pairs_skipped + 1
            next
        }
        
        df <- data.frame(
            expression = expr_vals,
            SNP = snp_vals,
            p_col = if(use_status_colors) p_colors else rep("darkred", length(common_samples)),
            p_pch = if(use_status_colors) p_shapes else rep(16, length(common_samples)),
            p_cex = if(use_status_colors) p_sizes else rep(0.7, length(common_samples)),
            stringsAsFactors = FALSE
        )
        
        # Remove rows with missing data
        df <- df[!is.na(df$SNP) & !is.na(df$expression), ]
        if (nrow(df) < 5) {
            pairs_skipped <- pairs_skipped + 1
            next
        }
 
        get_p <- function(formula, data) {
            tryCatch({
                fit <- lm(formula, data=data)
                formatC(summary(fit)$coefficients[2,4], format="e", digits=3)
            }, error = function(e) "NA")
        }
 
        par(mfrow=c(1,3), mar=c(5,4,4,2), oma=c(0,0,3,0))
        
        # Model 1: Additive
        add_counts <- table(factor(df$SNP, levels=0:2))
        add_labels <- paste0(c("Ref/Ref", "Het", "Hom Alt"), "\n(N=", add_counts, ")")
        df$SNP_f <- factor(df$SNP, levels=0:2, labels=add_labels)
        
        boxplot(expression ~ SNP_f, data=df, col="lightgreen", outline=FALSE, 
                ylab="Expression (Z-score)", xlab="Genotype",
                main=paste0("Additive Model\n(p = ", get_p(expression ~ SNP, df), ")"))
        points(jitter(as.numeric(df$SNP_f), amount=0.15), df$expression, 
                pch=df$p_pch, col=df$p_col, cex=df$p_cex)
        
        # Model 2: Dominant
        df$dom_val <- ifelse(df$SNP > 0, 1, 0)
        dom_counts <- table(factor(df$dom_val, levels=0:1))
        dom_labels <- paste0(c("Ref/Ref", "Any Alt"), "\n(N=", dom_counts, ")")
        df$dom_f <- factor(df$dom_val, levels=0:1, labels=dom_labels)
        
        boxplot(expression ~ dom_f, data=df, col="lightblue", outline=FALSE, 
                ylab="Expression (Z-score)", xlab="Genotype",
                main=paste0("Dominant Model\n(p = ", get_p(expression ~ dom_val, df), ")"))
        points(jitter(as.numeric(df$dom_f), amount=0.15), df$expression, 
                pch=df$p_pch, col=df$p_col, cex=df$p_cex)
 
        # Model 3: Recessive
        df$rec_val <- ifelse(df$SNP == 2, 1, 0)
        rec_counts <- table(factor(df$rec_val, levels=0:1))
        rec_labels <- paste0(c("Non-Hom Alt", "Hom Alt"), "\n(N=", rec_counts, ")")
        df$rec_f <- factor(df$rec_val, levels=0:1, labels=rec_labels)
        
        boxplot(expression ~ rec_f, data=df, col="lightpink", outline=FALSE, 
                ylab="Expression (Z-score)", xlab="Genotype",
                main=paste0("Recessive Model\n(p = ", get_p(expression ~ rec_val, df), ")"))
        points(jitter(as.numeric(df$rec_f), amount=0.15), df$expression, 
                pch=df$p_pch, col=df$p_col, cex=df$p_cex)
 
        mtext(paste("Gene:", G, "| SNP:", S), outer=TRUE, cex=1.1, font=2, line=0.5)
        plots_made <- plots_made + 1
    }
    
    dev.off()
    message(sprintf("# Plots made: %d, pairs skipped: %d", plots_made, pairs_skipped))
}
 
message(paste("# Processing", nrow(pairs), "gene-SNP pairs..."))
run_plotting(out_file_std, FALSE)
run_plotting(out_file_col, TRUE)
message("Done.")
