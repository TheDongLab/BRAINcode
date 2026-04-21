#!/usr/bin/env Rscript
###########################################
# _eQTL_boxplot.R - Fixed Duplication
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

message(sprintf("# Pairs to plot: %d", nrow(pairs)))
message(sprintf("# SNP matrix: %d SNPs x %d samples", nrow(snp_mat), ncol(snp_mat)-1))
message(sprintf("# Expression matrix: %d genes x %d samples", nrow(expr_mat), ncol(expr_mat)-1))

# Get the first column name (the row identifier)
cov_first_col <- names(cov_mat)[1]

# Find the 'is_als' row
als_row <- cov_mat[get(cov_first_col) == "is_als", ]

if (nrow(als_row) == 0) {
  message("# WARNING: 'is_als' row not found. Using all purple points.")
  als_vals <- rep(0, ncol(cov_mat) - 1)
} else {
  # Extract ALS status - use as.vector() to avoid duplication
  als_vals <- as.vector(as.numeric(als_row[, -1, with=FALSE]))
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
        
        # Get rows for this gene and SNP
        expr_row <- expr_mat[geneid == G]
        snp_row  <- snp_mat[snpid == S]
        
        if (nrow(expr_row) == 0 || nrow(snp_row) == 0) {
            pairs_skipped <- pairs_skipped + 1
            next
        }
        
        # CRITICAL FIX: Use as.vector() to extract without duplication
        # This converts the data.table row to a vector correctly
        expr_vals <- as.vector(as.numeric(expr_row[, -1, with=FALSE]))
        snp_vals  <- as.vector(as.numeric(snp_row[, -1, with=FALSE]))
        
        # Safety check
        if (length(expr_vals) == 0 || length(snp_vals) == 0) {
            pairs_skipped <- pairs_skipped + 1
            next
        }
        
        # Verify lengths match
        if (length(expr_vals) != length(snp_vals)) {
            message(sprintf("# WARNING: Length mismatch for %s x %s: expr=%d, snp=%d", 
                          G, S, length(expr_vals), length(snp_vals)))
            pairs_skipped <- pairs_skipped + 1
            next
        }
        
        df <- data.frame(
            expression = expr_vals,
            SNP = snp_vals,
            p_col = if(use_status_colors) p_colors else rep("darkred", length(expr_vals)),
            p_pch = if(use_status_colors) p_shapes else rep(16, length(expr_vals)),
            p_cex = if(use_status_colors) p_sizes else rep(0.7, length(expr_vals)),
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
