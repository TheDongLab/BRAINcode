#!/usr/bin/env Rscript
###########################################
# _eQTL_boxplot.R
# FIXED: Robust float-handling for ALS status coloring
###########################################

suppressPackageStartupMessages({
  library(data.table)
})

args <- commandArgs(TRUE)
GSfile    <- args[1]  # $TOP_PAIRS
snp_file  <- args[2]  # $SNP_FILE
expr_file <- args[3]  # $EXPR_FILE
cov_file  <- args[4]  # $COV_FILE
# args[5] is $SNP_LOC
out_dir   <- args[6]  # $OUTDIR
tissue    <- args[7]  # $TISSUE_DIR

out_file_std <- file.path(out_dir, paste0(tissue, "_top_eQTL_boxplots.pdf"))
out_file_col <- file.path(out_dir, paste0(tissue, "_top_eQTL_boxplots_colored.pdf"))

message("# Loading genomic matrices...")
pairs <- fread(GSfile, header=FALSE, col.names=c("geneid","snpid"))
snp_mat <- fread(snp_file, header=TRUE); setkey(snp_mat, snpid)
expr_mat <- fread(expr_file, header=TRUE); setkey(expr_mat, geneid)
cov_mat <- fread(cov_file, header=TRUE); setkey(cov_mat, id)

# --- ROBUST ALS STATUS EXTRACTION ---
# We force the row to numeric to handle '1.0' vs '1'. The '-1' drops the 'id' column so we only have the sample values
als_vals <- as.numeric(unlist(cov_mat[id == "is_als", -1, with=FALSE]))

# Check if we actually found the row
if(length(als_vals) == 0) {
    stop("ERROR: 'is_als' row not found in covariate file. Check your header labels.")
}

# Threshold check (handles floats perfectly)
# Red = ALS (> 0.5), DarkOrchid = Control (< 0.5)
dot_colors <- ifelse(als_vals > 0.5, "red", "darkorchid")
# ------------------------------------

run_plotting <- function(pdf_path, use_status_colors = FALSE) {
    pdf(pdf_path, width=12, height=5)
    
    for (i in seq_len(nrow(pairs))) {
        G <- pairs$geneid[i]; S <- pairs$snpid[i]
        expr_row <- expr_mat[geneid == G]
        snp_row <- snp_mat[snpid == S]
        
        if (nrow(expr_row) == 0 || nrow(snp_row) == 0) next
        
        df <- data.frame(
            expression = as.numeric(unlist(expr_row[, -1, with = FALSE])),
            SNP = as.numeric(unlist(snp_row[, -1, with = FALSE])),
            p_col = if(use_status_colors) dot_colors else rep("darkred", length(als_vals))
        )
        df <- df[!is.na(df$SNP), ]
        
        get_stats <- function(model_formula, data) {
            tryCatch({
                fit <- lm(model_formula, data=data)
                s <- summary(fit)$coefficients
                # Handle cases where a model might not be rank-sufficient
                if(nrow(s) < 2) return(list(p="NA"))
                list(p=signif(s[2,4], 3))
            }, error = function(e) list(p="Err"))
        }

        par(mfrow=c(1,3), mar=c(5,4,4,2), oma=c(0,0,3,0))
        
        plot_panel <- function(formula, data, title_prefix, box_color) {
            # Draw the boxplot
            boxplot(formula, data=data, col=box_color, outline=F, 
                    ylab="Expression (Z-score)", xlab="Genotype",
                    main=paste0(title_prefix, "\n(p = ", get_stats(formula, data)$p, ")"))
            
            # Overlay the individual points with our defined color vector
            stripchart(formula, data=data, vertical=T, method="jitter", add=T, 
                       pch=16, col=data$p_col, cex=0.7)
            
            # Add a legend only to the first panel if coloring is active
            if(use_status_colors && title_prefix == "Additive Model") {
                legend("topleft", legend=c("ALS", "Control"), 
                       col=c("red", "darkorchid"), pch=16, cex=0.8, bty="n")
            }
        }

        # Panel 1: Additive (0, 1, 2)
        add_counts <- table(factor(df$SNP, levels=0:2))
        add_labels <- paste0(c("Ref/Ref", "Het", "Hom Alt"), "\n(N=", add_counts, ")")
        df$SNP_f <- factor(df$SNP, levels=0:2, labels=add_labels)
        plot_panel(expression ~ SNP_f, df, "Additive Model", "lightgreen")

        # Panel 2: Dominant (0 vs 1+2)
        df$dom_val <- ifelse(df$SNP > 0, 1, 0)
        dom_counts <- table(factor(df$dom_val, levels=0:1))
        dom_labels <- paste0(c("Ref/Ref", "Any Alt"), "\n(N=", dom_counts, ")")
        df$dom <- factor(df$dom_val, levels=0:1, labels=dom_labels)
        plot_panel(expression ~ dom, df, "Dominant Model", "lightblue")

        # Panel 3: Recessive (0+1 vs 2)
        df$rec_val <- ifelse(df$SNP == 2, 1, 0)
        rec_counts <- table(factor(df$rec_val, levels=0:1))
        rec_labels <- paste0(c("Non-Hom Alt", "Hom Alt"), "\n(N=", rec_counts, ")")
        df$rec <- factor(df$rec_val, levels=0:1, labels=rec_labels)
        plot_panel(expression ~ rec, df, "Recessive Model", "lightpink")

        # Header for the whole page
        header_line <- paste("Gene:", G, "| SNP:", S)
        mtext(header_line, outer=TRUE, cex=1.2, font=2, line=0.5)
    }
    dev.off()
}

message("# Generating Standard Boxplots...")
run_plotting(out_file_std, FALSE)

message("# Generating Colored Boxplots (Red=ALS, Violet=Control)...")
run_plotting(out_file_col, TRUE)

message("DONE: Boxplots saved to results.")
