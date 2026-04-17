#!/usr/bin/env Rscript
###########################################
# _eQTL_boxplot.R
# FIXED: Added ALS vs Control coloring and tissue-specific naming
###########################################

suppressPackageStartupMessages({
  library(data.table)
})

args <- commandArgs(TRUE)
GSfile    <- args[1]  # $TOP_PAIRS
snp_file  <- args[2]  # $SNP_FILE
expr_file <- args[3]  # $EXPR_FILE
cov_file  <- args[4]  # We'll pass $COV_FILE here now
# args[5] is $SNP_LOC
out_dir   <- args[6]  # $OUTDIR
tissue    <- args[7]  # We'll pass $TISSUE_DIR as a new arg

# Define output file paths
out_file_std <- file.path(out_dir, paste0(tissue, "_top_eQTL_boxplots.pdf"))
out_file_col <- file.path(out_dir, paste0(tissue, "_top_eQTL_boxplots_colored.pdf"))

message("# Loading genomic matrices...")
pairs <- fread(GSfile, header=FALSE, col.names=c("geneid","snpid"))
snp_mat <- fread(snp_file, header=TRUE); setkey(snp_mat, snpid)
expr_mat <- fread(expr_file, header=TRUE); setkey(expr_mat, geneid)
cov_mat <- fread(cov_file, header=TRUE); setkey(cov_mat, id)

# Identify ALS status for coloring (is_als row)
als_status <- as.numeric(cov_mat[id == "is_als", -1, with=FALSE])
# Map 1 -> Red (ALS), 0 -> Violet (Control)
dot_colors <- ifelse(als_status == 1, "red", "darkorchid")

# Function to run the plotting loop
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
            col = if(use_status_colors) dot_colors else rep("darkred", length(als_status))
        )
        df <- df[!is.na(df$SNP), ]
        
        get_stats <- function(model_formula, data) {
            tryCatch({
                fit <- lm(model_formula, data=data)
                s <- summary(fit)$coefficients
                list(p=signif(s[2,4], 3), beta=signif(s[2,1], 2))
            }, error = function(e) list(p="Err", beta=0))
        }

        par(mfrow=c(1,3), mar=c(5,4,4,2), oma=c(0,0,3,0))
        
        get_n_labels <- function(vec, base_labels) {
            counts <- table(factor(vec, levels = 0:(length(base_labels)-1)))
            paste0(base_labels, "\n(N=", counts, ")")
        }

        # Sub-function for consistent panels
        plot_panel <- function(formula, data, labels, title_prefix, color) {
            boxplot(formula, data=data, col=color, outline=F, 
                    ylab="Expression (Z-score)", xlab="Genotype",
                    main=paste0(title_prefix, "\n(p = ", get_stats(formula, data)$p, ")"))
            stripchart(formula, data=data, vertical=T, method="jitter", add=T, 
                       pch=16, col=data$col, cex=0.7)
        }

        # Additive
        add_labels <- get_n_labels(df$SNP, c("Ref/Ref", "Het", "Hom Alt"))
        df$SNP_f <- factor(df$SNP, levels=0:2, labels=add_labels)
        plot_panel(expression ~ SNP_f, df, add_labels, "Additive Model", "lightgreen")

        # Dominant
        df$dom_val <- ifelse(df$SNP > 0, 1, 0)
        dom_labels <- get_n_labels(df$dom_val, c("Ref/Ref", "Any Alt"))
        df$dom <- factor(df$dom_val, levels=0:1, labels=dom_labels)
        plot_panel(expression ~ dom, df, dom_labels, "Dominant Model", "lightblue")

        # Recessive
        df$rec_val <- ifelse(df$SNP == 2, 1, 0)
        rec_labels <- get_n_labels(df$rec_val, c("Non-Hom Alt", "Hom Alt"))
        df$rec <- factor(df$rec_val, levels=0:1, labels=rec_labels)
        plot_panel(expression ~ rec, df, rec_labels, "Recessive Model", "lightpink")

        header_text <- paste("Gene:", G, "| SNP:", S)
        if(use_status_colors) header_text <- paste(header_text, "\n(Red = ALS, Violet = Control)")
        mtext(header_text, outer=TRUE, cex=1, font=2, line=0.5)
    }
    dev.off()
}

message("# Generating standard boxplots...")
run_plotting(out_file_std, use_status_colors = FALSE)

message("# Generating ALS-status colored boxplots...")
run_plotting(out_file_col, use_status_colors = TRUE)

message("Done. Files saved to results directory.")
