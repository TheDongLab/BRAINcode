#!/usr/bin/env Rscript
###########################################
# _eQTL_boxplot.R
# Fixed: Argument mapping and unlist coercion
###########################################

suppressPackageStartupMessages({
  library(data.table)
})

args <- commandArgs(TRUE)
# Mapping to run_eQTL.sh positions:
GSfile    <- args[1]  # $TOP_PAIRS
snp_file  <- args[2]  # $SNP_FILE
expr_file <- args[3]  # $EXPR_FILE
# args[4] is $BIM (ignored here)
# args[5] is $SNP_LOC (ignored here)
out_dir   <- args[6]  # $BOXPLOT_DIR

# Define output file path
out_file <- file.path(out_dir, "top_eQTL_boxplots.pdf")

# Load data
message("# Loading genomic matrices...")
pairs <- fread(GSfile, header=FALSE, col.names=c("geneid","snpid"))
snp_mat <- fread(snp_file, header=TRUE)
expr_mat <- fread(expr_file, header=TRUE)

# Start PDF
pdf(out_file, width=12, height=5)

for (i in seq_len(nrow(pairs))) {
    G <- pairs$geneid[i]; S <- pairs$snpid[i]
    
    # Efficient extraction using data.table keying logic
    expr_row <- expr_mat[gene_id == G]
    snp_row <- snp_mat[snpid == S]
    
    if (nrow(expr_row) == 0 || nrow(snp_row) == 0) next
    
    # Build data frame with explicit unlisting
    df <- data.frame(
        expression = as.numeric(unlist(expr_row[, -1, with = FALSE])),
        SNP = as.numeric(unlist(snp_row[, -1, with = FALSE]))
    )
    df <- df[!is.na(df$SNP), ]
    
    # Robust Stats Function
    get_stats <- function(model_formula, data) {
        if(var(data$expression, na.rm=TRUE) == 0) return(list(p="NS", beta=0))
        tryCatch({
            fit <- lm(model_formula, data=data)
            s <- summary(fit)$coefficients
            if(nrow(s) < 2) return(list(p="NA", beta=0))
            list(p=signif(s[2,4], 3), beta=signif(s[2,1], 2))
        }, error = function(e) list(p="Err", beta=0))
    }

    # Plotting Layout
    par(mfrow=c(1,3), mar=c(5,4,4,2), oma=c(0,0,2,0))
    
    # PANEL 1: Additive
    res <- get_stats(expression ~ SNP, df)
    df$SNP_f <- factor(df$SNP, levels=0:2, labels=c("0","1","2"))
    boxplot(expression ~ SNP_f, data=df, col="lightgreen", outline=F,
            ylab="Z-score", main=paste0("Additive (p=", res$p, ")"))
    stripchart(expression ~ SNP_f, data=df, vertical=T, method="jitter", add=T, pch=1, col="darkred")

    # PANEL 2: Dominant
    df$dom <- factor(ifelse(df$SNP > 0, 1, 0), levels=0:1, labels=c("0", "1+"))
    res_d <- get_stats(expression ~ dom, df)
    boxplot(expression ~ dom, data=df, col="lightblue", outline=F,
            main=paste0("Dom (p=", res_d$p, ")"))
    stripchart(expression ~ dom, data=df, vertical=T, method="jitter", add=T, pch=1, col="darkred")

    # PANEL 3: Recessive
    df$rec <- factor(ifelse(df$SNP == 2, 1, 0), levels=0:1, labels=c("<2", "2"))
    res_r <- get_stats(expression ~ rec, df)
    boxplot(expression ~ rec, data=df, col="lightpink", outline=F,
            main=paste0("Rec (p=", res_r$p, ")"))
    stripchart(expression ~ rec, data=df, vertical=T, method="jitter", add=T, pch=1, col="darkred")

    mtext(paste("Gene:", G, "| SNP:", S), outer=TRUE, cex=1.2, font=2)
}

dev.off()
message(paste("Done. Boxplots saved to:", out_file))
