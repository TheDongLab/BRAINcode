#!/usr/bin/env Rscript
###########################################
# _eQTL_boxplot.R
# Merges New layout with Z-score robustness
###########################################

suppressPackageStartupMessages({
  library(data.table)
})

args <- commandArgs(TRUE)
GSfile <- args[1]    # gene.snp.list (No header: geneid snpid)
snp_file <- args[2]  # snp_Cerebellum.txt
expr_file <- args[3] # expression_Zscores.txt
out_file <- args[4]  # output.pdf

# Load data once
message("# Loading genomic matrices...")
pairs <- fread(GSfile, header=FALSE, col.names=c("geneid","snpid"))
snp_mat <- fread(snp_file, header=TRUE)
expr_mat <- fread(expr_file, header=TRUE)

# Start single PDF
pdf(out_file, width=12, height=5)

for (i in seq_len(nrow(pairs))) {
    G <- pairs$geneid[i]; S <- pairs$snpid[i]
    
    # Efficient extraction
    expr_row <- expr_mat[gene_id == G]
    snp_row <- snp_mat[snpid == S]
    
    if (nrow(expr_row) == 0 || nrow(snp_row) == 0) next
    
    # Create the dataframe - Fixed the trailing comma and unlist logic
    df <- data.frame(
        expression = as.numeric(unlist(expr_row[, -1, with = FALSE])),
        SNP = as.numeric(unlist(snp_row[, -1, with = FALSE]))
    )
    
    # Remove NA genotypes
    df <- df[!is.na(df$SNP), ]
    
    # 2026 Robust Stats Function
    get_stats <- function(model_formula, data) {
        # Check for zero variance to avoid noise errors
        if(var(data$expression, na.rm=TRUE) == 0) return(list(p="NS", beta=0))
        tryCatch({
            fit <- lm(model_formula, data=data)
            s <- summary(fit)$coefficients
            # Check if the coefficient exists (enough levels in factor)
            if(nrow(s) < 2) return(list(p="NA", beta=0))
            list(p=signif(s[2,4], 3), beta=signif(s[2,1], 2))
        }, error = function(e) list(p="Err", beta=0))
    }

    # Layout Setup
    par(mfrow=c(1,3), mar=c(5,4,4,2), oma=c(0,0,2,0))
    
    # PANEL 1: Additive
    res <- get_stats(expression ~ SNP, df)
    df$SNP_f <- factor(df$SNP, levels=0:2, labels=c("Ref/Ref","Ref/Alt","Alt/Alt"))
    
    boxplot(expression ~ SNP_f, data=df, col="lightgreen", outline=F,
            ylab="Normalized Expression (Z-score)", 
            main=paste0("Additive (p=", res$p, ", b=", res$beta, ")"))
    stripchart(expression ~ SNP_f, data=df, vertical=T, method="jitter", add=T, pch=1, col="darkred")

    # PANEL 2: Dominant (Has any Alt)
    df$dom <- factor(ifelse(df$SNP > 0, 1, 0), levels=0:1, labels=c("No Alt", "With Alt"))
    res_d <- get_stats(expression ~ dom, df)
    
    boxplot(expression ~ dom, data=df, col="lightblue", outline=F,
            main=paste0("Dom (p=", res_d$p, ", b=", res_d$beta, ")"))
    stripchart(expression ~ dom, data=df, vertical=T, method="jitter", add=T, pch=1, col="darkred")

    # PANEL 3: Recessive (Two Alts)
    df$rec <- factor(ifelse(df$SNP == 2, 1, 0), levels=0:1, labels=c("Not Hom Alt", "Hom Alt"))
    res_r <- get_stats(expression ~ rec, df)
    
    boxplot(expression ~ rec, data=df, col="lightpink", outline=F,
            main=paste0("Rec (p=", res_r$p, ", b=", res_r$beta, ")"))
    stripchart(expression ~ rec, data=df, vertical=T, method="jitter", add=T, pch=1, col="darkred")

    # Add global title
    mtext(paste("Gene:", G, "| SNP:", S), outer=TRUE, cex=1.2, font=2)
}

dev.off()
