#!/usr/bin/env Rscript
###########################################
# _eQTL_boxplot.R
# FIXED: geneid variable mapping and data.table keys
###########################################

suppressPackageStartupMessages({
  library(data.table)
})

args <- commandArgs(TRUE)
GSfile    <- args[1]  # $TOP_PAIRS
snp_file  <- args[2]  # $SNP_FILE
expr_file <- args[3]  # $EXPR_FILE
out_dir   <- args[6]  # $BOXPLOT_DIR

out_file <- file.path(out_dir, "top_eQTL_boxplots.pdf")

message("# Loading genomic matrices...")
# Load and set keys for fast lookup
pairs <- fread(GSfile, header=FALSE, col.names=c("geneid","snpid"))
snp_mat <- fread(snp_file, header=TRUE); setkey(snp_mat, snpid)
expr_mat <- fread(expr_file, header=TRUE); setkey(expr_mat, geneid)

pdf(out_file, width=12, height=5)

for (i in seq_len(nrow(pairs))) {
    G <- pairs$geneid[i]; S <- pairs$snpid[i]
    
    # Use the 'geneid' column name we established in the Python prep script
    expr_row <- expr_mat[geneid == G]
    snp_row <- snp_mat[snpid == S]
    
    if (nrow(expr_row) == 0 || nrow(snp_row) == 0) next
    
    df <- data.frame(
        expression = as.numeric(unlist(expr_row[, -1, with = FALSE])),
        SNP = as.numeric(unlist(snp_row[, -1, with = FALSE]))
    )
    df <- df[!is.na(df$SNP), ]
    
    get_stats <- function(model_formula, data) {
        if(var(data$expression, na.rm=TRUE) == 0) return(list(p="NS", beta=0))
        tryCatch({
            fit <- lm(model_formula, data=data)
            s <- summary(fit)$coefficients
            if(nrow(s) < 2) return(list(p="NA", beta=0))
            list(p=signif(s[2,4], 3), beta=signif(s[2,1], 2))
        }, error = function(e) list(p="Err", beta=0))
    }

    par(mfrow=c(1,3), mar=c(5,4,4,2), oma=c(0,0,3,0))
    ylab_text <- "Normalized Expression (Z-score)"

    get_n_labels <- function(vec, base_labels) {
        counts <- table(factor(vec, levels = 0:(length(base_labels)-1)))
        paste0(base_labels, "\n(N=", counts, ")")
    }

    # PANEL 1: Additive
    res <- get_stats(expression ~ SNP, df)
    add_labels <- get_n_labels(df$SNP, c("Ref/Ref", "Het", "Hom Alt"))
    df$SNP_f <- factor(df$SNP, levels=0:2, labels=add_labels)
    
    boxplot(expression ~ SNP_f, data=df, col="lightgreen", outline=F,
            ylab=ylab_text, xlab="Genotype",
            main=paste0("Additive Model\n(p = ", res$p, ")"))
    stripchart(expression ~ SNP_f, data=df, vertical=T, method="jitter", add=T, pch=1, col="darkred")

    # PANEL 2: Dominant
    df$dom_val <- ifelse(df$SNP > 0, 1, 0)
    dom_labels <- get_n_labels(df$dom_val, c("Ref/Ref", "Any Alt"))
    df$dom <- factor(df$dom_val, levels=0:1, labels=dom_labels)
    res_d <- get_stats(expression ~ dom_val, df)
    
    boxplot(expression ~ dom, data=df, col="lightblue", outline=F,
            ylab=ylab_text, xlab="Genotype Grouping",
            main=paste0("Dominant Model\n(p = ", res_d$p, ")"))
    stripchart(expression ~ dom, data=df, vertical=T, method="jitter", add=T, pch=1, col="darkred")

    # PANEL 3: Recessive
    df$rec_val <- ifelse(df$SNP == 2, 1, 0)
    rec_labels <- get_n_labels(df$rec_val, c("Non-Hom Alt", "Hom Alt"))
    df$rec <- factor(df$rec_val, levels=0:1, labels=rec_labels)
    res_r <- get_stats(expression ~ rec_val, df)
    
    boxplot(expression ~ rec, data=df, col="lightpink", outline=F,
            ylab=ylab_text, xlab="Genotype Grouping",
            main=paste0("Recessive Model\n(p = ", res_r$p, ")"))
    stripchart(expression ~ rec, data=df, vertical=T, method="jitter", add=T, pch=1, col="darkred")

    mtext(paste("Gene:", G, "| SNP:", S), outer=TRUE, cex=1.3, font=2, line=0.5)
}

dev.off()
message(paste("Done. Boxplots saved to:", out_file))
