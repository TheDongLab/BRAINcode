#!/usr/bin/env Rscript
###########################################
# _eQTL_boxplot.R
# Purpose: Generate cleaned, statistically robust boxplots
#   - Filters for expression > 0.01 
#   - Dynamic fallback to t-test if a genotype bin is empty 
#   - Professional labeling ("Reference Alleles")
###########################################

suppressPackageStartupMessages(library(data.table))

args <- commandArgs(TRUE)
pairs_file <- args[1]; snp_file <- args[2]; expr_file <- args[3]
bim_file <- args[4]; snploc_file <- args[5]; out_dir <- args[6]

dir.create(out_dir, showWarnings=FALSE, recursive=TRUE)

# Load data
pairs  <- fread(pairs_file, header=FALSE, col.names=c("geneid","snpid"))
snploc <- fread(snploc_file, header=TRUE)
setnames(snploc, c("snpid","chr","pos"))
snp_mat  <- fread(snp_file, header=TRUE)[snpid %in% unique(pairs$snpid)]
expr_mat <- fread(expr_file, header=TRUE)[gene_id %in% unique(pairs$geneid)]

for (i in seq_len(nrow(pairs))) {
    G <- pairs$geneid[i]; S <- pairs$snpid[i]
    snp_row <- snp_mat[snpid == S]; expr_row <- expr_mat[gene_id == G]
    if (nrow(snp_row) == 0 || nrow(expr_row) == 0) next
    
    df <- data.frame(
        expression = as.numeric(unlist(expr_row[, -1, with=FALSE])),
        SNP        = as.numeric(unlist(snp_row[, -1, with=FALSE]))
    )
    df <- df[!is.na(df$SNP) & !is.na(df$expression), ]
    df$SNP_factor <- factor(df$SNP, levels=c(0, 1, 2), labels=c("Hom Ref", "Het", "Hom Alt"))

    # Diversity Check: At least 2 groups with N >= 3
    actual_counts <- table(df$SNP_factor)
    if (sum(actual_counts >= 3) < 2) next

    # Stats: Additive (Linear Reg) or Fallback (t-test)
    bins_with_data <- names(actual_counts[actual_counts > 0])
    if (length(bins_with_data) == 3) {
        p_val <- summary(lm(expression ~ SNP, data=df))$coefficients[2,4]
        test_label <- "Linear Reg"
    } else {
        top_two <- names(sort(actual_counts, decreasing=TRUE)[1:2])
        p_val <- t.test(df$expression[df$SNP_factor == top_two[1]], 
                        df$expression[df$SNP_factor == top_two[2]])$p.value
        test_label <- paste("t-test:", top_two[1], "vs", top_two[2])
    }

    # CALCULATE VERTICAL BUFFER
    y_range <- range(df$expression, na.rm=TRUE)
    y_buffer <- diff(y_range) * 0.15  # 15% room for comfort
    y_lims <- c(y_range[1] - y_buffer, y_range[2] + y_buffer)

    out_pdf <- file.path(out_dir, paste0("eQTLboxplot.", G, ".", gsub("[^A-Za-z0-9]","_", S), ".pdf"))
    pdf(out_pdf, width=7, height=6)
    par(mfrow=c(1,2), mar=c(6,4,4,1), oma=c(0,0,3,0))
    
    # PANEL 1: MAIN PLOT
    boxplot(expression ~ SNP_factor, data=df, ylim=y_lims, ylab="Normalized Expression", 
            xaxt='n', main=paste0(test_label, "\n(p=", signif(p_val, 3), ")"), 
            col=rgb(0.5, 0.8, 0.5, 0.5), outline=FALSE)
    
    # HORIZONTAL JITTER: Spreads dots left-to-right, keeps Y-axis perfect
    stripchart(expression ~ SNP_factor, data=df, vertical=TRUE, method="jitter", 
               jitter=0.2, pch=16, col=rgb(0.6, 0, 0, 0.4), cex=0.8, add=TRUE)
    
    axis(1, at=1:3, labels=levels(df$SNP_factor), cex.axis=0.7, font=2)
    mtext(paste0("N=", as.integer(actual_counts)), side=1, line=2, at=1:3, cex=0.7)

    # PANEL 2: DOMINANT MODEL
    df$has_alt <- factor(ifelse(df$SNP > 0, 1, 0), levels=0:1)
    dom_counts <- table(df$has_alt)
    boxplot(expression ~ has_alt, data=df, ylim=y_lims, ylab="Normalized Expression", 
            xaxt='n', main="Dominant Model", col=rgb(0.5, 0.7, 0.9, 0.5), outline=FALSE)
    
    stripchart(expression ~ has_alt, data=df, vertical=TRUE, method="jitter", 
               jitter=0.2, pch=16, col=rgb(0.6, 0, 0, 0.4), cex=0.8, add=TRUE)
    
    axis(1, at=1:2, labels=c("Reference Alleles", "Any Alt Alleles"), cex.axis=0.7, font=2)
    mtext(paste0("N=", as.integer(dom_counts)), side=1, line=2, at=1:2, cex=0.7)

    mtext(paste0("cis-eQTL: ", G, "\nSNP: ", (if(nrow(snploc[snpid==S])>0) paste0(snploc[snpid==S]$chr,":",snploc[snpid==S]$pos) else S)), 
          outer=TRUE, cex=1, font=2)
    dev.off()
}
