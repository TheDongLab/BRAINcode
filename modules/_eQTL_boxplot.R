#!/usr/bin/env Rscript
###########################################
# _eQTL_boxplot.R
# Purpose: Generate expression-vs-genotype boxplots for top
#           gene-SNP pairs from cis-eQTL results.
#
# Adapted from Xianjun Dong's original _eQTL_boxplot.R to work with:
#   - Unnamed PLINK SNP IDs (._AC format, no rsID)
#   - Allele info pulled from .bim file instead of SNP name parsing
#   - Expression matrix subsetted to tissue samples
#
# Usage:
#   Rscript _eQTL_boxplot.R \
#       top_pairs.txt \        # two-col: geneid snpid (no header)
#       snp_matrix.txt \       # snp_TISSUE.txt from prep
#       expression_matrix.txt \# expression_TISSUE.txt from prep
#       bim_file.bim \         # for allele lookup
#       snp_location.txt \     # for chr:pos lookup
#       output_dir
###########################################

suppressPackageStartupMessages(library(data.table))

args       <- commandArgs(TRUE)
pairs_file <- args[1]
snp_file   <- args[2]
expr_file  <- args[3]
bim_file   <- args[4]
snploc_file<- args[5]
out_dir    <- args[6]

dir.create(out_dir, showWarnings=FALSE, recursive=TRUE)

# Load metadata
pairs  <- fread(pairs_file, header=FALSE, col.names=c("geneid","snpid"))
snploc <- fread(snploc_file, header=TRUE)
setnames(snploc, c("snpid","chr","pos"))

# Load matrices (subsetted for speed)
snp_mat  <- fread(snp_file, header=TRUE)[snpid %in% unique(pairs$snpid)]
expr_mat <- fread(expr_file, header=TRUE)[gene_id %in% unique(pairs$geneid)]

extract_allele <- function(snp_name) { sub("^[^_]*_", "", snp_name) }

message("## Processing pairs...")
for (i in seq_len(nrow(pairs))) {
    G <- pairs$geneid[i]
    S <- pairs$snpid[i]

    snp_row <- snp_mat[snpid == S]
    expr_row <- expr_mat[gene_id == G]
    if (nrow(snp_row) == 0 || nrow(expr_row) == 0) next
    
    df <- data.frame(
        expression = as.numeric(expr_row[, -1, with=FALSE]),
        SNP        = as.numeric(snp_row[, -1, with=FALSE])
    )
    df <- df[!is.na(df$SNP) & !is.na(df$expression), ]
    df$SNP_factor <- factor(df$SNP, levels=2:0)

    # RELEVANT CHECK: Ensure at least 2 genotype groups have 3+ samples.
    # This prevents the "flat-line" p=0 plot you encountered[cite: 21, 22, 23].
    group_counts <- table(df$SNP_factor)
    if (sum(group_counts >= 3) < 2) {
        message(sprintf("   SKIP: %s x %s (Insufficent diversity: %s)", G, S, paste(group_counts, collapse=",")))
        next
    }

    alt_allele <- extract_allele(S)
    loc_row    <- snploc[snpid == S]
    pos_label  <- if (nrow(loc_row) > 0) paste0(loc_row$chr[1], ":", loc_row$pos[1]) else S
    
    out_pdf <- file.path(out_dir, paste0("eQTLboxplot.", G, ".", gsub("[^A-Za-z0-9]","_", S), ".pdf"))

    pdf(out_pdf, width=7, height=5)
    par(mfrow=c(1,2), mar=c(5,4,3,1), oma=c(0,0,3,0))

    # Additive Panel
    fit <- lm(expression ~ SNP, data=df)
    boxplot(expression ~ SNP_factor, data=df, ylab="Normalized Expression", xaxt='n',
            main=paste0("Additive (p=", signif(summary(fit)$coefficients[2,4], 3), ")"),
            col='lightgreen', outline=FALSE)
    stripchart(expression ~ SNP_factor, data=df, vertical=TRUE, method="jitter", pch=1, col="darkred", cex=0.7, add=TRUE)
    mtext(c("Hom Alt", "Het", "Hom Ref"), side=1, line=0.5, at=1:3, cex=0.7)
    mtext(paste0("N=", as.integer(group_counts)), side=1, line=1.5, at=1:3, cex=0.7)

    # Dominant Panel
    df$has_alt <- factor(ifelse(df$SNP > 0, 1, 0), levels=0:1)
    if (min(table(df$has_alt)) > 2) {
        boxplot(expression ~ has_alt, data=df, ylab="Normalized Expression", xaxt='n', col='lightblue', outline=FALSE)
        stripchart(expression ~ has_alt, data=df, vertical=TRUE, method="jitter", pch=1, col="darkred", cex=0.7, add=TRUE)
        mtext(c("Ref/Ref", "Any Alt"), side=1, line=0.5, at=1:2, cex=0.7)
    }

    mtext(paste0("cis-eQTL: ", G, "\nSNP: ", pos_label), outer=TRUE, cex=0.9, font=2)
    dev.off()
}
