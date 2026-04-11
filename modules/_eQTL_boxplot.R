#!/usr/bin/env Rscript
###########################################
# _eQTL_boxplot.R
# Purpose: Generate cleaned, statistically robust boxplots
#   - Filters for expression > 0.01 
#   - Dynamic fallback to t-test if a genotype bin is empty 
#   - Professional labeling ("Reference Alleles")
###########################################

suppressPackageStartupMessages(library(data.table))

args        <- commandArgs(TRUE)
pairs_file  <- args[1]
snp_file    <- args[2]
expr_file   <- args[3]
bim_file    <- args[4]
snploc_file <- args[5]
out_dir     <- args[6]

dir.create(out_dir, showWarnings=FALSE, recursive=TRUE)

# Load metadata
pairs  <- fread(pairs_file, header=FALSE, col.names=c("geneid","snpid"))
snploc <- fread(snploc_file, header=TRUE)
setnames(snploc, c("snpid","chr","pos"))

# Load matrices
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
    
    # 1. Prepare Data Frame
    df <- data.frame(
        expression = as.numeric(unlist(expr_row[, -1, with=FALSE])),
        SNP        = as.numeric(unlist(snp_row[, -1, with=FALSE]))
    )
    df <- df[!is.na(df$SNP) & !is.na(df$expression), ]
    df$SNP_factor <- factor(df$SNP, levels=c(0, 1, 2), labels=c("Hom Ref", "Het", "Hom Alt"))

    # 2. APPLY BIOLOGICAL FILTER: Expression > 0.01 in at least 1 group with N >= 3
    # This ensures we aren't plotting "noise" near zero 
    valid_groups <- df[df$expression > 0.01, ]
    group_counts_filtered <- table(valid_groups$SNP_factor)
    
    if (max(group_counts_filtered) < 3) {
        message(sprintf("   SKIP: %s x %s (Low expression/N: max group N=%d above 0.01)", G, S, max(group_counts_filtered)))
        next
    }

    # 3. DYNAMIC STATISTICAL TEST 
    # Check how many genotype bins actually have data (N > 0)
    actual_counts <- table(df$SNP_factor)
    bins_with_data <- names(actual_counts[actual_counts > 0])
    
    test_label <- ""
    p_val <- NA

    if (length(bins_with_data) == 3) {
        # Standard Additive Model (Linear Regression)
        fit <- lm(expression ~ SNP, data=df)
        p_val <- summary(fit)$coefficients[2,4]
        test_label <- "Linear Reg"
    } else if (length(bins_with_data) == 2) {
        # Fallback to Welch's t-test between the two existing groups
        group1 <- df$expression[df$SNP_factor == bins_with_data[1]]
        group2 <- df$expression[df$SNP_factor == bins_with_data[2]]
        t_test <- t.test(group1, group2)
        p_val <- t_test$p.value
        test_label <- paste("t-test:", bins_with_data[1], "vs", bins_with_data[2])
    } else {
        message(sprintf("   SKIP: %s x %s (Insufficient diversity for any test)", G, S))
        next
    }

    # 4. PLOTTING
    loc_row   <- snploc[snpid == S]
    pos_label <- if (nrow(loc_row) > 0) paste0(loc_row$chr[1], ":", loc_row$pos[1]) else S
    out_pdf   <- file.path(out_dir, paste0("eQTLboxplot.", G, ".", gsub("[^A-Za-z0-9]","_", S), ".pdf"))

    pdf(out_pdf, width=7, height=6)
    par(mfrow=c(1,2), mar=c(6,4,4,1), oma=c(0,0,3,0))

    # Panel A: Additive/Categorical
    boxplot(expression ~ SNP_factor, data=df, 
            ylab="Normalized Expression", xaxt='n',
            main=paste0(test_label, "\n(p=", signif(p_val, 3), ")"),
            col='lightgreen', outline=FALSE)
    stripchart(expression ~ SNP_factor, data=df, vertical=TRUE, method="jitter", 
               pch=1, col="darkred", cex=0.7, add=TRUE)
    
    # Custom labels with N counts [cite: 7, 9, 10]
    axis(1, at=1:3, labels=FALSE)
    mtext(levels(df$SNP_factor), side=1, line=0.5, at=1:3, cex=0.75, font=2)
    mtext(paste0("N=", as.integer(actual_counts)), side=1, line=1.8, at=1:3, cex=0.7)

    # Panel B: Dominant Model (Any Alt)
    df$has_alt <- factor(ifelse(df$SNP > 0, 1, 0), levels=0:1)
    dom_counts <- table(df$has_alt)
    
    # Use your new naming conventions
    boxplot(expression ~ has_alt, data=df, 
            ylab="Normalized Expression", xaxt='n', 
            main="Dominant Model", col='lightblue', outline=FALSE)
    stripchart(expression ~ has_alt, data=df, vertical=TRUE, method="jitter", 
               pch=1, col="darkred", cex=0.7, add=TRUE)
    
    axis(1, at=1:2, labels=FALSE)
    mtext(c("Reference Alleles", "Any Alt Alleles"), side=1, line=0.5, at=1:2, cex=0.75, font=2)
    mtext(paste0("N=", as.integer(dom_counts)), side=1, line=1.8, at=1:2, cex=0.7)

    # Global Title
    mtext(paste0("cis-eQTL: ", G, "\nSNP: ", pos_label), outer=TRUE, cex=1, font=2)
    dev.off()
}

message("## Done.")
