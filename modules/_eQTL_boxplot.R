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

message("## Reading gene-SNP pairs...")
pairs <- fread(pairs_file, header=FALSE, col.names=c("geneid","snpid"))
message(sprintf("   %d pairs to plot", nrow(pairs)))

message("## Reading SNP location file (for chr:pos labels)...")
snploc <- fread(snploc_file, header=TRUE)
setnames(snploc, c("snpid","chr","pos"))

message("## Reading .bim for allele information...")
# .bim columns: chr, snpid, cm, pos, allele1(ALT), allele2(REF)
bim <- fread(bim_file, header=FALSE,
             col.names=c("chr","bim_snpid","cm","pos","ALT","REF"))
# .raw SNP names are SNPID_ALLELE — build lookup from bim row index
# Since all SNPs are unnamed (.), we match by row position
# The .raw header SNP names follow the pattern ._ALLELE
# We extract the allele from the raw SNP name directly
extract_allele <- function(snp_name) {
    # ._AC -> ALT allele is "AC"
    sub("^[^_]*_", "", snp_name)
}

message("## Loading SNP matrix (this may take a moment)...")
# Read only the rows we need using fread with grep
snp_ids_needed <- unique(pairs$snpid)

# fread full header to get sample order, then subset rows
snp_header <- fread(snp_file, nrows=0)
all_samples <- colnames(snp_header)[-1]   # drop 'snpid' col

# Read SNP matrix subsetting to needed SNPs
message(sprintf("   Loading %d SNPs from matrix...", length(snp_ids_needed)))
snp_mat <- fread(snp_file, header=TRUE)
snp_mat <- snp_mat[snpid %in% snp_ids_needed]
message(sprintf("   Found %d of %d requested SNPs", nrow(snp_mat), length(snp_ids_needed)))

message("## Loading expression matrix...")
expr_mat <- fread(expr_file, header=TRUE)
gene_ids_needed <- unique(pairs$geneid)
expr_mat <- expr_mat[gene_id %in% gene_ids_needed]
message(sprintf("   Found %d of %d requested genes", nrow(expr_mat), length(gene_ids_needed)))

# Ensure sample columns match between SNP and expression matrices
snp_samples  <- colnames(snp_mat)[-1]
expr_samples <- colnames(expr_mat)[-1]
# Expression uses HRA IDs (underscore), SNP uses HDA IDs — they differ
# but are in the same column ORDER from prep alignment, so we use position
if (length(snp_samples) != length(expr_samples)) {
    stop("SNP and expression matrices have different numbers of samples. Re-run prep_eQTL.sh.")
}
n_samples <- length(snp_samples)
message(sprintf("   %d samples in both matrices", n_samples))

##############################################
# Generate one PDF per gene-SNP pair
##############################################
message("## Generating boxplots...")

for (i in seq_len(nrow(pairs))) {
    G <- pairs$geneid[i]
    S <- pairs$snpid[i]

    # Get genotype vector (numeric 0/1/2 = dosage of ALT allele)
    snp_row <- snp_mat[snpid == S]
    if (nrow(snp_row) == 0) {
        message(sprintf("   SKIP: SNP %s not found", S))
        next
    }
    geno <- as.numeric(snp_row[, -1, with=FALSE])

    # Get expression vector
    expr_row <- expr_mat[gene_id == G]
    if (nrow(expr_row) == 0) {
        message(sprintf("   SKIP: Gene %s not found", G))
        next
    }
    expr <- as.numeric(expr_row[, -1, with=FALSE])

    if (length(geno) != length(expr)) {
        message(sprintf("   SKIP: %s x %s — mismatched lengths", G, S))
        next
    }

    # Allele info from SNP name (._AC -> ALT=AC)
    alt_allele <- extract_allele(S)

    # Chr:pos label from snp_location
    loc_row <- snploc[snpid == S]
    pos_label <- if (nrow(loc_row) > 0)
        paste0(loc_row$chr[1], ":", loc_row$pos[1])
    else S

    df <- data.frame(
        expression = expr,
        SNP        = geno,
        stringsAsFactors = FALSE
    )
    df <- df[!is.na(df$SNP) & !is.na(df$expression), ]

    if (nrow(df) < 10) {
        message(sprintf("   SKIP: %s x %s — too few non-NA samples (%d)", G, S, nrow(df)))
        next
    }

    # Genotype factor: 0=HomRef, 1=Het, 2=HomAlt
    df$SNP_factor <- factor(df$SNP, levels=2:0)
    # REF is inferred as "not ALT" — label only shows ALT allele
    geno_labels  <- c("Hom Alt", "Het", "Hom Ref")
    allele_labels <- c(
        paste0(alt_allele, "/", alt_allele),
        paste0("Ref/", alt_allele),
        "Ref/Ref"
    )
    n_per_group <- table(factor(df$SNP, levels=2:0))

    out_pdf <- file.path(out_dir,
        paste0("eQTLboxplot.", G, ".", gsub("[^A-Za-z0-9]","_", S), ".pdf"))

    tryCatch({
        pdf(out_pdf, width=7, height=5)
        par(mfrow=c(1,2), mar=c(5,4,3,1), oma=c(0,0,3,0))

        # Panel 1: Additive model (0/1/2)
        bp <- boxplot(expression ~ SNP_factor, data=df,
                      ylab="Normalized Expression",
                      xaxt='n', main="", col='lightgreen',
                      outline=FALSE, outpch=NA)
        stripchart(expression ~ SNP_factor, data=df,
                   vertical=TRUE, method="jitter",
                   pch=1, col="darkred", cex=0.7, add=TRUE)

        fit  <- lm(expression ~ SNP, data=df)
        pval <- signif(summary(fit)$coefficients[2,4], 3)
        title(main=paste0("Additive (p=", pval, ")"), cex.main=0.9, line=0.5)
        mtext(geno_labels,  side=1, line=0.5, at=1:3, cex=0.7)
        mtext(allele_labels, side=1, line=1.5, at=1:3, cex=0.7)
        mtext(paste0("N=", as.integer(n_per_group)), side=1, line=2.5, at=1:3, cex=0.7)

        # Panel 2: Dominant model (0 vs 1+2 = any ALT)
        df$has_alt <- ifelse(df$SNP > 0, 1, 0)
        if (length(unique(df$has_alt)) > 1 && min(table(df$has_alt)) > 1) {
            df$has_alt_f <- factor(df$has_alt, levels=0:1)
            p_dom <- signif(t.test(expression ~ has_alt_f, df)$p.value, 3)
            bp2 <- boxplot(expression ~ has_alt_f, data=df,
                           ylab="Normalized Expression",
                           xaxt='n', main="", col='lightblue',
                           outline=FALSE, outpch=NA)
            stripchart(expression ~ has_alt_f, data=df,
                       vertical=TRUE, method="jitter",
                       pch=1, col="darkred", cex=0.7, add=TRUE)
            title(main=paste0("Dominant (p=", p_dom, ")"), cex.main=0.9, line=0.5)
            mtext(c("Ref/Ref", paste0("Any ", alt_allele)),
                  side=1, line=0.5, at=1:2, cex=0.7)
            mtext(paste0("N=", as.integer(table(df$has_alt_f))),
                  side=1, line=1.5, at=1:2, cex=0.7)
        } else {
            plot.new()
            text(0.5, 0.5, "Dominant model\nnot applicable\n(only one genotype group)",
                 cex=0.8, col="grey40")
        }

        mtext(paste0("cis-eQTL: ", G, "\nSNP: ", pos_label),
              outer=TRUE, cex=0.9, font=2)
        dev.off()
    }, error=function(e) {
        message(sprintf("   ERROR on %s x %s: %s", G, S, e$message))
        try(dev.off(), silent=TRUE)
    })

    if (i %% 10 == 0)
        message(sprintf("   %d / %d plots done", i, nrow(pairs)))
}

message(sprintf("## Boxplots written to: %s", out_dir))
message("## Done.")
