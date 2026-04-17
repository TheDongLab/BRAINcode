#!/usr/bin/env Rscript
###########################################
# _eQTL_boxplot.R - FINAL STABLE VERSION
###########################################

suppressPackageStartupMessages({
  library(data.table)
})

args <- commandArgs(TRUE)
GSfile    <- args[1]
snp_file  <- args[2]
expr_file <- args[3]
cov_file  <- args[4]
out_dir   <- args[6]
tissue    <- args[7]

out_file_std <- file.path(out_dir, paste0(tissue, "_top_eQTL_boxplots.pdf"))
out_file_col <- file.path(out_dir, paste0(tissue, "_top_eQTL_boxplots_colored.pdf"))

message("# Loading genomic matrices...")
pairs <- fread(GSfile, header=FALSE, col.names=c("geneid","snpid"))
snp_mat <- fread(snp_file, header=TRUE); setkey(snp_mat, snpid)
expr_mat <- fread(expr_file, header=TRUE); setkey(expr_mat, geneid)
cov_mat <- fread(cov_file, header=TRUE)

# --- FAILSAFE ALS STATUS EXTRACTION ---
als_idx <- which(cov_mat[[1]] == "is_als")
als_vals <- as.numeric(as.vector(cov_mat[als_idx, -1, with=FALSE]))

# Red (16) for ALS, Violet Hex (18) for Control
p_colors <- ifelse(als_vals > 0.5, "red", "#9932CC")
p_shapes <- ifelse(als_vals > 0.5, 16, 18)
p_sizes  <- ifelse(als_vals > 0.5, 0.7, 1.1) 

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
            p_col = if(use_status_colors) p_colors else rep("darkred", length(als_vals)),
            p_pch = if(use_status_colors) p_shapes else rep(16, length(als_vals)),
            p_cex = if(use_status_colors) p_sizes else rep(0.7, length(als_vals))
        )
        df <- df[!is.na(df$SNP), ]
        
        get_p <- function(formula, data) {
            tryCatch({
                fit <- lm(formula, data=data)
                signif(summary(fit)$coefficients[2,4], 3)
            }, error = function(e) "NA")
        }

        par(mfrow=c(1,3), mar=c(5,4,4,2), oma=c(0,0,3,0))
        
        # 1. ADDITIVE
        add_counts <- table(factor(df$SNP, levels=0:2))
        add_labels <- paste0(c("Ref/Ref", "Het", "Hom Alt"), "\n(N=", add_counts, ")")
        df$SNP_f <- factor(df$SNP, levels=0:2, labels=add_labels)
        
        boxplot(expression ~ SNP_f, data=df, col="lightgreen", outline=F, 
                ylab="Expression (Z-score)", xlab="Genotype",
                main=paste0("Additive Model\n(p = ", get_p(expression ~ SNP, df), ")"))
        points(jitter(as.numeric(df$SNP_f), amount=0.15), df$expression, 
               pch=df$p_pch, col=df$p_col, cex=df$p_cex)
        
        if(use_status_colors) {
            legend("topleft", legend=c("ALS", "Control"), 
                   col=c("red", "#9932CC"), pch=c(16, 18), 
                   pt.cex=c(1, 1.1), cex=0.8, bty="n")
        }

        # 2. DOMINANT
        df$dom_val <- ifelse(df$SNP > 0, 1, 0)
        dom_counts <- table(factor(df$dom_val, levels=0:1))
        dom_labels <- paste0(c("Ref/Ref", "Any Alt"), "\n(N=", dom_counts, ")")
        df$dom_f <- factor(df$dom_val, levels=0:1, labels=dom_labels)
        
        boxplot(expression ~ dom_f, data=df, col="lightblue", outline=F, 
                ylab="Expression (Z-score)", xlab="Genotype",
                main=paste0("Dominant Model\n(p = ", get_p(expression ~ dom_val, df), ")"))
        points(jitter(as.numeric(df$dom_f), amount=0.15), df$expression, 
               pch=df$p_pch, col=df$p_col, cex=df$p_cex)

        # 3. RECESSIVE
        df$rec_val <- ifelse(df$SNP == 2, 1, 0)
        rec_counts <- table(factor(df$rec_val, levels=0:1))
        rec_labels <- paste0(c("Non-Hom Alt", "Hom Alt"), "\n(N=", rec_counts, ")")
        df$rec_f <- factor(df$rec_val, levels=0:1, labels=rec_labels)
        
        boxplot(expression ~ rec_f, data=df, col="lightpink", outline=F, 
                ylab="Expression (Z-score)", xlab="Genotype",
                main=paste0("Recessive Model\n(p = ", get_p(expression ~ rec_val, df), ")"))
        points(jitter(as.numeric(df$rec_f), amount=0.15), df$expression, 
               pch=df$p_pch, col=df$p_col, cex=df$p_cex)

        mtext(paste("Gene:", G, "| SNP:", S), outer=TRUE, cex=1.2, font=2, line=0.5)
    }
    dev.off()
}

run_plotting(out_file_std, FALSE)
run_plotting(out_file_col, TRUE)
message("Done.")
