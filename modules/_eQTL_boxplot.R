#!/usr/bin/env Rscript
###########################################
# _eQTL_boxplot.R
# FIXED: Shape & Color differentiation for rare controls
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

# Logic: Red/Circle (16) for ALS, Violet/Diamond (18) for Control
# Using a slightly larger cex (size) for the diamonds to help them pop
p_colors <- ifelse(als_vals > 0.5, "red", "darkorchid")
p_shapes <- ifelse(als_vals > 0.5, 16, 18)
p_sizes  <- ifelse(als_vals > 0.5, 0.7, 1.2) 

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
        
        get_stats <- function(model_formula, data) {
            tryCatch({
                fit <- lm(model_formula, data=data); s <- summary(fit)$coefficients
                if(nrow(s) < 2) return(list(p="NA"))
                list(p=signif(s[2,4], 3))
            }, error = function(e) list(p="Err"))
        }

        par(mfrow=c(1,3), mar=c(5,4,4,2), oma=c(0,0,3,0))
        
        plot_panel <- function(formula, data, title_prefix, box_color) {
            boxplot(formula, data=data, col=box_color, outline=F, 
                    ylab="Expression (Z-score)", xlab="Genotype",
                    main=paste0(title_prefix, "\n(p = ", get_stats(formula, data)$p, ")"))
            
            # Using point-specific vectors for col, pch, and cex
            stripchart(formula, data=data, vertical=T, method="jitter", add=T, 
                       pch=data$p_pch, col=data$p_col, cex=data$p_cex)
            
            if(use_status_colors && title_prefix == "Additive Model") {
                legend("topleft", legend=c("ALS (Circle)", "Control (Diamond)"), 
                       col=c("red", "darkorchid"), pch=c(16, 18), 
                       pt.cex=c(1, 1.2), cex=0.8, bty="n")
            }
        }

        # Additive
        add_labels <- paste0(c("Ref/Ref", "Het", "Hom Alt"), "\n(N=", table(factor(df$SNP, levels=0:2)), ")")
        df$SNP_f <- factor(df$SNP, levels=0:2, labels=add_labels)
        plot_panel(expression ~ SNP_f, df, "Additive Model", "lightgreen")

        # Dominant
        df$dom_val <- ifelse(df$SNP > 0, 1, 0)
        dom_labels <- paste0(c("Ref/Ref", "Any Alt"), "\n(N=", table(factor(df$dom_val, levels=0:1)), ")")
        df$dom <- factor(df$dom_val, levels=0:1, labels=dom_labels)
        plot_panel(expression ~ dom, df, "Dominant Model", "lightblue")

        # Recessive
        df$rec_val <- ifelse(df$SNP == 2, 1, 0)
        rec_labels <- paste0(c("Non-Hom Alt", "Hom Alt"), "\n(N=", table(factor(df$rec_val, levels=0:1)), ")")
        df$rec <- factor(df$rec_val, levels=0:1, labels=rec_labels)
        plot_panel(expression ~ rec, df, "Recessive Model", "lightpink")

        mtext(paste("Gene:", G, "| SNP:", S), outer=TRUE, cex=1.2, font=2, line=0.5)
    }
    dev.off()
}

run_plotting(out_file_std, FALSE)
run_plotting(out_file_col, TRUE)
message("Done.")
