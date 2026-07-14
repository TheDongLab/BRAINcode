#!/usr/bin/env Rscript
###########################################
# _cQTL_boxplot_LINEAR_CROSS.R
# Interaction Layout: Split Boxplots + Regression Slopes (cQTL Edition)
###########################################

suppressPackageStartupMessages({
  library(data.table)
})

args <- commandArgs(TRUE)
if (length(args) < 7) {
  stop("Usage: Rscript _cQTL_boxplot_LINEAR_CROSS.R <pairs_file> <snp_file> <expr_file> <cov_file> <snp_loc> <out_dir> <tissue>")
}

GSfile    <- args[1]
snp_file  <- args[2]
expr_file <- args[3]
cov_file  <- args[4]
snp_loc   <- args[5]
out_dir   <- args[6]
tissue    <- args[7]

out_file_std <- file.path(out_dir, paste0(tissue, "_all_sig_Interaction_boxplots.pdf"))
out_file_col <- file.path(out_dir, paste0(tissue, "_all_sig_Interaction_boxplots_colored.pdf"))

# Read file structure (expecting headers/columns: circ_id and snpid)
pairs <- fread(GSfile, header=TRUE)

message("# Loading genomic matrices...")
snp_mat  <- fread(snp_file, header=TRUE)
expr_mat <- fread(expr_file, header=TRUE)
cov_mat  <- fread(cov_file, header=TRUE)

# Standardize headers for circular RNAs and SNPs
if("geneid" %in% names(pairs)) setnames(pairs, "geneid", "circ_id")
if("gene" %in% names(pairs)) setnames(pairs, "gene", "circ_id")
if("SNP" %in% names(pairs)) setnames(pairs, "SNP", "snpid")

setnames(expr_mat, names(expr_mat)[1], "circ_id")
setnames(snp_mat, names(snp_mat)[1], "snpid")

# Process the clinical interaction covariate status
cov_first_col <- names(cov_mat)[1]
als_row <- cov_mat[get(cov_first_col) == "is_als", ]
if (nrow(als_row) == 0) {
  als_vals <- rep(0, ncol(cov_mat) - 1)
} else {
  als_vals <- as.numeric(unlist(als_row[, -1, with=FALSE], use.names=FALSE))
}

# Jitter plotting definitions
p_colors <- ifelse(als_vals > 0.5, "red", "#9932CC")
p_shapes <- ifelse(als_vals > 0.5, 16, 18)
p_sizes  <- ifelse(als_vals > 0.5, 0.7, 1.1)

run_plotting <- function(pdf_path, use_status_colors = FALSE) {
    pdf(pdf_path, width=14, height=5.5)
    plots_made <- 0
    pairs_skipped <- 0
    
    for (i in seq_len(nrow(pairs))) {
        c_id   <- pairs$circ_id[i]
        S      <- pairs$snpid[i]
        
        expr_row <- expr_mat[circ_id == c_id]
        snp_row  <- snp_mat[snpid == S]
        
        if (nrow(expr_row) == 0 || nrow(snp_row) == 0) {
            pairs_skipped <- pairs_skipped + 1
            next
        }
        
        expr_vals <- as.numeric(unlist(expr_row[, -1, with=FALSE], use.names=FALSE))
        snp_vals  <- as.numeric(unlist(snp_row[, -1, with=FALSE], use.names=FALSE))
        
        if (length(expr_vals) == 0 || length(snp_vals) == 0 || length(expr_vals) != length(snp_vals)) {
            pairs_skipped <- pairs_skipped + 1
            next
        }
        
        df <- data.frame(
            expression = expr_vals,
            SNP = snp_vals,
            Status = ifelse(als_vals > 0.5, "ALS", "Control"),
            p_col = if(use_status_colors) p_colors else rep("darkred", length(expr_vals)),
            p_pch = if(use_status_colors) p_shapes else rep(16, length(expr_vals)),
            p_cex = if(use_status_colors) p_sizes else rep(0.7, length(expr_vals)),
            stringsAsFactors = FALSE
        )
        
        df <- df[!is.na(df$SNP) & !is.na(df$expression), ]
        if (nrow(df[df$Status == "ALS", ]) < 3 || nrow(df[df$Status == "Control", ]) < 3) {
            pairs_skipped <- pairs_skipped + 1
            next
        }

        # Pull p-value helper specifically for the Interaction term (Genotype:Status)
        get_interaction_p <- function(data) {
            tryCatch({
                fit <- lm(expression ~ SNP * factor(Status), data=data)
                p_val <- summary(fit)$coefficients["SNP:factor(Status)Control", 4]
                formatC(p_val, format="e", digits=3)
            }, error = function(e) "NA")
        }

        # 3-Panel Layout: Controls, ALS, and Overlay Line Trend
        par(mfrow=c(1,3), mar=c(5,4,4,1), oma=c(0,0,3,0))
        
        # 1. Controls Boxplot
        df_ctrl <- df[df$Status == "Control", ]
        ctrl_counts <- table(factor(df_ctrl$SNP, levels=0:2))
        ctrl_labels <- paste0(c("Ref/Ref", "Het", "Hom Alt"), "\n(N=", ctrl_counts, ")")
        df_ctrl$SNP_f <- factor(df_ctrl$SNP, levels=0:2, labels=ctrl_labels)
        
        boxplot(expression ~ SNP_f, data=df_ctrl, col="lightblue", outline=FALSE, 
                ylab="Expression (Z-score)", xlab="Genotype", main="Controls Only")
        points(jitter(as.numeric(df_ctrl$SNP_f), amount=0.15), df_ctrl$expression, 
               pch=df_ctrl$p_pch, col=df_ctrl$p_col, cex=df_ctrl$p_cex)
        
        # 2. ALS Boxplot
        df_als <- df[df$Status == "ALS", ]
        als_counts <- table(factor(df_als$SNP, levels=0:2))
        als_labels <- paste0(c("Ref/Ref", "Het", "Hom Alt"), "\n(N=", als_counts, ")")
        df_als$SNP_f <- factor(df_als$SNP, levels=0:2, labels=als_labels)
        
        boxplot(expression ~ SNP_f, data=df_als, col="lightpink", outline=FALSE, 
                ylab="Expression (Z-score)", xlab="Genotype", main="ALS Only")
        points(jitter(as.numeric(df_als$SNP_f), amount=0.15), df_als$expression, 
               pch=df_als$p_pch, col=df_als$p_col, cex=df_als$p_cex)
        
        # 3. Interaction Slope Visualization
        plot(df$SNP, df$expression, type="n", xaxt="n", xlab="Genotype", 
             ylab="Expression (Z-score)", main="Interaction Slope Context",
             sub=paste("Interaction p =", get_interaction_p(df)), col.sub="darkblue", font.sub=2)
        axis(1, at=0:2, labels=c("Ref/Ref", "Het", "Hom Alt"))
        
        # Add Jitter points for total breakdown visual
        points(jitter(df$SNP, amount=0.08), df$expression, pch=df$p_pch, col=df$p_col, cex=df$p_cex)
        
        # Add linear fit slopes per group
        tryCatch({
            abline(lm(expression ~ SNP, data=df_ctrl), col="#9932CC", lwd=3)
            abline(lm(expression ~ SNP, data=df_als), col="red", lwd=3)
        }, error = function(e) NULL)
        
        legend("topright", legend = c("ALS Slope", "Control Slope"), 
               col = c("red", "#9932CC"), lty=1, lwd=3, bty = "n", cex = 0.8)
        
        # Super-title label
        mtext(paste("circRNA:", c_id, " | SNP:", S), outer=TRUE, cex=1.1, font=2, line=0.5)
        
        plots_made <- plots_made + 1
    }
    dev.off()
    message(sprintf("# Interaction plots generated: %d, skipped: %d", plots_made, pairs_skipped))
}

message(paste("# Processing", nrow(pairs), "interaction circRNA-SNP pairs..."))
run_plotting(out_file_std, FALSE)
run_plotting(out_file_col, TRUE)
message("Done.")
