#!/usr/bin/env Rscript
###########################################
# _cQTL_boxplot.R
# Lightweight version adapted for circular RNA percentages
###########################################

suppressPackageStartupMessages({
  library(data.table)
})

args <- commandArgs(TRUE)
if (length(args) < 6) {
  stop("Usage: Rscript _cQTL_boxplot.R <pairs_file> <snp_file> <circ_file> <cov_file> <out_dir> <tissue>")
}

GSfile    <- args[1]
snp_file  <- args[2]
circ_file <- args[3]
cov_file  <- args[4]
out_dir   <- args[5]
tissue    <- args[6]

out_file_std <- file.path(out_dir, paste0(tissue, "_all_sig_cQTL_boxplots.pdf"))
out_file_col <- file.path(out_dir, paste0(tissue, "_all_sig_cQTL_boxplots_colored.pdf"))

# Read the pairs structure (expecting header containing circ_id and snpid)
pairs <- fread(GSfile, header=FALSE)
if (ncol(pairs) >= 2) {
  setnames(pairs, c("circ_id", "snpid"))
} else {
  stop("Error: Pairs file must contain at least two columns (circ_id and snpid).")
}

message("# Loading genomic matrices...")
snp_mat  <- fread(snp_file, header=TRUE)
circ_mat <- fread(circ_file, header=TRUE)
cov_mat  <- fread(cov_file, header=TRUE)

# Standardize internal matrix headers
setnames(circ_mat, names(circ_mat)[1], "circ_id")
setnames(snp_mat, names(snp_mat)[1], "snpid")

# Process covariates
cov_first_col <- names(cov_mat)[1]
als_row <- cov_mat[get(cov_first_col) == "is_als", ]
if (nrow(als_row) == 0) {
  als_vals <- rep(0, ncol(cov_mat) - 1)
} else {
  als_vals <- as.numeric(unlist(als_row[, -1, with=FALSE], use.names=FALSE))
}

p_colors <- ifelse(als_vals > 0.5, "red", "#9932CC")
p_shapes <- ifelse(als_vals > 0.5, 16, 18)
p_sizes  <- ifelse(als_vals > 0.5, 0.7, 1.1)

run_plotting <- function(pdf_path, use_status_colors = FALSE) {
    pdf(pdf_path, width=12, height=5)
    plots_made <- 0
    pairs_skipped <- 0
    
    for (i in seq_len(nrow(pairs))) {
        c_id <- pairs$circ_id[i]  # Used to pull data from matrix and title plot
        S    <- pairs$snpid[i]
        
        # Pull rows using the clean circRNA ID match
        circ_row <- circ_mat[circ_id == c_id]
        snp_row  <- snp_mat[snpid == S]
        
        if (nrow(circ_row) == 0 || nrow(snp_row) == 0) {
            pairs_skipped <- pairs_skipped + 1
            next
        }
        
        circ_vals <- as.numeric(unlist(circ_row[, -1, with=FALSE], use.names=FALSE))
        snp_vals  <- as.numeric(unlist(snp_row[, -1, with=FALSE], use.names=FALSE))
        
        if (length(circ_vals) == 0 || length(snp_vals) == 0 || length(circ_vals) != length(snp_vals)) {
            pairs_skipped <- pairs_skipped + 1
            next
        }
        
        df <- data.frame(
            circ_percent = circ_vals,
            SNP = snp_vals,
            p_col = if(use_status_colors) p_colors else rep("darkred", length(circ_vals)),
            p_pch = if(use_status_colors) p_shapes else rep(16, length(circ_vals)),
            p_cex = if(use_status_colors) p_sizes else rep(0.7, length(circ_vals)),
            stringsAsFactors = FALSE
        )
        
        df <- df[!is.na(df$SNP) & !is.na(df$circ_percent), ]
        if (nrow(df) < 5) {
            pairs_skipped <- pairs_skipped + 1
            next
        }

        get_p <- function(formula, data) {
            tryCatch({
                fit <- lm(formula, data=data)
                formatC(summary(fit)$coefficients[2,4], format="e", digits=3)
            }, error = function(e) "NA")
        }

        par(mfrow=c(1,3), mar=c(5,4,4,2), oma=c(0,0,3,0))
        
        # Additive Model
        add_counts <- table(factor(df$SNP, levels=0:2))
        add_labels <- paste0(c("Ref/Ref", "Het", "Hom Alt"), "\n(N=", add_counts, ")")
        df$SNP_f <- factor(df$SNP, levels=0:2, labels=add_labels)
        boxplot(circ_percent ~ SNP_f, data=df, col="lightgreen", outline=FALSE, 
                ylab="circRNA Percentage (%)", xlab="Genotype",
                main=paste0("Additive Model\n(p = ", get_p(circ_percent ~ SNP, df), ")"))
        points(jitter(as.numeric(df$SNP_f), amount=0.15), df$circ_percent, pch=df$p_pch, col=df$p_col, cex=df$p_cex)
        if (use_status_colors) legend("topleft", legend = c("ALS", "Control"), col = c("red", "#9932CC"), pch = c(16, 18), pt.cex = c(0.7, 1.1), bty = "n", cex = 0.8)
        
        # Dominant Model
        df$dom_val <- ifelse(df$SNP > 0, 1, 0)
        dom_counts <- table(factor(df$dom_val, levels=0:1))
        dom_labels <- paste0(c("Ref/Ref", "Any Alt"), "\n(N=", dom_counts, ")")
        df$dom_f <- factor(df$dom_val, levels=0:1, labels=dom_labels)
        boxplot(circ_percent ~ dom_f, data=df, col="lightblue", outline=FALSE, 
                ylab="circRNA Percentage (%)", xlab="Genotype",
                main=paste0("Dominant Model\n(p = ", get_p(circ_percent ~ dom_val, df), ")"))
        points(jitter(as.numeric(df$dom_f), amount=0.15), df$circ_percent, pch=df$p_pch, col=df$p_col, cex=df$p_cex)
        if (use_status_colors) legend("topleft", legend = c("ALS", "Control"), col = c("red", "#9932CC"), pch = c(16, 18), pt.cex = c(0.7, 1.1), bty = "n", cex = 0.8)

        # Recessive Model
        df$rec_val <- ifelse(df$SNP == 2, 1, 0)
        rec_counts <- table(factor(df$rec_val, levels=0:1))
        rec_labels <- paste0(c("Non-Hom Alt", "Hom Alt"), "\n(N=", rec_counts, ")")
        df$rec_f <- factor(df$rec_val, levels=0:1, labels=rec_labels)
        boxplot(circ_percent ~ rec_f, data=df, col="lightpink", outline=FALSE, 
                ylab="circRNA Percentage (%)", xlab="Genotype",
                main=paste0("Recessive Model\n(p = ", get_p(circ_percent ~ rec_val, df), ")"))
        points(jitter(as.numeric(df$rec_f), amount=0.15), df$circ_percent, pch=df$p_pch, col=df$p_col, cex=df$p_cex)
        if (use_status_colors) legend("topleft", legend = c("ALS", "Control"), col = c("red", "#9932CC"), pch = c(16, 18), pt.cex = c(0.7, 1.1), bty = "n", cex = 0.8)
        
        # Label the top of the page cleanly using the unique circRNA identifier
        mtext(paste("circRNA:", c_id, "| SNP:", S), outer=TRUE, cex=1.0, font=2, line=0.5)
        
        plots_made <- plots_made + 1
    }
    dev.off()
    message(sprintf("# Plots made: %d, pairs skipped: %d", plots_made, pairs_skipped))
}

message(paste("# Processing", nrow(pairs), "circ-SNP pairs..."))
run_plotting(out_file_std, FALSE)
run_plotting(out_file_col, TRUE)
message("Done.")
