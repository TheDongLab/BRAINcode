#!/usr/bin/env Rscript
###########################################
# _cQTL_postprocess_LINEAR_CROSS.R
# Post-processing script for interaction cQTLs
###########################################

library(data.table)

args <- commandArgs(TRUE)
if (length(args) < 4) {
    stop("Usage: Rscript _cQTL_postprocess_LINEAR_CROSS.R <cis_file> <snp_loc_file> <circ_loc_file> <out_prefix> [fdr_thresh] [top_n]")
}

cis_file      <- args[1]
snp_loc_file  <- args[2]
circ_loc_file <- args[3] # Replaced gene_loc_file
out_prefix    <- args[4]
fdr_thresh    <- ifelse(is.na(args[5]), 0.05, as.numeric(args[5]))
top_n_total   <- ifelse(is.na(args[6]), 1000000, as.integer(args[6]))

TELOMERE_DIST <- 500000 

message("## Reading cis-interaction-cQTL results...")
cqtl <- fread(cis_file, sep="\t", header=TRUE)
message("## Reading SNP location file...")
snploc <- fread(snp_loc_file, sep="\t", header=TRUE)
setnames(snploc, c("snpid","chr","pos"))
message("## Reading circRNA location file...")
circloc <- fread(circ_loc_file, sep="\t", header=TRUE)
setnames(circloc, c("circ_id","circ_chr","circ_start","circ_end"))

message("## Joining Coordinates...")
if("SNP" %in% names(cqtl)) setnames(cqtl, "SNP", "snpid")
if("gene" %in% names(cqtl)) setnames(cqtl, "gene", "circ_id")

cqtl <- merge(cqtl, snploc, by="snpid", all.x=TRUE)
cqtl <- merge(cqtl, circloc, by="circ_id", all.x=TRUE)

# Flag SNPs near telomeres using explicit assignment to prevent scoping issues
cqtl$telomeric_flag <- ifelse(!is.na(cqtl$pos) & cqtl$pos < TELOMERE_DIST, TRUE, FALSE)

# Filter for significant interaction effects (FDR-based)
cqtl_fdr <- cqtl[FDR < fdr_thresh]

# Fallback mechanism: If FDR filters out everything but raw p-values are strong
if (nrow(cqtl_fdr) == 0) {
  message("## No pairs passed strict FDR. Falling back to top hits by raw p-value...")
  
  # Sort by raw p-value and take pairs meeting a standard nominal threshold (e.g., p < 0.01)
  # Limit to top 250 pairs max so you don't generate thousands of junk plots
  cqtl_ordered <- cqtl[order(`p-value`)]
  cqtl_fdr <- head(cqtl_ordered[`p-value` < 0.01], 250)
}

# Lead interaction SNP selection (per circ_id)
cqtl_lead <- cqtl_fdr[order(circ_id, `p-value`, pos)]
cqtl_lead <- cqtl_lead[!duplicated(circ_id)]

message(sprintf("## Filtering for significant interaction hits (FDR < %s)...", fdr_thresh))
message(sprintf("## Total significant interaction pairs: %d", nrow(cqtl_fdr)))

# Subsetting and ordering for the boxplot generation step
cqtl_top <- cqtl_fdr[order(`p-value`), .(circ_id, snpid)]
cqtl_fdr_compact <- cqtl_fdr[, .(circ_id, snpid, beta, `t-stat`, `p-value`, FDR)]
message(sprintf("## Total interaction pairs to plot: %d", nrow(cqtl_top)))

# Save output files
fwrite(cqtl, file=paste0(out_prefix, ".full_annotated.txt"), sep="\t", quote=FALSE)
fwrite(cqtl_fdr_compact, file=paste0(out_prefix, ".FDR", fdr_thresh, ".txt"), sep="\t", quote=FALSE)
fwrite(cqtl_lead, file=paste0(out_prefix, ".lead_snps.txt"), sep="\t", quote=FALSE)
fwrite(cqtl_top, file=paste0(out_prefix, ".top_for_boxplot.txt"), sep="\t", quote=FALSE, col.names=FALSE)
message(sprintf("## Boxplot input file written with %d pairs", nrow(cqtl_top)))

########################################################################
# Meta-Analysis Logic (Adapted for Interaction Coefficients)
########################################################################
current_sex <- ifelse(grepl("Male", out_prefix), "Male", "Female")
other_sex   <- ifelse(current_sex == "Male", "Female", "Male")
other_sex_file <- gsub(paste0("/", current_sex, "/"), paste0("/", other_sex, "/"), paste0(out_prefix, ".full_annotated.txt"))

if (file.exists(other_sex_file)) {
    message("## Found opposite sex results. Running Interaction Meta-analysis...")
    
    combined_dir <- gsub(paste0("/", current_sex, "/results.*"), "/Combined/results", out_prefix)
    if(!dir.exists(combined_dir)) dir.create(combined_dir, showWarnings=FALSE, recursive=TRUE)
    
    cqtl_other <- fread(other_sex_file, sep="\t")
    
    meta <- merge(cqtl[, .(circ_id, snpid, chr, pos, beta_curr = beta, tstat_curr = `t-stat`, pval_curr = `p-value`)],
                  cqtl_other[, .(circ_id, snpid, beta_other = beta, tstat_other = `t-stat`, pval_other = `p-value`)],
                  by=c("circ_id", "snpid"))
    
    message(sprintf("## Performing Stouffer's Z-score interaction meta-analysis on %d shared pairs...", nrow(meta)))
    
    meta$z_curr  <- qnorm(meta$pval_curr / 2, lower.tail=FALSE) * sign(meta$tstat_curr)
    meta$z_other <- qnorm(meta$pval_other / 2, lower.tail=FALSE) * sign(meta$tstat_other)
    meta$z_meta  <- (meta$z_curr + meta$z_other) / sqrt(2)
    meta$p_meta  <- 2 * pnorm(abs(meta$z_meta), lower.tail=FALSE)
    meta$FDR_meta <- p.adjust(meta$p_meta, method="BH")
    
    tissue_name <- basename(dirname(dirname(dirname(dirname(out_prefix)))))
    
    meta_out <- file.path(combined_dir, paste0(tissue_name, "_Combined_Meta_Interaction_cQTL.txt"))
    fwrite(meta[order(p_meta)], file=meta_out, sep="\t", quote=FALSE)
    message(sprintf("## Meta-analysis saved to: %s", meta_out))
}

message("## Post-processing complete. Ready for plots.")
