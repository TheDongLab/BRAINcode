#!/usr/bin/env Rscript
###########################################
# _eQTL_postprocess.R - ALL SIGNIFICANT PAIRS
###########################################

library(data.table)

args <- commandArgs(TRUE)
cis_file      <- args[1]
snp_loc_file  <- args[2]
gene_loc_file <- args[3]
out_prefix    <- args[4]
fdr_thresh    <- ifelse(is.na(args[5]), 0.05, as.numeric(args[5]))
top_n_total   <- ifelse(is.na(args[6]), 1000000, as.integer(args[6]))

TELOMERE_DIST <- 500000 

message("## Reading cis-eQTL results...")
eqtl <- fread(cis_file, sep="\t", header=TRUE)
message("## Reading SNP location file...")
snploc <- fread(snp_loc_file, sep="\t", header=TRUE)
setnames(snploc, c("snpid","chr","pos"))
message("## Reading gene location file...")
geneloc <- fread(gene_loc_file, sep="\t", header=TRUE)
setnames(geneloc, c("geneid","gene_chr","gene_start","gene_end"))

message("## Joining Coordinates...")
if("SNP" %in% names(eqtl)) setnames(eqtl, "SNP", "snpid")
if("gene" %in% names(eqtl)) setnames(eqtl, "gene", "geneid")

eqtl <- merge(eqtl, snploc, by="snpid", all.x=TRUE)
eqtl <- merge(eqtl, geneloc, by="geneid", all.x=TRUE)

# Flag SNPs near telomeres (useful for filtering later if needed)
eqtl[, telomeric_flag := ifelse(!is.na(pos) & pos < TELOMERE_DIST, TRUE, FALSE)]

# Filter for significant hits only
eqtl_fdr <- eqtl[FDR < fdr_thresh]

# Lead SNP selection (per gene) for summary/reporting
eqtl_lead <- eqtl_fdr[order(geneid, `p-value`, pos)]
eqtl_lead <- eqtl_lead[!duplicated(geneid)]

########################################################################
# SELECT ALL SIGNIFICANT PAIRS, ORDERED BY P-VALUE
########################################################################
message(sprintf("## Filtering for significant hits (FDR < %s)...", fdr_thresh))
message(sprintf("## Total significant pairs: %d", nrow(eqtl_fdr)))

message(sprintf("## Selecting ALL significant pairs for boxplots (FDR < %s)...", fdr_thresh))

# Take all significant pairs, ordered by p-value (most to least significant)
eqtl_top <- eqtl_fdr[order(`p-value`), .(geneid, snpid)]

message(sprintf("## Total pairs to plot: %d", nrow(eqtl_top)))

# Save standard output files
fwrite(eqtl, file=paste0(out_prefix, ".full_annotated.txt"), sep="\t", quote=FALSE)
fwrite(eqtl_fdr, file=paste0(out_prefix, ".FDR", fdr_thresh, ".txt"), sep="\t", quote=FALSE)
fwrite(eqtl_lead, file=paste0(out_prefix, ".lead_snps.txt"), sep="\t", quote=FALSE)

# This file contains all significant pairs ordered by p-value, no header
fwrite(eqtl_top, file=paste0(out_prefix, ".top_for_boxplot.txt"), 
       sep="\t", quote=FALSE, col.names=FALSE)

message(sprintf("## Boxplot input file written with %d pairs (ordered by p-value)", nrow(eqtl_top)))

########################################################################
# Meta-Analysis Logic
########################################################################
# Detect if we are in a Sex-stratified folder
current_sex <- ifelse(grepl("Male", out_prefix), "Male", "Female")
other_sex   <- ifelse(current_sex == "Male", "Female", "Male")
other_sex_file <- gsub(paste0("/", current_sex, "/"), paste0("/", other_sex, "/"), paste0(out_prefix, ".full_annotated.txt"))

if (file.exists(other_sex_file)) {
    message("## Found opposite sex results. Running Meta-analysis...")
    
    # Path logic to find the 'Combined' output directory
    combined_dir <- gsub(paste0("/", current_sex, "/results.*"), "/Combined/results", out_prefix)
    if(!dir.exists(combined_dir)) dir.create(combined_dir, showWarnings=FALSE, recursive=TRUE)
    
    eqtl_other <- fread(other_sex_file, sep="\t")
    
    # Merge current and other sex results on Gene and SNP
    meta <- merge(eqtl[, .(geneid, snpid, chr, pos, beta_curr = beta, tstat_curr = `t-stat`, pval_curr = `p-value`)],
                  eqtl_other[, .(geneid, snpid, beta_other = beta, tstat_other = `t-stat`, pval_other = `p-value`)],
                  by=c("geneid", "snpid"))
    
    message(sprintf("## Performing Stouffer's Z-score meta-analysis on %d shared pairs...", nrow(meta)))
    
    # Stouffer's Z-score Method (assuming equal weights/sample sizes for simplicity)
    meta[, z_curr  := qnorm(pval_curr / 2, lower.tail=FALSE) * sign(tstat_curr)]
    meta[, z_other := qnorm(pval_other / 2, lower.tail=FALSE) * sign(tstat_other)]
    meta[, z_meta  := (z_curr + z_other) / sqrt(2)]
    meta[, p_meta  := 2 * pnorm(abs(z_meta), lower.tail=FALSE)]
    meta[, FDR_meta := p.adjust(p_meta, method="BH")]
    
    # Tissue name extraction based on your folder structure
    # Extracts the tissue name (e.g., Cerebellum) from the path
    tissue_name <- basename(dirname(dirname(dirname(dirname(out_prefix)))))
    
    meta_out <- file.path(combined_dir, paste0(tissue_name, "_Combined_Meta_eQTL.txt"))
    fwrite(meta[order(p_meta)], file=meta_out, sep="\t", quote=FALSE)
    message(sprintf("## Meta-analysis saved to: %s", meta_out))
}

message("## Post-processing complete. Ready for plots.")
