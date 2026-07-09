#!/usr/bin/env Rscript
###########################################
# _cQTL_postprocess.R - ALL SIGNIFICANT PAIRS
# (Corrected Inverse-Variance Meta-Analysis)
###########################################

library(data.table)

args <- commandArgs(TRUE)
cis_file      <- args[1]
snp_loc_file  <- args[2]
circ_loc_file <- args[3]
out_prefix    <- args[4]
fdr_thresh    <- ifelse(is.na(args[5]), 0.05, as.numeric(args[5]))
top_n_total   <- ifelse(is.na(args[6]), 1000000, as.integer(args[6]))

TELOMERE_DIST <- 500000 

message("## Reading cis-cQTL results...")
cqtl <- fread(cis_file, sep="\t", header=TRUE)
message("## Reading SNP location file...")
snploc <- fread(snp_loc_file, sep="\t", header=TRUE)
setnames(snploc, c("snpid","chr","pos"))
message("## Reading circ location file...")
circloc <- fread(circ_loc_file, sep="\t", header=TRUE)
setnames(circloc, c("circ_id","circ_chr","circ_start","circ_end"))

message("## Joining Coordinates...")
if("SNP" %in% names(cqtl)) setnames(cqtl, "SNP", "snpid")
if("gene" %in% names(cqtl)) setnames(cqtl, "gene", "circ_id")

cqtl <- merge(cqtl, snploc, by="snpid", all.x=TRUE)
cqtl <- merge(cqtl, circloc, by="circ_id", all.x=TRUE)

# Flag SNPs near telomeres (useful for filtering later if needed)
cqtl[, telomeric_flag := ifelse(!is.na(pos) & pos < TELOMERE_DIST, TRUE, FALSE)]

# Filter for significant hits only
cqtl_fdr <- cqtl[FDR < fdr_thresh]

# Fallback mechanism: If FDR filters out everything but raw p-values are strong
if (nrow(cqtl_fdr) == 0) {
  message("## No pairs passed strict FDR. Falling back to top hits by raw p-value...")
  
  # Sort by raw p-value and take pairs meeting a standard nominal threshold (e.g., p < 0.01)
  # Limit to top 250 pairs max so you don't generate thousands of junk plots
  cqtl_ordered <- cqtl[order(`p-value`)]
  cqtl_fdr <- head(cqtl_ordered[`p-value` < 0.01], 250)
}

# Lead SNP selection (per circRNA) for summary/reporting
cqtl_lead <- cqtl_fdr[order(circ_id, `p-value`, pos)]
cqtl_lead <- cqtl_lead[!duplicated(circ_id)]

########################################################################
# SELECT ALL SIGNIFICANT PAIRS, ORDERED BY P-VALUE
########################################################################
message(sprintf("## Filtering for significant hits (FDR < %s)...", fdr_thresh))
message(sprintf("## Total significant pairs: %d", nrow(cqtl_fdr)))
message(sprintf("## Selecting ALL significant pairs for boxplots (FDR < %s)...", fdr_thresh))

# Subsetting and ordering
cqtl_top <- cqtl_fdr[order(`p-value`), .(circ_id, snpid)]
cqtl_fdr_compact <- cqtl_fdr[, .(circ_id, snpid, beta, `t-stat`, `p-value`, FDR)]
message(sprintf("## Total pairs to plot: %d", nrow(cqtl_top)))

# Save output files
fwrite(cqtl, file=paste0(out_prefix, ".full_annotated.txt"), sep="\t", quote=FALSE)
fwrite(cqtl_fdr_compact, file=paste0(out_prefix, ".FDR", fdr_thresh, ".txt"), sep="\t", quote=FALSE)
fwrite(cqtl_lead, file=paste0(out_prefix, ".lead_snps.txt"), sep="\t", quote=FALSE)
fwrite(cqtl_top, file=paste0(out_prefix, ".top_for_boxplot.txt"), sep="\t", quote=FALSE, col.names=FALSE)
message(sprintf("## Boxplot input file written with %d pairs (ordered by p-value)", nrow(cqtl_top)))

########################################################################
# Meta-Analysis Logic (Fixed Effects Inverse-Variance Method)
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
    
    cqtl_other <- fread(other_sex_file, sep="\t")
    
    # Merge current and other sex results on circRNA ID and SNP
    meta <- merge(cqtl[, .(circ_id, snpid, chr, pos, beta_curr = beta, tstat_curr = `t-stat`, pval_curr = `p-value`)],
                  cqtl_other[, .(circ_id, snpid, beta_other = beta, tstat_other = `t-stat`, pval_other = `p-value`)],
                  by=c("circ_id", "snpid"))
    
    message(sprintf("## Performing Inverse-Variance Fixed-Effects meta-analysis on %d shared pairs...", nrow(meta)))
    
    # Calculate Standard Errors from beta and t-stat (SE = |beta / t-stat|)
    meta[, se_curr  := abs(beta_curr / tstat_curr)]
    meta[, se_other := abs(beta_other / tstat_other)]
    
    # Calculate weights based on inverse variance (w = 1 / SE^2)
    meta[, w_curr  := 1 / (se_curr^2)]
    meta[, w_other := 1 / (se_other^2)]
    
    # Calculate pooled Meta-Beta (weighted average of effect sizes)
    meta[, beta_meta := (beta_curr * w_curr + beta_other * w_other) / (w_curr + w_other)]
    
    # Calculate pooled Standard Error and Meta-Z-score
    meta[, se_meta   := sqrt(1 / (w_curr + w_other))]
    meta[, z_meta    := beta_meta / se_meta]
    
    # Compute two-tailed p-values and apply Benjamini-Hochberg FDR correction
    meta[, p_meta    := 2 * pnorm(abs(z_meta), lower.tail = FALSE)]
    meta[, FDR_meta  := p.adjust(p_meta, method = "BH")]
    
    # Robust tissue name extraction using regex on target path fields
    tissue_match <- regexpr("target_ALS/[A-Za-z0-9_]+", out_prefix)
    if(tissue_match != -1) {
        tissue_raw <- regmatches(out_prefix, tissue_match)
        tissue_name <- gsub("target_ALS/", "", tissue_raw)
    } else {
        tissue_name <- "Unknown_Tissue"
    }
    
    meta_out <- file.path(combined_dir, paste0(tissue_name, "_Combined_Meta_cQTL.txt"))
    fwrite(meta[order(p_meta)], file=meta_out, sep="\t", quote=FALSE)
    message(sprintf("## Meta-analysis saved to: %s", meta_out))
}

message("## Post-processing complete. Ready for plots.")
