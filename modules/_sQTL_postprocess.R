#!/usr/bin/env Rscript
###########################################
# _sQTL_postprocess.R - ALL SIGNIFICANT PAIRS
# Adapted from eQTL workflow v1.2
###########################################

library(data.table)

args <- commandArgs(TRUE)
cis_file      <- args[1]
snp_loc_file  <- args[2]
gene_loc_file <- args[3] # Splicing location map
out_prefix    <- args[4]
fdr_thresh    <- ifelse(is.na(args[5]), 0.05, as.numeric(args[5]))
top_n_total   <- ifelse(is.na(args[6]), 1000000, as.integer(args[6]))

TELOMERE_DIST <- 500000 

message("## Reading cis-sQTL results...")
sqtl <- fread(cis_file, sep="\t", header=TRUE)

message("## Reading SNP location file...")
snploc <- fread(snp_loc_file, sep="\t", header=TRUE)
setnames(snploc, c("snpid","chr","pos"))

message("## Reading splicing location file...")
geneloc <- fread(gene_loc_file, sep="\t", header=TRUE)
setnames(geneloc, c("junction_id","junction_chr","junction_start","junction_end"))

# Standardize feature IDs across data frames
if("SNP" %in% names(sqtl)) setnames(sqtl, "SNP", "snpid")
if("gene" %in% names(sqtl)) setnames(sqtl, "gene", "junction_id")

message("## Joining Coordinates...")
sqtl <- merge(sqtl, snploc, by="snpid", all.x=TRUE)
sqtl <- merge(sqtl, geneloc, by="junction_id", all.x=TRUE)

# Flag SNPs near telomeres
sqtl[, telomeric_flag := ifelse(!is.na(pos) & pos < TELOMERE_DIST, TRUE, FALSE)]

# Filter for significant hits only
sqtl_fdr <- sqtl[FDR < fdr_thresh]

# Lead SNP selection (per unique intron junction) for summary reporting
sqtl_lead <- sqtl_fdr[order(junction_id, `p-value`, pos)]
sqtl_lead <- sqtl_lead[!duplicated(junction_id)]

########################################################################
# SELECT ALL SIGNIFICANT PAIRS, ORDERED BY P-VALUE
########################################################################
message(sprintf("## Filtering for significant hits (FDR < %s)...", fdr_thresh))
message(sprintf("## Total significant pairs: %d", nrow(sqtl_fdr)))
message(sprintf("## Selecting ALL significant pairs for boxplots (FDR < %s)...", fdr_thresh))

# Take all significant pairs, ordered by p-value (most to least significant)
sqtl_top <- sqtl_fdr[order(`p-value`), .(junction_id, snpid)]

message(sprintf("## Total pairs to plot: %d", nrow(sqtl_top)))

# Save standard output files
fwrite(sqtl, file=paste0(out_prefix, ".full_annotated.txt"), sep="\t", quote=FALSE)
fwrite(sqtl_fdr, file=paste0(out_prefix, ".FDR", fdr_thresh, ".txt"), sep="\t", quote=FALSE)
fwrite(sqtl_lead, file=paste0(out_prefix, ".lead_snps.txt"), sep="\t", quote=FALSE)

# Boxplot driver map file (contains all significant pairs ordered by p-value, no header)
fwrite(sqtl_top, file=paste0(out_prefix, ".top_for_boxplot.txt"), 
       sep="\t", quote=FALSE, col.names=FALSE)

message(sprintf("## Boxplot input file written with %d pairs (ordered by p-value)", nrow(sqtl_top)))

########################################################################
# Meta-Analysis Logic
########################################################################
# Detect if we are running a Sex-stratified fold look
current_sex <- ifelse(grepl("Male", out_prefix), "Male", "Female")
other_sex   <- ifelse(current_sex == "Male", "Female", "Male")
other_sex_file <- gsub(paste0("/", current_sex, "/"), paste0("/", other_sex, "/"), paste0(out_prefix, ".full_annotated.txt"))

if (file.exists(other_sex_file)) {
    message("## Found opposite sex results. Running Meta-analysis...")
    
    combined_dir <- gsub(paste0("/", current_sex, "/results.*"), "/Combined/results", out_prefix)
    if(!dir.exists(combined_dir)) dir.create(combined_dir, showWarnings=FALSE, recursive=TRUE)
    
    sqtl_other <- fread(other_sex_file, sep="\t")
    
    # Standardize other file columns if necessary
    if("gene" %in% names(sqtl_other)) setnames(sqtl_other, "gene", "junction_id")
    if("SNP" %in% names(sqtl_other)) setnames(sqtl_other, "SNP", "snpid")
    
    # Merge current and other sex results on Junction and SNP
    meta <- merge(sqtl[, .(junction_id, snpid, chr, pos, beta_curr = beta, tstat_curr = `t-stat`, pval_curr = `p-value`)],
                  sqtl_other[, .(junction_id, snpid, beta_other = beta, tstat_other = `t-stat`, pval_other = `p-value`)],
                  by=c("junction_id", "snpid"))
    
    message(sprintf("## Performing Stouffer's Z-score meta-analysis on %d shared pairs...", nrow(meta)))
    
    # Stouffer's Z-score Method (assuming equal weights/sample sizes)
    meta[, z_curr  := qnorm(pval_curr / 2, lower.tail=FALSE) * sign(tstat_curr)]
    meta[, z_other := qnorm(pval_other / 2, lower.tail=FALSE) * sign(tstat_other)]
    meta[, z_meta  := (z_curr + z_other) / sqrt(2)]
    meta[, p_meta  := 2 * pnorm(abs(z_meta), lower.tail=FALSE)]
    meta[, FDR_meta := p.adjust(p_meta, method="BH")]
    
    # Robust tissue name extraction using regex on target path fields
    # Grabs names like Cerebellum, Lateral_Ventricle, etc.
    tissue_match <- regexpr("target_ALS/[A-Za-z0-9_]+", out_prefix)
    if(tissue_match != -1) {
        tissue_raw <- regmatches(out_prefix, tissue_match)
        tissue_name <- gsub("target_ALS/", "", tissue_raw)
    } else {
        tissue_name <- "Unknown_Tissue"
    }
    
    meta_out <- file.path(combined_dir, paste0(tissue_name, "_Combined_Meta_sQTL.txt"))
    fwrite(meta[order(p_meta)], file=meta_out, sep="\t", quote=FALSE)
    message(sprintf("## Meta-analysis saved to: %s", meta_out))
}

message("## Post-processing complete. Ready for plots.")
